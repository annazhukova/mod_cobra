from collections import defaultdict
import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
import libsbml

from mod_sbml.annotation.chebi.chebi_annotator import EQUIVALENT_RELATIONSHIPS
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from mod_cobra.mapping.metabolite_matcher import map_metabolites_compartments, get_model_data
from mod_cobra.mapping.model_merge_manager import join, merge
from mod_cobra.mapping.model_merger import simple_merge_models
from mod_cobra.efm.community_detector import detect_fm_communities, detect_reaction_community
from mod_cobra.efm.efm_pipeline import analyse_model_efm
from mod_cobra.html import describe
from mod_sbml.serialization.csv_manager import serialize_model_info, serialize_common_elements_to_csv
from mod_cobra.efm.serialization.coupled_reaction_sbml_manager import create_folded_sbml
from mod_sbml.onto import parse_simple
from sbml_vis.graph.color.color import get_n_colors
from mod_cobra import ZERO_THRESHOLD
from mod_cobra.fbva.serialization.fbva_serializer import serialize_fva, serialize_fluxes
from mod_cobra.fbva.serialization import format_r_id
from mod_cobra.sbml import r_ids2sbml
from mod_sbml.sbml.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, reverse_reaction, get_r_comps
from sbml_vis.mimoza_pipeline import process_sbml
from mod_cobra.fbva.fbva_analyser import analyse_by_fba, analyse_by_fva
from mod_cobra.fbva.serialization.fbva_sbml_manager import create_fva_model
from mod_sbml.utils.path_manager import create_dirs
from mod_sbml.serialization import get_sbml_r_formula
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, \
    select_metabolite_ids_by_term_ids
from mod_cobra.efm.serialization import coupled_reaction_group_serializer, community_serializer, efm_serializer

CHEBI_H = 'chebi:15378'
TREEEFM_PATH = "/home/anna/Applications/TreeEFM/tool/TreeEFMseq"

__author__ = 'anna'


def constraint_exchange_reactions(model, allowed_exchange_r_id2rev, cofactors=None, min_flux=0.01):
    if not cofactors:
        ub_ch_ids = get_ubiquitous_chebi_ids(add_common=True, add_cofactors=True)
        cofactors = select_metabolite_ids_by_term_ids(model, ub_ch_ids)

    for r in model.getListOfReactions():
        if r.id in allowed_exchange_r_id2rev:
            rev = allowed_exchange_r_id2rev[r.id]
            if rev:
                # reverse the reaction and set positive bounds,
                # as if both bounds are negative, the glp solver produces an error
                l_b, u_b = get_bounds(r)
                reverse_reaction(r)
                set_bounds(r, -u_b, -l_b)
                allowed_exchange_r_id2rev[r.id] = not rev
            r_l, r_u = get_bounds(r)
            set_bounds(r, max(min_flux, r_l), max(min_flux, r_u))
            continue
        rs, ps = set(get_reactants(r)), set(get_products(r))

        boundary_s_ids = {s_id for s_id in rs if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not rs:
            r_l, r_u = get_bounds(r)
            # if it's not only cofactors, constrain it
            set_bounds(r, min(r_l, 0), 0 if (ps - cofactors) else r_u)
            continue
        boundary_s_ids = {s_id for s_id in ps if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not ps:
            r_l, r_u = get_bounds(r)
            # if it's not only cofactors, constrain it
            set_bounds(r, 0 if (rs - cofactors) else r_l, max(r_u, 0))


def _prepare_dir(super_dir, dir_name, log):
    logging.info(log)
    cur_dir = os.path.join(super_dir, dir_name)
    create_dirs(cur_dir)
    return cur_dir


def analyse_model(sbml, out_r_id, out_rev, res_dir, in_m_id, out_m_id, in_r_id2rev=None, threshold=ZERO_THRESHOLD,
                  do_fva=True, do_fba=True, do_efm=True, max_efm_number=1000, mask_shift=4,
                  get_f_path=None, tree_efm_path=TREEEFM_PATH, model_name='', main_dir=None):
    # create directories to store results
    logging.info("Preparing directories...")
    create_dirs(res_dir, False)
    if not get_f_path:
        get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir))

    # if in_r_id2rev:
    #     logging.info("Constraining input reactions...")
    #     doc = libsbml.SBMLReader().readSBML(sbml)
    #     model = doc.getModel()
    #     constraint_exchange_reactions(model, allowed_exchange_r_id2rev=in_r_id2rev, cofactors=cofactors)
    #     sbml = os.path.join(res_dir, '%s_constrained.xml' % os.path.splitext(os.path.basename(sbml))[0])
    #     libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)

    # copy our model in the result directory
    if os.path.normpath(res_dir) != os.path.normpath(os.path.dirname(sbml)):
        shutil.copy(sbml, res_dir)
        sbml = os.path.join(res_dir, os.path.basename(sbml))

    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    logging.info("Serializing model info...")
    info_prefix = os.path.join(res_dir, 'model_info_')
    c_csv, m_csv, r_csv = serialize_model_info(model, info_prefix)

    r_string = lambda r, rev: '<b>%s</b>%s: %s' % (r.getId(), ' (reversed)' if rev else '',
                                                   get_sbml_r_formula(model, r, show_metabolite_ids=False))
    r = model.getReaction(out_r_id)
    description = describe('input_data.html',
                           model_name=(model.getName() if model.getName() else model.getId()),
                           sbml_filepath=get_f_path(sbml), c_csv=get_f_path(c_csv),
                           m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv), obj_rn=r_string(r, out_rev),
                           in_rn_len=len(in_r_id2rev) if in_r_id2rev else 0,
                           in_rns='; '.join(r_string(model.getReaction(r_id), rev)
                                            for (r_id, rev) in in_r_id2rev.iteritems()) if in_r_id2rev else '')

    r_id2mask = defaultdict(lambda: 0)
    layer2mask = {}
    main_layer = None
    objective_sense = 'minimize' if out_rev else 'maximize'
    vis_r_ids = set()
    cobra_model = None
    opt_val = None
    if do_fva:
        cur_dir = _prepare_dir(res_dir, 'fva', "Performing FVA...")
        cobra_model = create_cobra_model_from_sbml_file(sbml) if not cobra_model else cobra_model
        r_id2bounds, opt_val = analyse_by_fva(cobra_model=cobra_model, bm_r_id=out_r_id,
                                              objective_sense=objective_sense, threshold=threshold)
        ess_rn_num, var_rn_num, fva_file = 0, 0, None
        if opt_val:
            fva_file = os.path.join(cur_dir, 'fva.txt')
            ess_rn_num, var_rn_num = serialize_fva(cobra_model, r_id2bounds, fva_file, objective_sense, out_r_id)
            mask_shift = update_vis_layers(
                {r_id: (1 if u > 0 else -1) for (r_id, (l, u)) in r_id2bounds.iteritems() if l * u > 0},
                'FVA essential', r_id2mask, layer2mask, mask_shift, vis_r_ids)
            main_layer = 'FVA essential'
            fva_sbml = os.path.join(cur_dir, 'Model_FVA.xml')
            sbml = create_fva_model(sbml, r_id2bounds, fva_sbml)
            doc = libsbml.SBMLReader().readSBML(sbml)
            model = doc.getModel()
        description += describe('fva.html', optimal_value=opt_val, ess_r_num=ess_rn_num, var_r_num=var_rn_num,
                                description_filepath=get_f_path(fva_file),
                                sbml_filepath=get_f_path(sbml))

    if do_fba:
        cur_dir = _prepare_dir(res_dir, 'fba', "Performing FBA...")
        cobra_model = create_cobra_model_from_sbml_file(sbml) if not cobra_model else cobra_model
        r_id2val, opt_val = analyse_by_fba(cobra_model, bm_r_id=out_r_id, objective_sense=objective_sense,
                                           threshold=threshold)
        fba_file = None
        if opt_val:
            mask_shift = update_vis_layers(r_id2val, 'FBA', r_id2mask, layer2mask, mask_shift, vis_r_ids)
            main_layer = 'FBA'
            fba_file = os.path.join(cur_dir, 'fba.txt')
            serialize_fluxes(cobra_model, r_id2val, path=fba_file, objective_sense=objective_sense, out_r_id=out_r_id)
        description += describe('fba.html', optimal_value=opt_val, r_num=len(r_id2val),
                                description_filepath=get_f_path(fba_file))

    S, r_id2w = None, {}
    if do_efm:
        cur_dir = _prepare_dir(res_dir, 'efma', "Performing EFMA...")
        S, r_id2w = analyse_model_efm(sbml, out_r_id, out_rev, cur_dir, in_m_id, out_m_id, in_r_id2rev, tree_efm_path,
                                      threshold, max_efm_number)

        for serializer in (efm_serializer.serialize, coupled_reaction_group_serializer.serialize):
            description += \
                serializer(model=model, path=cur_dir, get_f_path=get_f_path, in_m_id=in_m_id, out_m_id=out_m_id,
                           out_r_id=out_r_id, S=S, model_name=model_name, main_dir=main_dir)
        logging.info('Serialized EFMA')

        if S.gr_id2r_id2c:
            clique_merged_sbml = os.path.join(cur_dir, 'Model_folded.xml')
            r_id2new_r_id = create_folded_sbml(S, sbml, clique_merged_sbml)
            sbml = clique_merged_sbml
            doc = libsbml.SBMLReader().readSBML(sbml)
            model = doc.getModel()

            vis_r_ids |= {cl_id for (r_id, cl_id) in r_id2new_r_id.iteritems() if r_id in vis_r_ids}
            for r_id, new_r_id in r_id2new_r_id.iteritems():
                if r_id in r_id2mask:
                    r_id2mask[new_r_id] |= r_id2mask[r_id]

    if not opt_val and not S.efm_id2i:
        description += describe('nothing_found.html')

    return S, r_id2w, sbml, vis_r_ids, description, mask_shift, r_id2mask, layer2mask, main_layer


def mean(values):
    return sum(values) / len(values)


def multimodel_pipeline(sbml2parameters, res_dir, do_fva=True, do_fba=True, do_efm=True, max_efm_number=1000):
    vis_dir = os.path.join(res_dir, 'visualization')
    create_dirs(vis_dir, True)
    get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir)) if f else None

    mask_shift = 4
    tab2html = {}
    layer2mask = {}
    invisible_layers = []
    model_id2id2mask = {}
    model_id2vis_r_ids = {}
    model_id2sbml = {}
    model_id2S = {}
    model_id2r_id2w = {}
    for model_id, (out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id, name) in sbml2parameters.iteritems():
        logging.info('Analysing %s...' % name)
        S, r_id2w, sub_sbml, r_ids, description, mask_shift, cur_id2mask, cur_layer2mask, main_layer = \
            analyse_model(sbml=model_id, out_r_id=out_r_id, out_rev=out_rev, in_r_id2rev=in_r_id2rev, in_m_id=in_m_id,
                          out_m_id=out_m_id, res_dir=os.path.join(res_dir, name), do_fva=do_fva, do_fba=do_fba,
                          do_efm=do_efm, mask_shift=mask_shift, get_f_path=get_f_path, max_efm_number=max_efm_number,
                          model_name=name, main_dir=vis_dir)

        tab2html['Analysis of %s' % name] = description, None
        layer2mask.update({'%s: %s' % (name, layer): mask for (layer, mask) in cur_layer2mask.iteritems()})
        invisible_layers.extend(['%s: %s' % (name, layer) for layer in cur_layer2mask.iterkeys() if layer != main_layer])
        model_id2id2mask[name] = cur_id2mask
        model_id2vis_r_ids[name] = r_ids
        model_id2sbml[name] = sub_sbml
        model_id2S[name] = S
        model_id2r_id2w[name] = r_id2w

    if len(model_id2sbml) > 1:
        mm_dir = os.path.join(res_dir, 'merged_model')
        create_dirs(mm_dir)

        logging.info('Going to merge models...')
        chebi = parse_simple(get_chebi())
        model_id2dfs = get_model_data(model_id2sbml, chebi=chebi)

        model_id2c_id_groups, model_id2m_id_groups, model_id2c_id2i = map_metabolites_compartments(model_id2dfs)
        logging.info('Mapped metabolites and compartments.')
        ignore_m_ids = get_ignored_metabolites(model_id2dfs, chebi)
        S = join(model_id2m_id_groups, model_id2S)
        ignore_m_ids |= {S.m_id2gr_id[m_id] for m_id in ignore_m_ids if m_id in S.m_id2gr_id}
        S = merge(S, ignore_m_ids)
        model_id2r_id_groups = get_r_id_groups(S)
        logging.info('Mapped reactions.')
        comp_csv, m_csv, r_csv = serialize_common_elements_to_csv(model_id2dfs, model_id2c_id_groups,
                                                                  model_id2m_id_groups, model_id2r_id_groups,
                                                                  os.path.join(mm_dir, 'Model_comparison_')) \
            if model_id2c_id_groups else (None, None, None)
        logging.info('Serialized the mappings.')

        tab2html['Model comparison'] = \
            describe('model_comparison.html',
                     c_num=len(model_id2c_id_groups), m_num=len(model_id2m_id_groups), r_num=len(model_id2r_id_groups),
                     c_csv=get_f_path(comp_csv), m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv)), None
        title = 'Combined model analysis'

        merged_sbml = os.path.join(res_dir, 'Merged_model.xml')
        model_id2id2id, common_ids, S = simple_merge_models(S, model_id2c_id2i, model_id2dfs, merged_sbml)
        r_id2w = defaultdict(lambda: 1)
        for model_id, r_id2w in model_id2r_id2w.iteritems():
            r_id2w.update({model_id2id2id[model_id][r_id]: w for (r_id, w) in r_id2w.iteritems()
                           if r_id in model_id2id2id[model_id]})

        id2mask = defaultdict(lambda: 0)
        vis_r_ids = set()
        for model_id in model_id2sbml.iterkeys():
            id2mask.update({model_id2id2id[model_id][s_id]: id2mask[model_id2id2id[model_id][s_id]] | mask
                            for (s_id, mask) in model_id2id2mask[model_id].iteritems() if s_id in model_id2id2id[model_id]})
            vis_r_ids |= {model_id2id2id[model_id][r_id] for r_id in model_id2vis_r_ids[model_id]
                          if r_id in model_id2id2id[model_id]}

        id2color, info = get_colors(common_ids, model_id2id2id, model_id2sbml)
        sbml = merged_sbml
    else:
        model_id, sbml = next(model_id2sbml.iteritems())
        vis_r_ids = model_id2vis_r_ids[model_id]
        id2mask = model_id2id2mask[model_id]
        id2color = None
        r_id2w = model_id2r_id2w[model_id]
        info = ''
        S = model_id2S[model_id].get_main_S()
        title = 'Model analysis'

    # Communities
    comm_dir = _prepare_dir(res_dir, 'communities', "Analysing communities...")
    id2cluster = detect_fm_communities(S, r_id2w)
    id2intersection = {cl_id: S.get_efm_intersection(cluster) for (cl_id, cluster) in id2cluster.iteritems()}
    id2imp_rns = {cl_id: detect_reaction_community(S, cluster, id2intersection[cl_id])
                  for (cl_id, cluster) in id2cluster.iteritems()}
    if id2cluster:
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        description = \
            community_serializer.serialize(model, S, id2cluster, id2intersection, id2imp_rns, comm_dir, get_f_path)
        if len(model_id2sbml) > 1:
            tab2html['Model comparison'] = tab2html['Model comparison'][0] + description, None
        else:
            tab2html['Pathway communities'] = description, None
        for cl_id, r_id2c in id2imp_rns.iteritems():
            mask_shift = update_vis_layers(r_id2c, 'Pathway community %s' % cl_id, id2mask, layer2mask, mask_shift,
                                           vis_r_ids)
            invisible_layers.append('Pathway community %s' % cl_id)

    visualize_model(sbml, vis_r_ids, id2mask, id2color, title, info, invisible_layers, layer2mask, res_dir, tab2html)


def get_colors(common_ids, model_id2id2id, model_id2sbml):
    colors = get_n_colors(len(model_id2sbml) + 1, 0.5, 0.8)
    mixed_color = None
    info = ''
    if common_ids:
        mixed_color = colors[0]
        r, g, b = mixed_color
        info = describe('color.html', r=r, g=g, b=b, name='common reactions/metabolites')
    model_id2color = dict(zip(model_id2sbml.iterkeys(), colors[1:]))
    id2color = {}
    for model_id, id2id in model_id2id2id.iteritems():
        color = model_id2color[model_id]
        r, g, b = color
        id2color.update({t_id: color if t_id not in common_ids else mixed_color for t_id in id2id.itervalues()})
        info += describe('color.html', r=r, g=g, b=b, name=model_id)
    return id2color, info


def get_r_id_groups(S_merged):
    model_id2r_id_groups = []
    for r_id2c in S_merged.gr_id2r_id2c.itervalues():
        model_id2r_ids = defaultdict(set)
        for (model_id, r_id) in r_id2c.iterkeys():
            model_id2r_ids[model_id].add(r_id)
        model_id2r_id_groups.append(model_id2r_ids)
    return model_id2r_id_groups


def get_ignored_metabolites(model_id2dfs, chebi):
    ignore_ch_ids = chebi.get_surroundings(chebi.get_term(CHEBI_H), relationships=EQUIVALENT_RELATIONSHIPS, radius=1)
    ignore_m_ids = set()
    for model_id, [df, _, _] in model_id2dfs.iteritems():
        for index, row in df.iterrows():
            m_id = row['Id']
            t_id = row["ChEBI"] if 'ChEBI' in row else None
            if t_id in ignore_ch_ids:
                ignore_m_ids.add((model_id, m_id))
    return ignore_m_ids


def visualize_model(sbml, vis_r_ids, id2mask, id2color, title, info, invisible_layers, layer2mask, res_dir, tab2html):
    combined_sbml = os.path.join(res_dir, 'Combined_model.xml')
    r_ids2sbml(vis_r_ids, sbml, combined_sbml, 'combined')
    doc = libsbml.SBMLReader().readSBML(combined_sbml)
    model = doc.getModel()
    id2mask = get_full_id2mask(id2mask, model)
    process_sbml(combined_sbml, verbose=True, path='visualization', generalize=False,
                 id2mask=id2mask, layer2mask=layer2mask, tab2html=tab2html, title=title,
                 id2color=id2color, tabs=None, info=info, invisible_layers=invisible_layers)


def update_vis_layers(r_id2val, layer, id2mask, layer2mask, mask_shift, vis_r_ids):
    if not r_id2val:
        return mask_shift
    l_mask = 1 << mask_shift
    layer2mask[layer] = l_mask
    for r_id in r_id2val.iterkeys():
        o_r_id = format_r_id(r_id)
        vis_r_ids |= {r_id, o_r_id}
        id2mask[r_id] |= l_mask
        id2mask[o_r_id] |= l_mask
    return mask_shift + 1


def get_full_id2mask(r_id2mask, model):
    id2mask = defaultdict(lambda: 0)
    for r_id, mask in r_id2mask.iteritems():
        r = model.getReaction(r_id)
        if r:
            id2mask[r_id] |= mask
            id2mask.update({c_id: id2mask[c_id] | mask for c_id in get_r_comps(r_id, model)})
            id2mask.update({s_id: id2mask[s_id] | mask for s_id in get_reactants(r)})
            id2mask.update({s_id: id2mask[s_id] | mask for s_id in get_products(r)})
    return id2mask


def process_n_fms(directory, id2fm, id2mask, layer2mask, mask_shift, vis_r_ids,
                  name, sorter, serializer, suffix, get_f_path, limit, update_vis=True):
    sbml_dir = os.path.join(directory, 'sbml')
    create_dirs(sbml_dir)

    fms = []
    for fm_id, fm in sorted(id2fm.iteritems(), key=sorter):
        r_id2coeff = fm.to_r_id2coeff()
        c_name = name[0].upper() + name[1:]
        fm_name = fm.id if fm.id else '%s %s' % (c_name, str(fm_id))
        fm_txt = os.path.join(sbml_dir, '%s.txt' % fm_name.replace(' ', '_'))
        serializer(fm_id, fm_txt)
        fms.append((c_name, fm_id, len(fm), suffix(fm_id), get_f_path(fm_txt)))

        if update_vis:
            mask_shift = update_vis_layers(r_id2coeff, fm_name, id2mask, layer2mask, mask_shift, vis_r_ids)

        limit -= 1
        if limit <= 0:
            break
    return mask_shift, fms

