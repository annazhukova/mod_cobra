from collections import defaultdict
import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
import libsbml
from metabolite_matcher import map_metabolites_compartments, get_model_data

from mod_cobra.constraint_based_analysis.efm.efm_pipeline import analyse_model_efm
from mod_cobra.constraint_based_analysis.html_descriptor import describe
from mod_sbml.serialization.csv_manager import serialize_common_part_to_csv, serialize_model_info, \
    serialize_common_metabolites_compartments_to_csv
from mod_cobra.constraint_based_analysis.efm.sbml.copled_reaction_sbml_manager import create_folded_sbml
from model_merge_manager import join, merge
from sbml_vis.graph.color.color import get_n_colors
from mod_cobra.constraint_based_analysis import ZERO_THRESHOLD
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import serialize_fva, \
    serialize_fluxes
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from mod_cobra.constraint_based_analysis.efm.efm_serialization_manager import r_ids2sbml
from mod_cobra.gibbs.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, reverse_reaction, get_r_comps
from mimoza_pipeline import process_sbml
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fba_analyser import analyse_by_fba
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva, \
    create_fva_model
from mod_cobra.constraint_based_analysis.efm.efm_analyser import TREEEFM_PATH
from mod_sbml.utils.path_manager import create_dirs
from mod_sbml.serialization.serialization_manager import get_sbml_r_formula
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, \
    select_metabolite_ids_by_term_ids
from model_comparison.model_merger import merge_models, simple_merge_models
from mod_cobra.constraint_based_analysis.efm.serialization import numpy_efm_serialization_manager, \
    numpy_coupled_reaction_group_serialization_manager

CHEBI_H = 'chebi:15378'

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
                  do_fva=True, do_fba=True, do_efm=True, max_efm_number=1000, cofactors=None, mask_shift=4,
                  get_f_path=None, tree_efm_path=TREEEFM_PATH):
    # create directories to store results
    logging.info("Preparing directories...")
    create_dirs(res_dir, False)
    if not get_f_path:
        get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir))

    if in_r_id2rev:
        logging.info("Constraining input reactions...")
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        constraint_exchange_reactions(model, allowed_exchange_r_id2rev=in_r_id2rev, cofactors=cofactors)
        sbml = os.path.join(res_dir, '%s_constrained.xml' % os.path.splitext(os.path.basename(sbml))[0])
        libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)

    # copy our model in the result directory
    if os.path.normpath(res_dir) != os.path.normpath(os.path.dirname(sbml)):
        shutil.copy(sbml, res_dir)
        sbml = os.path.join(res_dir, os.path.basename(sbml))

    logging.info("Serializing model info...")
    info_prefix = os.path.join(res_dir, 'model_info_')
    c_csv, m_csv, r_csv = serialize_model_info(sbml, info_prefix)
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    r_string = lambda r, rev: '<b>%s</b>%s: %s' % (r.getId(), ' (reversed)' if rev else '',
                                                   get_sbml_r_formula(model, r, show_metabolite_ids=False))
    r = model.getReaction(out_r_id)
    description = describe('input_data.html', model_name=model.getName() if model.getName() else model.getId(),
                           sbml_filepath=get_f_path(sbml), c_csv=get_f_path(c_csv),
                           m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv), obj_rn=r_string(r, out_rev),
                           in_rn_len=len(in_r_id2rev) if in_r_id2rev else 0,
                           in_rns=
                           '; '.join(r_string(model.getReaction(r_id), rev) for (r_id, rev) in in_r_id2rev.iteritems())
                           if in_r_id2rev else '')

    r_id2mask = defaultdict(lambda: 0)
    layer2mask = {}
    main_layer = None
    objective_sense = 'minimize' if out_rev else 'maximize'
    vis_r_ids = set()
    cobra_model = None
    if do_fva:
        cur_dir = _prepare_dir(res_dir, 'fva', "Performing FVA...")
        cobra_model = create_cobra_model_from_sbml_file(sbml) if not cobra_model else cobra_model
        r_id2bounds, opt_val = analyse_by_fva(cobra_model=cobra_model, bm_r_id=out_r_id,
                                              objective_sense=objective_sense, threshold=threshold)
        ess_rn_num, var_rn_num = 0, 0
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
                                description_filepath=get_f_path(fva_file) if opt_val else None,
                                sbml_filepath=get_f_path(fva_sbml) if opt_val else None)

    if do_fba:
        cur_dir = _prepare_dir(res_dir, 'fba', "Performing FBA...")
        cobra_model = create_cobra_model_from_sbml_file(sbml) if not cobra_model else cobra_model
        r_id2val, opt_val = analyse_by_fba(cobra_model, bm_r_id=out_r_id, objective_sense=objective_sense,
                                           threshold=threshold)
        if opt_val:
            mask_shift = update_vis_layers(r_id2val, 'FBA', r_id2mask, layer2mask, mask_shift, vis_r_ids)
            main_layer = 'FBA'
            fba_file = os.path.join(cur_dir, 'fba.txt')
            serialize_fluxes(cobra_model, r_id2val, path=fba_file, objective_sense=objective_sense, out_r_id=out_r_id)
        description += describe('fba.html', optimal_value=opt_val, r_num=len(r_id2val),
                                description_filepath=get_f_path(fba_file) if opt_val else None)

    S = None
    if do_efm:
        cur_dir = _prepare_dir(res_dir, 'efma', "Performing EFMA...")
        S_initial, S_folded, S_folded_wo_duplicates, S_merged = \
            analyse_model_efm(sbml, out_r_id, out_rev, cur_dir, in_m_id, out_m_id, in_r_id2rev,
                              threshold, tree_efm_path, max_efm_number)
        S = S_merged

        for serialize in (numpy_efm_serialization_manager.serialize,
                          numpy_coupled_reaction_group_serialization_manager.serialize):
            description += \
                serialize(model=model, path=res_dir, get_f_path=get_f_path, in_m_id=in_m_id, out_m_id=out_m_id,
                          out_r_id=out_r_id, S_initial=S_initial, S_coupled=S_folded,
                          S_no_duplicates=S_folded_wo_duplicates, S_merged=S_merged)

        if S_folded.gr_id2r_id2c or S_folded_wo_duplicates.gr_id2r_id2c:
            clique_merged_sbml = os.path.join(cur_dir, 'Model_folded.xml')
            r_id2new_r_id = create_folded_sbml(S_folded, S_folded_wo_duplicates, sbml, clique_merged_sbml)
            sbml = clique_merged_sbml

            vis_r_ids |= {cl_id for (r_id, cl_id) in r_id2new_r_id.iteritems() if r_id in vis_r_ids}
            for r_id, new_r_id in r_id2new_r_id.iteritems():
                if r_id in r_id2mask:
                    r_id2mask[new_r_id] |= r_id2mask[r_id]

        # all_fm_intersection, description, efm_num, fm_block, fm_num, id2fm, mask_shift, model, sbml = efm_analysis(
        #     cur_dir, description, efms, get_f_path, in_m_id, in_r_id2rev, layer2mask, mask_shift, max_efm_number, model,
        #     out_m_id, out_r_id, out_rev, r_id2mask, res_dir, sbml, threshold, vis_r_ids)
        # if efm_num:
        #     # Communities
        #     comm_dir = _prepare_dir(cur_dir, 'communities', "Analysing communities...")
        #     id2cluster = detect_fm_communities(id2fm, threshold=len(all_fm_intersection))
        #     id2intersection = {cl_id: reduce(lambda p1, p2: p1.intersection(p2),
        #                                      (id2fm[fm_id] for fm_id in cluster), id2fm[cluster[0]])
        #                        for (cl_id, cluster) in id2cluster.iteritems()}
        #     imp_fraction = 60
        #     r_ids = [r.getId() for r in model.getListOfReactions()]
        #     rev_r_ids = [r.getId() for r in model.getListOfReactions() if r.getReversible()]
        #     id2imp_rns = {cl_id: detect_reaction_community(r_ids, rev_r_ids, {fm_id: id2fm[fm_id] for fm_id in cluster},
        #                                                    imp_fraction, id2intersection[cl_id])
        #                   for (cl_id, cluster) in id2cluster.iteritems()}
        #     if id2cluster:
        #         comm_txt = os.path.join(comm_dir, 'communities.txt')
        #         serialize_communities(id2cluster, id2intersection, id2imp_rns, fm_num, all_fm_intersection, comm_txt)
        #         limit = len(id2cluster)
        #         mask_shift, fms \
        #             = process_n_fms(directory=comm_dir, id2fm=id2imp_rns, id2mask=r_id2mask,
        #                             layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
        #                             name='community', sorter=lambda (cl_id, cl): cl_id,
        #                             serializer=lambda cl_id, community_txt:
        #                             serialize_community(cl_id, id2intersection[cl_id], id2imp_rns[cl_id],
        #                                                 id2cluster[cl_id], model, community_txt),
        #                             suffix=lambda _: '', limit=limit, get_f_path=get_f_path)
        #         fm_block = describe('fm_block.html', element_num=limit, characteristics='(all)',
        #                             element_name='community', fms=fms, all=True)
        #     description += describe('communities.html', community_num=len(id2cluster),
        #                             description_filepath=get_f_path(comm_txt) if id2cluster else None,
        #                             selected_community_block=fm_block)

    if not vis_r_ids:
        description += describe('nothing_found.html')

    return S, sbml, vis_r_ids, description, mask_shift, r_id2mask, layer2mask, main_layer


def mean(values):
    return sum(values) / len(values)


def multimodel_pipeline(sbml2parameters, res_dir, do_fva=True, do_fba=True, do_efm=True, max_efm_number=1000):
    create_dirs(res_dir, True)
    get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir)) if f else None

    mask_shift = 4
    tab2html = {}
    layer2mask = {}
    invisible_layers = []
    model_id2id2mask = {}
    model_id2vis_r_ids = {}
    sbml2name = {}
    model_id2sbml = {}
    model_id2S = {}
    for model_id, (out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id, name) in sbml2parameters.iteritems():
        logging.info('Analysing %s...' % name)
        S, sub_sbml, r_ids, description, mask_shift, cur_id2mask, cur_layer2mask, main_layer = \
            analyse_model(sbml=model_id, out_r_id=out_r_id, out_rev=out_rev, in_r_id2rev=in_r_id2rev, in_m_id=in_m_id,
                          out_m_id=out_m_id, res_dir=os.path.join(res_dir, name), do_fva=do_fva, do_fba=do_fba,
                          do_efm=do_efm, mask_shift=mask_shift, get_f_path=get_f_path, max_efm_number=max_efm_number)

        tab2html['Analysis of %s' % name] = description, None
        layer2mask.update({'%s: %s' % (name, layer): mask for (layer, mask) in cur_layer2mask.iteritems()})
        invisible_layers.extend(['%s: %s' % (name, layer) for layer in cur_layer2mask.iterkeys() if layer != main_layer])
        model_id2id2mask[name] = cur_id2mask
        model_id2vis_r_ids[name] = r_ids
        sbml2name[sub_sbml] = name
        model_id2sbml[name] = sub_sbml
        model_id2S[name] = S

    if len(model_id2sbml) > 1:
        mm_dir = os.path.join(res_dir, 'merged_model')
        create_dirs(mm_dir)

        model_id2chebi_id2m_ids, model_id2c_id2data, model_id2m_id2data, model_id2r_id2data = \
            get_model_data(model_id2sbml)

        model_id2c_id_groups, model_id2m_id_groups, model_id2c_id2i = \
            map_metabolites_compartments(model_id2chebi_id2m_ids, model_id2c_id2data, model_id2m_id2data)
        model_id2r_id_groups = []
        ignore_m_ids = set()
        for model_id in model_id2S.iterkeys():
            if CHEBI_H in model_id2chebi_id2m_ids[model_id]:
                for m_id in model_id2chebi_id2m_ids[model_id][CHEBI_H]:
                    ignore_m_ids.add((model_id, m_id))
        S_merged = join(model_id2m_id_groups, model_id2S).remove_unused_metabolites()
        ignore_m_ids |= {S_merged.m_id2gr_id[m_id] for m_id in ignore_m_ids if m_id in S_merged.m_id2gr_id}
        S_merged = merge(S_merged, ignore_m_ids)
        for r_id2c in S_merged.gr_id2r_id2c.itervalues():
            model_id2r_ids = defaultdict(set)
            for (model_id, r_id) in r_id2c.iterkeys():
                model_id2r_ids[model_id].add(r_id)
            model_id2r_id_groups.append(model_id2r_ids)
        if model_id2c_id_groups or model_id2m_id_groups:
            comparison_prefix = os.path.join(mm_dir, 'Model_comparison_')
            comp_csv, m_csv, r_csv = \
                serialize_common_metabolites_compartments_to_csv(model_id2sbml, model_id2c_id_groups,
                                                                 model_id2m_id_groups, model_id2r_id_groups,
                                                                 comparison_prefix)
        else:
            comp_csv, m_csv, r_csv = None, None, None

        tab2html['Model comparison'] = \
            describe('model_comparison.html',
                     c_num=len(model_id2c_id_groups), m_num=len(model_id2m_id_groups), r_num=len(model_id2r_id_groups),
                     c_csv=get_f_path(comp_csv), m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv)), None
        title = 'Combined model analysis'

        merged_sbml = os.path.join(res_dir, 'Merged_model.xml')
        model_id2id2id, common_ids = \
            simple_merge_models(S_merged, model_id2c_id2i, model_id2c_id2data, model_id2m_id2data, model_id2r_id2data,
                                merged_sbml)

        id2mask = defaultdict(lambda: 0)
        vis_r_ids = set()
        for model_id in model_id2sbml.iterkeys():
            id2mask.update({model_id2id2id[model_id][s_id]: id2mask[model_id2id2id[model_id][s_id]] | mask
                            for (s_id, mask) in model_id2id2mask[model_id].iteritems() if s_id in model_id2id2id[model_id]})
            vis_r_ids |= {model_id2id2id[model_id][r_id] for r_id in model_id2vis_r_ids[model_id]
                          if r_id in model_id2id2id[model_id]}

        colors = get_n_colors(len(model_id2sbml) + 1, 0.5, 0.8)
        info = ''
        mixed_color = None
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

        sbml = merged_sbml
    else:
        model_id, sbml = next(model_id2sbml.iteritems())
        vis_r_ids = model_id2vis_r_ids[model_id]
        id2mask = model_id2id2mask[model_id]
        id2color = None
        info = ''
        title = 'Model analysis'

    visualize_model(sbml, vis_r_ids, id2mask, id2color, title, info, invisible_layers, layer2mask, res_dir, tab2html)


def visualize_model(sbml, vis_r_ids, id2mask, id2color, title, info, invisible_layers, layer2mask, res_dir, tab2html):
    combined_sbml = os.path.join(res_dir, 'Combined_model.xml')
    r_ids2sbml(vis_r_ids, sbml, combined_sbml, 'combined')
    doc = libsbml.SBMLReader().readSBML(combined_sbml)
    model = doc.getModel()
    id2mask = get_full_id2mask(id2mask, model)
    process_sbml(combined_sbml, verbose=True, path='visualization', generalize=False,
                 id2mask=id2mask, layer2mask=layer2mask, tab2html=tab2html, title=title,
                 id2color=id2color, tabs=None, info=info, invisible_layers=invisible_layers)


def join_models(model_id2sbml, get_f_path, res_dir, sbml2id2mask, sbml2name, sbml2vis_r_ids, sbml_files, tab2html):
    # merged_sbml = os.path.join(res_dir, 'Merged_model.xml')

    model_id2c_id_groups, model_id2m_id_groups = map_metabolites_compartments(model_id2sbml)


    # sbml2id2id, sbml2rev_r_ids, common_ids = merge_models(sbml_files, merged_sbml, sbml2name)

    if model_id2c_id_groups or model_id2m_id_groups:
        comparison_prefix = os.path.join(res_dir, 'Model_comparison_')
        (cc_num, cm_num), (comp_csv, m_csv) = \
            serialize_common_metabolites_compartments_to_csv(model_id2sbml,
                                                             model_id2c_id_groups, model_id2m_id_groups,
                                                             comparison_prefix)
    else:
        cc_num, cm_num = 0, 0
    tab2html['Model comparison'] = \
        describe('model_comparison.html', c_num=cc_num, m_num=cm_num, r_num=0,
                 c_csv=get_f_path(comp_csv) if cc_num else None,
                 m_csv=get_f_path(m_csv) if cm_num else None,
                 r_csv=None), None
    # id2mask = defaultdict(lambda: 0)
    # vis_r_ids = set()
    # for sbml in sbml_files:
    #     id2mask.update({sbml2id2id[sbml][s_id]: id2mask[sbml2id2id[sbml][s_id]] | mask
    #                     for (s_id, mask) in sbml2id2mask[sbml].iteritems() if s_id in sbml2id2id[sbml]})
    #     vis_r_ids |= {sbml2id2id[sbml][r_id] for r_id in sbml2vis_r_ids[sbml] if r_id in sbml2id2id[sbml]}
    #
    # colors = get_n_colors(len(sbml_files) + 1, 0.5, 0.8)
    # info = ''
    # if common_ids:
    #     mixed_color = colors[0]
    #     r, g, b = mixed_color
    #     info = describe('color.html', r=r, g=g, b=b, name='common reactions/metabolites')
    # sbml2color = dict(zip(sbml2id2id.iterkeys(), colors[1:]))
    # id2color = {}
    # for sbml, id2id in sbml2id2id.iteritems():
    #     color = sbml2color[sbml]
    #     r, g, b = color
    #     id2color.update({t_id: color if t_id not in common_ids else mixed_color for t_id in id2id.itervalues()})
    #     info += describe('color.html', r=r, g=g, b=b, name=sbml2name[sbml])
    # return id2color, id2mask, info, merged_sbml, vis_r_ids


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

