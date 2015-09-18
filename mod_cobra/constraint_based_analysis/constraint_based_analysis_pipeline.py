from collections import defaultdict
import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
import libsbml

from mod_cobra.constraint_based_analysis.efm.clique_detection import detect_cliques, clique2lumped_reaction
from mod_cobra.constraint_based_analysis.efm.community_detection import detect_fm_communities, detect_reaction_community
from mod_cobra.constraint_based_analysis.html_descriptor import describe
from mod_sbml.serialization.csv_manager import serialize_common_part_to_csv, serialize_model_info
from sbml_vis.graph.color.color import get_n_colors
from mod_cobra.constraint_based_analysis import ZERO_THRESHOLD
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import serialize_fva, \
    serialize_fluxes
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from mod_cobra.constraint_based_analysis.efm.efm_serialization_manager import r_ids2sbml, serialize_fm, \
    serialize_communities, serialize_community, serialize_cliques, serialize_clique, serialize_fms_txt
from mod_cobra.gibbs.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, reverse_reaction, get_r_comps
from mimoza_pipeline import process_sbml
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fba_analyser import analyse_by_fba
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva, \
    create_fva_model
from mod_cobra.constraint_based_analysis.efm.efm_analyser import get_efms, group_efms
from mod_sbml.utils.path_manager import create_dirs
from mod_sbml.serialization.serialization_manager import get_sbml_r_formula
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, \
    select_metabolite_ids_by_term_ids
from model_comparison.model_merger import merge_models

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


def analyse_model(sbml, out_r_id, out_rev, res_dir, in_r_id2rev=None, threshold=ZERO_THRESHOLD, do_fva=True,
                  do_fba=True, do_efm=True, efms=None, max_efm_number=1000,
                  cofactors=None, mask_shift=4, get_f_path=None):
    logging.info("Preparing directories...")
    # create directories to store results
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
    in_sbml = sbml
    doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = doc.getModel()
    r_string = lambda r, rev: '<b>%s</b>%s: %s' % (r.getId(), ' (reversed)' if rev else '',
                                                   get_sbml_r_formula(model, r, show_metabolite_ids=False))
    r = model.getReaction(out_r_id)
    description = describe('input_data.html', model_name=model.getName() if model.getName() else model.getId(),
                           sbml_filepath=get_f_path(in_sbml), c_csv=get_f_path(c_csv),
                           m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv), obj_rn=r_string(r, out_rev),
                           in_rn_len=len(in_r_id2rev) if in_r_id2rev else 0,
                           in_rns=
                           '; '.join(r_string(model.getReaction(r_id), rev) for (r_id, rev) in in_r_id2rev.iteritems())
                           if in_r_id2rev else '')

    id2mask = defaultdict(lambda: 0)
    layer2mask = {}
    main_layer=None

    objective_sense = 'minimize' if out_rev else 'maximize'
    vis_r_ids = set()
    if do_fva:
        logging.info("Performing FVA...")
        fva_dir = os.path.join(res_dir, 'fva')
        create_dirs(fva_dir)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        r_id2bounds, opt_val = analyse_by_fva(cobra_model=cobra_model, bm_r_id=out_r_id,
                                              objective_sense=objective_sense, threshold=threshold)
        ess_rn_num, var_rn_num = 0, 0
        if opt_val:
            fva_file = os.path.join(fva_dir, 'fva.txt')
            caption = '%s reaction %s (%s): %.4g\n==================================\n'\
                      % (objective_sense, out_r_id,
                         get_sbml_r_formula(model, model.getReaction(out_r_id), show_metabolite_ids=True,
                                            show_compartments=False), opt_val)
            ess_rn_num, var_rn_num = serialize_fva(cobra_model, r_id2bounds, fva_file,
                                                   title=caption)
            essential_r_id2rev = {r_id: u < 0 for (r_id, (l, u)) in r_id2bounds.iteritems() if l * u > 0}
            ess_r_ids = set(essential_r_id2rev.keys())
            update_vis_layers(ess_r_ids, 'FVA essential', id2mask, layer2mask, mask_shift, model, vis_r_ids)
            main_layer = 'FVA essential'
            mask_shift += 1

            fva_sbml = os.path.join(fva_dir, 'Model_FVA.xml')
            sbml = create_fva_model(sbml, r_id2bounds, fva_sbml)
        description += describe('fva.html', optimal_value=opt_val, ess_r_num=ess_rn_num, var_r_num=var_rn_num,
                                description_filepath=get_f_path(fva_file) if opt_val else None,
                                sbml_filepath=get_f_path(fva_sbml) if opt_val else None)

    if do_fba:
        logging.info("Performing FBA...")
        fba_dir = os.path.join(res_dir, 'fba')
        create_dirs(fba_dir)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        r_id2val, opt_val = analyse_by_fba(cobra_model, bm_r_id=out_r_id, objective_sense=objective_sense,
                                           threshold=threshold)
        if opt_val:
            fba_r_ids = set(r_id2val.keys())
            update_vis_layers(fba_r_ids, 'FBA', id2mask, layer2mask, mask_shift, model, vis_r_ids)
            main_layer = 'FBA'
            mask_shift += 1

            fba_file = os.path.join(fba_dir, 'fba.txt')
            caption = '%s reaction %s (%s): %.4g\n==================================\n'\
                      % (objective_sense, out_r_id,
                         get_sbml_r_formula(model, model.getReaction(out_r_id), show_metabolite_ids=True,
                                            show_compartments=False), opt_val)
            serialize_fluxes(cobra_model, r_id2val, path=fba_file, title=caption)
        description += describe('fba.html', optimal_value=opt_val, r_num=len(r_id2val),
                                description_filepath=get_f_path(fba_file) if opt_val else None)
    id2imp_rns = {}
    clique_merged_sbml = None
    if do_efm:
        logging.info("Performing EFMA...")
        efm_dir = os.path.join(res_dir, 'efms/')
        create_dirs(efm_dir)

        id2efm, efm_id2efficiency = get_efms(target_r_id=out_r_id, target_r_reversed=out_rev, r_id2rev=in_r_id2rev,
                                             sbml=sbml, directory=efm_dir, max_efm_number=max_efm_number,
                                             threshold=threshold, efms=efms)

        efm_num = len(id2efm)
        if efm_num:
            id2fm, fm_id2efficiency, fm_id2key, key2efm_ids = group_efms(id2efm, model, out_r_id, out_rev)
            fm_num = len(id2fm)

            all_efm_intersection = reduce(lambda p1, p2: p1.intersection(p2), id2efm.itervalues(),
                                          next(id2efm.itervalues()))
            all_efm_intersection_len = len(all_efm_intersection)

            if all_efm_intersection_len:
                update_vis_layers(set(all_efm_intersection.to_r_id2coeff(binary=True).iterkeys()),
                                  'EFM intersection', id2mask, layer2mask, mask_shift, model, vis_r_ids)
                if not main_layer:
                    main_layer = 'EFM intersection'
                mask_shift += 1

            logging.info('Detected %d EFMs with common part of length %d.' % (efm_num, all_efm_intersection_len))

            all_efm_file = os.path.join(efm_dir, 'efms.txt')
            serialize_fms_txt(id2fm, id2efm, fm_id2key, key2efm_ids, all_efm_file, efm_id2efficiency, fm_id2efficiency,
                              all_efm_intersection, model)

            # 3 most effective FMs
            limit = min(3, fm_num)
            mask_shift, fms \
                = process_n_fms(model, directory=efm_dir, id2fm=id2fm, id2mask=id2mask,
                                layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                name='FM', sorter=lambda (fm_id, _): -fm_id2efficiency[fm_id],
                                serializer=lambda fm_id, efm_txt:
                                serialize_fm(fm_id, id2fm[fm_id], model, efm_txt, fm_id2efficiency[fm_id]),
                                suffix=lambda _: '', get_f_path=get_f_path, limit=limit, update_vis=True)
            fm_block = describe('fm_block.html', element_num=limit, characteristics='most effective',
                                element_name='EFM', fms=fms)
            description += describe('efms.html', efm_num=efm_num, fm_num=fm_num,
                                    description_filepath=get_f_path(all_efm_file), selected_efm_block=fm_block)

            # Communities
            comm_dir = os.path.join(efm_dir, 'communities')
            create_dirs(comm_dir)
            id2cluster = detect_fm_communities(id2fm, threshold=all_efm_intersection_len)
            id2intersection = {cl_id: reduce(lambda p1, p2: p1.intersection(p2),
                                             (id2fm[fm_id] for fm_id in cluster), id2fm[cluster[0]])
                               for (cl_id, cluster) in id2cluster.iteritems()}
            imp_fraction = 60
            id2imp_rns = {cl_id: detect_reaction_community({fm_id: id2fm[fm_id] for fm_id in cluster},
                                                           imp_fraction, id2intersection[cl_id])
                          for (cl_id, cluster) in id2cluster.iteritems()}
            if id2cluster:
                comm_txt = os.path.join(comm_dir, 'communities.txt')
                serialize_communities(id2cluster, id2intersection, id2imp_rns, fm_num,
                                      all_efm_intersection, comm_txt)
                # All clusters
                limit = len(id2cluster)
                mask_shift, fms \
                    = process_n_fms(model, directory=comm_dir, id2fm=id2imp_rns, id2mask=id2mask,
                                    layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                    name='community', sorter=lambda (cl_id, cl): cl_id,
                                    serializer=lambda cl_id, community_txt:
                                    serialize_community(cl_id, id2intersection[cl_id], id2imp_rns[cl_id],
                                                        id2cluster[cl_id], model, community_txt),
                                    suffix=lambda _: '', limit=limit, get_f_path=get_f_path)
                fm_block = describe('fm_block.html', element_num=limit, characteristics='(all)',
                                    element_name='community', fms=fms, all=True)
            description += describe('communities.html', community_num=len(id2cluster),
                                    description_filepath=get_f_path(comm_txt) if id2cluster else None,
                                    selected_community_block=fm_block)

            # Cliques
            clique_dir = os.path.join(efm_dir, 'cliques')
            create_dirs(clique_dir)
            id2clique = detect_cliques(id2fm)

            if id2clique:
                cliques_txt = os.path.join(clique_dir, 'cliques.txt')
                serialize_cliques(id2clique, cliques_txt, all_efm_intersection)
                # 5 longest cliques
                limit = min(5, len(id2clique))
                mask_shift, fms \
                    = process_n_fms(model, directory=clique_dir, id2fm=id2clique, id2mask=id2mask,
                                    layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                    name='clique', sorter=lambda (_, cl): -len(cl),
                                    serializer=lambda cl_id, clique_txt:
                                    serialize_clique(cl_id, id2clique[cl_id], model, clique_txt),
                                    suffix=lambda _: '', limit=limit, get_f_path=get_f_path, update_vis=False)
                fm_block = describe('fm_block.html', element_num=limit, characteristics='longest',
                                    element_name='clique', fms=fms, all=limit == len(id2clique))

                clique_merged_sbml = os.path.join(clique_dir, 'Model_clique_merged.xml')
                r_id2cl_id = clique2lumped_reaction(id2clique, in_sbml, clique_merged_sbml)

                vis_r_ids |= {r_id2cl_id[r_id] for r_id in vis_r_ids if r_id in r_id2cl_id}
                for (id_, mask) in id2mask.items():
                    if id_ in r_id2cl_id:
                        id2mask[r_id2cl_id[id_]] |= mask

            description += describe('cliques.html', clique_num=len(id2clique),
                                    description_filepath=get_f_path(cliques_txt), selected_clique_block=fm_block)

    logging.info('Putting everything together...')

    if not vis_r_ids:
        description += describe('nothing_found.html')

    return sbml, clique_merged_sbml, vis_r_ids, description, mask_shift, id2mask, layer2mask, id2imp_rns, main_layer


def mean(values):
    return sum(values) / len(values)


def multimodel_pipeline(sbml2parameters, res_dir, do_fva=True, do_fba=True, do_efm=True, max_efm_number=1000):
    create_dirs(res_dir, True)
    mask_shift = 4
    tab2html = {}
    layer2mask = {}

    get_f_path = lambda f: os.path.join('..', '..', os.path.relpath(f, res_dir))

    invisible_layers = []
    sbml2cl_sbml = {}

    sbml2id2mask = {}
    sbml2vis_r_ids = {}
    sbml2name = {}
    for sbml, (out_r_id, out_rev, in_r_id2rev, name) in sbml2parameters.iteritems():
        logging.info('Analysing %s...' % name)
        sub_sbml, clique_merged_sbml, r_ids, description, mask_shift, cur_id2mask, cur_layer2mask, id2imp_rns, main_layer = \
            analyse_model(sbml, out_r_id, out_rev, os.path.join(res_dir, name), in_r_id2rev,
                          do_fva=do_fva, do_fba=do_fba, do_efm=do_efm, mask_shift=mask_shift,
                          get_f_path=get_f_path, max_efm_number=max_efm_number)
        if clique_merged_sbml:
            sbml2cl_sbml[sbml] = clique_merged_sbml
            sbml2name[clique_merged_sbml] = name
            sbml2id2mask[clique_merged_sbml] = cur_id2mask
            sbml2vis_r_ids[clique_merged_sbml] = r_ids
        tab2html['Analysis of %s' % name] = description, None
        layer2mask.update({'%s: %s' % (name, layer): mask for (layer, mask) in cur_layer2mask.iteritems()})
        invisible_layers.extend(['%s: %s' % (name, layer) for layer in cur_layer2mask.iterkeys() if layer != main_layer])
        sbml2id2mask[sbml] = cur_id2mask
        sbml2vis_r_ids[sbml] = r_ids
        sbml2name[sbml] = name

    mm_dir = os.path.join(res_dir, 'merged_model')
    create_dirs(mm_dir)
    visualize_merged_model(get_f_path, invisible_layers, layer2mask, mm_dir, sbml2id2mask, sbml2parameters.keys(),
                           sbml2name, sbml2vis_r_ids, tab2html)

    mm_dir = os.path.join(res_dir, 'merged_clique_model')
    create_dirs(mm_dir)
    visualize_merged_model(get_f_path, invisible_layers, layer2mask, mm_dir, sbml2id2mask, sbml2cl_sbml.values(),
                           sbml2name, sbml2vis_r_ids, tab2html)


def visualize_merged_model(get_f_path, invisible_layers, layer2mask, res_dir, sbml2id2mask, sbml_files, sbml2name,
                           sbml2vis_r_ids, tab2html):
    if not sbml_files:
        return
    info = ''
    id2color = None
    if len(sbml_files) > 1:
        merged_sbml = os.path.join(res_dir, 'Merged_model.xml')
        sbml2id2id, common_ids = merge_models(sbml_files, merged_sbml)
        if common_ids:
            logging.info('Merged all the models into %s.' % merged_sbml)
            comparison_prefix = os.path.join(res_dir, 'Model_comparison_')
            (cc_num, cm_num, cr_num), (comp_csv, m_csv, r_csv) = \
                serialize_common_part_to_csv(merged_sbml, sbml2id2id, common_ids,
                                             {sbml: sbml2name[sbml] for sbml in sbml_files},
                                             comparison_prefix)
        else:
            cc_num, cm_num, cr_num = 0, 0, 0
        tab2html['Model comparison'] = \
            describe('model_comparison.html', c_num=cc_num, m_num=cm_num, r_num=cr_num,
                     c_csv=get_f_path(comp_csv) if cc_num else None,
                     m_csv=get_f_path(m_csv) if cm_num else None,
                     r_csv=get_f_path(r_csv) if cr_num else None), None
        id2mask = defaultdict(lambda: 0)
        vis_r_ids = set()
        for sbml in sbml_files:
            id2mask.update({sbml2id2id[sbml][id_]: id2mask[sbml2id2id[sbml][id_]] | mask
                            for (id_, mask) in sbml2id2mask[sbml].iteritems() if id_ in sbml2id2id[sbml]})
            vis_r_ids |= {sbml2id2id[sbml][r_id] for r_id in sbml2vis_r_ids[sbml] if r_id in sbml2id2id[sbml]}
        colors = get_n_colors(len(sbml_files) + 1, 0.5, 0.8)
        if common_ids:
            mixed_color = colors[0]
            r, g, b = mixed_color
            info = describe('color.html', r=r, g=g, b=b, name='common reactions/metabolites')
        sbml2color = dict(zip(sbml2id2id.iterkeys(), colors[1:]))
        id2color = {}
        for sbml, id2id in sbml2id2id.iteritems():
            color = sbml2color[sbml]
            r, g, b = color
            id2color.update({t_id: color if t_id not in common_ids else mixed_color for t_id in id2id.itervalues()})
            info += describe('color.html', r=r, g=g, b=b, name=sbml2name[sbml])
            title='Combined model analysis'
    else:
        merged_sbml = next(iter(sbml_files))
        vis_r_ids = sbml2vis_r_ids[merged_sbml]
        id2mask = sbml2id2mask[merged_sbml]
        title='Model analysis'
    combined_sbml = os.path.join(res_dir, 'Combined_model.xml')
    r_ids2sbml(vis_r_ids, merged_sbml, combined_sbml, 'combined')
    process_sbml(combined_sbml, verbose=True, path='visualization', generalize=False,
                 id2mask=id2mask, layer2mask=layer2mask, tab2html=tab2html, title=title,
                 id2color=id2color, tabs=None, info=info, invisible_layers=invisible_layers)


def update_vis_layers(r_ids, layer, id2mask, layer2mask, mask_shift, model, vis_r_ids):
    r_ids |= {format_r_id(r_id, False) for r_id in r_ids}
    vis_r_ids |= r_ids
    l_mask = 1 << mask_shift
    layer2mask[layer] = l_mask
    for r_id in r_ids:
        r = model.getReaction(r_id)
        if r:
            id2mask[r_id] |= l_mask
            id2mask.update({c_id: id2mask[c_id] | l_mask for c_id in get_r_comps(r_id, model)})
            id2mask.update({s_id: id2mask[s_id] | l_mask for s_id in get_reactants(r)})
            id2mask.update({s_id: id2mask[s_id] | l_mask for s_id in get_products(r)})


def process_n_fms(model, directory, id2fm, id2mask, layer2mask, mask_shift, vis_r_ids,
                  name, sorter, serializer, suffix, get_f_path, limit, update_vis=True):
    sbml_dir = os.path.join(directory, 'sbml')
    create_dirs(sbml_dir)

    fms = []
    for fm_id, fm in sorted(id2fm.iteritems(), key=sorter):
        r_id2coeff = fm.to_r_id2coeff()
        c_name = name[0].upper() + name[1:]
        fm_name = '%s_%d_of_len_%d' % (c_name, fm_id, len(fm))
        fm_txt = os.path.join(sbml_dir, '%s.txt' % fm_name)
        serializer(fm_id, fm_txt)
        fms.append((c_name, fm_id, len(fm), suffix(fm_id), get_f_path(fm_txt)))

        if update_vis:
            update_vis_layers(set(r_id2coeff), '%s %d' % (c_name, fm_id), id2mask, layer2mask, mask_shift, model,
                              vis_r_ids)
            mask_shift += 1

        limit -= 1
        if limit <= 0:
            break
    return mask_shift, fms

