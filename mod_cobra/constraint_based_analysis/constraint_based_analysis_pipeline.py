from collections import defaultdict
import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
import libsbml

from community_detection import detect_efm_communities, detect_reaction_community
from html_descriptor import describe
from mod_sbml.serialization.csv_manager import serialize_common_part_to_csv, serialize_model_info
from graph.color.color import get_n_colors
from mod_cobra.constraint_based_analysis import ZERO_THRESHOLD
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import serialize_fva, \
    serialize_fluxes
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from mod_cobra.constraint_based_analysis.efm.efm_serialization_manager import r_ids2sbml, serialize_efms_txt, serialize_efm, \
    serialize_communities, serialize_community
from mod_cobra.gibbs.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, reverse_reaction, get_r_comps
from mimoza_pipeline import process_sbml
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fba_analyser import analyse_by_fba, \
    create_fba_model
from mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva, \
    create_fva_model, create_essential_r_ids_model
from mod_cobra.constraint_based_analysis.efm.efm_analyser import get_efms
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
                  do_fba=True, do_efm=True, efms=None, max_efm_number=1000, imp_rn_threshold=None,
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

    objective_sense = 'minimize' if out_rev else 'maximize'
    vis_r_ids = set()
    if do_fva:
        logging.info("Performing FVA...")
        fva_dir = os.path.join(res_dir, 'fva')
        create_dirs(fva_dir)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        r_id2bounds, opt_val = analyse_by_fva(cobra_model=cobra_model, bm_r_id=out_r_id,
                                              objective_sense=objective_sense, threshold=threshold)
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
            update_vis_layers(ess_r_ids, 'FVA: essential', id2mask, layer2mask, mask_shift, model, vis_r_ids)
            mask_shift += 1

            if opt_val and sbml:
                fva_sbml = os.path.join(fva_dir, 'Model_FVA.xml')
                sbml = create_fva_model(sbml, r_id2bounds, fva_sbml)
                ess_sbml = create_essential_r_ids_model(sbml, essential_r_id2rev, fva_dir)
                description += describe('fva.html', optimal_value=opt_val, ess_r_num=ess_rn_num, var_r_num=var_rn_num,
                                        description_filepath=get_f_path(fva_file), sbml_filepath=get_f_path(fva_sbml),
                                        ess_sbml_filepath=get_f_path(ess_sbml))

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
            mask_shift += 1

            fba_file = os.path.join(fba_dir, 'fba.txt')
            caption = '%s reaction %s (%s): %.4g\n==================================\n'\
                      % (objective_sense, out_r_id,
                         get_sbml_r_formula(model, model.getReaction(out_r_id), show_metabolite_ids=True,
                                            show_compartments=False), opt_val)
            serialize_fluxes(cobra_model, r_id2val, path=fba_file, title=caption)

            if opt_val and sbml:
                efm_sbml = os.path.join(fba_dir, 'Model_FBA.xml')
                create_fba_model(sbml, r_id2val, efm_sbml)
                description += describe('fba.html', optimal_value=opt_val, r_num=len(r_id2val),
                                        description_filepath=get_f_path(fba_file), sbml_filepath=get_f_path(efm_sbml))

    if do_efm:
        logging.info("Performing EFMA...")
        efm_dir = os.path.join(res_dir, 'efms/')
        create_dirs(efm_dir)

        id2efm = get_efms(target_r_id=out_r_id, target_r_reversed=out_rev, r_id2rev=in_r_id2rev, sbml=sbml,
                          directory=efm_dir, max_efm_number=max_efm_number, threshold=threshold, efms=efms,
                          rewrite=out_r_id)
        if id2efm:
            all_efm_intersection = reduce(lambda p1, p2: p1.intersection(p2), id2efm.itervalues(),
                                          next(id2efm.itervalues()))
            update_vis_layers(set(all_efm_intersection.to_r_id2coeff(binary=True).iterkeys()),
                              'EFM intersection', id2mask, layer2mask, mask_shift, model, vis_r_ids)
            mask_shift += 1

            logging.info('Detected %d EFMs with common part of length %d.' % (len(id2efm), len(all_efm_intersection)))
            all_efm_file = os.path.join(efm_dir, 'efms.txt')
            serialize_efms_txt(id2efm, all_efm_file, out_r_id, all_efm_intersection)
            # 3 shortest EFMs
            limit = min(3, len(id2efm))
            mask_shift, fms \
                = process_n_fms(sbml, model, directory=efm_dir, id2fm=id2efm, id2mask=id2mask,
                                layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                name='EFM', sorter=lambda (_, efm): len(efm),
                                serializer=lambda efm_id, efm_txt:
                                serialize_efm(efm_id, id2efm[efm_id], model, efm_txt, out_r_id),
                                suffix=lambda _: '', get_f_path=get_f_path, limit=limit, update_vis=False)
            fm_block = describe('fm_block.html', element_num=limit, characteristics='shortest', element_name='EFM',
                                fms=fms)
            description += describe('efms.html', efm_num=len(id2efm), description_filepath=get_f_path(all_efm_file),
                                    selected_efm_block=fm_block)

            # # Important reactions
            # imp_rn_threshold = calculate_imp_rn_threshold(len(id2efm))
            # rn_dir = os.path.join(efm_dir, 'important')
            # create_dirs(rn_dir)
            # r_id2efm_ids, important_r_ids = get_important_reactions(id2efm, imp_rn_threshold)
            # if important_r_ids:
            #     imp_rn_txt = os.path.join(rn_dir, 'r_list.txt')
            #     serialize_important_reactions(r_id2efm_ids, model, imp_rn_txt, imp_rn_threshold)
            #     imp_rn_sbml = os.path.join(rn_dir, 'Model_important.xml')
            #     r_ids2sbml(important_r_ids, sbml, imp_rn_sbml, suffix='important')
            #     update_vis_layers(important_r_ids, 'Important reactions', id2mask, layer2mask, mask_shift, model,
            #                       vis_r_ids)
            #     mask_shift += 1
            # description += describe('imp_reactions.html', r_num=len(important_r_ids), threshold=imp_rn_threshold,
            #                         description_filepath=get_f_path(imp_rn_txt) if important_r_ids else None,
            #                         sbml_filepath=get_f_path(imp_rn_sbml) if important_r_ids else None)

            # Communities
            comm_dir = os.path.join(efm_dir, 'communities')
            create_dirs(comm_dir)
            id2cluster = detect_efm_communities(id2efm, threshold=len(all_efm_intersection))
            id2intersection = {cl_id: reduce(lambda p1, p2: p1.intersection(p2),
                                             (id2efm[efm_id] for efm_id in cluster), id2efm[cluster[0]])
                               for (cl_id, cluster) in id2cluster.iteritems()}
            imp_fraction = 60
            # id2imp_rns = {cl_id: get_important_reaction_FM({efm_id: id2efm[efm_id] for efm_id in cluster},
            #                                                max(imp_fraction * len(cluster) / 100, 2))
            #               for (cl_id, cluster) in id2cluster.iteritems()}
            id2imp_rns = {cl_id: detect_reaction_community({efm_id: id2efm[efm_id] for efm_id in cluster},
                                                                     imp_fraction, id2intersection[cl_id])
                          for (cl_id, cluster) in id2cluster.iteritems()}
            if id2cluster:
                comm_txt = os.path.join(comm_dir, 'communities.txt')
                serialize_communities(id2cluster, id2intersection, id2imp_rns, len(id2efm),
                                      all_efm_intersection, comm_txt)
                # All clusters
                limit = len(id2cluster)
                mask_shift, fms \
                    = process_n_fms(sbml, model, directory=comm_dir, id2fm=id2imp_rns, id2mask=id2mask,
                                    layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                    name='community', sorter=lambda (cl_id, cl): cl_id,
                                    serializer=lambda cl_id, community_txt:
                                    serialize_community(cl_id, id2intersection[cl_id], id2imp_rns[cl_id],
                                                        len(id2imp_rns[cl_id]), model, community_txt),
                                    suffix=lambda _: '', limit=limit, get_f_path=get_f_path)
                fm_block = describe('fm_block.html', element_num=limit, characteristics='(all)',
                                    element_name='community', fms=fms, all=True)
            description += describe('communities.html', community_num=len(id2cluster),
                                    description_filepath=get_f_path(comm_txt) if id2cluster else None,
                                    selected_community_block=fm_block)

            # Reaction communities
            # comm_dir = os.path.join(efm_dir, 'reaction_communities')
            # create_dirs(comm_dir)
            # r_community_threshold = len(id2efm) / 2
            # id2r_community = detect_reaction_communities(id2efm, threshold=r_community_threshold, min_len=2)
            # if id2r_community:
            #     comm_txt = os.path.join(comm_dir, 'reaction_communities.txt')
            #     serialize_cliques(id2r_community, comm_txt, 2, all_efm_intersection, r_community_threshold)
            #     # 3 longest cliques
            #     limit = min(3, len(id2r_community))
            #     mask_shift, fms \
            #         = process_n_fms(sbml, model, directory=comm_dir, id2fm=id2r_community, id2mask=id2mask,
            #                         layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
            #                         name='clique', sorter=lambda (_, cl): -len(cl),
            #                         serializer=lambda cl_id, clique_txt:
            #                         serialize_clique(cl_id, id2r_community[cl_id], model, clique_txt),
            #                         suffix=lambda _: '', limit=limit, get_f_path=get_f_path)
            #     fm_block = describe('fm_block.html', element_num=limit, characteristics='longest', element_name='clique',
            #                         fms=fms)
            # description += describe('cliques.html', clique_num=len(id2r_community), min_len=2,
            #                         threshold=r_community_threshold, description_filepath=get_f_path(comm_txt),
            #                         selected_clique_block=fm_block)

            # # Cliques
            # min_clique_len = calculate_min_clique_len(len(all_efm_intersection))
            # clique_dir = os.path.join(efm_dir, 'cliques')
            # create_dirs(clique_dir)
            # id2clique = detect_cliques(id2efm, min_clique_size=min_clique_len, efm_num=imp_rn_threshold)
            #
            # if id2clique:
            #     cliques_txt = os.path.join(clique_dir, 'cliques.txt')
            #     serialize_cliques(id2clique, cliques_txt, min_clique_len, all_efm_intersection,
            #                       min_efm_num=imp_rn_threshold)
            #     # 3 longest cliques
            #     limit = min(3, len(id2clique))
            #     mask_shift, fms \
            #         = process_n_fms(sbml, model, directory=clique_dir, id2fm=id2clique, id2mask=id2mask,
            #                         layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
            #                         name='clique', sorter=lambda (_, cl): -len(cl),
            #                         serializer=lambda cl_id, clique_txt:
            #                         serialize_clique(cl_id, id2clique[cl_id], model, clique_txt),
            #                         suffix=lambda _: '', limit=limit, get_f_path=get_f_path)
            #     fm_block = describe('fm_block.html', element_num=limit, characteristics='longest', element_name='clique',
            #                         fms=fms)
            #
            # description += describe('cliques.html', clique_num=len(id2clique), min_len=min_clique_len,
            #                         threshold=imp_rn_threshold, description_filepath=get_f_path(cliques_txt),
            #                         selected_clique_block=fm_block)

    logging.info('Putting everything together...')

    if not vis_r_ids:
        description += describe('nothing_found.html')

    return sbml, vis_r_ids, description, mask_shift, id2mask, layer2mask


def mean(values):
    return sum(values) / len(values)


def multimodel_pipeline(sbml2parameters, res_dir, do_fva=True, do_fba=True, do_efm=True):
    create_dirs(res_dir, True)
    mask_shift = 4
    tab2html = {}
    layer2mask = {}

    get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir))
    merged_sbml = os.path.join(res_dir, 'Merged_model.xml')
    sbml2id2id, common_ids = merge_models(sbml2parameters.iterkeys(), merged_sbml)
    logging.info('Merged all the models into %s.' % merged_sbml)
    comparison_prefix = os.path.join(res_dir, 'Model_comparison_')
    (cc_num, cm_num, cr_num), (comp_csv, m_csv, r_csv) = \
        serialize_common_part_to_csv(merged_sbml, sbml2id2id, common_ids,
                                     {sbml: pars[-1] for (sbml, pars) in sbml2parameters.iteritems()}, comparison_prefix)

    tab2html['Model comparison'] = \
        describe('model_comparison.html', c_num=cc_num, m_num=cm_num, r_num=cr_num,
                 c_csv=get_f_path(comp_csv), m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv)), None

    id2mask = defaultdict(lambda: 0)
    vis_r_ids = set()
    for sbml, (out_r_id, out_rev, in_r_id2rev, name) in sbml2parameters.iteritems():
        logging.info('Analysing %s...' % name)
        sub_sbml, r_ids, description, mask_shift, cur_id2mask, cur_layer2mask = \
            analyse_model(sbml, out_r_id, out_rev, os.path.join(res_dir, name), in_r_id2rev,
                          do_fva=do_fva, do_fba=do_fba, do_efm=do_efm, max_efm_number=1000, mask_shift=mask_shift,
                          get_f_path=get_f_path)
        tab2html['Analysis of %s' % name] = description, None
        layer2mask.update({'%s: %s' % (name, layer): mask for (layer, mask) in cur_layer2mask.iteritems()})
        id2mask.update({sbml2id2id[sbml][id_]: id2mask[sbml2id2id[sbml][id_]] | mask
                        for (id_, mask) in cur_id2mask.iteritems() if id_ in sbml2id2id[sbml]})
        vis_r_ids |= {sbml2id2id[sbml][r_id] for r_id in r_ids if r_id in sbml2id2id[sbml]}

    combined_sbml = os.path.join(res_dir, 'Combined_model.xml')
    r_ids2sbml(vis_r_ids, merged_sbml, combined_sbml, 'combined')

    colors = get_n_colors(len(sbml2parameters) + 1, 0.5, 0.8)
    info = ''
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
        info += describe('color.html', r=r, g=g, b=b, name=sbml2parameters[sbml][-1])

    process_sbml(combined_sbml, verbose=True, path='visualization', generalize=False,
                 id2mask=id2mask, layer2mask=layer2mask, tab2html=tab2html, title='Combined model analysis',
                 id2color=id2color, tabs=None, info=info)


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


def process_n_fms(sbml, model, directory, id2fm, id2mask, layer2mask, mask_shift, vis_r_ids,
                  name, sorter, serializer, suffix, get_f_path, limit, update_vis=True):
    sbml_dir = os.path.join(directory, 'sbml')
    create_dirs(sbml_dir)

    fms = []
    for fm_id, fm in sorted(id2fm.iteritems(), key=sorter):
        r_id2coeff = fm.to_r_id2coeff()
        c_name = name[0].upper() + name[1:]
        fm_name = '%s_%d_of_len_%d' % (c_name, fm_id, len(fm))
        fm_sbml = os.path.join(sbml_dir, '%s.xml' % fm_name)
        r_ids2sbml(r_id2coeff.keys(), sbml, fm_sbml, fm_name)
        fm_txt = os.path.join(sbml_dir, '%s.txt' % fm_name)
        serializer(fm_id, fm_txt)
        fms.append((c_name, fm_id, len(fm), suffix(fm_id), get_f_path(fm_txt), get_f_path(fm_sbml)))

        if update_vis:
            update_vis_layers(set(r_id2coeff), '%s %d' % (c_name, fm_id), id2mask, layer2mask, mask_shift, model,
                              vis_r_ids)
            mask_shift += 1

        limit -= 1
        if limit <= 0:
            break
    return mask_shift, fms

