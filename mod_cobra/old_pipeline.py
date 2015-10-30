from collections import defaultdict
import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
import libsbml

from sbml_vis.graph.color.color import get_n_colors
from mod_cobra.fbva import MAXIMIZE, MINIMIZE
from mod_cobra.fbva.serialization import fva_serializer, fba_serializer
from mod_cobra.efm.community_detector import detect_communities_by_boundary_metabolites, detect_reaction_community
from mod_cobra.efm.efm_pipeline import analyse_model_efm
from mod_cobra.html import describe
from mod_cobra.efm.serialization.coupled_reaction_sbml_manager import create_folded_model
from mod_cobra.sbml.constraint_manager import constraint_exchange_reactions
from mod_cobra.sbml.serialization import model_serializer, mapping_serializer
from mod_cobra.sbml.mapping.mapping_pipeline import combine_models
from mod_cobra import ZERO_THRESHOLD
from mod_cobra.fbva.serialization import format_r_id
from mod_cobra.sbml import r_ids2sbml
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, get_r_comps, get_model_name
from sbml_vis.mimoza_pipeline import process_sbml
from mod_cobra.fbva.fbva_analyser import analyse_by_fba, analyse_by_fva
from mod_cobra.fbva.serialization.fbva_sbml_manager import create_fva_model
from mod_sbml.utils.path_manager import create_dirs
from mod_cobra.efm.serialization import coupled_reaction_group_serializer, community_serializer, efm_serializer

TREEEFM_PATH = "/home/anna/Applications/TreeEFM/tool/TreeEFMseq"

__author__ = 'anna'


def _prepare_dir(super_dir, dir_name, log, rewrite=True):
    logging.info(log)
    cur_dir = os.path.join(super_dir, dir_name)
    create_dirs(cur_dir, rewrite)
    return cur_dir


def analyse_model(sbml, out_r_id, out_rev, res_dir, in_m_id, out_m_id, in_r_id2rev=None, threshold=ZERO_THRESHOLD,
                  do_fva=True, do_fba=True, do_efm=True, max_efm_number=1000, mask_shift=4,
                  get_f_path=None, tree_efm_path=TREEEFM_PATH, main_dir=None, rewrite=True):
    model_name = get_model_name(sbml)
    logging.info('Analysing %s...' % model_name)

    # create directories to store results
    logging.info("Preparing directories...")
    res_dir = os.path.join(res_dir, ''.join(ch for ch in model_name.replace(' ', '_') if ch.isalnum() or '_' == ch))
    create_dirs(res_dir, False)
    if not get_f_path:
        get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir))

    if in_r_id2rev:
        sbml = constraint_exchange_reactions(sbml, forsed_r_id2rev=in_r_id2rev)

    # copy our model in the result directory
    if os.path.normpath(res_dir) != os.path.normpath(os.path.dirname(sbml)):
        shutil.copy(sbml, res_dir)
        sbml = os.path.join(res_dir, os.path.basename(sbml))

    description = model_serializer.serialize(sbml, out_r_id, out_rev, in_r_id2rev, res_dir, get_f_path)

    r_id2mask, layer2mask, vis_r_ids, main_layer = defaultdict(lambda: 0), {}, set(), None

    cobra_model, opt_val, objective_sense = None, None, MINIMIZE if out_rev else MAXIMIZE

    if do_fva:
        cur_dir = _prepare_dir(res_dir, 'fva', "Performing FVA...")
        cobra_model = create_cobra_model_from_sbml_file(sbml) if not cobra_model else cobra_model
        r_id2bounds, opt_val = analyse_by_fva(cobra_model=cobra_model, bm_r_id=out_r_id,
                                              objective_sense=objective_sense, threshold=threshold)
        if opt_val:
            mask_shift = update_vis_layers((r_id for (r_id, (l, u)) in r_id2bounds.iteritems() if l * u > 0),
                                           'FVA essential', r_id2mask, layer2mask, mask_shift, vis_r_ids)
            main_layer = 'FVA essential'
            fva_sbml = os.path.join(cur_dir, 'Model_FVA.xml')
            sbml = create_fva_model(sbml, r_id2bounds, fva_sbml)
        description += fva_serializer.serialize(cobra_model, opt_val, r_id2bounds, objective_sense, out_r_id,
                                                cur_dir, get_f_path, sbml)
    if do_fba:
        cur_dir = _prepare_dir(res_dir, 'fba', "Performing FBA...")
        cobra_model = create_cobra_model_from_sbml_file(sbml) if not cobra_model else cobra_model
        r_id2val, opt_val = analyse_by_fba(cobra_model, bm_r_id=out_r_id, objective_sense=objective_sense,
                                           threshold=threshold)
        if opt_val:
            mask_shift = update_vis_layers(r_id2val.iterkeys(), 'FBA', r_id2mask, layer2mask, mask_shift, vis_r_ids)
            main_layer = 'FBA'
        description += fba_serializer.serialize(cobra_model, opt_val, r_id2val, objective_sense, out_r_id, cur_dir,
                                                get_f_path)

    S, r_id2w = None, {}
    if do_efm:
        cur_dir = _prepare_dir(res_dir, 'efma', "Performing EFMA...", rewrite=rewrite)

        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()

        S = analyse_model_efm(model, out_r_id, out_rev, cur_dir, in_m_id, out_m_id, in_r_id2rev, tree_efm_path,
                              threshold, max_efm_number, rewrite=rewrite)

        for serializer in (efm_serializer.serialize, coupled_reaction_group_serializer.serialize):
            description += \
                serializer(model=model, path=cur_dir, get_f_path=get_f_path, in_m_id=in_m_id, out_m_id=out_m_id,
                           out_r_id=out_r_id, S=S, model_name=model_name, main_dir=main_dir)

        if S.gr_id2r_id2c:
            clique_merged_sbml = os.path.join(cur_dir, 'Model_folded.xml')
            r_id2new_r_id = create_folded_model(S, sbml, clique_merged_sbml)
            sbml = clique_merged_sbml

            vis_r_ids |= {cl_id for (r_id, cl_id) in r_id2new_r_id.iteritems() if r_id in vis_r_ids}
            for r_id, new_r_id in r_id2new_r_id.iteritems():
                if r_id in r_id2mask:
                    r_id2mask[new_r_id] |= r_id2mask[r_id]

    if not opt_val and (not S or not S.efm_id2i):
        description += describe('nothing_found.html')

    return model_name, S, sbml, vis_r_ids, description, mask_shift, r_id2mask, layer2mask, main_layer


def multimodel_pipeline(sbml2parameters, res_dir, do_fva=True, do_fba=True, do_efm=True, max_efm_number=1000,
                        rewrite=True):
    vis_dir = os.path.join(res_dir, 'visualization')
    create_dirs(vis_dir, rewrite)
    get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir)) if f else None

    mask_shift, tab2html, layer2mask, invisible_layers = 4, {}, {}, []
    model_id2id2mask, model_id2vis_r_ids = {}, {}
    model_id2sbml, model_id2S = {}, {}

    for model_id, (out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id) in sbml2parameters.iteritems():
        name, model_id2S[name], model_id2sbml[name], model_id2vis_r_ids[name], description, mask_shift, \
        model_id2id2mask[name], cur_layer2mask, main_layer = \
            analyse_model(sbml=model_id, out_r_id=out_r_id, out_rev=out_rev, in_r_id2rev=in_r_id2rev, in_m_id=in_m_id,
                          out_m_id=out_m_id, res_dir=res_dir, do_fva=do_fva, do_fba=do_fba,
                          do_efm=do_efm, mask_shift=mask_shift, get_f_path=get_f_path, max_efm_number=max_efm_number,
                          main_dir=vis_dir, rewrite=rewrite)

        tab2html['Analysis of %s' % name] = description, None
        layer2mask.update({'%s: %s' % (name, layer): mask for (layer, mask) in cur_layer2mask.iteritems()})
        invisible_layers.extend(['%s: %s' % (name, layer)
                                 for layer in cur_layer2mask.iterkeys() if layer != main_layer])

    if len(model_id2sbml) > 1:
        mm_dir = os.path.join(res_dir, 'merged_model')
        create_dirs(mm_dir)

        sbml, S, model_id2id2id, common_ids, model_id2dfs, mappings = combine_models(model_id2sbml, model_id2S, mm_dir)

        tab2html['Model comparison'] = mapping_serializer.serialize(model_id2dfs, mappings, mm_dir, get_f_path), None
        title = 'Combined model analysis'

        id2mask = defaultdict(lambda: 0)
        vis_r_ids = set()
        for model_id in model_id2sbml.iterkeys():
            id2mask.update({model_id2id2id[model_id][s_id]: id2mask[model_id2id2id[model_id][s_id]] | mask
                            for (s_id, mask) in model_id2id2mask[model_id].iteritems()
                            if s_id in model_id2id2id[model_id]})
            vis_r_ids |= {model_id2id2id[model_id][r_id] for r_id in model_id2vis_r_ids[model_id]
                          if r_id in model_id2id2id[model_id]}

        id2color, info = get_colors(common_ids, model_id2id2id, model_id2sbml.keys())
    else:
        model_id, sbml = next(model_id2sbml.iteritems())
        S = model_id2S[model_id].get_main_S()
        vis_r_ids, id2mask, id2color = model_id2vis_r_ids[model_id], model_id2id2mask[model_id], None
        info, title = '', 'Model analysis'

    # Communities
    comm_dir = _prepare_dir(res_dir, 'communities', "Analysing communities...")
    id2cluster = detect_communities_by_boundary_metabolites(S)
    id2intersection = {cl_id: S.get_efm_intersection(cluster) for (cl_id, cluster) in id2cluster.iteritems()}
    # id2imp_rns = {cl_id: {} for (cl_id, cluster) in id2cluster.iteritems()}
    id2bounds = {cl_id: S.get_boundary_metabolite_distribution(cluster) for (cl_id, cluster) in id2cluster.iteritems()}
    id2imp_rns = {cl_id: detect_reaction_community(S, cluster, id2intersection[cl_id])
                  for (cl_id, cluster) in id2cluster.iteritems()}
    if id2cluster:
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        description = \
            community_serializer.serialize(model, S, id2cluster, id2intersection, id2bounds, id2imp_rns,
                                           comm_dir, get_f_path)
        if len(model_id2sbml) > 1:
            tab2html['Model comparison'] = tab2html['Model comparison'][0] + description, None
        else:
            tab2html['Pathway communities'] = description, None
        for cl_id, r_id2c in id2intersection.iteritems():
            mask_shift = update_vis_layers(r_id2c.iterkeys(), 'Pathway community %s' % cl_id, id2mask, layer2mask,
                                           mask_shift, vis_r_ids)
            invisible_layers.append('Pathway community %s' % cl_id)

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


def update_vis_layers(r_ids, layer, id2mask, layer2mask, mask_shift, vis_r_ids):
    if not r_ids:
        return mask_shift
    l_mask = 1 << mask_shift
    layer2mask[layer] = l_mask
    for r_id in r_ids:
        o_r_id = format_r_id(r_id, remove_prefix=False)
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


def get_colors(common_ids, model_id2id2id, model_ids):
    colors = get_n_colors(len(model_ids) + 1, 0.5, 0.8)
    mixed_color = None
    description = ''
    if common_ids:
        mixed_color = colors[0]
        r, g, b = mixed_color
        description = describe('color.html', r=r, g=g, b=b, name='common reactions/metabolites')
    model_id2color = dict(zip(model_ids, colors[1:]))
    id2color = {}
    for model_id, id2id in model_id2id2id.iteritems():
        color = model_id2color[model_id]
        r, g, b = color
        id2color.update({t_id: color if t_id not in common_ids else mixed_color for t_id in id2id.itervalues()})
        description += describe('color.html', r=r, g=g, b=b, name=model_id)
    return id2color, description


