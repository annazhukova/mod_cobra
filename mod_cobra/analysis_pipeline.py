import io
import logging
import os
from shutil import copytree

from jinja2 import Environment, PackageLoader
import libsbml
import shutil

from mod_sbml.annotation.chebi.chebi_annotator import get_species_id2chebi_id
from mod_sbml.annotation.chebi.chebi_serializer import CHEBI
from mod_sbml.annotation.kegg.pathway_manager import get_name2pw
from mod_sbml.annotation.pts.pts_serializer import get_pts
from mod_cobra.efm.community_detector import detect_communities_by_inputs_of_type, \
    detect_communities_by_boundary_metabolites
from mod_cobra.efm.serialization import community_serializer
from mod_sbml.onto import parse_simple, parse
from mod_cobra.sbml.mapping.pathway_mapper import get_pathways
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, select_metabolite_ids_by_term_ids
from mod_sbml.annotation.annotator import annotate
from mod_cobra.efm.efm_pipeline import analyse_model_efm
from mod_cobra.efm.serialization import efm_serializer, coupled_reaction_group_serializer
from mod_cobra.efm.serialization.coupled_reaction_sbml_manager import create_folded_model
from mod_cobra.html import describe
from mod_cobra.sbml.constraint_manager import constraint_exchange_reactions
from mod_sbml.sbml.sbml_manager import get_model_name
from mod_cobra.sbml.serialization import model_serializer
from mod_sbml.utils.path_manager import create_dirs
from mod_cobra.sbml.mapping.mapping_pipeline import combine_models
from mod_cobra.sbml.serialization import mapping_serializer
import sbml_vis

__author__ = 'anna'


def multimodel_pipeline(sbml2parameters, res_dir, treeefm_path, max_efm_number=1000, rewrite=True, org=None):
    create_dirs(res_dir, rewrite)
    get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir)) if f else None
    tab2html = {}

    model_id2sbml, model_id2S, model_id2efm_id2pws = {}, {}, {}

    name2pw = get_name2pw()
    pts = parse_simple(get_pts())
    root_ids = {t.get_id() for t in pts.get_roots()}

    chebi = parse(CHEBI)
    ub_ch_ids = get_ubiquitous_chebi_ids(add_common=True, add_cofactors=True, chebi=chebi)

    efm_id2pws = {}

    model_id2cofactors = {}
    modeld_id2m_id2chebi_id = {}

    for sbml, (r_id2rev, r_id2rev_banned) in sbml2parameters.iteritems():
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()

        model_name = get_model_name(model=model)
        short_model_name = model_name
        if len(model_name) > 12:
            short_model_name = model_name[:10].strip('-_ ')
            if len(short_model_name) == 10:
                short_model_name += '...'
        safe_m_name = ''.join(ch for ch in short_model_name.replace(' ', '_') if ch.isalnum() or '_' == ch)
        logging.info('Analysing %s...' % model_name)

        # create directories to store results
        logging.info("Preparing directories...")
        m_dir = os.path.join(res_dir, safe_m_name)
        create_dirs(m_dir, rewrite)

        # exchange_rs = get_exchange_reactions(model)
        # csv = '%s/%s.exchanges.csv' % (m_dir, safe_m_name)
        # df2csv(reactions2df(model, exchange_rs), csv)

        cofactors = select_metabolite_ids_by_term_ids(model, ub_ch_ids)

        if r_id2rev:
            constraint_exchange_reactions(model, forsed_r_id2rev=r_id2rev, prohibited_r_id2rev=r_id2rev_banned,
                                          cofactors=cofactors if not r_id2rev_banned else None)

        logging.info("Annotating the model...")
        annotate(model, org=org, reactions=False, pathways=False, chebi=chebi)
        m_id2ch_id = get_species_id2chebi_id(model)

        # copy our model in the result directory
        sbml = os.path.join(m_dir, '%s.constrained.xml' % safe_m_name)
        libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)

        description = model_serializer.serialize(sbml, model, model_name, r_id2rev, m_dir, get_f_path)

        pw2rs = get_pathways(model, pts, name2pw, root_ids)

        logging.info("Performing EFMA...")
        efma_dir = os.path.join(m_dir, 'efma')
        create_dirs(efma_dir, rewrite)

        S, efm_id2pws = analyse_model_efm(model, efma_dir, r_id2rev, tree_efm_path=treeefm_path,
                                          max_efm_number=max_efm_number, rewrite=rewrite, pw2rs=pw2rs)

        for serializer in (efm_serializer.serialize, coupled_reaction_group_serializer.serialize):
            description += \
                serializer(model=model, path=efma_dir, get_f_path=get_f_path, S=S, model_name=model_name)

        if S.gr_id2r_id2c:
            sbml = os.path.join(efma_dir, '%s.folded.xml' % safe_m_name)
            create_folded_model(S, model)
            libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)

        if not S or not S.efm_id2i:
            description += describe('nothing_found.html')

        model_id2sbml[safe_m_name] = sbml
        model_id2S[safe_m_name] = S
        model_id2efm_id2pws[safe_m_name] = efm_id2pws
        model_id2cofactors[safe_m_name] = cofactors
        modeld_id2m_id2chebi_id[safe_m_name] = m_id2ch_id

        tab2html['Analysis of %s' % short_model_name.decode('utf-8')] = description, None

    cofactors = set()
    m_id2ch_id = {}
    if len(model_id2sbml) > 1:
        mm_dir = os.path.join(res_dir, 'merged_model')
        create_dirs(mm_dir)

        sbml, S, model_id2id2id, common_ids, model_id2dfs, mappings = combine_models(model_id2sbml, model_id2S, mm_dir)

        for model_id in model_id2sbml.iterkeys():
            efm_id2pws.update({model_id2id2id[model_id][efm_id]: pws
                               for (efm_id, pws) in model_id2efm_id2pws[model_id].iteritems()
                               if efm_id in model_id2id2id[model_id]})
            cofactors |= {model_id2id2id[model_id][m_id] for m_id in model_id2cofactors[model_id]
                          if m_id in model_id2id2id[model_id]}
            m_id2ch_id.update({model_id2id2id[model_id][m_id]: ch_id
                               for (m_id, ch_id) in modeld_id2m_id2chebi_id[model_id].iteritems()
                               if m_id in model_id2id2id[model_id]})

        tab2html['Model comparison'] = mapping_serializer.serialize(model_id2dfs, mappings, mm_dir, get_f_path), None
        title = 'Combined model analysis'
    else:
        model_id, sbml = next(model_id2sbml.iteritems())
        efm_id2pws = model_id2efm_id2pws[model_id]
        cofactors = model_id2cofactors[model_id]
        m_id2ch_id = modeld_id2m_id2chebi_id[model_id]
        S = model_id2S[model_id].get_main_S()
        info, title, id2color = '', 'Model analysis', None


    # Communities
    logging.info("Analysing communities...")
    comm_dir = os.path.join(res_dir, 'communities')
    create_dirs(comm_dir, rewrite)

    # id2cluster = detect_communities_by_inputs_of_type(S, 'AMINO ACID', m_id2ch_id, chebi)
    id2cluster = detect_communities_by_boundary_metabolites(S, cofactors=cofactors, threshold=50)

    if id2cluster:
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        description = \
            community_serializer.serialize(model, S, id2cluster, comm_dir, get_f_path, m_id2ch_id, chebi)
        if len(model_id2sbml) > 1:
            tab2html['Model comparison'] = tab2html['Model comparison'][0] + description, None
        else:
            tab2html['EFM communities'] = description, None

    serialize(res_dir, tab2html, title)


def serialize(res_dir, tab2html, title):
    vis_dir = os.path.join(res_dir, 'visualisation')
    lib_dir = os.path.join(vis_dir, 'lib')
    shutil.rmtree(lib_dir, True)
    copytree(os.path.join(os.path.dirname(os.path.abspath(sbml_vis.__file__)), '..', 'lib'),
             lib_dir)
    env = Environment(loader=PackageLoader('sbml_vis.html', 'templates'))
    template = env.get_template('tabbed_page.html')
    page = template.render(css_list=[], js_list=[], title=title, tab2html=tab2html)
    with io.open(os.path.join(res_dir, 'visualisation', 'index.html'), 'w+', encoding='utf-8') as f:
        f.write(page)