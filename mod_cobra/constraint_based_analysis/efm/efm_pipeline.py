import logging
import os
import libsbml
import sys
import shutil
from SBMLTestCase import create_test_sbml, TEST_SBML, DATA_DIR
from System import System, FoldedSystem
from chebi.chebi_annotator import get_species_to_chebi
from chebi.chebi_serializer import get_chebi
from compartment.compartment_manager import separate_boundary_species
from csv_manager import serialize_model_info

from mod_cobra.constraint_based_analysis import ZERO_THRESHOLD
from mod_cobra.constraint_based_analysis.efm.efm_analyser import TREEEFM_PATH
from mod_cobra.constraint_based_analysis.efm.efm_manager import compute_efms
from mod_cobra.constraint_based_analysis.efm.numpy_efm_manager import get_element2id_mapping, get_stoichiometric_matrix, \
    get_efm_matrix, get_yield, get_control_efficiency, get_len, get_coupled_reactions, lump_coupled_reactions, \
    get_reaction_duplicates, remove_reaction_duplicates, get_efm_duplicates, remove_efm_duplicates, \
    get_efm_groups_based_on_boundary_metabolites, merge_efm_groups, remove_invalid_efms
from mod_cobra.constraint_based_analysis.efm.numpy_efm_serialization_manager import serialize_efms, serialize_r_cliques
from models.model_list import MODEL_DIR
from models.rr.model_rr import RR_SER_PRODUCTION, RR_MODEL_DIR, RR_MODEL, RR_GLN_B, RR_SER_B, RR_GLN_EXCHANGE
from mod_sbml.utils.path_manager import create_dirs
from mod_cobra.gibbs.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import reverse_reaction, get_reactants, get_products
from models.serine.serine_sr import SR_SER_PRODUCTION, SR_GLN_B, SR_SER_B, SR_GLN_EXCHANGE
from models.smith_robinson.model_SmithRobinson import SR_MODEL_SER, SR_MODEL_DIR
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, select_metabolite_ids_by_term_ids
from onto import parse_simple
from recon.model_recon import RECON_MODEL_BOUNDARIES
from serine.serine_recon import RECON_SER_PRODUCTION, RECON_GLN_B, RECON_SER_B, RECON_GLN_EXCHANGE

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


def analyse_model(sbml, out_r_id, out_rev, res_dir, in_m_id, out_m_id, in_r_id2rev=None, threshold=ZERO_THRESHOLD,
                  tree_efm_path=TREEEFM_PATH, max_efm_number=1000, serializers=(serialize_efms, serialize_r_cliques)):
    create_dirs(res_dir, True)

    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    # s_id2chebi = get_species_to_chebi(model, parse_simple(get_chebi()))
    # separate_boundary_species(model, s_id2chebi)
    if in_r_id2rev:
        constraint_exchange_reactions(model, allowed_exchange_r_id2rev=in_r_id2rev)
    sbml = os.path.join(res_dir, '%s_constrained.xml' % os.path.splitext(os.path.basename(sbml))[0])
    libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)
    serialize_model_info(sbml, os.path.join(res_dir, 'model_info_'))

    analyse_model_efm(sbml, out_r_id, out_rev, res_dir, in_m_id, out_m_id, in_r_id2rev, threshold,
                      tree_efm_path, max_efm_number, serializers)


def analyse_model_efm(sbml, out_r_id, out_rev, res_dir, in_m_id, out_m_id, in_r_id2rev=None, threshold=ZERO_THRESHOLD,
                      tree_efm_path=TREEEFM_PATH, max_efm_number=1000, serializers=None):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    # Step 1: compute EFMs
    S_initial = step1(sbml, model, res_dir, out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id, tree_efm_path, max_efm_number,
              threshold)

    # # Step 2: lump coupled reactions
    S_coupled, S_no_duplicates = step2(S_initial)

    # # Step 3: merge pathways based on boundary metabolites
    S_merged = step3(model, S_no_duplicates)

    if serializers:
        for serialize in serializers:
            serialize(model=model, path=res_dir, in_m_id=in_m_id, out_m_id=out_m_id, out_r_id=out_r_id,
                      S_initial=S_initial, S_coupled=S_coupled, S_no_duplicates=S_no_duplicates, S_merged=S_merged)


def step3(model, S):
    boundary_m_indices = \
        [S.m_id2i[m.getId()] for m in model.getListOfSpecies() if m.getBoundaryCondition()]
    boundary_efm_id_groups = \
        get_efm_groups_based_on_boundary_metabolites(S.N[boundary_m_indices, :],
                                                     S.V, S.efm_id2i)
    # boundary_efm_id_groups = get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i)
    V_m, m_efm_id2i, efm_id2merged_id = merge_efm_groups(S.V, boundary_efm_id_groups,
                                                         S.efm_id2i)
    S_merged = FoldedSystem(st_matrix=S.st_matrix, V=V_m, efm_id2i=m_efm_id2i,
                            efm_id2gr_id=efm_id2merged_id, boundary_m_ids=S.boundary_m_ids)
    return S_merged


def step2(S):
    coupled_r_id_groups = get_coupled_reactions(S.V, S.r_id2i)
    N_c, V_c, c_r_id2i, r_id2lr_id, lr_id2r_id2c = lump_coupled_reactions(S.N, S.V, coupled_r_id_groups, S.r_id2i)
    S_coupled = FoldedSystem(N=N_c, V=V_c, m_id2i=S.m_id2i, r_id2i=c_r_id2i, efm_id2i=S.efm_id2i,
                             r_id2gr_id=r_id2lr_id, gr_id2r_id2c=lr_id2r_id2c, boundary_m_ids=S.boundary_m_ids)

    duplicated_r_id_groups = get_reaction_duplicates(S_coupled.N, S_coupled.r_id2i)
    N_d, V_d, d_r_id2i, r_id2gr_id, gr_id2r_id2c = \
        remove_reaction_duplicates(S_coupled.N, S_coupled.V, duplicated_r_id_groups, S_coupled.r_id2i)
    duplicated_efm_id_groups = get_efm_duplicates(V_d, S_coupled.efm_id2i)
    V_dd, d_efm_id2i, efm_id2gr_id = remove_efm_duplicates(V_d, duplicated_efm_id_groups, S_coupled.efm_id2i)
    S_no_duplicates = FoldedSystem(N=N_d, V=V_dd, m_id2i=S_coupled.m_id2i, r_id2i=d_r_id2i, efm_id2i=d_efm_id2i,
                                   r_id2gr_id=r_id2gr_id, gr_id2r_id2c=gr_id2r_id2c, efm_id2gr_id=efm_id2gr_id,
                                   boundary_m_ids=S_coupled.boundary_m_ids)

    return S_coupled, S_no_duplicates


def step1(sbml, model, res_dir, out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id,
          tree_efm_path, max_efm_number, threshold):
    m_id2i, r_id2i = get_element2id_mapping(model)
    N = get_stoichiometric_matrix(model, m_id2i, r_id2i)
    logging.info('Calculated N of shape %s' % 'x'.join((str(it) for it in N.shape)))
    r_id2coefficient_list = compute_efms(sbml, res_dir, max_efm_number, out_r_id, out_rev, r_id2rev=in_r_id2rev,
                                         tree_efm_path=tree_efm_path, threshold=threshold)
    V = get_efm_matrix(r_id2coefficient_list, r_id2i)
    # inner_m_ids = [m.getId() for m in model.getListOfSpecies() if not m.getBoundaryCondition()]
    # V = remove_invalid_efms(N[tuple(m_id2i[m_id] for m_id in inner_m_ids), :], V)
    logging.info('Calculated V of shape %s' % 'x'.join((str(it) for it in V.shape)))
    v_yield = lambda v_i: -get_yield(N, V[:, v_i], m_id2i[out_m_id], m_id2i[in_m_id]) \
        if out_m_id and in_m_id and out_m_id in m_id2i and in_m_id in m_id2i else None
    i2efficiency = {i: (v_yield(i), -get_control_efficiency(V[:, i], r_id2i[out_r_id]), get_len(V[:, i]))
                    for i in xrange(0, V.shape[1])}
    efm_id2i = dict(zip(xrange(0, len(i2efficiency)), sorted(i2efficiency.iterkeys(), key=lambda i: i2efficiency[i])))
    return System(N=N, V=V, m_id2i=m_id2i, r_id2i=r_id2i, efm_id2i=efm_id2i,
                  boundary_m_ids=[m.getId() for m in model.getListOfSpecies() if m.getBoundaryCondition()])


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    res_dir = os.path.join(MODEL_DIR, 'test_efm')

    # analyse_model(sbml=RECON_MODEL_BOUNDARIES, out_r_id=RECON_SER_PRODUCTION, out_rev=False, res_dir=res_dir,
    #               in_m_id=RECON_GLN_B, out_m_id=RECON_SER_B, in_r_id2rev={RECON_GLN_EXCHANGE: True},
    #               max_efm_number=1000)

    analyse_model(sbml=SR_MODEL_SER, out_r_id=SR_SER_PRODUCTION, out_rev=False, res_dir=res_dir,
                  in_m_id=SR_GLN_B, out_m_id=SR_SER_B, in_r_id2rev={SR_GLN_EXCHANGE: False}, max_efm_number=1000)

    # analyse_model(sbml=RR_MODEL, out_r_id=RR_SER_PRODUCTION, out_rev=False, res_dir=res_dir,
    #               in_m_id=RR_GLN_B, out_m_id=RR_SER_B, in_r_id2rev={RR_GLN_EXCHANGE: False}, max_efm_number=1000)


    # create_test_sbml()
    # analyse_model(sbml=TEST_SBML, out_r_id='r3', out_rev=False, res_dir=res_dir,
    #               in_m_id='m1_b', out_m_id='m2_b', in_r_id2rev={'r1': False})



if "__main__" == __name__:
    sys.exit(main())
