import logging
from functools import reduce

from mod_cobra.efm.System import System
from mod_cobra import ZERO_THRESHOLD
from mod_cobra.efm.tree_efm_manager import compute_efms
from mod_cobra.efm.numpy_efm_manager import get_element2id_mapping, model2stoichiometric_matrix, \
    get_efm_matrix
from mod_cobra.sbml.mapping.pathway_mapper import r_ids2pws

__author__ = 'anna'


def analyse_model_efm(model, res_dir, r_id2rev, tree_efm_path, threshold=ZERO_THRESHOLD, max_efm_number=1000,
                      rewrite=True, pw2rs=None):
    system, efm_id2pws = get_initial_system(model, res_dir, r_id2rev, tree_efm_path, max_efm_number, threshold,
                                            rewrite=rewrite, pw2rs=pw2rs)
    logging.info('Computed initial system')
    if not system.efm_id2i:
        return system

    system.lump_coupled_reactions()
    system.remove_unused_metabolites()
    logging.info('Folded pathways')

    system.remove_reaction_duplicates()
    system.remove_efm_duplicates()
    logging.info('Removed duplicates')

    for gr_id, efm_ids in system.gr_id2efm_ids.items():
        efm_id2pws[gr_id] = reduce(lambda s1, s2: s1 | s2, (efm_id2pws[efm_id] for efm_id in efm_ids), set())

    # system.merge_efm_groups()
    # logging.info('Merged pathways')

    return system, efm_id2pws


def get_initial_system(model, res_dir, r_id2rev, tree_efm_path, max_efm_number, threshold, rewrite=True,
                       pw2rs=None):
    m_id2i, r_id2i = get_element2id_mapping(model)
    boundary_m_ids = [m.getId() for m in model.getListOfSpecies() if m.getBoundaryCondition()]
    N = model2stoichiometric_matrix(model, m_id2i, r_id2i)
    logging.info('Calculated N of shape %s' % 'x'.join((str(it) for it in N.shape)))
    r_id2coefficient_list = compute_efms(model, res_dir, max_efm_number, r_id2rev=r_id2rev,
                                         tree_efm_path=tree_efm_path, threshold=threshold, rewrite=rewrite)

    efm_id2pws = {'efm_%d' % i: r_ids2pws(set(r_id2coeff.keys()), pw2rs)
                  for i, r_id2coeff in enumerate(r_id2coefficient_list)}
    V = get_efm_matrix(r_id2coefficient_list, r_id2i)
    # inner_m_ids = [m.getId() for m in model.getListOfSpecies() if not m.getBoundaryCondition()]
    # rev = [model.getReaction(r_id).getReversible() for (r_id, i) in
    #        sorted(r_id2i.items(), key=lambda (_, i): i)]
    # N_in = N[[m_id2i[m_id] for m_id in inner_m_ids], :]
    # V = sampler(N_in, rev, K=None, debug=True)
    # V = V[:, [i for i in range(0, V.shape[1])
    #           if V[r_id2i[out_r_id], i] * (-1 if out_rev else 1) > 0 and
    #           next((False for (in_r_id, in_rev) in in_r_id2rev.items() if
    #                 V[r_id2i[in_r_id], i] * (-1 if in_rev else 1) <= 0), True)]]
    # V = remove_invalid_efms(N[tuple(m_id2i[m_id] for m_id in inner_m_ids), :], V)

    logging.info('Calculated V of shape %s' % 'x'.join((str(it) for it in V.shape)))
    efm_id2i = {'efm_%d' % i: i for i in range(0, V.shape[1])}
    S = System(N=N, V=V, m_id2i=m_id2i, r_id2i=r_id2i, efm_id2i=efm_id2i, boundary_m_ids=boundary_m_ids)
    return S.get_used_system(), efm_id2pws
