import logging

import libsbml

from mod_cobra.efm.System import System
from mod_cobra import ZERO_THRESHOLD
from mod_cobra.efm.tree_efm_manager import compute_efms
from mod_cobra.efm.numpy_efm_manager import get_element2id_mapping, model2stoichiometric_matrix, \
    get_efm_matrix, get_yield, get_control_efficiency, get_len

__author__ = 'anna'


def analyse_model_efm(sbml, out_r_id, out_rev, res_dir, in_m_id, out_m_id, in_r_id2rev, tree_efm_path,
                      threshold=ZERO_THRESHOLD, max_efm_number=1000):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    system_initial = get_initial_system(sbml, model, res_dir, out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id,
                                        tree_efm_path, max_efm_number, threshold)
    system_folded = system_initial.lump_coupled_reactions()
    system_folded_no_duplicates = system_folded.remove_reaction_duplicates().remove_efm_duplicates()
    system_merged = system_folded_no_duplicates.merge_efm_groups()

    return system_initial, system_folded, system_folded_no_duplicates, system_merged


def get_initial_system(sbml, model, res_dir, out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id,
                       tree_efm_path, max_efm_number, threshold):
    m_id2i, r_id2i = get_element2id_mapping(model)
    N = model2stoichiometric_matrix(model, m_id2i, r_id2i)
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
    efm_id2i = dict(zip(('efm_%d' % it for it in xrange(0, len(i2efficiency))),
                        sorted(i2efficiency.iterkeys(), key=lambda i: i2efficiency[i])))
    return System(N=N, V=V, m_id2i=m_id2i, r_id2i=r_id2i, efm_id2i=efm_id2i,
                  boundary_m_ids=[m.getId() for m in model.getListOfSpecies() if m.getBoundaryCondition()])
