from collections import defaultdict
import logging

import libsbml
from emsampler import sampler

from mod_cobra.efm.System import System, FoldedSystem
from mod_cobra import ZERO_THRESHOLD
from mod_cobra.efm.tree_efm_manager import compute_efms
from mod_cobra.efm.numpy_efm_manager import get_element2id_mapping, model2stoichiometric_matrix, \
    get_efm_matrix, get_yield, get_control_efficiency, get_len, remove_invalid_efms

__author__ = 'anna'


def analyse_model_efm(sbml, out_r_id, out_rev, res_dir, in_m_id, out_m_id, in_r_id2rev, tree_efm_path,
                      threshold=ZERO_THRESHOLD, max_efm_number=1000):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    system_initial = get_initial_system(sbml, model, res_dir, out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id,
                                        tree_efm_path, max_efm_number, threshold)
    if not system_initial.efm_id2i:
        return system_initial, FoldedSystem(st_matrix=system_initial.st_matrix, pws=system_initial.pws), \
               FoldedSystem(st_matrix=system_initial.st_matrix, pws=system_initial.pws), \
               FoldedSystem(st_matrix=system_initial.st_matrix, pws=system_initial.pws), defaultdict(lambda: 1)

    logging.info('Computed initial system')
    system_folded = system_initial.lump_coupled_reactions().remove_unused_metabolites()
    r_id2w = defaultdict(lambda: 1)
    r_id2w.update({gr_r_id: len(r_id2c) for (gr_r_id, r_id2c) in system_folded.gr_id2r_id2c.iteritems()})
    logging.info('Folded pathways')
    system_folded_no_duplicates = \
        system_folded.remove_reaction_duplicates().remove_efm_duplicates()
    r_id2w.update({gr_r_id: max(r_id2w[r_id] for r_id in r_id2c.iterkeys())
                   for (gr_r_id, r_id2c) in system_folded_no_duplicates.gr_id2r_id2c.iteritems()})
    logging.info('Removed duplicates')
    system_merged = system_folded_no_duplicates.merge_efm_groups()
    logging.info('Merged pathways')
    return system_initial, system_folded, system_folded_no_duplicates, system_merged, r_id2w


def get_initial_system(sbml, model, res_dir, out_r_id, out_rev, in_r_id2rev, in_m_id, out_m_id,
                       tree_efm_path, max_efm_number, threshold):
    m_id2i, r_id2i = get_element2id_mapping(model)
    boundary_m_ids = [m.getId() for m in model.getListOfSpecies() if m.getBoundaryCondition()]
    N = model2stoichiometric_matrix(model, m_id2i, r_id2i)
    logging.info('Calculated N of shape %s' % 'x'.join((str(it) for it in N.shape)))
    r_id2coefficient_list = compute_efms(sbml, res_dir, max_efm_number, out_r_id, out_rev, r_id2rev=in_r_id2rev,
                                         tree_efm_path=tree_efm_path, threshold=threshold)
    V = get_efm_matrix(r_id2coefficient_list, r_id2i)
    # inner_m_ids = [m.getId() for m in model.getListOfSpecies() if not m.getBoundaryCondition()]
    # rev = [model.getReaction(r_id).getReversible() for (r_id, i) in
    #        sorted(r_id2i.iteritems(), key=lambda (_, i): i)]
    # N_in = N[[m_id2i[m_id] for m_id in inner_m_ids], :]
    # V = sampler(N_in, rev, K=None, debug=True)
    # V = V[:, [i for i in xrange(0, V.shape[1])
    #           if V[r_id2i[out_r_id], i] * (-1 if out_rev else 1) > 0 and
    #           next((False for (in_r_id, in_rev) in in_r_id2rev.iteritems() if
    #                 V[r_id2i[in_r_id], i] * (-1 if in_rev else 1) <= 0), True)]]
    # V = remove_invalid_efms(N[tuple(m_id2i[m_id] for m_id in inner_m_ids), :], V)
    logging.info('Calculated V of shape %s' % 'x'.join((str(it) for it in V.shape)))
    v_yield = lambda v_i: -get_yield(N, V[:, v_i], m_id2i[out_m_id], m_id2i[in_m_id]) \
        if out_m_id and in_m_id and out_m_id in m_id2i and in_m_id in m_id2i else None
    i2efficiency = {i: (v_yield(i), -get_control_efficiency(V[:, i], r_id2i[out_r_id]), get_len(V[:, i]))
                    for i in xrange(0, V.shape[1])}
    efm_id2i = dict(zip(('efm_%d' % it for it in xrange(0, len(i2efficiency))),
                        sorted(i2efficiency.iterkeys(), key=lambda i: i2efficiency[i])))
    return System(N=N, V=V, m_id2i=m_id2i, r_id2i=r_id2i, efm_id2i=efm_id2i,
                  boundary_m_ids=boundary_m_ids)
