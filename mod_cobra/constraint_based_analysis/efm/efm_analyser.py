from collections import defaultdict
import logging

from mod_cobra.constraint_based_analysis.efm.control_effective_flux_calculator import get_fm_efficiency
from mod_cobra.constraint_based_analysis import ZERO_THRESHOLD
from mod_cobra.constraint_based_analysis.efm.efm_manager import compute_efms
from mod_sbml.sbml.submodel_manager import compress_reaction_participants

TREEEFM_PATH = "/home/anna/Applications/TreeEFM/tool/TreeEFMseq"

__author__ = 'anna'


def get_efms(target_r_id, target_r_reversed, r_id2rev, sbml, directory, max_efm_number=1000, threshold=ZERO_THRESHOLD,
             efms=None, r_ids=None, rewrite=True, tree_efm_path=TREEEFM_PATH):
    if not efms:
        if tree_efm_path:
            efms, r_ids, rev_r_ids = compute_efms(sbml, directory, max_efm_number, target_r_id, target_r_reversed,
                                                  tree_efm_path, r_id2rev, threshold=threshold, r_ids=r_ids,
                                                  rewrite=rewrite)
            if not efms:
                logging.info('Found no EFMs of interest.')
                return None
        else:
            raise ValueError('No EFMs found :(. Probably, you forgot to specify TreeEFM path?')

    efm2efficiency = {efm: get_fm_efficiency(efm, target_r_id, target_r_reversed) for efm in efms}
    id2efm = dict(zip((str(it) for it in xrange(0, len(efms))), sorted(efms, key=lambda efm: -efm2efficiency[efm])))

    return id2efm, {efm_id: efm2efficiency[efm] for (efm_id, efm) in id2efm.iteritems()}


def group_efms(id2efm, model, target_r_id, target_r_rev):
    key2efm_ids = defaultdict(list)
    for efm_id, efm in id2efm.iteritems():
        m_id2stoichiometry = tuple(sorted(compress_reaction_participants(model, efm.to_r_id2coeff()).iteritems()))
        key2efm_ids[m_id2stoichiometry].append(efm_id)
    fm2key = {id2efm[efm_ids[0]].join([id2efm[efm_id] for efm_id in efm_ids[1:]]): key
              for (key, efm_ids) in key2efm_ids.iteritems()}

    fm2efficiency = {fm: get_fm_efficiency(fm, target_r_id, target_r_rev) for fm in fm2key.iterkeys()}
    id2fm = dict(zip((str(it) for it in xrange(0, len(fm2key))),
                     sorted(fm2key.iterkeys(), key=lambda fm: (-fm2efficiency[fm], min(key2efm_ids[fm2key[fm]])))))

    return id2fm, {fm_id: fm2efficiency[fm] for (fm_id, fm) in id2fm.iteritems()}, \
           {fm_id: fm2key[fm] for (fm_id, fm) in id2fm.iteritems()}, key2efm_ids


def calculate_imp_rn_threshold(total_len, imp_rn_threshold=None):
    if not imp_rn_threshold:
        imp_rn_threshold = max(2, int(total_len * 0.3))
    return imp_rn_threshold


def calculate_min_pattern_len(avg_efm_len, essential_rn_number, min_efm_len, min_pattern_len):
    if not min_pattern_len:
        min_pattern_len = min(int(avg_efm_len * 0.3), min_efm_len)
        if essential_rn_number:
            min_pattern_len = min(max(essential_rn_number + 3, min_pattern_len), avg_efm_len)
    return min_pattern_len


def calculate_min_clique_len(intersection_len, min_clique_len=None):
    if not min_clique_len:
        min_clique_len = intersection_len + 1
    return min_clique_len


