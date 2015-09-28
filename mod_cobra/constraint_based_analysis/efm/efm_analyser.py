from collections import defaultdict
import logging

from mod_cobra.constraint_based_analysis.efm.control_effective_flux_calculator import get_fm_efficiency, get_fm_yield
from mod_cobra.constraint_based_analysis import ZERO_THRESHOLD
from mod_cobra.constraint_based_analysis.efm.efm_manager import compute_efms
from mod_sbml.sbml.submodel_manager import compress_reaction_participants

TREEEFM_PATH = "/home/anna/Applications/TreeEFM/tool/TreeEFMseq"

__author__ = 'anna'


def get_efms(target_r_id, target_r_reversed, r_id2rev, sbml, directory, max_efm_number=1000, threshold=ZERO_THRESHOLD,
             efms=None, rewrite=True, tree_efm_path=TREEEFM_PATH):
    if not efms:
        if tree_efm_path:
            efms, r_ids, rev_r_ids = compute_efms(sbml, directory, max_efm_number, target_r_id, target_r_reversed,
                                                  tree_efm_path, r_id2rev, threshold=threshold, rewrite=rewrite)
            if not efms:
                logging.info('Found no EFMs of interest.')
                return None
        else:
            raise ValueError('No EFMs found :(. Probably, you forgot to specify TreeEFM path?')

    efm2efficiency = {efm: get_fm_efficiency(efm, target_r_id, target_r_reversed) for efm in efms}
    id2efm = {}
    id2efficiency = {}

    i = 0
    for efm in sorted(efms, key=lambda efm: -efm2efficiency[efm]):
        efm_id = 'EFM %d' % i
        efm.id = efm_id
        id2efm[efm_id] = efm
        id2efficiency[efm_id] = efm2efficiency[efm]
        i += 1

    return id2efm, id2efficiency


def group_efms(id2efm, model, in_m_id, out_m_id):
    key2efm_ids = defaultdict(list)
    for efm_id, efm in id2efm.iteritems():
        r_id2st, p_id2st = compress_reaction_participants(model, efm.to_r_id2coeff())
        r_id2st, p_id2st = tuple(sorted(r_id2st.iteritems())), tuple(sorted(p_id2st.iteritems()))
        key2efm_ids[r_id2st, p_id2st].append(efm_id)
    fm2key = {id2efm[efm_ids[0]].join([id2efm[efm_id] for efm_id in efm_ids[1:]]): key
              for (key, efm_ids) in key2efm_ids.iteritems()}
    id2fm = {}
    id2key = {}
    i = 0
    for fm in sorted(fm2key.iterkeys(), reverse=True,
                     key=lambda fm: (get_fm_yield(fm2key[fm], in_m_id, out_m_id), max(key2efm_ids[fm2key[fm]]))):
        if not fm.id:
            fm.id = 'Pathway %d' % i
            i += 1
        id2fm[fm.id] = fm
        id2key[fm.id] = fm2key[fm]
    return id2fm, id2key, key2efm_ids


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


