from collections import defaultdict
import logging
import os
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from mod_sbml.onto import parse_simple
from mod_cobra.sbml.mapping.metabolite_matcher import get_model_data, map_metabolites_compartments
from mod_sbml.sbml.ubiquitous_manager import get_proton_ch_ids
from mod_cobra.sbml.mapping.model_merger import simple_merge_models, join, merge

__author__ = 'anna'


def combine_models(model_id2sbml, model_id2S, path):
    logging.info('Going to merge models...')
    chebi = parse_simple(get_chebi())
    model_id2dfs = get_model_data(model_id2sbml)
    model_id2c_id_groups, model_id2m_id_groups, model_id2c_id2i = \
        map_metabolites_compartments(model_id2dfs, chebi=chebi)
    logging.info('Mapped metabolites and compartments.')
    ignore_m_ids = get_ignored_metabolites(model_id2dfs, get_proton_ch_ids())
    S = join(model_id2m_id_groups, model_id2S)
    ignore_m_ids |= {S.m_id2gr_id[m_id] for m_id in ignore_m_ids if m_id in S.m_id2gr_id}
    merge(S, ignore_m_ids)
    model_id2r_id_groups = get_r_id_groups(S)
    logging.info('Mapped reactions.')
    sbml = os.path.join(path, 'Merged_model.xml')
    model_id2id2id, common_ids, S = simple_merge_models(S, model_id2c_id2i, model_id2dfs, sbml)
    return sbml, S, model_id2id2id, common_ids, model_id2dfs, \
           (model_id2c_id_groups, model_id2m_id_groups, model_id2r_id_groups)


def get_r_id_groups(S_merged):
    model_id2r_id_groups = []
    for r_id2c in S_merged.gr_id2r_id2c.itervalues():
        model_id2r_ids = defaultdict(set)
        for (model_id, r_id) in r_id2c.iterkeys():
            model_id2r_ids[model_id].add(r_id)
        model_id2r_id_groups.append(model_id2r_ids)
    return model_id2r_id_groups


def get_ignored_metabolites(model_id2dfs, ignored_ch_ids):
    ignored_m_ids = set()
    for model_id, [df, _, _] in model_id2dfs.iteritems():
        for index, row in df.iterrows():
            m_id = row['Id']
            t_id = row["ChEBI"] if 'ChEBI' in row else None
            if t_id in ignored_ch_ids:
                ignored_m_ids.add((model_id, m_id))
    return ignored_m_ids
