from collections import defaultdict

import libsbml

from mod_sbml.serialization.csv_manager import metabolites2df, reactions2df, compartments2df
from mod_sbml.sbml.compartment.compartment_positioner import get_go_term, GO_CYTOSOL, GO_CYTOPLASM
from mod_sbml.annotation. chebi.chebi_annotator import get_species_to_chebi
from mod_sbml.onto import parse_simple
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from mod_sbml.annotation.gene_ontology.go_serializer import get_go

__author__ = 'anna'


def map_comps(model_id2dfs):
    key2comp = defaultdict(lambda: defaultdict(set))
    for model_id, [_, _, df] in model_id2dfs.iteritems():
        for index, row in df.iterrows():
            c_id, c_name, t_id = row['Id'], row['Name'], row['GO']
            key = t_id if t_id else c_name
            key = key.strip().lower() if key else ''
            if key:
                key2comp[key][model_id].add(c_id)
    if GO_CYTOPLASM in key2comp and GO_CYTOSOL in key2comp \
            and not set(it[0] for it in key2comp[GO_CYTOSOL]) & set(it[0] for it in key2comp[GO_CYTOPLASM]):
        key2comp[GO_CYTOPLASM] |= key2comp[GO_CYTOSOL]
        del key2comp[GO_CYTOSOL]
    return [comps for comps in key2comp.itervalues() if len(comps) > 1]


def map_metabolites(model_id2dfs):
    key2ms = defaultdict(lambda: defaultdict(set))
    for model_id, [df, _, _] in model_id2dfs.iteritems():
        for index, row in df.iterrows():
            m_id, name, formula, kegg, chebi = row['Id'], row['Name'], row['Formula'], row['KEGG'], row["ChEBI"]
            key = chebi if chebi else (kegg if kegg else (formula if formula else name))
            key = key.strip().lower() if key else None
            if key:
                key2ms[key][model_id].add(m_id)
    return [model_id2m_ids for model_id2m_ids in key2ms.itervalues() if len(model_id2m_ids) > 1]


def get_model_data(model_id2sbml, chebi=None, go=None):
    if not chebi:
        chebi = parse_simple(get_chebi())
    if not go:
        go = parse_simple(get_go())

    def get_t_id(comp):
        term = get_go_term(comp, go)
        return term.get_id() if term else None

    model_id2dfs = {}

    for model_id, sbml in model_id2sbml.iteritems():
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        m_id2chebi_id = get_species_to_chebi(model, chebi)
        model_id2dfs[model_id] = \
            [metabolites2df(model, m_id2chebi_id), reactions2df(model),
             compartments2df(model, get_t_id)]

    return model_id2dfs


def map_metabolites_compartments(model_id2dfs):

    model_id2m_ids_groups = map_metabolites(model_id2dfs)
    model_id2c_ids_groups = map_comps(model_id2dfs)

    model_id2c_id2i = defaultdict(dict)
    for i, model_id2c_ids in enumerate(model_id2c_ids_groups):
        for model_id, c_ids in model_id2c_ids.iteritems():
            model_id2c_id2i[model_id].update({c_id: i for c_id in c_ids})

    model_id2m_ids_same_comp_groups = []
    for model_id2m_ids in model_id2m_ids_groups:
        i2m_ids = defaultdict(list)
        for model_id, m_ids in model_id2m_ids.iteritems():
            for m_id in m_ids:
                df, _, _ = model_id2dfs[model_id]
                c_id = df.at[m_id, 'Compartment']
                if c_id in model_id2c_id2i[model_id]:
                    i2m_ids[model_id2c_id2i[model_id][c_id]].append((model_id, m_id))
        for m_ids in i2m_ids.itervalues():
            model_id2m_ids = defaultdict(set)
            if len(m_ids) > 1:
                for model_id, m_id in m_ids:
                    model_id2m_ids[model_id].add(m_id)
            if model_id2m_ids:
                model_id2m_ids_same_comp_groups.append(model_id2m_ids)
    return model_id2c_ids_groups, model_id2m_ids_same_comp_groups, model_id2c_id2i