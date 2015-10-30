from collections import defaultdict, Counter
from itertools import chain

import libsbml

from mod_sbml.serialization.csv_manager import metabolites2df, reactions2df, compartments2df
from mod_sbml.sbml.compartment.compartment_positioner import GO_CYTOSOL, GO_CYTOPLASM
from mod_sbml.annotation.chebi.chebi_annotator import EQUIVALENT_RELATIONSHIPS
from mod_sbml.onto import parse_simple
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi

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
    # if there are models with GO_CYTOPLASM and models with GO_CYTOSOL but none with both,
    # merge GO_CYTOSOL and GO_CYTOPLASM
    if GO_CYTOPLASM in key2comp and GO_CYTOSOL in key2comp \
            and not set(key2comp[GO_CYTOSOL].iterkeys()) & set(key2comp[GO_CYTOPLASM].iterkeys()):
        key2comp[GO_CYTOPLASM].update(key2comp[GO_CYTOSOL])
        del key2comp[GO_CYTOSOL]
    return [comps for comps in key2comp.itervalues() if len(comps) > 1]


def map_chebi_ids(model_id2dfs, chebi):
    initial_chebi_ids = set()
    chebi_id2model_ids = defaultdict(set)
    chebi_id2ancestors = {}
    ancestor2chebi_ids = defaultdict(set)
    for model_id, [df, _, _] in model_id2dfs.iteritems():
        for chebi_id in df["ChEBI"].unique():
            chebi_id2model_ids[chebi_id].add(model_id)
            if chebi_id not in initial_chebi_ids:
                initial_chebi_ids.add(chebi_id)
                ancestors = {chebi_id}
                term = chebi.get_term(chebi_id)
                if term:
                    ancestors |= {t.get_id()
                                  for t in chebi.get_generalized_ancestors(term, direct=False,
                                                                           relationships=EQUIVALENT_RELATIONSHIPS,
                                                                           depth=4)}
                    ancestors |= {t.get_id() for t in chebi.get_equivalents(term=term,
                                                                            relationships=EQUIVALENT_RELATIONSHIPS)}
                chebi_id2ancestors[chebi_id] = ancestors
                for ancestor in ancestors:
                    ancestor2chebi_ids[ancestor].add(chebi_id)
    ancestor2model_id2count = defaultdict(Counter)
    for ancestor, chebi_ids in ancestor2chebi_ids.iteritems():
        for chebi_id in chebi_ids:
            ancestor2model_id2count[ancestor].update({model_id: 1 for model_id in chebi_id2model_ids[chebi_id]})
    # keep only terms that cover at most one initial chebi term per model
    ancestor2count = {ancestor: len(model_id2count)
                      for (ancestor, model_id2count) in ancestor2model_id2count.iteritems()
                      if {1} == set(model_id2count.itervalues())}

    available_ancestors = set(ancestor2count.iterkeys())
    chebi_id2chebi_id = {}
    for chebi_id in initial_chebi_ids:
        ancestors = chebi_id2ancestors[chebi_id] & available_ancestors
        if ancestors:
            the_ancestor = max(ancestors, key=lambda it: ancestor2count[it])
            ancestors -= {the_ancestor}
            available_ancestors -= ancestors
            chebi_id2chebi_id[chebi_id] = the_ancestor

    return chebi_id2chebi_id


def map_metabolites(model_id2dfs, chebi):
    chebi_id2chebi_id = map_chebi_ids(model_id2dfs, chebi)
    key2ms = defaultdict(lambda: defaultdict(set))
    for model_id, [df, _, _] in model_id2dfs.iteritems():
        for index, row in df.iterrows():
            m_id, name, formula, kegg, chebi_id = row['Id'], row['Name'], row['Formula'], row['KEGG'], row["ChEBI"]
            if chebi_id and chebi_id in chebi_id2chebi_id:
                chebi_id = chebi_id2chebi_id[chebi_id]
            key = chebi_id if chebi_id else (kegg if kegg else (formula if formula else name))
            key = key.strip().lower() if key else None
            if key:
                key2ms[key][model_id].add(m_id)
    return [model_id2m_ids for model_id2m_ids in key2ms.itervalues() if len(model_id2m_ids) > 1]


def get_model_data(model_id2sbml):
    model_id2dfs = {}

    for model_id, sbml in model_id2sbml.iteritems():
        if isinstance(sbml, libsbml.Model):
            model = sbml
        else:
            doc = libsbml.SBMLReader().readSBML(sbml)
            model = doc.getModel()
        model_id2dfs[model_id] = [metabolites2df(model), reactions2df(model), compartments2df(model)]

    return model_id2dfs


def map_metabolites_compartments(model_id2dfs, chebi=None):
    if not chebi:
        chebi = parse_simple(get_chebi())

    model_id2m_ids_groups = map_metabolites(model_id2dfs, chebi)
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
    return model_id2c_ids_groups, [it for it in model_id2m_ids_same_comp_groups if len(it) > 1], model_id2c_id2i
