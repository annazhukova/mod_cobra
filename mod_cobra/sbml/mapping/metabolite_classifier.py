from collections import defaultdict
from mod_sbml.annotation.chebi.chebi_annotator import EQUIVALENT_RELATIONSHIPS

from mod_sbml.annotation.chebi.chebi_serializer import CHEBI
from mod_sbml.onto import parse

__author__ = 'anna'


CHEBI_AMINO_ACID = 'chebi:33709'
CHEBI_LIPID = 'chebi:18059'
CHEBI_SALT = 'chebi:24866'
CHEBI_VITAMIN = 'chebi:33229'
CHEBI_HEME = 'chebi:30413'


def classify_m_ids(m_ids, m_id2chebi_id,
                   chebi_classes=(CHEBI_AMINO_ACID, CHEBI_LIPID, CHEBI_SALT, CHEBI_HEME, CHEBI_VITAMIN), chebi=None):
    if not chebi:
        chebi = parse(CHEBI)
    terms = set()
    chebi_id2m_ids = defaultdict(set)
    for m_id in m_ids:
        if m_id in m_id2chebi_id:
            term = chebi.get_term(m_id2chebi_id[m_id])
            if term:
                terms.add(term)
                chebi_id2m_ids[m_id2chebi_id[m_id]].add(m_id)
    t2chebi_ids = {t: {it.get_id() for it in chebi.get_sub_tree(t, relationships=EQUIVALENT_RELATIONSHIPS) & terms}
                   for t in (chebi.get_term(t_id) for t_id in chebi_classes)}
    return {t.get_name(): reduce(lambda s1, s2: s1 | s2, (chebi_id2m_ids[ch_id] for ch_id in ch_ids), set())
            for (t, ch_ids) in t2chebi_ids.iteritems() if ch_ids}
