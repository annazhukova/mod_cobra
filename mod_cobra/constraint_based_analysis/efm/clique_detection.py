from collections import Counter
import logging

import libsbml
from networkx import find_cliques, Graph

from mod_cobra.constraint_based_analysis.efm.EFM import EFM
from mod_sbml.sbml.submodel_manager import remove_unused_species, create_lumped_reaction

__author__ = 'anna'


def detect_cliques(id2fm, min_clique_size=2):
    """
    The method takes the found FMs and constructs the reaction graph in the following way:
    nodes are marked with reaction ids (or -r_id for the reversed versions of reversible reactions),
    there exists an edge between nodes r_i and r_j iff they are related, i.e.
    there exists at least efm_num of EFMs that contain both reaction r_i and reaction r_j.
    The method then detects the maximal cliques of size greater or equal to min_clique_size.

    :param id2fm: dictionary that maps FM identifiers (int) to the FMs.
    :param min_clique_size: int, minimal size of a clique for it to be considered.

    :return: id2clique
    """
    logging.info("Going to rank reactions by FM number.")
    r_id_pair2count = Counter()
    r_id2count = Counter()

    for fm_id, fm in id2fm.iteritems():
        r_id2coeff = fm.to_r_id2coeff(binary=False)
        r_ids = list(r_id2coeff.iterkeys())
        i = 0
        for r_id in r_ids:
            i += 1
            r_id2count.update({r_id: 1})
            for r_id2 in r_ids[i:]:
                r_0, r_1 = sorted([r_id, r_id2])
                r_id_pair2count.update(
                    {(r_0, r_1, 1 if r_id2coeff[r_0] > 0 else -1, r_id2coeff[r_1] / r_id2coeff[r_0]): 1})
    gr = Graph()
    r_id2direction = {}
    r_id_pair2ratio = {}
    for (r_0, r_1, direction, ratio), count in r_id_pair2count.iteritems():
        if count == r_id2count[r_0] and count == r_id2count[r_1]:
            gr.add_edge(r_0, r_1)
            r_id2direction[r_0] = direction
            r_id_pair2ratio[(r_0, r_1)] = ratio
    sample_fm = next(id2fm.itervalues())
    r_ids, rev_r_ids, int_size = sample_fm.r_ids, sample_fm.rev_r_ids, sample_fm.int_size

    def clique2r_id2coeff(clique):
        clique = sorted(clique)
        r_0 = clique[0]
        clique = clique[1:]
        r_id2coeff = {r_0: r_id2direction[r_0]}
        r_id2coeff.update({r_1: r_id_pair2ratio[(r_0, r_1)] * r_id2direction[r_0] for r_1 in clique})
        return r_id2coeff

    cliques = [EFM(r_ids=r_ids, rev_r_ids=rev_r_ids, int_size=int_size,
                   r_id2coeff=clique2r_id2coeff(clique)) for clique in
               (clique for clique in find_cliques(gr) if len(clique) >= min_clique_size)]
    id2clique = dict(zip(xrange(0, len(cliques)), cliques))
    return id2clique


def clique2lumped_reaction(id2clique, in_sbml, out_sbml):
    doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = doc.getModel()
    r_id2cl_id = {}
    for cl_id, clique in id2clique.iteritems():
        r_id2coeff = clique.to_r_id2coeff()
        new_r_id = create_lumped_reaction(r_id2coeff, model, id_prefix='rl_%d' % cl_id)
        r_id2cl_id.update({r_id: new_r_id for r_id in r_id2coeff.iterkeys()})
    remove_unused_species(model)
    model.setId('%s_merged_cliques' % model.getId())
    model.setName('%s_merged_cliques' % model.getName())
    libsbml.SBMLWriter().writeSBMLToFile(doc, out_sbml)
    return r_id2cl_id






