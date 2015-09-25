from collections import Counter, defaultdict
import logging

import libsbml

from networkx import find_cliques, Graph
from mod_sbml.utils.misc import invert_map

from mod_cobra.constraint_based_analysis.efm.EFM import EFM, TYPE_PATTERN
from mod_sbml.sbml.submodel_manager import remove_unused_species, compress_reaction_participants
from mod_sbml.sbml.sbml_manager import create_reaction

__author__ = 'anna'


def detect_cliques(id2fm, model, min_clique_size=2):
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
    r_id2fm_ids = defaultdict(set)

    symbol2rc = lambda r: (r[1:], -1) if '-' == r[0] else (r, 1)
    rc2symbol = lambda r, c: ('-%s' % r) if c < 0 else r

    for fm_id, fm in id2fm.iteritems():
        r_id2coeff = fm.to_r_id2coeff(binary=False)
        r_ids = list(r_id2coeff.iterkeys())
        i = 0

        update_r = lambda r: rc2symbol(r, r_id2coeff[r])
        for rr in r_ids:
            i += 1

            r_0 = update_r(rr)
            r_id2fm_ids[r_0].add(fm_id)

            for r in r_ids[i:]:
                r_1 = update_r(r)
                r_pair = (r_0, r_1) if r_0 <= r_1 else (r_1, r_0)
                ratio = r_id2coeff[rr] / r_id2coeff[r] if r_0 <= r_1 else r_id2coeff[r] / r_id2coeff[rr]

                r_id_pair2count.update(
                    {(r_pair, abs(ratio)): 1})

    gr = Graph()
    r_id_pair2ratio = {}
    for ((r_0, r_1), ratio), count in r_id_pair2count.iteritems():
        if count == len(r_id2fm_ids[r_0]) and count == len(r_id2fm_ids[r_1]) and count > 1:
            gr.add_edge(r_0, r_1)
            r_id_pair2ratio[(r_0, r_1)] = ratio

    def clique2r_id2coeff(clique):
        clique = sorted(clique)
        rr = clique[0]
        clique = clique[1:]
        r_0, c_0 = symbol2rc(rr)
        r_id2coeff = {r_0: c_0}
        for r in clique:
            r_1, c_1 = symbol2rc(r)
            ratio = r_id_pair2ratio[(rr, r)]
            r_id2coeff[r_1] = ratio * c_1
        return r_id2coeff

    clique2key = {}
    key2cliques = defaultdict(list)
    for clique in (clique for clique in find_cliques(gr) if len(clique) >= min_clique_size):
        r_id2coeff = clique2r_id2coeff(clique)
        r_id2st, p_id2st = compress_reaction_participants(model, r_id2coeff)
        key = tuple(sorted(r_id2st.iteritems())), tuple(sorted(p_id2st.iteritems()))
        clique = EFM(r_id2coeff=r_id2coeff, type=TYPE_PATTERN)
        clique2key[clique] = key
        key2cliques[key].append(clique)
    id2key = dict(zip(xrange(0, len(key2cliques)), sorted(key2cliques.iterkeys(), key=lambda k: -len(key2cliques[k]))))
    key2id = {k: k_id for (k_id, k) in id2key.iteritems()}
    id2clique = dict(zip(xrange(0, len(clique2key)), sorted(clique2key.iterkeys(),
                                                            key=lambda cl: (key2id[clique2key[cl]], -len(cl)))))
    clique2id = {cl: cl_id for (cl_id, cl) in id2clique.iteritems()}
    return id2clique, {clique_id: r_id2fm_ids[rc2symbol(*clique.get_sample_r_id_coefficient_pair())]
                       for (clique_id, clique) in id2clique.iteritems()}, \
           id2key, {key: [clique2id[cl] for cl in cliques] for (key, cliques) in key2cliques.iteritems()}


def clique2lumped_reaction(id2clique, id2key, key2cl_ids, in_sbml, out_sbml):
    doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = doc.getModel()
    cl_id2r_ids = {}

    for cl_id, clique in id2clique.iteritems():
        r_id2coeff = clique.to_r_id2coeff()
        cl_id2r_ids[cl_id] = [(r_id, 1 if coeff > 0 else -1) for (r_id, coeff) in r_id2coeff.iteritems()]

    r_id2new_r_id = {}
    cl_id2new_r_id = {}
    r_id2direction = defaultdict(lambda: [False, False])
    for k_id, key in id2key.iteritems():
        cl_ids = key2cl_ids[key]
        name = ('lumped reaction groups of type %d' % k_id) if len(cl_ids) > 1 else \
            'lumped reaction group %d' % next(iter(cl_ids))
        id_ = ('rg_t_%d' % k_id) if len(cl_ids) > 1 else 'rg_%d' % next(iter(cl_ids))
        new_r_id = create_reaction(model, dict(key[0]), dict(key[1]), reversible=False,
                                   id_=id_, name=name).getId()
        for cl_id in cl_ids:
            cl_id2new_r_id[cl_id] = new_r_id
            for r_id, coeff in cl_id2r_ids[cl_id]:
                r_id2new_r_id[r_id, coeff] = new_r_id
                r_id2direction[r_id][0 if coeff > 0 else 1] = True

    for r_id, [d1, d2] in r_id2direction.iteritems():
        if d1 and (d2 or not model.getReaction(r_id).getReversible()):
            model.removeReaction(r_id)

    remove_unused_species(model)
    model.setId('%s_merged_cliques' % model.getId())
    model.setName('%s_merged_cliques' % model.getName())
    libsbml.SBMLWriter().writeSBMLToFile(doc, out_sbml)
    return r_id2new_r_id, cl_id2new_r_id


def fold_efms(id2efm, id2clique, cl_id2new_r_id):
    efm_id2folded_efm = {efm_id: efm.fold_cliques(id2clique, cl_id2new_r_id) for (efm_id, efm) in id2efm.iteritems()}
    folded_efm2ids = invert_map(efm_id2folded_efm)
    id2folded_efm = {}
    folded_efm_id2efm_ids = {}
    i = 0
    for efm_id, folded_efm in efm_id2folded_efm.iteritems():
        if efm_id != folded_efm.id:
            folded_efm.id = 'Folded EFM %d' % i
            i += 1
        id2folded_efm[folded_efm.id] = folded_efm
        folded_efm_id2efm_ids[folded_efm.id] = folded_efm2ids[folded_efm]
    return id2folded_efm, folded_efm_id2efm_ids






