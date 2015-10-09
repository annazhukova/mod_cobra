from collections import Counter
import logging

from igraph import Graph
from louvain import find_partition

__author__ = 'anna'


def detect_fm_communities(S, r_id2w):
    g = Graph(directed=False)
    g.add_vertices(len(S.efm_id2i))

    fm_ids = list(S.efm_id2i.iterkeys())
    fm_id_pair2count = {}
    i = 1
    for fm_id_1 in fm_ids:
        for fm_id_2 in fm_ids[i:]:
            fm_id_pair2count[(fm_id_1, fm_id_2)] = \
                sum(r_id2w[r_id] for r_id in S.pws.get_efm_intersection(efm_ids=(fm_id_1, fm_id_2)))
        i += 1

    avg_intersection = sum(fm_id_pair2count.itervalues()) / len(fm_id_pair2count)
    logging.info("Going to detect communities in FMs, using %d as a threshold." % avg_intersection)

    for (fm_id_1, fm_id_2), w in fm_id_pair2count.iteritems():
        if w >= avg_intersection:
            g.add_edge(S.efm_id2i[fm_id_1], S.efm_id2i[fm_id_2], weight=w - avg_intersection)

    partition = find_partition(graph=g, method='Modularity', weight='weight')
    i2efm_id = {i: efm_id for (efm_id, i) in S.efm_id2i.iteritems()}
    return dict(zip(('pc_%d' % it for it in xrange(0, len(partition))),
                    ([i2efm_id[i] for i in cluster] for cluster in partition)))


def detect_reaction_community(S, efm_ids, selected_fm):

    g = Graph(directed=False)
    # the id of the reversed version of a reaction r_id reversed would be len(S.r_id2i) + S.r_id2i[r_id]
    g.add_vertices(2 * len(S.r_id2i))

    r_id_pair2count = Counter()
    for efm_id in efm_ids:
        r_id2coeff = S.pws.get_r_id2coeff(efm_id)
        r_ids = [S.r_id2i[r_id] + (0 if coeff > 0 else len(S.r_id2i)) for (r_id, coeff) in r_id2coeff.iteritems()]
        i = 0
        for r_id in r_ids:
            i += 1
            for r_id2 in r_ids[i:]:
                r_id_pair2count.update({tuple(sorted([r_id, r_id2])): 1})

    avg_intersection = sum(r_id_pair2count.itervalues()) / len(r_id_pair2count)
    logging.info("Going to detect the largest reaction community, using %d as a threshold." % avg_intersection)

    for (r_i1, r_i2), w in r_id_pair2count.iteritems():
        if w >= avg_intersection:
            g.add_edge(r_i1, r_i2, weight=w - avg_intersection)

    partition = find_partition(graph=g, method='Modularity', weight='weight')

    selected_fm_indices = {S.r_id2i[r_id] + (0 if coeff > 0 else len(S.r_id2i))
                           for (r_id, coeff) in selected_fm.iteritems()}
    cluster = next(it for it in partition if set(it) & selected_fm_indices)
    i2r_id = {i: (r_id, 1) for (r_id, i) in S.r_id2i.iteritems()}
    i2r_id.update({len(S.r_id2i) + i: (r_id, -1) for (r_id, i) in S.r_id2i.iteritems()})
    return dict(i2r_id[i] for i in cluster)
