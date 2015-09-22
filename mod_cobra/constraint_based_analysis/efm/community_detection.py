from collections import Counter
import logging

from igraph import Graph
from louvain import find_partition

from mod_cobra.constraint_based_analysis.efm.EFM import EFM

__author__ = 'anna'


def detect_fm_communities(id2fm, threshold):
    logging.info("Going to detect communities in FMs, using %d as a threshold." % threshold)

    fm_ids = list(id2fm.iterkeys())
    i = 1
    g = Graph(directed=False)
    for fm_id in fm_ids:
        g.add_vertex(fm_id)
    for fm_id_1 in fm_ids:
        for fm_id_2 in fm_ids[i:]:
            w = id2fm[fm_id_1].intersection_len(id2fm[fm_id_2]) - threshold
            if w > 0:
                g.add_edge(fm_id_1, fm_id_2, weight=w)
        i += 1
    partition = find_partition(graph=g, method='Modularity', weight='weight')
    return dict(zip(xrange(0, len(partition)), ([g.vs[v_id]['name'] for v_id in cluster] for cluster in partition)))


def detect_reaction_communities(id2fm, threshold, min_len=2):
    logging.info("Going to detect reaction communities, using %d as a threshold." % threshold)

    sample_fm = next(id2fm.itervalues())
    g = Graph(directed=False)
    for r_id in sample_fm.r_ids:
        g.add_vertex(r_id)
    for r_id in sample_fm.rev_r_ids:
        g.add_vertex('-%s' % r_id)

    r_id_pair2count = Counter()
    for fm_id, fm in id2fm.iteritems():
        r_id2coeff = fm.to_r_id2coeff(binary=True)
        r_ids = [r_id if coeff > 0 else '-%s' % r_id for (r_id, coeff) in r_id2coeff.iteritems()]
        i = 0
        for r_id in r_ids:
            i += 1
            for r_id2 in r_ids[i:]:
                r_id_pair2count.update({tuple(sorted([r_id, r_id2])): 1})
    for (r_id1, r_id2), count in r_id_pair2count.iteritems():
        if count >= threshold:
            g.add_edge(r_id1, r_id2, weight=count - threshold)

    partition = find_partition(graph=g, method='Modularity', weight='weight')
    clusters = [EFM(r_ids=sample_fm.r_ids, rev_r_ids=sample_fm.rev_r_ids, int_size=sample_fm.int_size,
                    r_id2coeff=
                    {r_id if r_id.find('-') != 0 else r_id[1:]: 1 if r_id.find('-') != 0 else -1 for r_id in
                     (g.vs[v_id]['name'] for v_id in cluster)})
                for cluster in partition if len(cluster) >= min_len]
    return dict(zip(xrange(0, len(partition)), clusters))


def detect_reaction_community(id2efm, threshold_percent, selected_fm):
    logging.info("Going to detect the largest reaction community, using %d%% as a threshold." % threshold_percent)

    g = Graph(directed=False)
    for r_id in selected_fm.r_ids:
        g.add_vertex(r_id)
    for r_id in selected_fm.rev_r_ids:
        g.add_vertex('-%s' % r_id)

    threshold = threshold_percent * len(id2efm) / 100

    r_id_pair2count = Counter()
    for efm in id2efm.itervalues():
        r_id2coeff = efm.to_r_id2coeff(binary=True)
        r_ids = [r_id if coeff > 0 else '-%s' % r_id for (r_id, coeff) in r_id2coeff.iteritems()]
        i = 0
        for r_id in r_ids:
            i += 1
            for r_id2 in r_ids[i:]:
                r_id_pair2count.update({tuple(sorted([r_id, r_id2])): 1})
    for (r_id1, r_id2), count in r_id_pair2count.iteritems():
        if count >= threshold:
            g.add_edge(r_id1, r_id2, weight=count - threshold)

    partition = find_partition(graph=g, method='Modularity', weight='weight')
    max_fm = cluster2fm(g, max(partition, key=len), selected_fm)
    if not max_fm.intersection_len(selected_fm):
        max_fm = next((fm for fm in (cluster2fm(g, cluster, selected_fm) for cluster in partition)
                       if fm.intersection_len(selected_fm)), max_fm)
    return max_fm


def cluster2fm(g, cluster, sample_efm):
    return EFM(r_ids=sample_efm.r_ids, rev_r_ids=sample_efm.rev_r_ids, int_size=sample_efm.int_size,
               r_id2coeff={r_id if r_id.find('-') != 0 else r_id[1:]: 1 if r_id.find('-') != 0 else -1 for r_id in
                           (g.vs[v_id]['name'] for v_id in cluster)})
