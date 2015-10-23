import logging

from igraph import Graph
from louvain import find_partition

__author__ = 'anna'


def detect_fm_communities(S, r_id2w):
    if len(S.efm_id2i) == 1:
        return {'pc_0': [next(S.efm_id2i.iterkeys())]}

    g = Graph(directed=False)
    g.add_vertices(len(S.efm_id2i))

    fm_ids = list(S.efm_ids)
    fm_id2rs = {fm_id: set(S.pws.get_r_id2coeff(fm_id, binary=True).iteritems()) for fm_id in S.efm_ids}
    fm_id_pair2count = {}
    i = 1
    for fm_id_1 in fm_ids:
        rs = fm_id2rs[fm_id_1]
        for fm_id_2 in fm_ids[i:]:
            fm_id_pair2count[(fm_id_1, fm_id_2)] = sum(r_id2w[r_id] for r_id, _ in rs & fm_id2rs[fm_id_2])
        i += 1

    avg_intersection = sum(fm_id_pair2count.itervalues()) / len(fm_id_pair2count)
    logging.info("Going to detect communities in FMs, using %d as a threshold." % avg_intersection)

    for (fm_id_1, fm_id_2), w in fm_id_pair2count.iteritems():
        if w >= avg_intersection:
            g.add_edge(S.efm_id2i[fm_id_1], S.efm_id2i[fm_id_2], weight=w - avg_intersection)

    logging.info("Looking for a partition...")
    partition = find_partition(graph=g, method='Modularity', weight='weight')
    logging.info("Found a partition...")
    i2efm_id = {i: efm_id for (efm_id, i) in S.efm_id2i.iteritems()}
    return dict(zip(('pc_%d' % it for it in xrange(0, len(partition))),
                    ([i2efm_id[i] for i in cluster] for cluster in partition)))


def detect_pathway_communities(S):
    if len(S.efm_id2i) == 1:
        return {'pc_0': [next(S.efm_id2i.iterkeys())]}

    g = Graph(directed=False)
    g.add_vertices(len(S.efm_id2i))

    fm_ids = list(S.efm_id2i.iterkeys())
    fm_id_pair2count = {}

    def get_key(fm_id):
        r2st, p2st = S.get_boundary_inputs_outputs(fm_id)
        return set(r2st.iterkeys()), set(p2st.iterkeys())

    fm_id2boundary_ms = {fm_id: get_key(fm_id) for fm_id in S.efm_id2i.iterkeys()}
    i = 1
    for fm_id_1 in fm_ids:
        rs_1, ps_1 = fm_id2boundary_ms[fm_id_1]
        for fm_id_2 in fm_ids[i:]:
            rs_2, ps_2 = fm_id2boundary_ms[fm_id_2]
            fm_id_pair2count[(fm_id_1, fm_id_2)] = len(rs_1 & rs_2) + len(ps_1 & ps_2)
        i += 1

    avg_intersection = sum(fm_id_pair2count.itervalues()) / len(fm_id_pair2count)

    logging.info("Building the FM graph, using %d as a threshold." % avg_intersection)
    edges = []
    weights = []
    for (fm_id_1, fm_id_2), w in fm_id_pair2count.iteritems():
        if w > avg_intersection:
            edges.append((S.efm_id2i[fm_id_1], S.efm_id2i[fm_id_2]))
            weights.append(w - avg_intersection)

    if not edges:
        return {}

    g.add_edges(edges)
    g.es['weight'] = weights

    logging.info("Looking for a partition...")
    partition = find_partition(graph=g, method='Modularity', weight='weight')
    logging.info("Found a partition...")
    i2efm_id = {i: efm_id for (efm_id, i) in S.efm_id2i.iteritems()}
    return dict(zip(('pc_%d' % it for it in xrange(0, len(partition))),
                    ([i2efm_id[i] for i in cluster] for cluster in partition if len(cluster) > 1)))


def detect_reaction_community(S, efm_ids, selected_fm):
    if len(efm_ids) == 1:
        return None

    r_ids = list(S.r_ids)

    g = Graph(directed=False)
    # the id of the reversed version of a reaction r_id reversed would be len(S.r_id2i) + S.r_id2i[r_id]
    g.add_vertices(2 * len(r_ids))

    r_i2fms = {}
    for r_id in r_ids:
        pos_count = S.pws.get_efm_ids(r_id, efm_ids=efm_ids, reversed=False)
        neg_count = S.pws.get_efm_ids(r_id, efm_ids=efm_ids, reversed=True)
        if pos_count:
            r_i2fms[S.r_id2i[r_id]] = pos_count
        if neg_count:
            r_i2fms[S.r_id2i[r_id] + len(r_ids)] = neg_count

    r_i_pair2count = {}
    i = 0
    total_count = 0
    for r_id_1 in r_ids:
        pos_efms_1 = r_i2fms[S.r_id2i[r_id_1]] if S.r_id2i[r_id_1] in r_i2fms else None
        neg_efms_1 = r_i2fms[S.r_id2i[r_id_1] + len(r_ids)] if (S.r_id2i[r_id_1] + len(r_ids)) in r_i2fms else None
        i += 1
        for r_id_2 in r_ids[i:]:
            pos_efms_2 = r_i2fms[S.r_id2i[r_id_2]] if S.r_id2i[r_id_2] in r_i2fms else None
            neg_efms_2 = r_i2fms[S.r_id2i[r_id_2] + len(r_ids)] if (S.r_id2i[r_id_2] + len(r_ids)) in r_i2fms else None
            if pos_efms_1 and pos_efms_2:
                total_count += 1
                value = len(pos_efms_1 & pos_efms_2)
                if value:
                    r_i_pair2count[S.r_id2i[r_id_1], S.r_id2i[r_id_2]] = value
            if pos_efms_1 and neg_efms_2:
                total_count += 1
                value = len(pos_efms_1 & neg_efms_2)
                if value:
                    r_i_pair2count[S.r_id2i[r_id_1], S.r_id2i[r_id_2] + len(r_ids)] = value
            if neg_efms_1 and pos_efms_2:
                total_count += 1
                value = len(neg_efms_1 & pos_efms_2)
                if value:
                    r_i_pair2count[S.r_id2i[r_id_1] + len(r_ids), S.r_id2i[r_id_2]] = value
            if neg_efms_1 and neg_efms_2:
                total_count += 1
                value = len(neg_efms_1 & neg_efms_2)
                if value:
                    r_i_pair2count[S.r_id2i[r_id_1] + len(r_ids), S.r_id2i[r_id_2] + len(r_ids)] = value
    max_intersection = max(r_i_pair2count.itervalues())
    if max_intersection == 0:
        return {}
    avg_intersection = sum(r_i_pair2count.itervalues()) / total_count

    logging.info("Building the reaction graph, using %d as a threshold." % avg_intersection)
    edges = []
    weights = []
    for (r_i1, r_i2), w in r_i_pair2count.iteritems():
        if w > avg_intersection:
            edges.append((r_i1, r_i2))
            weights.append(w - avg_intersection)

    if not edges:
        return {}

    g.add_edges(edges)
    g.es['weight'] = weights

    logging.info("Going to detect the largest reaction community...")
    partition = find_partition(graph=g, method='Modularity', weight='weight')

    if selected_fm:
        selected_fm_indices = {S.r_id2i[r_id] + (0 if c > 0 else len(S.r_id2i)) for (r_id, c) in selected_fm.iteritems()}
        cluster = next(it for it in partition if set(it) & selected_fm_indices)
    else:
        cluster = max(partition, key=lambda it: len(it))
    i2r_id = {i: (r_id, 1) for (r_id, i) in S.r_id2i.iteritems()}
    i2r_id.update({len(S.r_id2i) + i: (r_id, -1) for (r_id, i) in S.r_id2i.iteritems()})
    return dict(i2r_id[i] for i in cluster)
