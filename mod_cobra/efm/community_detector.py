from collections import defaultdict
import logging

from igraph import Graph
from louvain import find_partition
from tarjan import tarjan_iter
from mod_cobra.sbml.mapping.metabolite_classifier import classify_m_ids

__author__ = 'anna'


def detect_communities_by_inputs_of_type(S, superclas, m_id2ch_id, chebi, threshold=75):
    if len(S.efm_id2i) == 1:
        return {'pc_0': [next(S.efm_id2i.keys())]}

    fm_ids = list(S.efm_id2i.keys())

    g = Graph(directed=False)
    g.add_vertices(len(fm_ids))

    def get_key(fm_id):
        r2st, _ = S.get_boundary_inputs_outputs(fm_id)
        class2m_ids = classify_m_ids(r2st.keys(), m_id2ch_id, chebi=chebi)
        return set(class2m_ids[superclas]) if superclas in class2m_ids else set()

    fm_id2boundary_ms = {fm_id: get_key(fm_id) for fm_id in fm_ids}
    edges = []
    weights = []
    i = 1
    for fm_id_1 in fm_ids:
        rs1 = fm_id2boundary_ms[fm_id_1]
        len1 = len(rs1)
        for fm_id_2 in fm_ids[i:]:
            rs2 = fm_id2boundary_ms[fm_id_2]
            len2 = len(rs2)
            common_len = len(rs1 & rs2)
            local_total = max(len1, len2)
            if not local_total or 100 * common_len / local_total > threshold:
                edges.append((S.efm_id2i[fm_id_1], S.efm_id2i[fm_id_2]))
                weights.append((100 - 100 * common_len / local_total) if local_total else 100)
        i += 1

    if not edges:
        return {}

    g.add_edges(edges)
    g.es['weight'] = weights

    logging.info("Looking for a partition...")
    partition = find_partition(graph=g, method='Modularity', weight='weight')
    logging.info("Found a partition...")
    i2efm_id = {i: efm_id for (efm_id, i) in S.efm_id2i.items()}
    return dict(zip(('pc_%d' % it for it in range(0, len(partition))),
                    ([i2efm_id[i] for i in cluster] for cluster in partition if len(cluster) > 1)))

    # clustering = defaultdict(list)
    # for efm_id in fm_ids:
    #     clustering[get_key(efm_id)].append(efm_id)
    #
    # return dict(zip(('pc_%d' % it for it in range(0, len(clustering))),
    #                 (cluster for cluster in clustering if len(cluster) > 1)))


def detect_communities_by_reactions(S, r_id2w):
    if len(S.efm_id2i) == 1:
        return {'pc_0': [next(S.efm_id2i.keys())]}

    g = Graph(directed=False)
    g.add_vertices(len(S.efm_id2i))

    fm_ids = list(S.efm_ids)
    fm_id2rs = {fm_id: set(S.pws.get_r_id2coeff(fm_id, binary=True).items()) for fm_id in S.efm_ids}
    fm_id_pair2count = {}
    i = 1
    for fm_id_1 in fm_ids:
        rs = fm_id2rs[fm_id_1]
        for fm_id_2 in fm_ids[i:]:
            fm_id_pair2count[(fm_id_1, fm_id_2)] = sum(r_id2w[r_id] for r_id, _ in rs & fm_id2rs[fm_id_2])
        i += 1

    avg_intersection = sum(fm_id_pair2count.values()) / len(fm_id_pair2count)
    logging.info("Going to detect communities in FMs, using %d as a threshold." % avg_intersection)

    for (fm_id_1, fm_id_2), w in fm_id_pair2count.items():
        if w >= avg_intersection:
            g.add_edge(S.efm_id2i[fm_id_1], S.efm_id2i[fm_id_2], weight=w - avg_intersection)

    logging.info("Looking for a partition...")
    partition = find_partition(graph=g, method='Modularity', weight='weight')
    logging.info("Found a partition...")
    i2efm_id = {i: efm_id for (efm_id, i) in S.efm_id2i.items()}
    return dict(zip(('pc_%d' % it for it in range(0, len(partition))),
                    ([i2efm_id[i] for i in cluster] for cluster in partition)))


def detect_communities_by_boundary_metabolites(S, cofactors=None, threshold=50):
    if len(S.efm_id2i) == 1:
        return {'pc_0': [next(S.efm_id2i.keys())]}
    if not cofactors:
        cofactors = set()

    g = Graph(directed=False)
    g.add_vertices(len(S.efm_id2i))

    fm_ids = list(S.efm_id2i.keys())

    def get_key(fm_id):
        r2st, p2st = S.get_boundary_inputs_outputs(fm_id)
        return set(r2st.keys()) - cofactors, set(p2st.keys()) - cofactors

    fm_id2boundary_ms = {fm_id: get_key(fm_id) for fm_id in fm_ids}
    edges = []
    weights = []
    i = 1
    for fm_id_1 in fm_ids:
        rs_1, ps_1 = fm_id2boundary_ms[fm_id_1]
        len_1 = len(rs_1) + len(ps_1)
        for fm_id_2 in fm_ids[i:]:
            rs_2, ps_2 = fm_id2boundary_ms[fm_id_2]
            len_2 = len(rs_2) + len(ps_2)
            common_len = len(rs_1 & rs_2) + len(ps_1 & ps_2)
            local_total = max(len_1, len_2)
            if not local_total or 100 * common_len / local_total > threshold:
                edges.append((S.efm_id2i[fm_id_1], S.efm_id2i[fm_id_2]))
                weights.append((100 - 100 * common_len / local_total) if local_total else 100)
        i += 1

    if not edges:
        return {}

    g.add_edges(edges)
    g.es['weight'] = weights

    logging.info("Looking for a partition...")
    partition = find_partition(graph=g, method='Modularity', weight='weight')
    logging.info("Found a partition...")
    i2efm_id = {i: efm_id for (efm_id, i) in S.efm_id2i.items()}
    return dict(zip(('pc_%d' % it for it in range(0, len(partition))),
                    ([i2efm_id[i] for i in cluster] for cluster in partition if len(cluster) > 1)))


def detect_communities_by_subsystems(efm_ids, id2pws, threshold=50):
    # g = Graph(directed=False)
    # g.add_vertices(len(efm_ids))

    efm_ids = sorted(efm_ids)

    # i2efm_id = dict(enumerate(efm_ids))
    # efm_id2i = {efm_id: i for (i, efm_id) in i2efm_id.items()}

    # edges = []
    # weights = []
    i = 1
    graph = defaultdict(list)
    for efm_id_1 in efm_ids:
        pws1 = id2pws[efm_id_1]
        if not pws1:
            continue
        # len_1 = len(pws1)
        for efm_id_2 in efm_ids[i:]:
            pws2 = id2pws[efm_id_2]
            # len_2 = len(pws2)
            # common_len = len(pws1 & pws2)
            # local_total = max(len_1, len_2)
            # if 100 * common_len / local_total > threshold:
            #     edges.append((efm_id2i[efm_id_1], efm_id2i[efm_id_2]))
            #     weights.append(100 - 100 * common_len / local_total)
            if len(pws1 & pws2):
                graph[efm_id_1].append(efm_id_2)
                graph[efm_id_2].append(efm_id_1)
                # edges.append((efm_id2i[efm_id_1], efm_id2i[efm_id_2]))
                # edges.append((efm_id2i[efm_id_2], efm_id2i[efm_id_1]))
                # weights.append(common_len)
        i += 1

    # if not edges:
    #     return {}
    #
    # g.add_edges(edges)
    # # g.es['weight'] = weights

    logging.info("Looking for a partition...")
    # partition = find_partition(graph=g, method='Modularity')
    partition = [cluster for cluster in tarjan_iter(graph) if len(cluster) > 1]
    logging.info("Found a partition...")
    return dict(zip(('sub_pc_%d' % it for it in range(0, len(partition))), partition))


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
    max_intersection = max(r_i_pair2count.values())
    if max_intersection == 0:
        return {}
    avg_intersection = sum(r_i_pair2count.values()) / total_count

    logging.info("Building the reaction graph, using %d as a threshold." % avg_intersection)
    edges = []
    weights = []
    for (r_i1, r_i2), w in r_i_pair2count.items():
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
        selected_fm_indices = {S.r_id2i[r_id] + (0 if c > 0 else len(S.r_id2i)) for (r_id, c) in selected_fm.items()}
        cluster = next(it for it in partition if set(it) & selected_fm_indices)
    else:
        cluster = max(partition, key=lambda it: len(it))
    i2r_id = {i: (r_id, 1) for (r_id, i) in S.r_id2i.items()}
    i2r_id.update({len(S.r_id2i) + i: (r_id, -1) for (r_id, i) in S.r_id2i.items()})
    return dict(i2r_id[i] for i in cluster)
