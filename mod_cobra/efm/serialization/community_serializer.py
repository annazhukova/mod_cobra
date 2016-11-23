from collections import Counter
import logging
import os
from functools import reduce

from mod_cobra.html import describe
from mod_cobra.efm.serialization import THIN_DELIMITER, \
    write_detailed_r_id2c
from mod_sbml.serialization import format_m_name
from mod_cobra.sbml.mapping.metabolite_classifier import classify_m_ids

__author__ = 'anna'


def serialize_communities(model, S, id2cluster, path, m_id2ch_id, chebi):
    comm_txt = os.path.join(path, 'Pathway_communities.txt')
    with open(comm_txt, 'w+') as f:
        f.write('Analysed %d EFMs, found %d following communities:\n\n' % (len(S.efm_id2i), len(id2cluster)))
        for clu_id in sorted(id2cluster.keys()):
            f.write(THIN_DELIMITER)
            cluster = id2cluster[clu_id]
            # intersection = S.get_efm_intersection(cluster)
            # hidden_r_ids = set(intersection.keys())
            # imp_rns = detect_reaction_community(S, cluster, intersection)
            r_id2count, p_id2count = Counter(), Counter()
            for fm_id in cluster:
                r2st, p2st = S.get_boundary_inputs_outputs(fm_id)
                r_id2count.update({it: 1 for it in r2st.keys()})
                p_id2count.update({it: 1 for it in p2st.keys()})
            class2rs = classify_m_ids(r_id2count.keys(), m_id2ch_id, chebi=chebi)
            class2ps = classify_m_ids(p_id2count.keys(), m_id2ch_id, chebi=chebi)
            # pws = count_pathways(cluster, efm_id2pws)
            f.write('%s of %d following pathways:\n\n\t%s\n\n'
                    % (clu_id, len(cluster), ', '.join(sorted(cluster))))
            # if intersection:
            #     f.write('Reactions contained by all of these EFMs:\n\n\t%s\n\n'
            #             % (r_id2c_to_string(intersection, binary=True, get_key=get_key)))
            # if imp_rns and imp_rns != intersection:
            #     if intersection:
            #         f.write('The largest reaction community contains also:\n\n\t %s\n\n'
            #                 % (r_id2c_to_string(imp_rns, binary=True, hidden_r_ids=hidden_r_ids)))
            #     else:
            #         f.write('The largest reaction community contains:\n\n\t %s\n\n'
            #                 % (r_id2c_to_string(imp_rns, binary=True, hidden_r_ids=hidden_r_ids)))
            def format_class2ms(class2m_ids, m_id2count):
                format_m_ids = lambda m_ids: \
                    ', '.join(('%s [%d]' % (format_m_name(model.getSpecies(m_id), model, False, False), m_id2count[m_id])
                               for m_id in sorted(m_ids, key=lambda m_id: (-m_id2count[m_id], m_id))))
                result = '\n\t'.join('%sS: %s' % (name.upper(), format_m_ids(m_ids)) for (name, m_ids) in class2m_ids.items())
                others = set(m_id2count.keys()) - reduce(lambda s1, s2: s1 | s2, class2m_ids.values(), set())
                if others:
                    result += ('\n\tOTHER: ' if result else '\tOTHER: ') + format_m_ids(others)
                return result

            if r_id2count:
                f.write('Input metabolites [numbers of EFMs where they are present] are:\n\n\t%s\n\n'
                        % format_class2ms(class2rs, r_id2count))
            if p_id2count:
                f.write('Output metabolites [numbers of EFMs where they are present] are:\n\n\t%s\n\n'
                        % format_class2ms(class2ps, p_id2count))
            # if pws:
            #     f.write('Subsystems [numbers of EFMs where they are present] are:\n\n\t%s\n\n'
            #             % ', '.join(('%s [%d]' % (pw, count)
            #                          for (pw, count) in sorted(pws.items(), key=lambda (pw, c): (-c, pw)))))
    return comm_txt


def serialize_n_longest_communities(n, id2cluster, id2intersection, id2imp_rns, model, path):
    limit = min(n if n is not None and n >= 0 else len(id2cluster), len(id2cluster))
    rg_data = [(cl_id, len(id2cluster[cl_id]),
                serialize_community(cl_id, id2cluster[cl_id], id2intersection[cl_id], id2imp_rns[cl_id], model, path))
               for cl_id in sorted(id2cluster.keys(), key=lambda cl_id: (-len(id2cluster[cl_id]), cl_id))[0: limit]]
    return limit, rg_data


def serialize_community(cl_id, cluster, intersection, imp_r_id2c, model, path):
    cluster_txt = os.path.join(path, 'EFM_community_%s.txt' % cl_id)
    with open(cluster_txt, 'w+') as f:
        f.write('Community %s of %d following EFMs:\n\n\t%s\n\n'
                % (cl_id, len(cluster), ', '.join(sorted(cluster))))
        f.write(THIN_DELIMITER)
        if intersection:
            f.write('Reactions that participate in all EFMs of this community (%d):\n\n' % len(intersection))
            f.write(THIN_DELIMITER)
            write_detailed_r_id2c(model, intersection, f)
        if imp_r_id2c:
            f.write(THIN_DELIMITER)
            if intersection:
                f.write('Other reactions that participate in the same reaction community (%d):\n\n' % len(imp_r_id2c))
            else:
                f.write('Reactions that participate in the largest reaction community (%d):\n\n' % len(imp_r_id2c))
            f.write(THIN_DELIMITER)
            write_detailed_r_id2c(model, imp_r_id2c, f)
    return cluster_txt


def serialize(model, S, id2cluster, path, get_f_path, m_id2ch_id, chebi):
    logging.info('Serializing communities...')
    communities_txt = serialize_communities(model, S, id2cluster, path, m_id2ch_id, chebi=chebi)
    # limit, community_data = \
    #     serialize_n_longest_communities(len(id2cluster), id2cluster, id2intersection, id2imp_rns, model, path)
    # fm_block = describe('fm_block.html', element_num=limit, characteristics='longest',
    #                     element_name='pathway community',
    #                     fms=[(cl_id, cl_len, get_f_path(cl_txt)) for (cl_id, cl_len, cl_txt) in community_data],
    #                     all=limit == len(id2cluster))
    return describe('communities.html', community_num=len(id2cluster), description_filepath=get_f_path(communities_txt),
                    selected_community_block=None)
