import logging
import os

from mod_cobra.html import describe
from mod_cobra.efm.serialization import THICK_DELIMITER, r_id2c_to_string, THIN_DELIMITER, \
    write_detailed_r_id2c

__author__ = 'anna'


def serialize_communities(S, id2cluster, id2intersection, id2imp_rns, path):
    comm_txt = os.path.join(path, 'Pathway_communities.txt')
    with open(comm_txt, 'w+') as f:
        f.write('Analysed %d pathways. ' % len(S.efm_id2i))
        all_fm_intersection = S.get_efm_intersection()

        get_key = lambda r_id: (1 if r_id in all_fm_intersection else 0, (None, 0),
                                (1 if r_id in all_fm_intersection else 0, 1), r_id)

        if all_fm_intersection:
            f.write('All pathways contain following %d reactions:\n\t%s.\n\n'
                    % (len(all_fm_intersection),
                       r_id2c_to_string(all_fm_intersection, binary=True, get_key=get_key)))
            f.write(THICK_DELIMITER)

        f.write('Found %d following communities.\n\n' % len(id2cluster))
        for clu_id in sorted(id2cluster.iterkeys()):
            f.write(THIN_DELIMITER)
            cluster = id2cluster[clu_id]
            intersection = id2intersection[clu_id]
            hidden_r_ids = set(intersection.iterkeys())
            imp_rns = id2imp_rns[clu_id]
            f.write('%s of %d following pathways:\n\n\t%s\n\n'
                    % (clu_id, len(cluster), ', '.join(sorted(cluster))))
            if intersection:
                f.write('all of which contain following reactions:\n\n\t%s\n\n'
                        % (r_id2c_to_string(intersection, binary=True, get_key=get_key)))
            if imp_rns and imp_rns != intersection:
                if intersection:
                    f.write('which form a reaction community, which also contains following reactions:\n\n\t %s\n\n'
                            % (r_id2c_to_string(imp_rns, binary=True, hidden_r_ids=hidden_r_ids)))
                else:
                    f.write('its largest reaction community contains following reactions:\n\n\t %s\n\n'
                            % (r_id2c_to_string(imp_rns, binary=True, hidden_r_ids=hidden_r_ids)))
    return comm_txt


def serialize_n_longest_communities(n, id2cluster, id2intersection, id2imp_rns, model, path):
    limit = min(n if n is not None and n >= 0 else len(id2cluster), len(id2cluster))
    rg_data = [(cl_id, len(id2cluster[cl_id]),
                serialize_community(cl_id, id2cluster[cl_id], id2intersection[cl_id], id2imp_rns[cl_id], model, path))
               for cl_id in sorted(id2cluster.iterkeys(), key=lambda cl_id: (-len(id2cluster[cl_id]), cl_id))[0: limit]]
    return limit, rg_data


def serialize_community(cl_id, cluster, intersection, imp_r_id2c, model, path):
    cluster_txt = os.path.join(path, 'Pathway_community_%s.txt' % cl_id)
    with open(cluster_txt, 'w+') as f:
        f.write('Community %s of %d following pathways:\n\n\t%s\n\n'
                % (cl_id, len(cluster), ', '.join(sorted(cluster))))
        f.write(THIN_DELIMITER)
        if intersection:
            f.write('Reactions that participate in all pathways of this community (%d):\n\n' % len(intersection))
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


def serialize(model, S, id2cluster, id2intersection, id2imp_rns, path, get_f_path):
    logging.info('Serializaing reaction communities...')
    communities_txt = serialize_communities(S, id2cluster, id2intersection, id2imp_rns, path)
    limit, community_data = \
        serialize_n_longest_communities(len(id2cluster), id2cluster, id2intersection, id2imp_rns, model, path)
    fm_block = describe('fm_block.html', element_num=limit, characteristics='longest',
                        element_name='pathway community',
                        fms=[(cl_id, cl_len, get_f_path(cl_txt)) for (cl_id, cl_len, cl_txt) in community_data],
                        all=limit == len(id2cluster))
    return describe('communities.html', community_num=len(id2cluster), description_filepath=get_f_path(communities_txt),
                    selected_community_block=fm_block)
