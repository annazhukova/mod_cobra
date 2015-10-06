from mod_cobra.constraint_based_analysis.efm.serialization import THICK_DELIMITER, r_id2c_to_string, THIN_DELIMITER

__author__ = 'anna'


def serialize_communities(S, id2cluster, id2intersection, id2imp_rns, path):
    with open(path, 'w+') as f:
        f.write('Analysed %d pathways. ' % len(S.efm_id2i))
        all_fm_intersection = S.pws.get_efm_intersection()

        get_key = lambda r_id: (1 if r_id in all_fm_intersection else 0, (None, 0),
                                (1 if r_id in all_fm_intersection else 0, 1), r_id)

        if all_fm_intersection:
            f.write('All pathways contain following %d reactions:\n\t%s.\n\n'
                    % (len(all_fm_intersection),
                       r_id2c_to_string(all_fm_intersection, binary=True, get_key=get_key)))
            f.write(THICK_DELIMITER)
        else:
            f.write('\n\n')
        f.write(THICK_DELIMITER)

        f.write('Found %d communities\n\n' % len(id2cluster))
        for clu_id in sorted(id2cluster.iterkeys()):
            f.write(THIN_DELIMITER)
            cluster = id2cluster[clu_id]
            intersection = id2intersection[clu_id]
            hidden_r_ids = set(intersection.iterkeys())
            imp_rns = id2imp_rns[clu_id]
            f.write('Community %d contains %d following pathways:\n\n\t%s\n\n'
                    % (clu_id, len(cluster), ', '.join(str(it) for it in sorted(cluster))))
            f.write('all of which contain following reactions:\n\n\t%s\n\n'
                    % (r_id2c_to_string(intersection, binary=True, get_key=get_key)))
            f.write('which form a reaction community, which also contains following reactions:\n\n\t %s\n\n'
                    % (r_id2c_to_string(imp_rns, binary=True, hidden_r_ids=hidden_r_ids)))