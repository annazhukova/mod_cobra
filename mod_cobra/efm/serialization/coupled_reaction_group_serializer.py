import os

from mod_cobra.html import describe
from mod_cobra.efm.serialization import THIN_DELIMITER, THICK_DELIMITER, write_inputs_outputs, \
    r_id2c_to_string, write_detailed_r_id2c

__author__ = 'anna'


def serialize_coupled_reaction_groups(model, path, **kwargs):
    S_coupled, S_no_duplicates = kwargs['S_coupled'], kwargs['S_no_duplicates']
    lr_ids = set(S_coupled.gr_id2r_id2c.iterkeys())
    merged_lr_ids = set(S_no_duplicates.r_id2gr_id.iterkeys()) & lr_ids
    non_merged_lr_ids = lr_ids - merged_lr_ids
    mlr_ids = set(S_no_duplicates.r_id2gr_id[lr_id] for lr_id in merged_lr_ids)
    mr_ids = set(S_no_duplicates.gr_id2r_id2c.iterkeys()) - mlr_ids

    def write_coupled_rn_group(lr_id, inverted=False, tab=''):
        r_id2c = S_coupled.gr_id2r_id2c[lr_id]
        f.write('%sCoupled reaction group %s of length %d%s:\n\n\t\t%s\n\n'
                % (tab, lr_id, len(r_id2c), ' (Inputs and Outputs are inverted)' if inverted else '',
                   r_id2c_to_string(r_id2c)))
        efm_ids = sorted(S_coupled.pws.get_efm_ids_by_r_id(lr_id))
        f.write('%sFound in %d EFMs: %s.\n\n' % (tab, len(efm_ids), ', '.join((str(efm_id) for efm_id in efm_ids))))

    rgroups_txt = os.path.join(path, 'coupled_reactions.txt')
    with open(rgroups_txt, 'w+') as f:
        f.write('Found %d coupled reaction groups of %d types (Types are based on inputs and outputs).\n\n'
                % (len(lr_ids), len(mlr_ids) + len(non_merged_lr_ids)))
        f.write(THICK_DELIMITER)

        # groups that include coupled reactions
        for mr_id in sorted(mlr_ids, key=lambda mr_id: (-len(S_no_duplicates.gr_id2r_id2c[mr_id]), mr_id)):
            r_id2c = S_no_duplicates.gr_id2r_id2c[mr_id]
            cur_r_ids = set(r_id2c.iterkeys())
            cur_lr_ids = cur_r_ids & lr_ids
            cur_non_lr_ids = cur_r_ids - lr_ids

            write_inputs_outputs(f, model, S_no_duplicates.st_matrix.get_inputs_outputs(mr_id))

            f.write('Type %s common to %d reaction group%s%s:\n\n'
                    % (mr_id, len(cur_lr_ids), 's' if len(cur_lr_ids) != 1 else '',
                       (' and %d reaction%s' % (len(cur_non_lr_ids), 's' if len(cur_non_lr_ids) > 1 else ''))
                       if cur_non_lr_ids else ''))
            f.write(THIN_DELIMITER)
            for lr_id in sorted(cur_lr_ids, key=lambda lr_id: (-len(S_coupled.gr_id2r_id2c[lr_id]), lr_id)):
                write_coupled_rn_group(lr_id, inverted=r_id2c[lr_id] < 0, tab='\t')
                f.write('\t%s' % THIN_DELIMITER)

            if cur_non_lr_ids:
                f.write('\tReaction%s %s also ha%s the same structure.\n\n'
                        % ('s' if len(cur_non_lr_ids) != 1 else '',
                           ', '.join((('-%s' % r_id) if r_id2c[r_id] < 0 else r_id for r_id in sorted(cur_non_lr_ids))),
                           's' if len(cur_non_lr_ids) == 1 else 've'))
            f.write(THICK_DELIMITER)

        # coupled reaction that do not belong to any group
        for lr_id in sorted(non_merged_lr_ids, key=lambda lr_id: (-len(S_coupled.gr_id2r_id2c[lr_id]), lr_id)):
            write_inputs_outputs(f, model, S_no_duplicates.st_matrix.get_inputs_outputs(lr_id))
            write_coupled_rn_group(lr_id)
            f.write(THICK_DELIMITER)

        # non-coupled reaction groups
        f.write(THICK_DELIMITER)
        for mr_id in sorted(mr_ids, key=lambda mr_id: (-len(S_no_duplicates.gr_id2r_id2c[mr_id]), mr_id)):
            write_inputs_outputs(f, model, S_no_duplicates.st_matrix.get_inputs_outputs(mr_id))
            r_id2c = S_no_duplicates.gr_id2r_id2c[mr_id]
            f.write('Type %s is common to reaction%s %s.\n\n'
                    % (mr_id, 's' if len(r_id2c) != 1 else '',
                       ', '.join((('-%s' % r_id) if r_id2c[r_id] < 0 else r_id for r_id in sorted(r_id2c.iterkeys())))))
            f.write(THICK_DELIMITER)
    return rgroups_txt, len(lr_ids), len(mlr_ids) + len(non_merged_lr_ids)


def serialize_n_longest_coupled_reaction_groups(model, path, **kwargs):
    S = kwargs['S_coupled']
    n = kwargs['n'] if 'n' in kwargs else 5
    limit = min(n if n is not None and n >= 0 else len(S.gr_id2r_id2c), len(S.gr_id2r_id2c))
    rg_data = [(cl_id, len(S.gr_id2r_id2c[cl_id]), serialize_coupled_reaction_group(cl_id, S, model, path))
               for cl_id in sorted(S.gr_id2r_id2c.iterkeys(),
                                   key=lambda cl_id: (-len(S.gr_id2r_id2c[cl_id]), cl_id))[0: limit]]
    return limit, rg_data


def serialize_coupled_reaction_group(cl_id, S, model, path):
    rg_txt = os.path.join(path, 'Coupled_reactions_%s.txt' % cl_id)
    with open(rg_txt, 'w+') as f:
        r_id2c = S.gr_id2r_id2c[cl_id]
        efm_ids = sorted(S.pws.get_efm_ids_by_r_id(cl_id))
        f.write('Reaction group %s of length %d found in %d EFMs:\n\n\t%s\n\n'
                % (cl_id, len(r_id2c), len(efm_ids), ', '.join((str(efm_id) for efm_id in efm_ids))))
        f.write(THIN_DELIMITER)
        write_detailed_r_id2c(model, r_id2c, f)
    return rg_txt


def serialize(model, path, get_f_path, **kwargs):
    rgroups_txt, rg_num, type_num = serialize_coupled_reaction_groups(model, path, **kwargs)
    limit, rg_data = serialize_n_longest_coupled_reaction_groups(model, path, n=5, **kwargs)
    fm_block = describe('fm_block.html', element_num=limit, characteristics='longest',
                        element_name='reaction group',
                        fms=[(rg_id, rg_len, get_f_path(rg_txt)) for (rg_id, rg_len, rg_txt) in rg_data],
                        all=limit == rg_num)
    return describe('cliques.html', clique_num=rg_num, description_filepath=get_f_path(rgroups_txt),
                    selected_clique_block=fm_block, type_num=type_num)

