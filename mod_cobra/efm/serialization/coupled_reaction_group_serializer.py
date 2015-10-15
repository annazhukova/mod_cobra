import os

from mod_cobra.html import describe
from mod_cobra.efm.serialization import THIN_DELIMITER, THICK_DELIMITER, write_inputs_outputs, \
    r_id2c_to_string, write_detailed_r_id2c

__author__ = 'anna'


def serialize_coupled_reaction_groups(model, path, **kwargs):
    S = kwargs['S']
    coupled_r_ids = S.coupled_rs
    r_types = S.r_types

    initial_efm_ids = [efm_id for efm_id in S.efm_id2i.iterkeys() if efm_id not in S.gr_id2efm_ids]

    def write_coupled_rn_group(cr_id, inverted=False, tab=''):
        r_id2c = S.gr_id2r_id2c[cr_id]
        f.write('%sCoupled reaction group %s of length %d%s:\n\n\t\t%s\n\n'
                % (tab, cr_id, len(r_id2c), ' (Inputs and Outputs are inverted)' if inverted else '',
                   r_id2c_to_string(r_id2c)))
        efm_ids = sorted(S.get_efm_ids_by_r_id(cr_id, efm_ids=initial_efm_ids))
        f.write('%sFound in %d EFMs: %s.\n\n' % (tab, len(efm_ids), ', '.join((str(efm_id) for efm_id in efm_ids))))

    rgroups_txt = os.path.join(path, 'coupled_reactions.txt')
    with open(rgroups_txt, 'w+') as f:
        f.write('Found %d coupled reaction groups.\n\n' % len(coupled_r_ids))
        f.write(THICK_DELIMITER)

        # types that include coupled reactions
        for rtype_id in sorted((rtype_id for rtype_id in r_types
                                if set(S.gr_id2r_id2c[rtype_id].iterkeys()) & coupled_r_ids),
                               key=lambda rtype_id: (-len(S.gr_id2r_id2c[rtype_id]), rtype_id)):
            r_id2c = S.gr_id2r_id2c[rtype_id]
            cur_r_ids = set(r_id2c.iterkeys())
            cur_coupled_ids = cur_r_ids & coupled_r_ids
            cur_simple_ids = cur_r_ids - coupled_r_ids

            write_inputs_outputs(f, model, S.st_matrix.get_inputs_outputs(rtype_id))

            f.write('Type %s common to %d reaction group%s%s:\n\n'
                    % (rtype_id, len(cur_coupled_ids), 's' if len(cur_coupled_ids) != 1 else '',
                       (' and %d reaction%s' % (len(cur_simple_ids), 's' if len(cur_simple_ids) > 1 else ''))
                       if cur_simple_ids else ''))
            f.write(THIN_DELIMITER)
            for cr_id in sorted(cur_coupled_ids, key=lambda lr_id: (-len(S.gr_id2r_id2c[lr_id]), lr_id)):
                write_coupled_rn_group(cr_id, inverted=r_id2c[cr_id] < 0, tab='\t')
                f.write('\t%s' % THIN_DELIMITER)

            if cur_simple_ids:
                f.write('\tReaction%s %s also ha%s the same structure.\n\n'
                        % ('s' if len(cur_simple_ids) != 1 else '',
                           ', '.join((('-%s' % r_id) if r_id2c[r_id] < 0 else r_id for r_id in sorted(cur_simple_ids))),
                           's' if len(cur_simple_ids) == 1 else 've'))
            f.write(THICK_DELIMITER)

        # coupled reaction that form an unique type
        for cr_id in sorted(S.coupled_rs & S.r_ids, key=lambda lr_id: (-len(S.gr_id2r_id2c[lr_id]), lr_id)):
            write_inputs_outputs(f, model, S.st_matrix.get_inputs_outputs(cr_id))
            write_coupled_rn_group(cr_id)
            f.write(THICK_DELIMITER)

        # non-coupled reaction groups
        f.write(THICK_DELIMITER)
        for rtype_id in sorted((rtype_id for rtype_id in r_types
                                if not set(S.gr_id2r_id2c[rtype_id].iterkeys()) & coupled_r_ids),
                               key=lambda rtype_id: (-len(S.gr_id2r_id2c[rtype_id]), rtype_id)):
            write_inputs_outputs(f, model, S.st_matrix.get_inputs_outputs(rtype_id))
            r_id2c = S.gr_id2r_id2c[rtype_id]
            f.write('Type %s is common to reaction%s %s.\n\n'
                    % (rtype_id, 's' if len(r_id2c) != 1 else '',
                       ', '.join((('-%s' % r_id) if r_id2c[r_id] < 0 else r_id for r_id in sorted(r_id2c.iterkeys())))))
            f.write(THICK_DELIMITER)
    return rgroups_txt, len(coupled_r_ids)


def serialize_n_longest_coupled_reaction_groups(model, path, **kwargs):
    S = kwargs['S']
    n = kwargs['n'] if 'n' in kwargs else 5
    initial_efm_ids = [efm_id for efm_id in S.efm_id2i.iterkeys() if efm_id not in S.gr_id2efm_ids]
    limit = min(n if n is not None and n >= 0 else len(S.coupled_rs), len(S.coupled_rs))
    rg_data = [(cl_id, len(S.gr_id2r_id2c[cl_id]),
                serialize_coupled_reaction_group(initial_efm_ids, cl_id, S, model, path))
               for cl_id in sorted(S.coupled_rs, key=lambda cl_id: (-len(S.gr_id2r_id2c[cl_id]), cl_id))[0: limit]]
    return limit, rg_data


def serialize_coupled_reaction_group(initial_efm_ids, cl_id, S, model, path):
    rg_txt = os.path.join(path, 'Coupled_reactions_%s.txt' % cl_id)
    with open(rg_txt, 'w+') as f:
        r_id2c = S.gr_id2r_id2c[cl_id]
        efm_ids = sorted(S.get_efm_ids_by_r_id(cl_id, efm_ids=initial_efm_ids))
        f.write('Reaction group %s of length %d found in %d EFMs:\n\n\t%s\n\n'
                % (cl_id, len(r_id2c), len(efm_ids), ', '.join((str(efm_id) for efm_id in efm_ids))))
        f.write(THIN_DELIMITER)
        write_detailed_r_id2c(model, r_id2c, f)
    return rg_txt


def serialize(model, path, get_f_path, **kwargs):
    rgroups_txt, rg_num = serialize_coupled_reaction_groups(model, path, **kwargs)
    limit, rg_data = serialize_n_longest_coupled_reaction_groups(model, path, n=5, **kwargs)
    fm_block = describe('fm_block.html', element_num=limit, characteristics='longest',
                        element_name='reaction group',
                        fms=[(rg_id, rg_len, get_f_path(rg_txt)) for (rg_id, rg_len, rg_txt) in rg_data],
                        all=limit == rg_num)
    return describe('cliques.html', clique_num=rg_num, description_filepath=get_f_path(rgroups_txt),
                    selected_clique_block=fm_block)

