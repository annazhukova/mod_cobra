import os
import numpy as np

from numpy_efm_manager import get_yield, get_len, get_r_id2coeff, get_control_efficiency, get_efm_intersection

THICK_DELIMITER = '==============================================================\n\n'
THIN_DELIMITER = '------------------------------------------\n\n'
TINY_DELIMITER = '-------------\n\n'

SIMPLE_PATTERN_SORTER = lambda p_id: -p_id

__author__ = 'anna'


def coeff_to_binary(coeff):
    return 1 if coeff > 0 else -1


def coeff_to_string(coeff, binary=False):
    if binary or abs(coeff) == 1:
        return '' if coeff > 0 else '-'
    return '%g ' % coeff


def r_id2coeff_to_string(r_id2coeff, binary=False, get_key=lambda r_id: (0, (None, 0), (None, 0), r_id),
                         hidden_r_ids=None):
    result = []
    clique_started = None
    started = False
    for r_id in sorted(r_id2coeff.iterkeys(), key=lambda r_id: get_key(r_id), reverse=True):
        coeff = r_id2coeff[r_id]
        _, (group, gc), (clique, cc), _ = get_key(r_id)
        if clique != clique_started:
            if clique_started:
                result.append(')')
                clique_started = None
                started = True

            if clique and (not hidden_r_ids or r_id not in hidden_r_ids):
                clique_started = clique
                gr_c = coeff_to_string(((coeff / cc) / gc) if gc else 0, binary)
                cl_c = coeff_to_string(coeff / cc if not gc else cc, binary)
                result.append('%s%s: ('
                              % ('\t' if started else '',
                                 ('%s%s(%s%s)' % (gr_c, group, cl_c, clique)) if group else '%s%s' % (cl_c, clique)))
                started = False
        if not hidden_r_ids or r_id not in hidden_r_ids:
            coeff = coeff_to_string(cc if cc else coeff, binary)
            result.append('%s%s%s' % ('\t' if started else '', coeff, r_id))
            started = True

    if clique_started:
        result.append(')')

    return ''.join(result)


def write_metabolites(m_id2st, model, f, prefix=''):
    f.write('%s%s;\n\n' % (prefix, ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                             for (m_id, st) in m_id2st)))


def write_inputs_outputs(f, model, (r_id2st, p_id2st)):
    write_metabolites(r_id2st, model, f, prefix='Inputs: ')
    write_metabolites(p_id2st, model, f, prefix='Outputs: ')


def serialize_r_cliques(model, path, **kwargs):
    S_initial, S_coupled, S_no_duplicates = kwargs['S_initial'], kwargs['S_coupled'], kwargs['S_no_duplicates']
    lr_ids = set(S_coupled.gr_id2r_id2c.iterkeys())
    merged_lr_ids = set(S_no_duplicates.r_id2gr_id.iterkeys()) & lr_ids
    non_merged_lr_ids = lr_ids - merged_lr_ids
    mlr_ids = set(S_no_duplicates.r_id2gr_id[lr_id] for lr_id in merged_lr_ids)
    mr_ids = set(S_no_duplicates.gr_id2r_id2c.iterkeys()) - mlr_ids

    def write_coupled_rn_group(lr_id, inverted=False, tab=''):
        r_id2c = S_coupled.gr_id2r_id2c[lr_id]
        f.write('%sCoupled reaction group %s of length %d%s:\n\n\t\t%s\n\n'
                % (tab, lr_id, len(r_id2c), ' (Inputs and Outputs are inverted)' if inverted else '',
                   r_id2coeff_to_string(r_id2c)))
        sample_r_index = S_initial.r_id2i[next(r_id2c.iterkeys())]
        f.write('%sFound in %d EFMs: %s.\n\n'
                % (tab, np.count_nonzero(S_initial.V[sample_r_index, :]),
                   ', '.join((str(efm_id) for (efm_id, i) in sorted(S_initial.efm_id2i.iteritems(),
                                                                    key=lambda (efm_id, _): efm_id)
                              if S_initial.V[sample_r_index, i]))))

    with open(os.path.join(path, 'coupled_reactions.txt'), 'w+') as f:
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


def serialize_efms(model, path, **kwargs):
    S_initial, S_coupled, S_no_duplicates, S_merged = \
        kwargs['S_initial'], kwargs['S_coupled'], kwargs['S_no_duplicates'], kwargs['S_merged']
    out_m_id, in_m_id, out_r_id = kwargs['out_m_id'], kwargs['in_m_id'], kwargs['out_r_id']

    all_fm_intersection = get_efm_intersection(S_initial.V, S_initial.r_id2i)
    all_fm_intersection_folded = get_efm_intersection(S_no_duplicates.V, S_no_duplicates.r_id2i)

    def get_key(r_id):
        if r_id in S_coupled.r_id2gr_id:
            gr_id = S_coupled.r_id2gr_id[r_id]
            gc = S_coupled.gr_id2r_id2c[gr_id][r_id]
        elif r_id in S_no_duplicates.r_id2gr_id:
            gr_id = S_no_duplicates.r_id2gr_id[r_id]
            gc = S_no_duplicates.gr_id2r_id2c[gr_id][r_id]
        else:
            gr_id = None
            gc = 0
        if gr_id and gr_id in S_no_duplicates.r_id2gr_id:
            mr_id = S_no_duplicates.r_id2gr_id[gr_id]
            mc = S_no_duplicates.gr_id2r_id2c[mr_id][gr_id]
        else:
            mr_id = None
            mc = 0

        in_intersection = 1 if r_id in all_fm_intersection or r_id in all_fm_intersection_folded else 0
        return in_intersection, (mr_id, mc), (gr_id, gc), r_id

    with open(os.path.join(path, 'pathways.txt'), 'w+') as f:
        f.write('Found %d EFMs, folded into %d EFMs, grouped into %d pathways:\n\n'
                % (len(S_initial.efm_id2i), len(S_no_duplicates.efm_id2i), len(S_merged.efm_id2i)))
        f.write(THICK_DELIMITER)

        if all_fm_intersection and len(all_fm_intersection):
            f.write('All EFMs contain following %d reactions:\n\t%s.\n\n'
                    % (len(all_fm_intersection),
                       r_id2coeff_to_string(all_fm_intersection, binary=True, get_key=get_key)))
            f.write(THICK_DELIMITER)

        for fm_id in sorted(S_merged.efm_id2i.iterkeys()):
            write_pathway(model=model, pathway_id=fm_id,
                          S=S_initial, S_folded=S_no_duplicates, S_merged=S_merged,
                          out_m_id=out_m_id, in_m_id=in_m_id, out_r_id=out_r_id, get_key=get_key, f=f)
            f.write(THICK_DELIMITER)


def write_pathway(model, pathway_id, S, S_folded, S_merged, out_m_id, in_m_id, out_r_id, get_key, f):
    write_inputs_outputs(f, model, S_merged.get_boundary_inputs_outputs(pathway_id))
    v = S_merged.V[:, S_merged.efm_id2i[pathway_id]]
    fm_yield = get_yield(S_merged.N, v, S_merged.m_id2i[out_m_id], S_merged.m_id2i[in_m_id])
    if fm_yield is not None:
        f.write('Yield: %g;\n\n' % fm_yield)

    if pathway_id not in S_merged.gr_id2efm_ids:
        write_folded_efm(f_efm_id=pathway_id, S=S, S_folded=S_folded, out_r_id=out_r_id, get_key=get_key, f=f)
        return

    f.write('Pathway %s of length %d:\n\n\t%s\n\n'
            % (pathway_id, get_len(v), r_id2coeff_to_string(get_r_id2coeff(v, S_merged.r_id2i), get_key=get_key)))
    folded_efm_ids = S_merged.gr_id2efm_ids[pathway_id]
    f.write('\tContains %d elements:\n\n' % len(folded_efm_ids))
    started = False
    for f_efm_id in folded_efm_ids:
        if started:
            f.write('\t\t%s' % THIN_DELIMITER)
        else:
            started = True
        write_folded_efm(f_efm_id=f_efm_id, S=S, S_folded=S_folded,
                         out_r_id=out_r_id, get_key=get_key, f=f, tab='\t\t')


def write_folded_efm(f_efm_id, S, S_folded, out_r_id, get_key, f, tab=''):
    if f_efm_id not in S_folded.gr_id2efm_ids:
        write_efm(efm_id=f_efm_id, S=S, out_r_id=out_r_id, get_key=get_key, f=f, tab=tab)
        return

    v = S_folded.V[:, S_folded.efm_id2i[f_efm_id]]

    f.write('%sFolded EFM %s of length %d:\n\n%s\t%s\n\n'
            % (tab, f_efm_id, get_len(v), tab, r_id2coeff_to_string(get_r_id2coeff(v, S_folded.r_id2i),
                                                                    get_key=get_key)))

    efm_ids = S_folded.gr_id2efm_ids[f_efm_id]
    f.write('%s\tUnfolds into %s' % (tab, ('the following %d EFMs:\n\n' % len(efm_ids)) if len(efm_ids) != 1 else ''))

    started = False
    for efm_id in efm_ids:
        if started:
            f.write('%s\t\t%s' % (tab, TINY_DELIMITER))
        else:
            started = True
        write_efm(efm_id=efm_id, S=S, out_r_id=out_r_id, get_key=get_key, f=f, tab='%s\t\t' % tab,
                  no_first_tab=len(efm_ids) == 1)


def write_efm(efm_id, S, out_r_id, get_key, f, tab='', no_first_tab=False):
    v = S.V[:, S.efm_id2i[efm_id]]
    efficiency = get_control_efficiency(v, S.r_id2i[out_r_id])
    f.write('%sEFM %s of length %d of efficiency %g:\n\n%s\t%s\n\n'
            % ('' if no_first_tab else tab, efm_id, get_len(v), efficiency, tab,
               r_id2coeff_to_string(get_r_id2coeff(v, S.r_id2i), get_key=get_key)))
