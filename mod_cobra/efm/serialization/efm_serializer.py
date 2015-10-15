import os
from mod_cobra import round_value

from mod_cobra.html import describe
from mod_cobra.efm.serialization import r_id2c_to_string, write_inputs_outputs, THIN_DELIMITER, \
    THICK_DELIMITER, TINY_DELIMITER, write_detailed_r_id2c

__author__ = 'anna'


def serialize_n_most_efficient_efms(model, path, **kwargs):
    S = kwargs['S_initial']
    out_r_id = kwargs['out_r_id']
    in_m_id = kwargs['in_m_id']
    out_m_id = kwargs['out_m_id']
    n = kwargs['n'] if 'n' in kwargs else 3
    limit = min(n if n is not None and n >= 0 else len(S.efm_id2i), len(S.efm_id2i))
    efm_data = [(efm_id, S.pws.get_len(efm_id), serialize_efm(S, efm_id, in_m_id, out_m_id, out_r_id, model, path))
                for efm_id in sorted(S.efm_id2i.iterkeys(),
                                     key=lambda efm_id:
                                     (-S.pws.get_control_efficiency(efm_id, out_r_id), efm_id))[0: limit]]
    return limit, efm_data


def serialize_efm(S, efm_id, in_m_id, out_m_id, out_r_id, model, path):
    efm_txt = os.path.join(path, '%s.txt' % efm_id)
    with open(efm_txt, 'w+') as f:
        r_id2c = S.pws.get_r_id2coeff(efm_id)
        f.write('EFM %s of length %d, of yield %g, of efficiency %g\n\n'
                % (efm_id, len(r_id2c), round_value(S.get_yield(efm_id, in_m_id, out_m_id)),
                   round_value(S.pws.get_control_efficiency(efm_id, out_r_id))))
        write_inputs_outputs(f, model, S.get_boundary_inputs_outputs(efm_id), in_m_id)
        f.write(THIN_DELIMITER)
        write_detailed_r_id2c(model, r_id2c, f)
    return efm_txt


def serialize_efms(model, path, **kwargs):
    S_initial, S_coupled, S_no_duplicates, S_merged = \
        kwargs['S_initial'], kwargs['S_coupled'], kwargs['S_no_duplicates'], kwargs['S_merged']
    out_m_id, in_m_id, out_r_id = kwargs['out_m_id'], kwargs['in_m_id'], kwargs['out_r_id']

    all_fm_intersection = S_initial.pws.get_efm_intersection()
    all_fm_intersection_folded = S_no_duplicates.pws.get_efm_intersection()

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

    pathways_txt = os.path.join(path, 'pathways.txt')
    with open(pathways_txt, 'w+') as f:
        f.write('Found %d EFMs, folded into %d EFMs, grouped into %d pathways:\n\n'
                % (len(S_initial.efm_id2i), len(S_no_duplicates.efm_id2i), len(S_merged.efm_id2i)))
        f.write(THICK_DELIMITER)

        if all_fm_intersection and len(all_fm_intersection):
            f.write('All EFMs contain following %d reactions:\n\t%s.\n\n'
                    % (len(all_fm_intersection),
                       r_id2c_to_string(all_fm_intersection, binary=True, get_key=get_key)))
            f.write(THICK_DELIMITER)

        for fm_id in sorted(S_merged.efm_id2i.iterkeys()):
            write_pathway(model=model, pathway_id=fm_id,
                          S=S_initial, S_folded=S_no_duplicates, S_merged=S_merged,
                          out_m_id=out_m_id, in_m_id=in_m_id, out_r_id=out_r_id, get_key=get_key, f=f)
            f.write(THICK_DELIMITER)
    return pathways_txt, (len(S_initial.efm_id2i), len(S_no_duplicates.efm_id2i), len(S_merged.efm_id2i))


def write_pathway(model, pathway_id, S, S_folded, S_merged, out_m_id, in_m_id, out_r_id, get_key, f):
    write_inputs_outputs(f, model, S_merged.get_boundary_inputs_outputs(pathway_id), in_m_id)
    fm_yield = S_merged.get_yield(pathway_id, in_m_id, out_m_id)
    if fm_yield is not None:
        f.write('Yield: %g;\n\n' % round_value(fm_yield))

    if pathway_id not in S_merged.gr_id2efm_ids:
        write_folded_efm(f_efm_id=pathway_id, S=S, S_folded=S_folded, out_r_id=out_r_id, in_m_id=in_m_id,
                         out_m_id=out_m_id, get_key=get_key, f=f)
        return

    f.write('%s of length %d:\n\n\t%s\n\n'
            % (pathway_id, S_merged.pws.get_len(pathway_id),
               r_id2c_to_string(S_merged.pws.get_r_id2coeff(pathway_id), get_key=get_key)))
    folded_efm_ids = S_merged.gr_id2efm_ids[pathway_id]
    f.write('\tContains %d elements:\n\n' % len(folded_efm_ids))
    started = False
    for f_efm_id in folded_efm_ids:
        if started:
            f.write('\t\t%s' % THIN_DELIMITER)
        else:
            started = True
        write_folded_efm(f_efm_id=f_efm_id, S=S, S_folded=S_folded, out_r_id=out_r_id, in_m_id=in_m_id,
                         out_m_id=out_m_id, get_key=get_key, f=f, tab='\t\t')


def write_folded_efm(f_efm_id, S, S_folded, out_r_id, in_m_id, out_m_id, get_key, f, tab=''):
    if f_efm_id not in S_folded.gr_id2efm_ids:
        write_efm(efm_id=f_efm_id, S=S, out_r_id=out_r_id, in_m_id=in_m_id, out_m_id=out_m_id, get_key=get_key, f=f,
                  tab=tab)
        return

    f.write('%s%s of length %d:\n\n%s\t%s\n\n'
            % (tab, f_efm_id, S_folded.pws.get_len(f_efm_id), tab,
               r_id2c_to_string(S_folded.pws.get_r_id2coeff(f_efm_id), get_key=get_key)))

    efm_ids = S_folded.gr_id2efm_ids[f_efm_id]
    f.write('%s\tUnfolds into %s' % (tab, ('the following %d EFMs:\n\n' % len(efm_ids)) if len(efm_ids) != 1 else ''))

    started = False
    for efm_id in efm_ids:
        if started:
            f.write('%s\t\t%s' % (tab, TINY_DELIMITER))
        else:
            started = True
        write_efm(efm_id=efm_id, S=S, out_r_id=out_r_id, in_m_id=in_m_id, out_m_id=out_m_id, get_key=get_key, f=f,
                  tab='%s\t\t' % tab, no_first_tab=len(efm_ids) == 1)


def write_efm(efm_id, S, out_r_id, in_m_id, out_m_id, get_key, f, tab='', no_first_tab=False):
    f.write('%s%s of length %d, of efficiency %g, of yield %g:\n\n%s\t%s\n\n'
            % ('' if no_first_tab else tab, efm_id, S.pws.get_len(efm_id),
               round_value(S.pws.get_control_efficiency(efm_id, out_r_id)),
               round_value(S.get_yield(efm_id, in_m_id, out_m_id)), tab,
               r_id2c_to_string(S.pws.get_r_id2coeff(efm_id), get_key=get_key)))


def serialize(model, path, get_f_path, **kwargs):
    pathways_txt, (efm_num, f_efm_num, pw_num) = serialize_efms(model, path, **kwargs)
    limit, efm_data = serialize_n_most_efficient_efms(model, path, n=3, **kwargs)
    fm_block = describe('fm_block.html', element_num=limit, characteristics='most effective',
                        element_name='EFM',
                        fms=[(efm_id, efm_len, get_f_path(efm_txt)) for (efm_id, efm_len, efm_txt) in efm_data],
                        all=limit == efm_num)
    return describe('efms.html', efm_num=efm_num, fm_num=pw_num, description_filepath=get_f_path(pathways_txt),
                    selected_efm_block=fm_block)
