import os
from mod_cobra import round_value

from mod_cobra.html import describe
from mod_cobra.efm.serialization import r_id2c_to_string, write_inputs_outputs, THIN_DELIMITER, \
    THICK_DELIMITER, TINY_DELIMITER, write_detailed_r_id2c
from mod_cobra.efm.serialization import pathway_matrix_serializer

__author__ = 'anna'


def serialize_n_shortest_efms(model, path, S, n=3):
    initial_efm_ids = {efm_id for efm_id in S.efm_id2i.iterkeys() if efm_id not in S.gr_id2efm_ids}
    initial_r_ids = {r_id for r_id in S.r_id2i.iterkeys() if r_id not in S.gr_id2r_id2c}
    limit = min(n if n is not None and n >= 0 else len(S.efm_id2i), len(S.efm_id2i))
    efm_data = [(efm_id, S.get_len(efm_id, r_ids=initial_r_ids),
                 serialize_efm(S, efm_id, model, path))
                for efm_id in sorted(initial_efm_ids,
                                     key=lambda efm_id: (S.get_len(efm_id, r_ids=initial_r_ids), efm_id))[0: limit]]
    return limit, efm_data


def serialize_efm(S, efm_id, model, path):
    efm_txt = os.path.join(path, '%s.txt' % efm_id)
    initial_r_ids = {r_id for r_id in S.r_id2i.iterkeys() if r_id not in S.gr_id2r_id2c}
    r_id2c = S.pws.get_r_id2coeff(efm_id, r_ids=initial_r_ids)
    with open(efm_txt, 'w+') as f:
        f.write('EFM %s of length %d\n\n' % (efm_id, len(r_id2c)))
        write_inputs_outputs(f, model, S.get_boundary_inputs_outputs(efm_id))
        f.write(THIN_DELIMITER)
        write_detailed_r_id2c(model, r_id2c, f)
    return efm_txt


def serialize_efms(model, path, S):
    initial_efm_ids = {efm_id for efm_id in S.efm_id2i.iterkeys() if efm_id not in S.gr_id2efm_ids}
    initial_r_ids = {r_id for r_id in S.r_id2i.iterkeys() if r_id not in S.gr_id2r_id2c}
    all_fm_intersection = S.get_efm_intersection(initial_efm_ids, r_ids=initial_r_ids)
    all_fm_intersection_folded = S.get_efm_intersection(initial_efm_ids)

    def get_key(r_id):
        mr_id, mc = None, 0
        gr_id, gc = None, 0
        if r_id in S.r_id2gr_id:
            gr_id = S.r_id2gr_id[r_id]
            gc = S.gr_id2r_id2c[gr_id][r_id]
            if gr_id in S.r_id2gr_id:
                mr_id = S.r_id2gr_id[gr_id]
                mc = S.gr_id2r_id2c[mr_id][gr_id]
        in_intersection = 1 if r_id in all_fm_intersection or r_id in all_fm_intersection_folded \
                               or gr_id and gr_id in all_fm_intersection_folded \
                               or mr_id and mr_id in all_fm_intersection_folded else 0
        return in_intersection, (mr_id, mc), (gr_id, gc), r_id

    pathways_txt = os.path.join(path, 'pathways.txt')
    with open(pathways_txt, 'w+') as f:
        f.write('Found %d EFMs, grouped into %d behaviours:\n\n' % (len(initial_efm_ids), len(S.efm_ids)))
        f.write(THICK_DELIMITER)

        if all_fm_intersection and len(all_fm_intersection):
            f.write('All EFMs contain following %d reactions:\n\t%s.\n\n'
                    % (len(all_fm_intersection),
                       r_id2c_to_string(all_fm_intersection, binary=True, get_key=get_key)))
            f.write(THICK_DELIMITER)

        for fm_id in sorted(S.efm_ids):
            write_pathway(model=model, pathway_id=fm_id, S=S, get_key=get_key, f=f)
            f.write(THICK_DELIMITER)
    return pathways_txt, (len(initial_efm_ids), len(S.efm_ids))


def write_pathway(model, pathway_id, S, get_key, f):
    write_inputs_outputs(f, model, S.get_boundary_inputs_outputs(pathway_id))

    if pathway_id not in S.pathways:
        write_folded_efm(f_efm_id=pathway_id, S=S, get_key=get_key, f=f)
        return

    f.write('%s of length %d:\n\n\t%s\n\n'
            % (pathway_id, S.get_len(pathway_id),
               r_id2c_to_string(S.pws.get_r_id2coeff(pathway_id, r_ids=S.r_ids), get_key=get_key)))
    folded_efm_ids = S.gr_id2efm_ids[pathway_id]
    f.write('\tContains %d elements:\n\n' % len(folded_efm_ids))
    started = False
    for f_efm_id in folded_efm_ids:
        if started:
            f.write('\t\t%s' % THIN_DELIMITER)
        else:
            started = True
        write_folded_efm(f_efm_id=f_efm_id, S=S, get_key=get_key, f=f, tab='\t\t')


def write_folded_efm(f_efm_id, S, get_key, f, tab=''):
    if f_efm_id not in S.folded_efms:
        write_efm(efm_id=f_efm_id, S=S, get_key=get_key, f=f, tab=tab)
        return

    f.write('%s%s of length %d:\n\n%s\t%s\n\n'
            % (tab, f_efm_id, S.get_len(f_efm_id), tab,
               r_id2c_to_string(S.pws.get_r_id2coeff(f_efm_id, r_ids=S.r_ids), get_key=get_key)))

    efm_ids = S.gr_id2efm_ids[f_efm_id]
    f.write('%s\tContains %s' % (tab, ('the following %d EFMs:\n\n' % len(efm_ids)) if len(efm_ids) != 1 else ''))

    started = False
    for efm_id in efm_ids:
        if started:
            f.write('%s\t\t%s' % (tab, TINY_DELIMITER))
        else:
            started = True
        write_efm(efm_id=efm_id, S=S, get_key=get_key, f=f, tab='%s\t\t' % tab, no_first_tab=len(efm_ids) == 1)


def write_efm(efm_id, S, get_key, f, tab='', no_first_tab=False):
    initial_r_ids = {r_id for r_id in S.r_id2i.iterkeys() if r_id not in S.gr_id2r_id2c}
    f.write('%s%s of length %d:\n\n%s\t%s\n\n'
            % ('' if no_first_tab else tab, efm_id, S.get_len(efm_id, r_ids=initial_r_ids),
               tab, r_id2c_to_string(S.pws.get_r_id2coeff(efm_id, r_ids=initial_r_ids), get_key=get_key)))


def serialize(model, path, get_f_path, **kwargs):
    model_name, S = kwargs['model_name'], kwargs['S']
    pathways_txt, (efm_num, pw_num) = serialize_efms(model, path, S)
    limit, efm_data = serialize_n_shortest_efms(model, path, S, n=3)
    fm_block = describe('fm_block.html', element_num=limit, characteristics='shortest',
                        element_name='EFM',
                        fms=[(efm_id, efm_len, get_f_path(efm_txt)) for (efm_id, efm_len, efm_txt) in efm_data],
                        all=limit == efm_num)
    # pw_matrix_html = pathway_matrix_serializer.serialize(S, model, in_m_id, out_m_id, model_name, main_dir)
    return describe('efms.html', efm_num=efm_num, fm_num=pw_num, description_filepath=get_f_path(pathways_txt),
                    selected_efm_block=fm_block, pw_matrix_html=None)
