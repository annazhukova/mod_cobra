import libsbml

from mod_cobra.constraint_based_analysis.efm.control_effective_flux_calculator import get_fm_yield
from mod_sbml.utils.misc import invert_map
from mod_cobra.constraint_based_analysis.efm.EFM import EFM, TYPE_EFM, TYPE_FOLDED_EFM, TYPE_PATHWAY
from mod_sbml.sbml.submodel_manager import submodel
from mod_sbml.serialization.serialization_manager import get_sbml_r_formula
from numpy_efm_manager import get_boundary_metabolites, get_yield, get_len, get_r_id2coeff, get_control_efficiency

THICK_DELIMITER = '==============================================================\n\n'
THIN_DELIMITER = '------------------------------------------\n\n'
TINY_DELIMITER = '-------------\n\n'

SIMPLE_PATTERN_SORTER = lambda p_id: -p_id

__author__ = 'anna'


def coeff_to_binary(coeff):
    return 1 if coeff > 0 else -1


def r_id2coeff_to_string(r_id2coeff, binary=False, get_key=lambda r_id, c: r_id, get_group=lambda r_id, c: None,
                         hidden_r_ids=None):
    if binary:
        return '\t'.join(('%s%s' % ('-' if c < 0 else '', r_id) for (r_id, c) in r_id2coeff.iteritems()))
    return '\t'.join(('%g %s' % (c, r_id) for (r_id, c) in r_id2coeff.iteritems()))
    # result = []
    # clique_started = None
    # started = False
    # for r_id in sorted(r_id2coeff.iterkeys(),
    #                    key=lambda r_id: get_key(r_id, coeff_to_binary(r_id2coeff[r_id]))):
    #     cl_id = get_group(r_id, coeff_to_binary(r_id2coeff[r_id]))
    #     if cl_id != clique_started:
    #         if clique_started:
    #             result.append(')')
    #             clique_started = None
    #         if cl_id and (not hidden_r_ids or r_id not in hidden_r_ids):
    #             clique_started = cl_id
    #             result.append('%s(' % ('\t' if started else ''))
    #             started = False
    #
    #     if not hidden_r_ids or r_id not in hidden_r_ids:
    #         if binary:
    #             result.append('%s%s%s' % ('\t' if started else '', '-' if r_id2coeff[r_id] < 0 else '', r_id))
    #         else:
    #             result.append('%s%g %s' % ('\t' if started else '', r_id2coeff[r_id], r_id))
    #         started = True
    #
    # if clique_started:
    #     result.append(')')
    #
    # return ''.join(result)


def serialize_efms(model, m_id2i, N_m, V, V_c, V_m, r_id2i, c_r_id2i, efm_id2i, gr_efm_id2i, m_efm_id2i,
                   efm_id2gr_id, efm_id2merged_id, out_m_id, in_m_id, out_r_id, path):
    merged_id2efm_ids = invert_map(efm_id2merged_id)
    grouped_id2efm_ids = invert_map(efm_id2gr_id)

    get_key = lambda r_id, c: r_id
    get_group = lambda r_id, c: None

    with open(path, 'w+') as f:
        f.write('Found %d EFMs, folded into %d EFMs, grouped into %d pathways:\n\n'
                % (len(efm_id2i), len(gr_efm_id2i), len(m_efm_id2i)))
        f.write(THICK_DELIMITER)

        # def get_key(r_id, c):
        #     key = [0] if r_id in all_fm_intersection or r_id in all_fm_intersection_folded else []
        #     if (r_id, c) in r_id2new_r_id:
        #         key.append(r_id2new_r_id[(r_id, c)])
        #     if (r_id, c) in r_id2cl_id:
        #         key.append(r_id2cl_id[(r_id, c)])
        #     key.append(r_id)
        #     return tuple(key)
        #
        # get_group = lambda r_id, c: r_id2cl_id[(r_id, c)] if (r_id, c) in r_id2cl_id else None
        #
        # if all_fm_intersection and len(all_fm_intersection):
        #     f.write('All EFMs contain following %d reactions:\n\t%s,\n\n'
        #             % (len(all_fm_intersection),
        #                all_fm_intersection.to_string(binary=True, get_key=get_key, get_group=get_group)))
        #     if len(all_fm_intersection) > len(all_fm_intersection_folded):
        #         f.write('or (with reaction groups folded):\n\t%s.\n\n'
        #                 % all_fm_intersection_folded.to_string(binary=True, get_key=get_key))
        #     f.write(THICK_DELIMITER)

        for fm_id in m_efm_id2i.iterkeys():
            write_pathway(model=model, fm_id=fm_id, m_fm_id2i=m_efm_id2i, g_fm_id2i=gr_efm_id2i, fm_id2i=efm_id2i,
                          N_m=N_m, V_m=V_m, V_c=V_c, V=V,
                          m_id2i=m_id2i, c_r_id2i=c_r_id2i, r_id2i=r_id2i,
                          out_m_id=out_m_id, in_m_id=in_m_id, out_r_id=out_r_id,
                          merged_id2efm_ids=merged_id2efm_ids, gr_id2efm_ids=grouped_id2efm_ids,
                          get_key=get_key, get_group=get_group, f=f)
            f.write(THICK_DELIMITER)


def write_pathway(model, fm_id, m_fm_id2i, g_fm_id2i, fm_id2i, N_m, V_m, V_c, V,
                  m_id2i, c_r_id2i, r_id2i, out_m_id, in_m_id, out_r_id,
                  merged_id2efm_ids, gr_id2efm_ids, get_key, get_group, f):
    v = V_m[:, m_fm_id2i[fm_id]]
    boundary_m_ids = {m.getId() for m in model.getListOfSpecies() if m.getBoundaryCondition()}
    boundary_m_indices = [m_id2i[m_id] for m_id in boundary_m_ids]
    bm_id2i = dict(zip(boundary_m_ids, xrange(0, len(boundary_m_ids))))
    b_ms = get_boundary_metabolites(N_m[boundary_m_indices, :], v)
    # bm_id2i = m_id2i
    # b_ms = get_boundary_metabolites(N_m, v)
    r_id2st = {(m_id, -b_ms[bm_id2i[m_id]]) for m_id in boundary_m_ids if b_ms[bm_id2i[m_id]] < 0}
    p_id2st = {(m_id, b_ms[bm_id2i[m_id]]) for m_id in boundary_m_ids if b_ms[bm_id2i[m_id]] > 0}
    f.write('Inputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                        for (m_id, st) in r_id2st))
    f.write('Outputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                           for (m_id, st) in p_id2st))
    fm_yield = get_yield(N_m, v, m_id2i[out_m_id], m_id2i[in_m_id])
    if fm_yield is not None:
        f.write('Yield: %g;\n\n' % fm_yield)

    if fm_id not in merged_id2efm_ids:
        write_folded_efm(fm_id, g_fm_id2i, fm_id2i, V_c, V, c_r_id2i, r_id2i, gr_id2efm_ids, out_r_id, get_key, get_group, f)
        return

    f.write('Pathway %s of length %d:\n\n\t%s\n\n'
            % (fm_id, get_len(v), r_id2coeff_to_string(get_r_id2coeff(v, c_r_id2i),
                                                       get_key=get_key, get_group=get_group)))
    folded_efm_ids = merged_id2efm_ids[fm_id]
    f.write('\tContains %d elements:\n\n' % len(folded_efm_ids))
    started = False
    for f_efm_id in folded_efm_ids:
        if started:
            f.write('\t\t%s' % THIN_DELIMITER)
        else:
            started = True
        write_folded_efm(f_efm_id, g_fm_id2i, fm_id2i,
                         V_c, V, c_r_id2i, r_id2i, gr_id2efm_ids, out_r_id, get_key, get_group, f, tab='\t\t')


def write_folded_efm(fm_id, g_fm_id2i, fm_id2i, V_g, V,
                     c_r_id2i, r_id2i, gr_id2efm_ids, out_r_id, get_key, get_group, f, tab=''):
    if fm_id not in gr_id2efm_ids:
        write_efm(fm_id, fm_id2i, V, r_id2i, out_r_id, get_key, get_group, f, tab=tab)
        return

    v = V_g[:, g_fm_id2i[fm_id]]

    f.write('%sFolded EFM %s of length %d:\n\n%s\t%s\n\n'
            % (tab, fm_id, get_len(v), tab, r_id2coeff_to_string(get_r_id2coeff(v, c_r_id2i),
                                                                 get_key=get_key, get_group=get_group)))

    efm_ids = gr_id2efm_ids[fm_id]
    f.write('%s\tUnfolds into %s' % (tab, ('the following %d EFMs:\n\n' % len(efm_ids)) if len(efm_ids) != 1 else ''))

    started = False
    for efm_id in efm_ids:
        if started:
            f.write('%s\t\t%s' % (tab, TINY_DELIMITER))
        else:
            started = True
        write_efm(efm_id, fm_id2i, V, r_id2i, out_r_id, get_key, get_group, f, tab='%s\t\t' % tab,
                  no_first_tab=len(efm_ids) == 1)


def write_efm(efm_id, efm_id2i, V, r_id2i, out_r_id, get_key, get_group, f, tab='', no_first_tab=False):
    v = V[:, efm_id2i[efm_id]]
    efficiency = get_control_efficiency(v, r_id2i[out_r_id])
    f.write('%sEFM %s of length %d of efficiency %g:\n\n%s\t%s\n\n'
            % ('' if no_first_tab else tab, efm_id, get_len(v), efficiency, tab,
               r_id2coeff_to_string(get_r_id2coeff(v, r_id2i), get_key=get_key, get_group=get_group)))

