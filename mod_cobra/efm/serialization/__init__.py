from mod_sbml.utils.misc import invert_map
from mod_sbml.serialization import get_sbml_r_formula

__author__ = 'anna'

THICK_DELIMITER = '==============================================================\n\n'
THIN_DELIMITER = '------------------------------------------\n\n'
TINY_DELIMITER = '-------------\n\n'


def coefficient_to_string(c, binary=False):
    if binary or abs(c) == 1:
        return '' if c > 0 else '-'
    return '%g ' % c


def r_id2c_to_string(r_id2c, binary=False, get_key=lambda r_id: (0, ('', 0), ('', 0), r_id), hidden_r_ids=None):
    result = []
    clique_started = None
    started = False
    for r_id in sorted(r_id2c.keys(), key=lambda r_id: get_key(r_id), reverse=True):
        coefficient = r_id2c[r_id]
        _, (group, gc), (clique, cc), _ = get_key(r_id)
        if clique != clique_started:
            if clique_started:
                result.append(')')
                clique_started = None
                started = True

            if clique and (not hidden_r_ids or r_id not in hidden_r_ids):
                clique_started = clique
                gr_c = coefficient_to_string(((coefficient / cc) / gc) if gc else 0, binary)
                cl_c = coefficient_to_string(coefficient / cc if not gc else cc, binary)
                result.append('%s%s: ('
                              % ('\t' if started else '',
                                 ('%s%s(%s%s)' % (gr_c, group, cl_c, clique)) if group else '%s%s' % (cl_c, clique)))
                started = False
        if not hidden_r_ids or r_id not in hidden_r_ids:
            coefficient = coefficient_to_string(cc if cc else coefficient, binary)
            result.append('%s%s%s' % ('\t' if started else '', coefficient, r_id))
            started = True

    if clique_started:
        result.append(')')

    return ''.join(result)


def write_detailed_r_id2c(model, r_id2c, f):
    c2r_ids = invert_map(r_id2c)
    for c, r_ids in sorted(c2r_ids.items(), key=lambda it: (-abs(it[0]), -it[0])):
        for r_id in sorted(r_ids):
            f.write('%g\t%s:\t%s\n'
                    % (c, r_id, get_sbml_r_formula(model, model.getReaction(r_id), show_compartments=False,
                                                   show_metabolite_ids=True)))
        f.write('\n')


def write_metabolites(m_id2st, model, f, prefix=''):
    f.write('%s%s;\n\n' % (prefix, ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                             for (m_id, st) in sorted(m_id2st.items(), key=lambda it: it[0]))))


def write_inputs_outputs(f, model, inputs, in_m_id=None):
    (r_id2st, p_id2st) = inputs
    if in_m_id and in_m_id in r_id2st:
        c = r_id2st[in_m_id]
        if c:
            r_id2st = {m_id: 1.0 * st / c for (m_id, st) in r_id2st.items()}
            p_id2st = {m_id: 1.0 * st / c for (m_id, st) in p_id2st.items()}
    write_metabolites(r_id2st, model, f, prefix='Inputs: ')
    write_metabolites(p_id2st, model, f, prefix='Outputs: ')
