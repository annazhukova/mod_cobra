from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
import libsbml
from constraint_based_analysis import ZERO_THRESHOLD, round_value

from sbml_manager import get_reactants, get_products, create_compartment, create_species, create_reaction

__author__ = 'anna'


def get_element2id_mapping(model):
    s_id2i = dict(zip((s.getId() for s in sorted(model.getListOfSpecies(),
                                                 key=lambda s: 1 if s.getBoundaryCondition() else 0)),
                      xrange(0, model.getNumSpecies())))
    r_id2i = dict(zip((r.getId() for r in model.getListOfReactions()), xrange(0, model.getNumReactions())))
    return s_id2i, r_id2i


def get_efm_matrix(r_id2coefficient_list, r_id2i):
    rows, cols, data = [], [], []
    for efm_id, r_id2coef in zip(xrange(0, len(r_id2coefficient_list)), r_id2coefficient_list):
        cols.extend([efm_id] * len(r_id2coef))
        rows.extend((r_id2i[r_id] for r_id in r_id2coef.iterkeys()))
        data.extend(r_id2coef.itervalues())
    V = csr_matrix((np.array(data), (np.array(rows), np.array(cols))),
                   shape=(len(r_id2i), len(r_id2coefficient_list))).toarray()
    replace_zeros(V)
    return V


def remove_invalid_efms(N, V, threshold=ZERO_THRESHOLD):
    bm = np.dot(N, V)
    return V[:, [i for i in xrange(0, bm.shape[1]) if not next((it for it in bm[:, i] if abs(it) > threshold), False)]]
    # return np.array([v for v in V.T if np.count_nonzero(np.dot(N, v.T)) == 0]).T


def get_control_efficiency(v, r_index):
    return 1.0 * v[r_index] / sum((abs(c) for c in v))


def get_len(v):
    return np.count_nonzero(v)


def get_yield(N, v, out_m_index, in_m_index):
    in_value = np.dot(N[in_m_index, :], v)
    if not in_value:
        return float('inf')
    return np.dot(N[out_m_index, :], v) / (-1.0 * in_value)


def model2stoichiometric_matrix(model, s_id2i, r_id2i):
    rows, cols, data = [], [], []
    for r in model.getListOfReactions():
        cols.extend([r_id2i[r.getId()]] * (r.getNumReactants() + r.getNumProducts()))
        for s_id, st in get_reactants(r, stoichiometry=True):
            rows.append(s_id2i[s_id])
            data.append(-st)
        for s_id, st in get_products(r, stoichiometry=True):
            rows.append(s_id2i[s_id])
            data.append(st)
    return csr_matrix((np.array(data), (np.array(rows), np.array(cols))), shape=(len(s_id2i), len(r_id2i))).toarray()


def get_coupled_reactions(V, r_id2i):
    return get_equivalent_rows(V, r_id2i, axis=0)


def replace_zeros(array, threshold=ZERO_THRESHOLD):
    np.place(array, abs(array) < threshold, 0)


def lump_coupled_reactions(N, V, coupled_r_id_groups, r_id2i):
    coupled_r_ids = reduce(lambda s1, s2: s1 | set(s2), coupled_r_id_groups, set())
    ordered_r_ids = [r_id for r_id in sorted(r_id2i.iterkeys(), key=lambda r_id: r_id2i[r_id])
                     if r_id not in coupled_r_ids]
    new_r_id2i = dict(zip(ordered_r_ids, xrange(0, len(ordered_r_ids))))
    indices = tuple((r_id2i[r_id] for r_id in ordered_r_ids))
    N_new = N[:, indices]
    V_new = V[indices, :]
    r_id2lr_id = {}
    lr_id2r_id2c = {}
    for r_ids, i in zip(coupled_r_id_groups, xrange(len(ordered_r_ids), len(ordered_r_ids) + len(coupled_r_id_groups))):
        lumped_r_id = 'lr_%d' % i
        r_id2lr_id.update({r_id: lumped_r_id for r_id in r_ids})
        new_r_id2i[lumped_r_id] = i

        sample_r_index = r_id2i[r_ids[0]]
        sample_efm_index = np.nonzero(V[sample_r_index, :])[0][0]
        sample_efm = V[:, sample_efm_index]
        c = abs(sample_efm[sample_r_index]) * 1.0
        sign = 1 if sample_efm[sample_r_index] > 0 else -1
        lr_id2r_id2c[lumped_r_id] = {r_id: sample_efm[r_id2i[r_id]] / c for r_id in r_ids}
        lambdas = np.array([[lr_id2r_id2c[lumped_r_id][r_id]] for r_id in r_ids])
        r_new = np.dot(N[:, tuple((r_id2i[r_id] for r_id in r_ids))], lambdas)
        replace_zeros(r_new)
        N_new = np.concatenate((N_new, r_new), axis=1)
        V_new = np.concatenate((V_new, [V[sample_r_index, :] * sign]), axis=0)
    return N_new, V_new, new_r_id2i, r_id2lr_id, lr_id2r_id2c


def get_reaction_duplicates(N, r_id2i, m_indices=None):
    return get_equivalent_rows(N[m_indices, :] if m_indices else N, r_id2i, axis=1)


def remove_reaction_duplicates(N, V, r_id_groups, r_id2i):
    grouped_r_ids = reduce(lambda s1, s2: s1 | set(s2), r_id_groups, set())
    ordered_r_ids = [r_id for r_id in r_id2i.iterkeys() if r_id not in grouped_r_ids]
    r_id2gr_id = {}
    new_r_id2i = dict(zip(ordered_r_ids, xrange(0, len(ordered_r_ids))))
    indices = tuple((r_id2i[r_id] for r_id in ordered_r_ids))
    N_new = N[:, indices]
    V_new = V[indices, :]
    gr_id2r_id2c = {}
    for r_ids, i in zip(r_id_groups, xrange(len(ordered_r_ids), len(ordered_r_ids) + len(r_id_groups))):
        grouped_r_id = 'gr_%d' % i
        r_id2gr_id.update({r_id: grouped_r_id for r_id in r_ids})
        new_r_id2i[grouped_r_id] = i

        sample_r_index = r_id2i[r_ids[0]]

        sample_m_index = np.nonzero(N[:, sample_r_index])[0][0]
        sample_m = N[sample_m_index, :]
        st = sample_m[sample_r_index] * 1.0
        indices = tuple((r_id2i[r_id] for r_id in r_ids))
        gr_id2r_id2c[grouped_r_id] = {r_id: sample_m[r_id2i[r_id]] / st for r_id in r_ids}

        N_new = np.concatenate((N_new, N[:, (sample_r_index, )]), axis=1)
        v_new = np.dot(np.array([gr_id2r_id2c[grouped_r_id][r_id] for r_id in r_ids]), V[indices, :])
        replace_zeros(v_new)
        V_new = np.vstack((V_new, v_new))
    return N_new, V_new, new_r_id2i, r_id2gr_id, gr_id2r_id2c


def get_efm_duplicates(V, efm_id2i):
    return get_equivalent_rows(V, efm_id2i, axis=1, consider_sign=True)


def remove_efm_duplicates(V, efm_id_groups, efm_id2i):
    grouped_efm_ids = reduce(lambda s1, s2: s1 | set(s2), efm_id_groups, set())
    ordered_efm_ids = [efm_id for efm_id in efm_id2i.iterkeys() if efm_id not in grouped_efm_ids]
    efm_id2gr_id = {}
    new_efm_id2i = dict(zip(ordered_efm_ids, xrange(0, len(ordered_efm_ids))))
    indices = tuple((efm_id2i[efm_id] for efm_id in ordered_efm_ids))
    V_new = V[:, indices]
    grouped_efm_id = max(efm_id2i.iterkeys())
    for efm_ids, i in zip(efm_id_groups, xrange(len(ordered_efm_ids), len(ordered_efm_ids) + len(efm_id_groups))):
        grouped_efm_id += 1
        efm_id2gr_id.update({efm_id: grouped_efm_id for efm_id in efm_ids})
        new_efm_id2i[grouped_efm_id] = i

        sample_efm_index = efm_id2i[efm_ids[0]]
        V_new = np.concatenate((V_new, V[:, (sample_efm_index, )]), axis=1)
    return V_new, new_efm_id2i, efm_id2gr_id


def get_boundary_metabolites(N, v):
    ms = np.dot(N, v)
    replace_zeros(ms)
    divider = 1.0 * min((abs(it) for it in ms if it))
    ms /= divider
    return np.array([round_value(it) for it in ms])


def get_efm_intersection(V, r_id2i):
    result = {}
    for r_id, i in r_id2i.iteritems():
        r = V[i, :]
        replace_zeros(r)
        if r[r > 0].shape[0] == r.shape[0]:
            result[r_id] = 1
        elif r[r < 0].shape[0] == r.shape[0]:
            result[r_id] = -1
    return result


def get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i):
    key2efm_ids = defaultdict(list)
    for efm_id, i in efm_id2i.iteritems():
        key2efm_ids[tuple(get_boundary_metabolites(N, V[:, i]))].append(efm_id)
    return [e_ids for e_ids in key2efm_ids.itervalues() if len(e_ids) > 1]


def merge_efm_groups(V, efm_id_groups, efm_id2i):
    grouped_efm_ids = reduce(lambda s1, s2: s1 | set(s2), efm_id_groups, set())
    ordered_efm_ids = [efm_id for efm_id in efm_id2i.iterkeys() if efm_id not in grouped_efm_ids]
    efm_id2gr_id = {}
    new_efm_id2i = dict(zip(ordered_efm_ids, xrange(0, len(ordered_efm_ids))))
    indices = tuple((efm_id2i[efm_id] for efm_id in ordered_efm_ids))
    V_new = V[:, indices]
    grouped_efm_id = max(efm_id2i.iterkeys())
    for efm_ids, i in zip(efm_id_groups, xrange(len(ordered_efm_ids), len(ordered_efm_ids) + len(efm_id_groups))):
        grouped_efm_id += 1
        efm_id2gr_id.update({efm_id: grouped_efm_id for efm_id in efm_ids})
        new_efm_id2i[grouped_efm_id] = i

        indices = tuple((efm_id2i[efm_id] for efm_id in efm_ids))
        lambdas = np.array([[1] for _ in efm_ids])
        v_new = np.dot(V[:, indices], lambdas)
        replace_zeros(v_new)
        V_new = np.concatenate((V_new, v_new), axis=1)
    return V_new, new_efm_id2i, efm_id2gr_id


# -----------------------Utility-method--------------------------#
def get_equivalent_rows(M, e_id2i, axis=0, consider_sign=False, threshold=ZERO_THRESHOLD):
    """
    :param M: matrix
    :param e_id2i: dict {element_id: its_index}
    :param axis: 0 to look for equivalent rows, 1 for columns
    :return: list of groups of ids of equivalent elements
    """
    key2e_ids = defaultdict(list)
    for e_id, i in e_id2i.iteritems():
        e_i = M[i, :] if 0 == axis else M[:, i]
        replace_zeros(e_i, threshold)
        if np.count_nonzero(e_i) <= 1:
            continue
        coeff = e_i[np.nonzero(e_i)[0][0]] * 1.0
        key2e_ids[(tuple(e_i / coeff), (coeff > 0) if consider_sign else True)].append(e_id)
    return [e_ids for e_ids in key2e_ids.itervalues() if len(e_ids) > 1]

