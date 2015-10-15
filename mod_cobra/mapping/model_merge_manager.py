from numpy import hstack, zeros, vstack
from scipy.linalg import block_diag
from mod_cobra.efm.System import System

__author__ = 'anna'


def join(model_id2m_ids_groups, model_id2S):
    common_m_ids = set()
    for model_id2m_ids in model_id2m_ids_groups:
        for model_id, m_ids in model_id2m_ids.iteritems():
            common_m_ids |= {(model_id, m_id) for m_id in m_ids}

    N = block_diag(*[S.N[sorted(S.m_id2i[m_id] for m_id in S.m_ids if (model_id, m_id) not in common_m_ids), :]
                     [:, sorted(S.r_id2i[r_id] for r_id in S.r_ids)]
                     for (model_id, S) in model_id2S.iteritems()])
    m_id2i, r_id2i, efm_id2i = {}, {}, {}
    r_shift, m_shift, efm_shift = 0, 0, 0
    boundary_m_ids = []
    for model_id, S in model_id2S.iteritems():
        r_id2i.update({(model_id, r_id): i + r_shift
                       for (i, r_id) in enumerate(sorted(S.r_ids, key=lambda r_id: S.r_id2i[r_id]))})
        r_shift += len(S.r_ids)
        efm_id2i.update({(model_id, efm_id): i + efm_shift
                         for (i, efm_id) in enumerate(sorted(S.efm_ids, key=lambda efm_id: S.efm_id2i[efm_id]))})
        efm_shift += len(S.efm_ids)
        m_ids = {m_id for m_id in S.m_ids if (model_id, m_id) not in common_m_ids}
        m_id2i.update({(model_id, m_id): i + m_shift
                       for (i, m_id) in enumerate(sorted(m_ids, key=lambda m_id: S.m_id2i[m_id]))})
        boundary_m_ids += [(model_id, m_id) for m_id in set(S.boundary_m_ids) & m_ids]
        m_shift += len(m_ids)
    m_id2new_m_id = dict()
    for model_id2m_ids in model_id2m_ids_groups:
        model_id2m_row = {model_id: zeros(shape=(1, len(S.r_ids))) for model_id, S in model_id2S.iteritems()}
        new_m_id = tuple(((model_id, tuple(m_ids)) for (model_id, m_ids) in model_id2m_ids.iteritems()))
        m_id2i[new_m_id] = m_shift
        if next((True for (model_id, m_ids) in model_id2m_ids.iteritems()
                 if m_ids & set(model_id2S[model_id].boundary_m_ids)), False):
            boundary_m_ids.append(new_m_id)
        for model_id, m_ids in model_id2m_ids.iteritems():
            S = model_id2S[model_id]
            r_is = sorted(S.r_id2i[r_id] for r_id in S.r_ids)
            for m_id in m_ids:
                m_id2new_m_id[(model_id, m_id)] = new_m_id
                model_id2m_row[model_id] += S.N[S.m_id2i[m_id], r_is]
        row = hstack(tuple(model_id2m_row[model_id] for model_id in model_id2S.iterkeys()))
        N = vstack((N, row))
        m_shift += 1
    V = block_diag(*[S.V[sorted(S.r_id2i[r_id] for r_id in S.r_ids), :]
                     [:, sorted(S.efm_id2i[efm_id] for efm_id in S.efm_ids)]
                     for (model_id, S) in model_id2S.iteritems()])
    return System(N=N, m_id2i=m_id2i, r_id2i=r_id2i, V=V, efm_id2i=efm_id2i, m_id2gr_id=m_id2new_m_id,
                  boundary_m_ids=boundary_m_ids)


def merge(S, ignored_ids=None):
    return S.remove_reaction_duplicates(ignored_ids)
