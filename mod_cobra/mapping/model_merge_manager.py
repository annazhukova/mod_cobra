from numpy import hstack, zeros, vstack
from scipy.linalg import block_diag

from mod_cobra.efm.System import FoldedSystem

__author__ = 'anna'


def join(model_id2m_ids_groups, model_id2S):
    common_m_ids = set()
    for model_id2m_ids in model_id2m_ids_groups:
        for model_id, m_ids in model_id2m_ids.iteritems():
            common_m_ids |= {(model_id, m_id) for m_id in m_ids}
    N = block_diag(*[S.N[sorted(i for (m_id, i) in S.m_id2i.iteritems() if (model_id, m_id) not in common_m_ids), :]
                     for (model_id, S) in model_id2S.iteritems()])
    m_id2i, r_id2i, efm_id2i = {}, {}, {}
    r_shift, m_shift, efm_shift = 0, 0, 0
    boundary_m_ids = []
    for model_id, S in model_id2S.iteritems():
        r_id2i.update({(model_id, r_id): i + r_shift for (r_id, i) in S.r_id2i.iteritems()})
        r_shift += len(S.r_id2i)
        efm_id2i.update({(model_id, efm_id): i + efm_shift for (efm_id, i) in S.efm_id2i.iteritems()})
        efm_shift += len(S.efm_id2i)
        m_ids = [(model_id, m_id) for (m_id, i) in sorted(S.m_id2i.iteritems(), key=lambda (_, i): i)
                 if (model_id, m_id) not in common_m_ids]
        m_id2i.update(dict(zip(m_ids, xrange(m_shift, len(m_ids) + m_shift))))
        boundary_m_ids += [(model_id, m_id) for m_id in S.boundary_m_ids if (model_id, m_id) not in common_m_ids]
        m_shift += len(m_ids)
    m_id2new_m_id = dict()
    for model_id2m_ids in model_id2m_ids_groups:
        model_id2m_row = {model_id: zeros(shape=(1, model_id2S[model_id].N.shape[1]))
                          for model_id in model_id2S.iterkeys()}
        new_m_id = tuple(((model_id, tuple(m_ids)) for (model_id, m_ids) in model_id2m_ids.iteritems()))
        m_id2i[new_m_id] = m_shift
        if next((True for (model_id, m_ids) in model_id2m_ids.iteritems()
                 if m_ids & set(model_id2S[model_id].boundary_m_ids)), False):
            boundary_m_ids.append(new_m_id)
        for model_id, m_ids in model_id2m_ids.iteritems():
            for m_id in m_ids:
                m_id2new_m_id[(model_id, m_id)] = new_m_id
                model_id2m_row[model_id] += model_id2S[model_id].N[model_id2S[model_id].m_id2i[m_id], :]
        row = hstack(tuple(model_id2m_row[model_id] for model_id in model_id2S.iterkeys()))
        N = vstack((N, row))
        m_shift += 1
    V = block_diag(*[S.V for (model_id, S) in model_id2S.iteritems()])
    set_join = lambda s1, s2: s1 | set(s2)
    return FoldedSystem(N=N, m_id2i=m_id2i, r_id2i=r_id2i, V=V, efm_id2i=efm_id2i, m_id2gr_id=m_id2new_m_id,
                        boundary_m_ids=
                        sorted(reduce(set_join, (st_m.boundary_m_ids for st_m in model_id2S.itervalues()), set())))


def merge(S, ignored_ids=None):
    return S.remove_reaction_duplicates(ignored_ids)
