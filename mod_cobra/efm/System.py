from collections import Counter
import numpy as np

from mod_cobra.efm import coefficient_to_binary
from mod_sbml.utils.misc import invert_map
from mod_cobra.efm.numpy_efm_manager import get_boundary_metabolites, get_coupled_reactions, get_reaction_duplicates, \
    get_efm_duplicates, get_efm_groups_based_on_boundary_metabolites, get_len, get_yield, get_control_efficiency, \
    get_efm_intersection, get_support, get_2_efm_intersection, get_unique_id, replace_zeros

__author__ = 'anna'


class StoichiometricMatrix(object):
    def __init__(self, N, m_id2i, r_id2i, boundary_m_ids):
        self.N = N
        self.m_id2i = m_id2i
        self.r_id2i = r_id2i
        self.boundary_m_ids = boundary_m_ids

    def get_inputs_outputs(self, r_id):
        ms = self.N[:, self.r_id2i[r_id]]
        return {m_id: -ms[i] for (m_id, i) in self.m_id2i.iteritems() if ms[i] < 0},\
               {m_id: ms[i] for (m_id, i) in self.m_id2i.iteritems() if ms[i] > 0}


class PathwaySet(object):
    def __init__(self, V, r_id2i, efm_id2i):
        self.V = V
        self.r_id2i = r_id2i
        self.efm_id2i = efm_id2i
        self.support = None
        self.i2r_id = None

    def get_r_id2coeff(self, efm_id, binary=False, r_ids=None):
        if not r_ids:
            r_ids = self.r_id2i.iterkeys()
        if binary:
            return {r_id: coefficient_to_binary(self.V[self.r_id2i[r_id], self.efm_id2i[efm_id]]) for r_id in r_ids
                    if self.V[self.r_id2i[r_id], self.efm_id2i[efm_id]]}
        return {r_id: self.V[self.r_id2i[r_id], self.efm_id2i[efm_id]] for r_id in r_ids
                if self.V[self.r_id2i[r_id], self.efm_id2i[efm_id]]}

    def get_efm_ids(self, r_id, reversed=False, efm_ids=None):
        if not efm_ids:
            efm_ids = self.efm_id2i.iterkeys()
        check = lambda it: (it < 0) if reversed else (it > 0)
        return {efm_id for efm_id in efm_ids if check(self.V[self.r_id2i[r_id], self.efm_id2i[efm_id]])}

    def get_control_efficiency(self, efm_id, r_id):
        v = self.V[:, self.efm_id2i[efm_id]]
        return get_control_efficiency(v, self.r_id2i[r_id])

    def get_support_V(self):
        if self.support is not None:
            return self.support
        return get_support(self.V)


class System(object):
    def __init__(self, st_matrix=None, pws=None,
                 N=None, V=None, m_id2i=None, r_id2i=None, efm_id2i=None, boundary_m_ids=None,
                 r_ids=None, m_ids=None, efm_ids=None,
                 r_id2gr_id=None, gr_id2r_id2c=None, efm_id2gr_id=None, m_id2gr_id=None,):
        if not r_id2i:
            r_id2i = st_matrix.r_id2i if st_matrix else (pws.r_id2i if pws else {})
        self.st_matrix = st_matrix
        if N is not None and m_id2i is not None and r_id2i is not None:
            self.st_matrix = StoichiometricMatrix(N, m_id2i, r_id2i, boundary_m_ids)
        self.pws = pws
        if V is not None and r_id2i is not None:
            self.pws = PathwaySet(V, r_id2i, efm_id2i)

        self.r_ids = set(r_ids) if r_ids else set(self.r_id2i.iterkeys())
        self.m_ids = set(m_ids) if m_ids else set(self.m_id2i.iterkeys())
        self.efm_ids = set(efm_ids) if efm_ids else set(self.efm_id2i.iterkeys())

        self.r_id2gr_id = r_id2gr_id if r_id2gr_id else {}
        self.gr_id2r_id2c = gr_id2r_id2c if gr_id2r_id2c else {}
        self.efm_id2gr_id = efm_id2gr_id if efm_id2gr_id else {}
        self.gr_id2efm_ids = invert_map(self.efm_id2gr_id, list)
        self.m_id2gr_id = m_id2gr_id if m_id2gr_id else {}

        self.coupled_rs = set()
        self.r_types = set()
        self.folded_efms = set()
        self.pathways = set()

    @property
    def V(self):
        return self.pws.V

    @V.setter
    def V(self, V):
        self.pws.V = V

    @property
    def efm_id2i(self):
        return self.pws.efm_id2i

    @property
    def r_id2i(self):
        return self.st_matrix.r_id2i

    @property
    def N(self):
        return self.st_matrix.N

    @N.setter
    def N(self, N):
        self.st_matrix.N = N

    @property
    def m_id2i(self):
        return self.st_matrix.m_id2i

    @property
    def boundary_m_ids(self):
        return self.st_matrix.boundary_m_ids

    def get_efm_intersection(self, efm_ids=None, r_ids=None):
        V = self.pws.get_support_V()
        if efm_ids and len(efm_ids) == 1:
            return self.pws.get_r_id2coeff(next(iter(efm_ids)), r_ids)
        r_ids = sorted(r_ids if r_ids else self.r_ids)
        i2r_id = dict(enumerate(r_ids))
        V = V[[self.r_id2i[r_id] for r_id in r_ids], :]
        if efm_ids and len(efm_ids) == 2:
            return get_2_efm_intersection(V, self.efm_id2i[efm_ids[0]], self.efm_id2i[efm_ids[1]], i2r_id)
        if efm_ids:
            V = V[:, [self.efm_id2i[efm_id] for efm_id in efm_ids]]
        return get_efm_intersection(V, i2r_id)

    def get_efm_ids_by_r_id(self, r_id, efm_ids=None):
        if not efm_ids:
            efm_ids = self.efm_id2i.iterkeys()
        return [efm_id for efm_id in efm_ids if self.V[self.r_id2i[r_id], self.efm_id2i[efm_id]]]

    def get_len(self, efm_id, r_ids=None):
        return get_len(self.V[[self.r_id2i[r_id] for r_id in (r_ids if r_ids else self.r_ids)], self.efm_id2i[efm_id]])

    def get_yield(self, efm_id, in_m_id, out_m_id):
        r_is = [self.r_id2i[r_id] for r_id in self.r_ids]
        return get_yield(self.N[:, r_is], self.V[r_is, self.efm_id2i[efm_id]], self.m_id2i[out_m_id],
                         self.m_id2i[in_m_id])

    def get_boundary_inputs_outputs(self, pathway_id):
        r_is = [self.r_id2i[r_id] for r_id in self.r_ids]
        v = self.V[r_is, self.efm_id2i[pathway_id]]
        bm_ids = sorted(set(self.boundary_m_ids) & self.m_ids)
        b_ms = get_boundary_metabolites(self.N[[self.m_id2i[m_id] for m_id in bm_ids], :][:, r_is], v)
        bm_id2i = dict(zip(bm_ids, xrange(0, len(bm_ids))))
        # b_ms = get_boundary_metabolites(self.N, v)
        r_id2st = {m_id: -b_ms[bm_id2i[m_id]] for m_id in bm_ids if b_ms[bm_id2i[m_id]] < 0}
        p_id2st = {m_id: b_ms[bm_id2i[m_id]] for m_id in bm_ids if b_ms[bm_id2i[m_id]] > 0}
        return r_id2st, p_id2st

    def lump_coupled_reactions(self):
        coupled_r_id_groups = get_coupled_reactions(self.V, self.r_id2i)
        lr_i = 0
        for r_ids in coupled_r_id_groups:
            lumped_r_id, lr_i = get_unique_id(self.r_id2i, 'coupled_r', lr_i)
            self.coupled_rs.add(lumped_r_id)
            self.r_id2gr_id.update({r_id: lumped_r_id for r_id in r_ids})
            self.r_id2i[lumped_r_id] = self.V.shape[0]

            sample_r_index = self.r_id2i[r_ids[0]]
            sample_efm_index = np.nonzero(self.V[sample_r_index, :])[0][0]
            sample_efm = self.V[:, sample_efm_index]
            c = abs(sample_efm[sample_r_index]) * 1.0
            sign = coefficient_to_binary(sample_efm[sample_r_index])
            self.gr_id2r_id2c[lumped_r_id] = {r_id: sample_efm[self.r_id2i[r_id]] / c for r_id in r_ids}
            lambdas = np.array([[self.gr_id2r_id2c[lumped_r_id][r_id]] for r_id in r_ids])
            r_new = np.dot(self.N[:, tuple((self.r_id2i[r_id] for r_id in r_ids))], lambdas)
            replace_zeros(r_new)
            self.N = np.concatenate((self.N, r_new), axis=1)
            self.V = np.concatenate((self.V, [self.V[sample_r_index, :] * sign]), axis=0)

            self.r_ids -= set(r_ids)
            self.r_ids.add(lumped_r_id)

    def merge_efm_groups(self):
        r_id2i = self.get_main_r_id2i()
        m_id2i = self.get_main_m_id2i(boundary=True)
        efm_id2i = self.get_main_efm_id2i()
        V = self.get_main_V(r_id2i, efm_id2i)
        N = self.get_main_N(r_id2i, m_id2i)
        efm_id_groups = \
            get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i)

        gr_i = 0
        for efm_ids in efm_id_groups:
            grouped_efm_id, gr_i = get_unique_id(self.efm_id2i, 'pathway', gr_i)
            self.pathways.add(grouped_efm_id)
            self.efm_id2gr_id.update({efm_id: grouped_efm_id for efm_id in efm_ids})
            self.gr_id2efm_ids[grouped_efm_id] = efm_ids
            self.efm_id2i[grouped_efm_id] = self.V.shape[1]

            indices = tuple((self.efm_id2i[efm_id] for efm_id in efm_ids))
            lambdas = np.array([[1] for _ in efm_ids])
            v_new = np.dot(self.V[:, indices], lambdas)
            replace_zeros(v_new)
            self.V = np.concatenate((self.V, v_new), axis=1)

            self.efm_ids -= set(efm_ids)
            self.efm_ids.add(grouped_efm_id)

    def remove_efm_duplicates(self):
        efm_id2i = self.get_main_efm_id2i()
        V = self.get_main_V(efm_id2i=efm_id2i)
        efm_id_groups = get_efm_duplicates(V, efm_id2i)

        gr_i = 0
        for efm_ids in efm_id_groups:
            grouped_efm_id, gr_i = get_unique_id(self.efm_id2i, 'efm_type', gr_i)
            self.folded_efms.add(grouped_efm_id)
            self.efm_id2gr_id.update({efm_id: grouped_efm_id for efm_id in efm_ids})
            self.gr_id2efm_ids[grouped_efm_id] = efm_ids
            self.efm_id2i[grouped_efm_id] = self.V.shape[1]

            sample_efm_index = self.efm_id2i[efm_ids[0]]
            self.V = np.concatenate((self.V, self.V[:, (sample_efm_index, )]), axis=1)

            self.efm_ids -= set(efm_ids)
            self.efm_ids.add(grouped_efm_id)

    def remove_reaction_duplicates(self, ignored_m_ids=None):
        r_id2i = self.get_main_r_id2i()
        m_id2i = self.get_main_m_id2i(ignored_m_ids=ignored_m_ids)
        N = self.get_main_N(r_id2i=r_id2i, m_id2i=m_id2i)
        r_id_groups = get_reaction_duplicates(N, r_id2i)

        gr_i = 0
        for r_ids in r_id_groups:
            grouped_r_id, gr_i = get_unique_id(self.r_id2i, 'r_type', gr_i)
            self.r_types.add(grouped_r_id)
            self.r_id2gr_id.update({r_id: grouped_r_id for r_id in r_ids})
            self.r_id2i[grouped_r_id] = self.V.shape[0]

            sample_r_index = self.r_id2i[r_ids[0]]
            sample_m_index = np.nonzero(self.N[:, sample_r_index])[0][0]
            sample_m = self.N[sample_m_index, :]
            st = sample_m[sample_r_index] * 1.0
            indices = tuple((self.r_id2i[r_id] for r_id in r_ids))
            self.gr_id2r_id2c[grouped_r_id] = {r_id: sample_m[self.r_id2i[r_id]] / st for r_id in r_ids}

            self.N = np.concatenate((self.N, self.N[:, (sample_r_index, )]), axis=1)
            v_new = np.dot(np.array([self.gr_id2r_id2c[grouped_r_id][r_id] for r_id in r_ids]), self.V[indices, :])
            replace_zeros(v_new)
            self.V = np.vstack((self.V, v_new))

            self.r_ids -= set(r_ids)
            self.r_ids.add(grouped_r_id)

    def remove_unused_metabolites(self):
        self.m_ids = {m_id for (m_id, i) in self.m_id2i.iteritems()
                      if get_len(self.N[i, [self.r_id2i[r_id] for r_id in self.r_ids]])}

    def get_used_system(self):
        m_ids = {m_id for (m_id, i) in self.m_id2i.iteritems() if get_len(self.N[i, :])}
        r_ids = {r_id for (r_id, i) in self.r_id2i.iteritems() if get_len(self.V[i, :])}
        r_id2i = {r_id: i for (i, r_id) in enumerate(r_ids)}
        m_id2i = {m_id: i for (i, m_id) in enumerate(m_ids)}
        efm_id2i = self.efm_id2i
        V = self.get_main_V(r_id2i=r_id2i, efm_id2i=efm_id2i)
        N = self.get_main_N(r_id2i=r_id2i, m_id2i=m_id2i)
        return System(N=N, V=V, r_id2i=r_id2i, m_id2i=m_id2i, efm_id2i=efm_id2i,
                      boundary_m_ids=sorted(set(self.boundary_m_ids) & m_ids))

    def get_main_m_id2i(self, ignored_m_ids=None, boundary=False):
        m_ids = set(self.m_ids)
        if boundary:
            m_ids &= set(self.boundary_m_ids)
        if ignored_m_ids:
            m_ids -= ignored_m_ids
        if len(m_ids) == len(self.m_id2i):
            return self.m_id2i
        return {m_id: i for (i, m_id) in enumerate(m_ids)}

    def get_main_bm_ids(self):
        return set(self.boundary_m_ids) & self.m_ids

    def get_main_r_id2i(self):
        if len(self.r_ids) == len(self.r_id2i):
            return self.r_id2i
        return {r_id: i for (i, r_id) in enumerate(self.r_ids)}

    def get_main_efm_id2i(self):
        if len(self.efm_ids) == len(self.efm_id2i):
            return self.efm_id2i
        return {efm_id: i for (i, efm_id) in enumerate(self.efm_ids)}

    def get_main_V(self, r_id2i=None, efm_id2i=None):
        if not r_id2i:
            r_id2i = self.get_main_r_id2i()
        if not efm_id2i:
            efm_id2i = self.get_main_efm_id2i()
        r_is = [self.r_id2i[r_id] for r_id in sorted(r_id2i.iterkeys(), key=lambda r_id: r_id2i[r_id])]
        efm_is = [self.efm_id2i[efm_id] for efm_id in sorted(efm_id2i.iterkeys(), key=lambda efm_id: efm_id2i[efm_id])]
        return self.V[r_is, :][:, efm_is]

    def get_main_N(self, r_id2i=None, m_id2i=None):
        if not r_id2i:
            r_id2i = self.get_main_r_id2i()
        if not m_id2i:
            m_id2i = self.get_main_m_id2i()
        r_is = [self.r_id2i[r_id] for r_id in sorted(r_id2i.iterkeys(), key=lambda r_id: r_id2i[r_id])]
        m_is = [self.m_id2i[m_id] for m_id in sorted(m_id2i.iterkeys(), key=lambda m_id: m_id2i[m_id])]
        return self.N[:, r_is][m_is, :]

    def get_main_S(self):
        r_id2i = self.get_main_r_id2i()
        efm_id2i = self.get_main_efm_id2i()
        m_id2i = self.get_main_m_id2i()
        return System(N=self.get_main_N(r_id2i, m_id2i), V=self.get_main_V(r_id2i, efm_id2i),
                      m_id2i=m_id2i,  r_id2i=r_id2i, efm_id2i=efm_id2i,
                      boundary_m_ids=self.get_main_bm_ids())
