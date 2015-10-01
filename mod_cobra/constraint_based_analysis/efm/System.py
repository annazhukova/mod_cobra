from misc import invert_map
from numpy_efm_manager import get_boundary_metabolites, get_coupled_reactions, lump_coupled_reactions, \
    get_reaction_duplicates, remove_reaction_duplicates, get_efm_duplicates, remove_efm_duplicates, \
    get_efm_groups_based_on_boundary_metabolites, merge_efm_groups, get_len, get_yield, get_control_efficiency, \
    get_efm_intersection

__author__ = 'anna'


class StoichiometricMatrix(object):
    def __init__(self, N, m_id2i, r_id2i, boundary_indices):
        self.N = N
        self.m_id2i = m_id2i
        self.r_id2i = r_id2i
        self.boundary_m_ids = boundary_indices

    def get_inputs_outputs(self, r_id):
        ms = self.N[:, self.r_id2i[r_id]]
        return {(m_id, -ms[i]) for (m_id, i) in self.m_id2i.iteritems() if ms[i] < 0},\
               {(m_id, ms[i]) for (m_id, i) in self.m_id2i.iteritems() if ms[i] > 0}

    def get_reaction_duplicates(self):
        return get_reaction_duplicates(self.N, self.r_id2i)


class PathwaySet(object):
    def __init__(self, V, r_id2i, efm_id2i):
        self.V = V
        self.r_id2i = r_id2i
        self.efm_id2i = efm_id2i

    def get_coupled_reactions(self):
        return get_coupled_reactions(self.V, self.r_id2i)

    def get_efm_duplicates(self):
        return get_efm_duplicates(self.V, self.efm_id2i)

    def get_r_id2coeff(self, efm_id):
        return {r_id: self.V[i, self.efm_id2i[efm_id]] for (r_id, i) in self.r_id2i.iteritems()
                if self.V[i, self.efm_id2i[efm_id]]}

    def get_len(self, efm_id):
        return get_len(self.V[:, self.efm_id2i[efm_id]])

    def get_control_efficiency(self, efm_id, r_id):
        v = self.V[:, self.efm_id2i[efm_id]]
        return get_control_efficiency(v, self.r_id2i[r_id])

    def get_efm_intersection(self):
        return get_efm_intersection(self.V, self.r_id2i)


class System(object):
    def __init__(self, st_matrix=None, pws=None,
                 N=None, V=None, m_id2i=None, r_id2i=None, efm_id2i=None, boundary_m_ids=None):
        if not r_id2i:
            r_id2i = st_matrix.r_id2i if st_matrix else (pws.r_id2i if pws else None)
        self.st_matrix = st_matrix
        if N is not None and m_id2i and r_id2i:
            self.st_matrix = StoichiometricMatrix(N, m_id2i, r_id2i, boundary_m_ids)
        self.pws = pws
        if V is not None and r_id2i and efm_id2i:
            self.pws = PathwaySet(V, r_id2i, efm_id2i)

    @property
    def V(self):
        return self.pws.V

    @property
    def efm_id2i(self):
        return self.pws.efm_id2i

    @property
    def r_id2i(self):
        return self.st_matrix.r_id2i

    @property
    def N(self):
        return self.st_matrix.N

    @property
    def m_id2i(self):
        return self.st_matrix.m_id2i

    @property
    def boundary_m_ids(self):
        return self.st_matrix.boundary_m_ids

    def get_yield(self, efm_id, in_m_id, out_m_id):
        v = self.V[:, self.efm_id2i[efm_id]]
        return get_yield(self.N, v, self.m_id2i[out_m_id], self.m_id2i[in_m_id])

    def get_boundary_inputs_outputs(self, pathway_id):
        v = self.V[:, self.efm_id2i[pathway_id]]
        b_ms = get_boundary_metabolites(self.N[[self.m_id2i[m_id] for m_id in self.boundary_m_ids], :], v)
        bm_id2i = dict(zip(self.boundary_m_ids, xrange(0, len(self.boundary_m_ids))))
        # b_ms = get_boundary_metabolites(self.N, v)
        r_id2st = {(m_id, -b_ms[bm_id2i[m_id]]) for m_id in self.boundary_m_ids if b_ms[bm_id2i[m_id]] < 0}
        p_id2st = {(m_id, b_ms[bm_id2i[m_id]]) for m_id in self.boundary_m_ids if b_ms[bm_id2i[m_id]] > 0}
        return r_id2st, p_id2st

    def get_efm_groups_based_on_boundary_metabolites(self):
        return get_efm_groups_based_on_boundary_metabolites(
            self.N[[self.m_id2i[m_id] for m_id in self.boundary_m_ids], :], self.V, self.efm_id2i)

    def lump_coupled_reactions(self):
        coupled_r_id_groups = self.pws.get_coupled_reactions()
        N_c, V_c, c_r_id2i, r_id2lr_id, lr_id2r_id2c = \
            lump_coupled_reactions(self.N, self.V, coupled_r_id_groups, self.r_id2i)
        return FoldedSystem(N=N_c, V=V_c, m_id2i=self.m_id2i, r_id2i=c_r_id2i, efm_id2i=self.efm_id2i,
                            r_id2gr_id=r_id2lr_id, gr_id2r_id2c=lr_id2r_id2c, boundary_m_ids=self.boundary_m_ids)

    def remove_reaction_duplicates(self):
        duplicated_r_id_groups = self.st_matrix.get_reaction_duplicates()
        N_d, V_d, d_r_id2i, r_id2gr_id, gr_id2r_id2c = \
            remove_reaction_duplicates(self.N, self.V, duplicated_r_id_groups, self.r_id2i)
        return FoldedSystem(N=N_d, V=V_d, m_id2i=self.m_id2i, r_id2i=d_r_id2i, efm_id2i=self.efm_id2i,
                            r_id2gr_id=r_id2gr_id, gr_id2r_id2c=gr_id2r_id2c, boundary_m_ids=self.boundary_m_ids)

    def merge_efm_groups(self):
        boundary_efm_id_groups = self.get_efm_groups_based_on_boundary_metabolites()
        V_m, m_efm_id2i, efm_id2merged_id = merge_efm_groups(self.V, boundary_efm_id_groups, self.efm_id2i)
        return FoldedSystem(st_matrix=self.st_matrix, V=V_m, efm_id2i=m_efm_id2i, efm_id2gr_id=efm_id2merged_id)


class FoldedSystem(System):
    def __init__(self, r_id2gr_id=None, gr_id2r_id2c=None, efm_id2gr_id=None, st_matrix=None, pws=None,
                 N=None, V=None, m_id2i=None, r_id2i=None, efm_id2i=None, boundary_m_ids=None):
        System.__init__(self, st_matrix=st_matrix, pws=pws, N=N, V=V, m_id2i=m_id2i, r_id2i=r_id2i, efm_id2i=efm_id2i,
                        boundary_m_ids=boundary_m_ids)
        self.r_id2gr_id = r_id2gr_id
        self.gr_id2r_id2c = gr_id2r_id2c
        self.efm_id2gr_id = efm_id2gr_id if efm_id2gr_id else {}
        self.gr_id2efm_ids = invert_map(self.efm_id2gr_id)

    def remove_efm_duplicates(self):
        duplicated_efm_id_groups = self.pws.get_efm_duplicates()
        V_dd, d_efm_id2i, efm_id2gr_id = remove_efm_duplicates(self.V, duplicated_efm_id_groups, self.efm_id2i)
        return FoldedSystem(st_matrix=self.st_matrix, V=V_dd, efm_id2i=d_efm_id2i, efm_id2gr_id=efm_id2gr_id,
                            r_id2gr_id=self.r_id2gr_id, gr_id2r_id2c=self.gr_id2r_id2c)
