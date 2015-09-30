import os
import unittest

import libsbml

from SBMLTestCase import create_test_sbml, TEST_SBML
from numpy_efm_manager import get_element2id_mapping, get_stoichiometric_matrix, get_efm_matrix, get_control_efficiency, \
    get_yield, get_coupled_reactions, lump_coupled_reactions, get_reaction_duplicates, remove_reaction_duplicates, \
    get_efm_duplicates, remove_efm_duplicates, get_boundary_metabolites, get_efm_groups_based_on_boundary_metabolites, \
    merge_efm_groups, get_len, remove_invalid_efms
import numpy as np

__author__ = 'anna'


class NumpyEfmTestCase(unittest.TestCase):

    def setUp(self):
        """
                         r2-> m2 <-r3-> m2_b
                      <-      ^
        m1_b <-r1-> m1         r6
                      <-      |
                         r4-> m3 --r5-> m3_b


        Expect to find 2 EFMs that include r3:
            10 r1	10 r2	10 r3
            10 r1	10 r3	10 r4	10 r6
        """
        create_test_sbml()
        self.doc = libsbml.SBMLReader().readSBML(TEST_SBML)
        self.model =self.doc.getModel()
        self.s_id2i, self.r_id2i = get_element2id_mapping(self.model)

    def tearDown(self):
        if os.path.exists(TEST_SBML):
            os.remove(TEST_SBML)

    def test_st_matrix_product(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        self.assertEqual(1, N[self.s_id2i['m2'], self.r_id2i['r6']],
                         'N[m2, r6] was supposed to be 1, got %g instead.' % N[self.s_id2i['m2'], self.r_id2i['r6']])

    def test_st_matrix_reactant(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        self.assertEqual(-1, N[self.s_id2i['m1'], self.r_id2i['r2']],
                         'N[m1, r2] was supposed to be -1, got %g instead.' % N[self.s_id2i['m2'], self.r_id2i['r6']])

    def test_st_matrix_unrelated_metabolite(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        self.assertEqual(0, N[self.s_id2i['m2'], self.r_id2i['r1']],
                         'N[m2, r1] was supposed to be -1, got %g instead.' % N[self.s_id2i['m2'], self.r_id2i['r6']])

    def test_efm_conversion_r_coefficient(self):
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        self.assertEqual(10, V[self.r_id2i['r1'], 0],
                         'V[r1, efm0] was supposed to be 10, got %g instead.' % V[self.r_id2i['r1'], 0])

    def test_efm_conversion_absent_r_coefficient(self):
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        self.assertEqual(0, V[self.r_id2i['r4'], 0],
                         'V[r4, efm0] was supposed to be 0, got %g instead.' % V[self.r_id2i['r1'], 0])

    def test_control_efficiency(self):
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        ce = get_control_efficiency(V[:, 0], r_index=self.r_id2i['r3'])
        self.assertEqual(1/3.0, ce,
                         'Control efficiency of r3 in EFM 0 was supposed to be 1/3, got %g instead.' % ce)

    def test_yield(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        y = get_yield(N, V[:, 0], out_m_index=self.s_id2i['m2_b'], in_m_index=self.s_id2i['m1_b'])
        self.assertEqual(1, y,
                         'Yield of m2_b with respect to m1_b in EFM 0 was supposed to be 1, got %g instead.' % y)

    def test_yield_no_in_metabolite(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        y = get_yield(N, V[:, 0], out_m_index=self.s_id2i['m2_b'], in_m_index=self.s_id2i['m3_b'])
        self.assertEqual(float('inf'), y,
                         'Yield of m2_b with respect to m3_b in EFM 0 was supposed to be inf, got %g instead.' % y)

    def test_yield_no_out_metabolite(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        y = get_yield(N, V[:, 0], out_m_index=self.s_id2i['m3_b'], in_m_index=self.s_id2i['m1_b'])
        self.assertEqual(0, y,
                         'Yield of m3_b with respect to m1_b in EFM 0 was supposed to be 0, got %g instead.' % y)

    def test_coupled_reactions_len(self):
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        c_r_id_lists = get_coupled_reactions(V, self.r_id2i)
        self.assertEqual(1, len(c_r_id_lists),
                         'Was supposed to find one coupled reaction group, got %g instead.' % len(c_r_id_lists))

    def test_coupled_reactions_content(self):
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        coupled_r_id_groups = get_coupled_reactions(V, self.r_id2i)
        self.assertIn(['r1', 'r3'], coupled_r_id_groups, 'Was supposed to find [r1, r3] coupled reaction group.')

    def test_lump_coupled_reactions_N_1(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 1, 'r3': 1, 'r4': 1, 'r6': 1}], self.r_id2i)
        coupled_r_id_groups = get_coupled_reactions(V, self.r_id2i)
        N_new, V_new, new_r_id2i, r_id2lr_id = lump_coupled_reactions(N, V, coupled_r_id_groups, self.r_id2i)
        lr_id = next(r_id2lr_id.itervalues())
        st = N_new[self.s_id2i['m1_b'], new_r_id2i[lr_id]]
        self.assertEqual(-1, st, 'The lumped reaction was supposed to consume 1 m1_b, consumes %g instead' % -st)

    def test_lump_coupled_reactions_N_2(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 1, 'r3': 1, 'r4': 1, 'r6': 1}], self.r_id2i)
        coupled_r_id_groups = get_coupled_reactions(V, self.r_id2i)
        N_new, V_new, new_r_id2i, r_id2lr_id = lump_coupled_reactions(N, V, coupled_r_id_groups, self.r_id2i)
        lr_id = next(r_id2lr_id.itervalues())
        st = N_new[self.s_id2i['m2_b'], new_r_id2i[lr_id]]
        self.assertEqual(1, st, 'The lumped reaction was supposed to produce 1 m2_b, consumes %g instead' % -st)

    def test_lump_coupled_reactions_V(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 1, 'r3': 1, 'r4': 1, 'r6': 1}], self.r_id2i)
        coupled_r_id_groups = get_coupled_reactions(V, self.r_id2i)
        N_new, V_new, new_r_id2i, r_id2lr_id = lump_coupled_reactions(N, V, coupled_r_id_groups, self.r_id2i)
        lr_id = next(r_id2lr_id.itervalues())
        r_i = V_new[new_r_id2i[lr_id], :]
        self.assertListEqual([1, 0.1], list(r_i),
                             'The lumped reaction was supposed to participate in EFMs with coefficients [1, 0.1] not %s'
                             % r_i)

    def test_reaction_duplicates_len(self):
        N = np.array([[1, 2, 0, -1],
                      [0, 0, 1, 0],
                      [-2, -4, -1, 2]])
        r_id2i = {'r1': 0, 'r2': 1, 'r3': 2, 'r4': 3}
        r_duplicates = get_reaction_duplicates(N, r_id2i)
        self.assertEqual(1, len(r_duplicates),
                         'Was supposed to find one group of duplicates, found %d' % len(r_duplicates))

    def test_reaction_duplicates_content(self):
        N = np.array([[1, 2, 0, -1],
                      [0, 0, 1, 0],
                      [-2, -4, -1, 2]])
        r_id2i = {'r1': 0, 'r2': 1, 'r3': 2, 'r4': 3}
        r_duplicates = get_reaction_duplicates(N, r_id2i)
        self.assertSetEqual({'r1', 'r2', 'r4'}, set(r_duplicates[0]),
                            'Was supposed to find [r1, r2, r4] as duplicates, found %s' % r_duplicates[0])

    def test_remove_reaction_duplicates_N_shape(self):
        N = np.array([[1, 2, 0, -1],
                      [0, 0, 1, 0],
                      [-2, -4, -1, 2]])
        V = np.array([[-1, 0],
                      [1, 1],
                      [0, -1],
                      [2, -1]])
        r_id2i = {'r1': 0, 'r2': 1, 'r3': 2, 'r4': 3}
        r_duplicates = get_reaction_duplicates(N, r_id2i)
        N_new, V_new, new_r_id2i, r_id2gr_id, r_id2lambda = remove_reaction_duplicates(N, V, r_duplicates, r_id2i)
        self.assertTupleEqual((3, 2), N_new.shape,
                              'Was supposed to get 3x2 N, got %s' % [str(it) for it in N_new.shape])

    def test_remove_reaction_duplicates_N_content_1(self):
        N = np.array([[1, 2, 0, -1],
                      [0, 0, 1, 0],
                      [-2, -4, -1, 2]])
        V = np.array([[-1, 0],
                      [1, 1],
                      [0, -1],
                      [2, -1]])
        r_id2i = {'r1': 0, 'r2': 1, 'r3': 2, 'r4': 3}
        r_duplicates = get_reaction_duplicates(N, r_id2i)
        N_new, V_new, new_r_id2i, r_id2gr_id, r_id2lambda = remove_reaction_duplicates(N, V, r_duplicates, r_id2i)
        m2_st = N[1, new_r_id2i[r_id2gr_id['r1']]]
        self.assertEqual(0, m2_st,
                         'Was supposed to find no metabolite m2 in the grouped reaction, found %g' % m2_st)

    def test_remove_reaction_duplicates_N_content_2(self):
        N = np.array([[1, 2, 0, -1],
                      [0, 0, 1, 0],
                      [-2, -4, -1, 2]])
        V = np.array([[-1, 0],
                      [1, 1],
                      [0, -1],
                      [2, -1]])
        r_id2i = {'r1': 0, 'r2': 1, 'r3': 2, 'r4': 3}
        r_duplicates = get_reaction_duplicates(N, r_id2i)
        N_new, V_new, new_r_id2i, r_id2gr_id, r_id2lambda = remove_reaction_duplicates(N, V, r_duplicates, r_id2i)
        new_r_index = new_r_id2i[r_id2gr_id['r1']]
        m1_st = N[0, new_r_index]
        m3_st = N[2, new_r_index]
        self.assertEqual(-2, m3_st / m1_st,
                         'Was supposed to get -2 as a proportion between m3 and m1 in the new reaction, found %g'
                         % (m3_st / m1_st))

    def test_remove_reaction_duplicates_V_shape(self):
        N = np.array([[1, 2, 0, -1],
                      [0, 0, 1, 0],
                      [-2, -4, -1, 2]])
        V = np.array([[-1, 0],
                      [1, 1],
                      [0, -1],
                      [2, -1]])
        r_id2i = {'r1': 0, 'r2': 1, 'r3': 2, 'r4': 3}
        r_duplicates = get_reaction_duplicates(N, r_id2i)
        N_new, V_new, new_r_id2i, r_id2gr_id, r_id2lambda = remove_reaction_duplicates(N, V, r_duplicates, r_id2i)
        self.assertTupleEqual((2, 2), V_new.shape,
                              'Was supposed to get 2x2 V, got %s' % [str(it) for it in V_new.shape])

    def test_remove_reaction_duplicates_V_content(self):
        N = np.array([[1, 2, 0, -1],
                      [0, 0, 1, 0],
                      [-2, -4, -1, 2]])
        V = np.array([[-1, 0],
                      [1, 1],
                      [0, -1],
                      [2, -1]])
        r_id2i = {'r1': 0, 'r2': 1, 'r3': 2, 'r4': 3}
        r_duplicates = get_reaction_duplicates(N, r_id2i)
        N_new, V_new, new_r_id2i, r_id2gr_id, r_id2lambda = remove_reaction_duplicates(N, V, r_duplicates, r_id2i)
        new_r_index = new_r_id2i[r_id2gr_id['r1']]
        efm1_c = V_new[new_r_index, 0]
        efm2_c = V_new[new_r_index, 1]
        ratio = efm2_c / efm1_c
        self.assertEqual(-3, ratio,
                         'Was supposed to get -3 as a ratio between the coefficients of the grouped reaction in EFMs 2 and 1, got %g'
                         % ratio)

    def test_efm_duplicates_len(self):
            V = np.array([[1, 2, -1, -1],
                          [2, 4, -2, 0],
                          [-1, -2, 1, 2]])
            efm_id2i = {0: 0, 1: 1, 2: 2, 3: 3}
            efm_duplicates = get_efm_duplicates(V, efm_id2i)
            self.assertEqual(1, len(efm_duplicates),
                             'Was supposed to find one group of duplicates, found %d' % len(efm_duplicates))

    def test_efm_duplicates_content(self):
        V = np.array([[1, 2, -1, -1],
                      [2, 4, -2, 0],
                      [-1, -2, 1, 2]])
        efm_id2i = {0: 0, 1: 1, 2: 2, 3: 3}
        efm_duplicates = get_efm_duplicates(V, efm_id2i)
        self.assertSetEqual({0, 1}, set(efm_duplicates[0]),
                            'Was supposed to find [0, 1] as duplicates, found %s' % efm_duplicates[0])

    def test_remove_efm_duplicates_V_shape(self):
        V = np.array([[1, 2, -1, -1],
                      [2, 4, -2, 0],
                      [-1, -2, 1, 2]])
        efm_id2i = {0: 0, 1: 1, 2: 2, 3: 3}
        efm_duplicates = get_efm_duplicates(V, efm_id2i)
        V_new, new_efm_id2i, efm_id2gr_id = remove_efm_duplicates(V, efm_duplicates, efm_id2i)
        self.assertTupleEqual((3, 3), V_new.shape,
                              'Was supposed to get 3x3 V, got %s' % [str(it) for it in V_new.shape])

    def test_remove_efm_duplicates_V_content(self):
        V = np.array([[1, 2, -1, -1],
                      [2, 4, -2, 0],
                      [-1, -2, 1, 2]])
        efm_id2i = {0: 0, 1: 1, 2: 2, 3: 3}
        efm_duplicates = get_efm_duplicates(V, efm_id2i)
        V_new, new_efm_id2i, efm_id2gr_id = remove_efm_duplicates(V, efm_duplicates, efm_id2i)
        new_efm_index = new_efm_id2i[efm_id2gr_id[0]]
        r1_c = V_new[0, new_efm_index]
        r2_c = V_new[1, new_efm_index]
        ratio = r2_c / r1_c
        self.assertEqual(2, ratio,
                         'Was supposed to get 2 as a ratio between the coefficients of r1 and r2 in the grouped EFMs, got %g'
                         % ratio)

    def test_boundary_metabolites(self):
        N = get_stoichiometric_matrix(self.model, self.s_id2i, self.r_id2i)
        V = get_efm_matrix([{'r1': 10, 'r2': 10, 'r3': 10}, {'r1': 10, 'r3': 10, 'r4': 10, 'r6': 10}], self.r_id2i)
        m = get_boundary_metabolites(N, V[:, 0])
        result = {(m[self.s_id2i[m_id]], m_id) for m_id in self.s_id2i.iterkeys() if m[self.s_id2i[m_id]]}
        self.assertSetEqual({(-1, 'm1_b'), (1, 'm2_b')},
                            result,
                            'Was supposed to get -1 m1_b 1 m2_b as boundary metabolites, got %s' % result)


    def test_remove_invalid_efms_shape(self):
        N = np.array([[-1, 0, -1, 0],
                      [1, -1, 0, 0],
                      [0, 1, 1, 1],
                      [0, 0, 0, -1]])
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [-1, 10, 0],
                      [0, 10, 0]])
        V_new = remove_invalid_efms(N, V)
        self.assertTupleEqual((4, 1), V_new.shape,
                              'Was supposed to get 4x1 V, got %s' % [str(it) for it in V_new.shape])

    def test_remove_invalid_efms_content(self):
        N = np.array([[-1, 0, -1, 0],
                      [1, -1, 0, 0],
                      [0, 1, 1, 1],
                      [0, 0, 0, -1]])
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [-1, 10, 0],
                      [0, 10, 0]])
        v_new = remove_invalid_efms(N, V)[:, 0]
        self.assertListEqual(list(V[:,0]), list(v_new), 'EFM 0 was supposed to be valid')

    def test_get_efm_groups_based_on_b_ms_size(self):
        """
           -r1-> m2 -r2->
        m1                m3
           -r3-> m4 -r4->
        :return:
        """
        N = np.array([[-1, 0, -1, 0],
                      [1, -1, 0, 0],
                      [0, 1, 0, 1],
                      [0, 0, 1, -1]])
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [0, 10, 0],
                      [0, 10, 0]])
        efm_id2i = {0: 0, 1: 1, 2: 2}
        groups = get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i)
        self.assertEqual(1, len(groups), 'Was supposed to find 1 group, got %d' % len(groups))

    def test_get_efm_groups_based_on_b_ms_content(self):
        """
           -r1-> m2 -r2->
        m1                m3
           -r3-> m4 -r4->
        :return:
        """
        N = np.array([[-1, 0, -1, 0],
                      [1, -1, 0, 0],
                      [0, 1, 0, 1],
                      [0, 0, 1, -1]])
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [0, 10, 0],
                      [0, 10, 0]])
        efm_id2i = {0: 0, 1: 1, 2: 2}
        groups = get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i)
        self.assertSetEqual({0, 1}, set(groups[0]),
                            'Was supposed to find group {0, 1}, got %s' % {str(it) for it in groups[0]})

    def test_merge_efm_groups_V_shape(self):
        N = np.array([[-1, 0, -1, 0],
                      [1, -1, 0, 0],
                      [0, 1, 0, 1],
                      [0, 0, 1, -1]])
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [0, 10, 0],
                      [0, 10, 0]])
        efm_id2i = {0: 0, 1: 1, 2: 2}
        groups = get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i)
        V_new, new_efm_id2i, efm_id2gr_id = merge_efm_groups(V, groups, efm_id2i)
        self.assertTupleEqual((4, 2), V_new.shape,
                              'Was supposed to get 4x2 V, got %s' % [str(it) for it in V_new.shape])

    def test_merge_efm_groups_V_content_1(self):
        N = np.array([[-1, 0, -1, 0],
                      [1, -1, 0, 0],
                      [0, 1, 0, 1],
                      [0, 0, 1, -1]])
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [0, 10, 0],
                      [0, 10, 0]])
        efm_id2i = {0: 0, 1: 1, 2: 2}
        groups = get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i)
        V_new, new_efm_id2i, efm_id2gr_id = merge_efm_groups(V, groups, efm_id2i)
        new_efm_index = new_efm_id2i[efm_id2gr_id[0]]
        r1 = V_new[0, new_efm_index] / V_new[1, new_efm_index]
        self.assertEqual(1, r1, 'Ratio r1/r2 in the new EFM was supposed to be 1, got %d' % r1)

    def test_merge_efm_groups_V_content_2(self):
        N = np.array([[-1, 0, -1, 0],
                      [1, -1, 0, 0],
                      [0, 1, 0, 1],
                      [0, 0, 1, -1]])
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [0, 10, 0],
                      [0, 10, 0]])
        efm_id2i = {0: 0, 1: 1, 2: 2}
        groups = get_efm_groups_based_on_boundary_metabolites(N, V, efm_id2i)
        V_new, new_efm_id2i, efm_id2gr_id = merge_efm_groups(V, groups, efm_id2i)
        new_efm_index = new_efm_id2i[efm_id2gr_id[0]]
        r2 = V_new[2, new_efm_index] / V_new[3, new_efm_index]
        self.assertEqual(1, r2, 'Ratio r3/r4 in the new EFM was supposed to be 1, got %d' % r2)

    def test_len(self):
        V = np.array([[1, 0, 1],
                      [1, 0, 0],
                      [0, 10, 0],
                      [0, 10, 0]])
        v = V[:, 0]
        v_len = get_len(v)
        self.assertEqual(2, v_len, 'EFM length was supposed to be 2, got %d' % v_len)