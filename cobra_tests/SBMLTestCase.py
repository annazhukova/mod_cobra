import os
import unittest

import libsbml

from metatool_manager import convert_metabolite, convert_reaction
import cobra_tests

__author__ = 'anna'

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(cobra_tests.__file__)), 'data')
TEST_SBML = os.path.join(DATA_DIR, 'test.xml')


def create_test_sbml():
    """
                     r2-> m2 <-r3-> m2_b
                  <-      ^
    m1_b <-r1-> m1         r6
                  <-      |
                     r4-> m3 --r5-> m3_b
    """
    document = libsbml.SBMLDocument(2, 4)
    model = document.createModel()
    model.setId('test_model')
    c = model.createCompartment()
    c_id = 'c'
    c.setId(c_id)
    b = model.createCompartment()
    b_id = 'b'
    b.setId(b_id)
    for m_id in ('m1', 'm2', 'm3'):
        convert_metabolite(m_id, model, False, c_id='c')
    for m_id in ('m1_b', 'm2_b', 'm3_b'):
        convert_metabolite(m_id, model, True, c_id='b')
    for (r_id, r_m_id2st, p_m_id2st, rev) in (('r1', {'m1_b': 1}, {'m1': 1}, True), ('r2', {'m1': 1}, {'m2': 1}, True),
                                              ('r3', {'m2': 1}, {'m2_b': 1}, True), ('r4', {'m1': 1}, {'m3': 1}, True),
                                              ('r5', {'m3': 1}, {'m3_b': 1}, False),
                                              ('r6', {'m3': 1}, {'m2': 1}, False)):
        convert_reaction(model, r_id, rev, r_m_id2st, p_m_id2st)
    libsbml.SBMLWriter().writeSBMLToFile(document, TEST_SBML)


class SBMLTestCase(unittest.TestCase):
    def setUp(self):
        create_test_sbml()

    def tearDown(self):
        if os.path.exists(TEST_SBML):
            os.remove(TEST_SBML)

    def test_species_num(self):
        doc = libsbml.SBMLReader().readSBML(TEST_SBML)
        model = doc.getModel()
        num_sps = model.getNumSpecies()
        self.assertEqual(6, num_sps, "Number of species was supposed to be 6, got %d" % num_sps)

    def test_reactions_num(self):
        doc = libsbml.SBMLReader().readSBML(TEST_SBML)
        model = doc.getModel()
        num_rs = model.getNumReactions()
        self.assertEqual(6, num_rs, "Number of reactions was supposed to be 6, got %d" % num_rs)

    def test_reactants_num(self):
        doc = libsbml.SBMLReader().readSBML(TEST_SBML)
        model = doc.getModel()
        r = model.getReaction('r6')
        self.assertEqual(1, r.getNumReactants(), "Number of reactants of r6 was supposed to be 1, got %d"
                         % r.getNumReactants())

    def test_products_num(self):
        doc = libsbml.SBMLReader().readSBML(TEST_SBML)
        model = doc.getModel()
        r = model.getReaction('r6')
        self.assertEqual(1, r.getNumProducts(), "Number of products of r6 was supposed to be 1, got %d"
                         % r.getNumProducts())

    def test_comps_num(self):
        doc = libsbml.SBMLReader().readSBML(TEST_SBML)
        model = doc.getModel()
        num_comps = model.getNumCompartments()
        self.assertEqual(2, num_comps, "Number of compartments was supposed to be 2, got %d" % num_comps)
