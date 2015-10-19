import sys

import libsbml
from mod_sbml.annotation.chebi.chebi_annotator import get_species_to_chebi
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi

from examples.models import RECON_MODEL_BOUNDARY, RECON_MODEL
from mod_sbml.sbml.compartment.compartment_manager import create_boundary_species_in_boundary_reactions
from mod_sbml.onto import parse_simple

__author__ = 'anna'


def main():
    add_boundary_metabolites(RECON_MODEL, RECON_MODEL_BOUNDARY)


def add_boundary_metabolites(in_sbml, out_sbml):
    """
    Creates a boundary compartment with the id 'Boundary',
    and transforms each input/output reaction '*_e <-> ' into '*_e <-> *_b'.
    :param in_sbml: path to the SBML file with the original model
    :param out_sbml: path where to store the resulting SBML file
    """
    input_doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = input_doc.getModel()
    m_id2chebi_id = get_species_to_chebi(model, parse_simple(get_chebi()))
    create_boundary_species_in_boundary_reactions(model, m_id2chebi_id)
    libsbml.SBMLWriter().writeSBMLToFile(input_doc, out_sbml)


if "__main__" == __name__:
    sys.exit(main())
