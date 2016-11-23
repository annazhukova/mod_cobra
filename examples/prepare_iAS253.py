import logging
import re
import sys

import libsbml

from mod_sbml.annotation.annotator import annotate
from mod_sbml.annotation.gene_ontology.go_serializer import get_go
from mod_sbml.annotation.gene_ontology.go_annotator import annotate_compartments
from mod_sbml.annotation.kegg.pathway_manager import ORG_HUMAN
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from models import SR_MODEL, SR_MODEL_ANNOTATED, RECON_MODEL, SR_MODEL_SER
from mod_cobra.sbml.mapping.metabolite_matcher import get_model_data, map_metabolites_compartments
from mod_sbml.sbml.compartment.compartment_manager import BOUNDARY_C_ID
from mod_sbml.annotation.chebi.chebi_annotator import get_chebi_id, annotate_metabolites, CHEBI_PREFIX
from mod_sbml.annotation.kegg.kegg_annotator import get_kegg_m_id, get_kegg_r_id, KEGG_COMPOUND_PREFIX, KEGG_REACTION_PREFIX
from mod_sbml.sbml.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import create_reaction, create_species, get_reactants, get_products, get_r_comps, generate_unique_id
from mod_sbml.annotation.rdf_annotation_helper import add_annotation
from mod_sbml.onto import parse_simple

__author__ = 'anna'


def main():
    convert_annotations(SR_MODEL, SR_MODEL_ANNOTATED)
    separate_boundary_species(SR_MODEL_ANNOTATED, SR_MODEL_ANNOTATED)
    create_serine_sythesis(SR_MODEL_ANNOTATED, RECON_MODEL, SR_MODEL_SER)
    input_doc = libsbml.SBMLReader().readSBML(SR_MODEL_SER)
    model = input_doc.getModel()
    annotate(model, pw_threshold=0.5, org=ORG_HUMAN)
    libsbml.SBMLWriter().writeSBMLToFile(input_doc, SR_MODEL_SER)


def convert_annotations(in_sbml, out_sbml):
    """
    Converts the initial iAS253 model by Smith and Robinson (BMC Syst Biol 2011; 5: 102),
    which uses KEGG ids followed by compartment (e.g., C00001Cyto) as identifiers of the elements,
    into a model with KEGG reaction and metabolite annotations.
    :param in_sbml: path to the SBML file with the original iAS253 model
    :param out_sbml: path where to store the annotated SBML file
    """
    input_doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = input_doc.getModel()
    for species in model.getListOfSpecies():
        _id = species.getId()
        s_ids = re.findall(r'C\d{4,7}', _id)
        if s_ids:
            add_annotation(species, libsbml.BQB_IS, s_ids.pop(), KEGG_COMPOUND_PREFIX)
        if species.getName().find(_id) != -1:
            species.setName(species.getName().replace(_id, "").strip())
    for reaction in model.getListOfReactions():
        _id = reaction.getId()
        r_ids = re.findall(r'R\d{4,10}', _id)
        if r_ids:
            add_annotation(reaction, libsbml.BQB_IS, r_ids.pop(), KEGG_REACTION_PREFIX)
        if reaction.getName().find(_id) != -1:
            reaction.setName(reaction.getName().replace(_id, "").strip())
    libsbml.SBMLWriter().writeSBMLToFile(input_doc, out_sbml)


def separate_boundary_species(in_sbml, out_sbml):
    """
    Creates a boundary compartment with the id 'Boundary' and moves the boundary species (with the ids '*_b') there.
    :param in_sbml: path to the SBML file with the original model
    :param out_sbml: path where to store the modified SBML file
    :return: void
    """
    input_doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = input_doc.getModel()
    boundary_comp = model.createCompartment()
    id_ = generate_unique_id(model, BOUNDARY_C_ID)
    if libsbml.LIBSBML_OPERATION_SUCCESS != boundary_comp.setId(id_):
        logging.error("boundary compartment  %s creation error" % id_)
    boundary_comp.setName("Boundary")
    for s in model.getListOfSpecies():
        if s.getId().find("_b") != -1:
            s.setCompartment(boundary_comp.getId())
    libsbml.SBMLWriter().writeSBMLToFile(input_doc, out_sbml)


def create_serine_sythesis(sr_in_sbml, recon_sbml, sr_out_sbml):
    sr_doc = libsbml.SBMLReader().readSBML(sr_in_sbml)
    sr_model = sr_doc.getModel()

    r_doc = libsbml.SBMLReader().readSBML(recon_sbml)
    r_model = r_doc.getModel()

    # NADPH --> NADH; NADP --> NAD
    metabolite_id_mapping_recon2sr = {'M_h_c': 'C00080Cyto', 'M_nadph_c': 'C00004Cyto', 'M_nadp_c': 'C00003Cyto',
                                      'M_glu_L_c': 'C00025Cyto'}
    add_recon_reactions(sr_model, r_model, {'R_PGCD', 'R_PSERT', 'R_PSP_L', "R_PEPCK", "R_ME2", "R_NDPK1", "R_GHMT2r"},
                        metabolite_id_mapping_recon2sr,
                        recon_m_id2boundary_condition_to_add={'M_gtp_c': False, 'M_pser_L_c': False, 'M_3php_c': False,
                                                              'M_mlthf_c': True, 'M_thf_c': True})
    # Reversible transporters
    for tr_id in ['Transporter30', 'Transporter1', 'Transporter9', 'Transporter17']:
        r = sr_model.getReaction(tr_id)
        make_reversible(r)

    # Reversible glutamate dehydrogenase
    for r_id in ['R00243MM', 'R00248MM']:
        r = sr_model.getReaction(r_id)
        make_reversible(r)

    # allow for proton input
    r = sr_model.getReaction('Boundary60')  # H+(C00080_b) <=> H+(C00080Cyto)
    l_b, u_b = get_bounds(r)
    if u_b == 0:
        u_b = -l_b if l_b else 1000
    set_bounds(r, l_b, u_b)

    # allow for NH3 transport Mito -> Cyto
    r = sr_model.getReaction('Transporter85')  # NH3(C00014MM) <=> NH3(C00014Cyto)
    l_b, u_b = get_bounds(r)
    if u_b == 0:
        u_b = -l_b if l_b else 1000
    set_bounds(r, l_b, u_b)

    # create NH3 output: NH3_cyto -> NH3_b
    nh3_cyto = sr_model.getSpecies('C00014Cyto')
    chebi_id = get_chebi_id(nh3_cyto)
    kegg_id = get_kegg_m_id(nh3_cyto)
    nh3_b = create_species(sr_model, compartment_id=BOUNDARY_C_ID, name=nh3_cyto.getName(), bound=True, id_='C00014_b')
    if kegg_id:
        add_annotation(nh3_b, libsbml.BQB_IS, kegg_id.upper(), KEGG_COMPOUND_PREFIX)
    if chebi_id:
        add_annotation(nh3_b, libsbml.BQB_IS, chebi_id.upper(), CHEBI_PREFIX)
    r = create_reaction(sr_model, {nh3_b.getId(): 1}, {nh3_cyto.getId(): 1}, "NH3", True, "BoundaryNH3")
    set_bounds(r, -1000, 0)

    # prohibit ATP exit
    r = sr_model.getReaction('Boundary59')  # ATP
    if r:
        l_b, u_b = get_bounds(r)
        set_bounds(r, max(0, l_b), u_b)

    # constraint pyruvate kinase reactions
    for r in sr_model.getListOfReactions():
        rs_names, ps_names = {sr_model.getSpecies(m_id).getName() for m_id in get_reactants(r)}, \
                             {sr_model.getSpecies(m_id).getName() for m_id in get_products(r)}
        if 'Phosphoenolpyruvate' in rs_names and 'Pyruvate' in ps_names and \
                next((m for m in rs_names if m[-2:] == 'DP'), None):
            l_b, u_b = get_bounds(r)
            set_bounds(r, max(0, l_b), u_b)
        elif 'Phosphoenolpyruvate' in ps_names and 'Pyruvate' in rs_names and \
                next((m for m in ps_names if m[-2:] == 'DP'), None):
            l_b, u_b = get_bounds(r)
            set_bounds(r, l_b, min(0, u_b))

    libsbml.SBMLWriter().writeSBMLToFile(sr_doc, sr_out_sbml)


def add_recon_reactions(sr_model, r_model, recon_r_ids_to_add, metabolite_id_mapping_recon2sr=None,
                        recon_m_id2boundary_condition_to_add=None):
    RECON2_ID = 'Recon2'
    SR_ID = 'SR'
    model_id2sbml = {SR_ID: sr_model, RECON2_ID: r_model}
    go = parse_simple(get_go())
    annotate_compartments(sr_model, go)
    annotate_compartments(r_model, go)

    chebi = parse_simple(get_chebi())
    annotate_metabolites(sr_model, chebi)
    annotate_metabolites(r_model, chebi)

    model_id2dfs = get_model_data(model_id2sbml)
    model_id2c_ids_groups, model_id2m_ids_groups, model_id2c_id2i = \
        map_metabolites_compartments(model_id2dfs, chebi=chebi)

    c_id_recon2sr = {}
    c_id_sr2recon = {}
    m_id_recon2sr = {}
    for model_id2m_ids in model_id2m_ids_groups:
        m_id_recon2sr.update({it: next(iter(model_id2m_ids[SR_ID])) for it in model_id2m_ids[RECON2_ID]})
    for model_id2c_ids in model_id2c_ids_groups:
        c_id_recon2sr.update({it: next(iter(model_id2c_ids[SR_ID])) for it in model_id2c_ids[RECON2_ID]})
        c_id_sr2recon.update({it: next(iter(model_id2c_ids[RECON2_ID])) for it in model_id2c_ids[SR_ID]})

    if metabolite_id_mapping_recon2sr:
        m_id_recon2sr.update(metabolite_id_mapping_recon2sr)

    if recon_m_id2boundary_condition_to_add:
        for recon_m_id in recon_m_id2boundary_condition_to_add.keys():
            df = model_id2dfs[RECON2_ID][0]
            row = df.loc[recon_m_id]
            kegg_id = row['KEGG']
            chebi_id = row['ChEBI']
            name = row['Name']
            c_id = c_id_recon2sr[row['Compartment']]
            sr_m = create_species(name=name, compartment_id=c_id,
                                     model=sr_model, bound=recon_m_id2boundary_condition_to_add[recon_m_id],
                                     id_=(recon_m_id if not kegg_id else kegg_id.upper()) + c_id)
            if kegg_id:
                add_annotation(sr_m, libsbml.BQB_IS, kegg_id.upper(), KEGG_COMPOUND_PREFIX)
            if chebi_id:
                add_annotation(sr_m, libsbml.BQB_IS, chebi_id.upper(), CHEBI_PREFIX)
            m_id_recon2sr[recon_m_id] = sr_m.getId()

    for r_id in recon_r_ids_to_add:
        r = r_model.getReaction(r_id)
        r_id2st, p_id2st = {}, {}
        for r_m_id, st in get_reactants(r, stoichiometry=True):
            if r_m_id not in m_id_recon2sr:
                raise ValueError("Could not map %s (%s) to any metabolite in SR"
                                 % (r_model.getSpecies(r_m_id).getName(), r_m_id))
            r_id2st[m_id_recon2sr[r_m_id]] = st
        for r_m_id, st in get_products(r, stoichiometry=True):
            if r_m_id not in m_id_recon2sr:
                raise ValueError("Could not map %s (%s) to any metabolite in SR"
                                 % (r_model.getSpecies(r_m_id).getName(), r_m_id))
            p_id2st[m_id_recon2sr[r_m_id]] = st
        kegg_id = get_kegg_r_id(r)
        new_r = create_reaction(sr_model, r_id2st, p_id2st, r.getName(), reversible=r.getReversible(),
                                id_=r.getId() if not kegg_id
                                else (kegg_id.upper() + "_".join(get_r_comps(r_id, r_model))))
        if kegg_id:
            add_annotation(new_r, libsbml.BQB_IS, kegg_id.upper(), KEGG_REACTION_PREFIX)


def make_reversible(r):
    l_b, u_b = get_bounds(r)
    if u_b == 0:
        u_b = -l_b if l_b else 1000
    if l_b == 0:
        l_b = -u_b if u_b else -1000
    set_bounds(r, l_b, u_b)
    r.setReversible(True)


if "__main__" == __name__:
    sys.exit(main())
