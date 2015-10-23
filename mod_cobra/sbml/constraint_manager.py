import logging
import os
import libsbml
from mod_sbml.sbml.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import reverse_reaction, get_reactants, get_products
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, select_metabolite_ids_by_term_ids

__author__ = 'anna'


def constraint_exchange_reactions(sbml, allowed_exchange_r_id2rev, path, min_flux=0.01):
    logging.info("Constraining input reactions...")
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    ub_ch_ids = get_ubiquitous_chebi_ids(add_common=True, add_cofactors=True)
    cofactors = select_metabolite_ids_by_term_ids(model, ub_ch_ids)

    for r in model.getListOfReactions():
        if r.id in allowed_exchange_r_id2rev:
            rev = allowed_exchange_r_id2rev[r.id]
            if rev:
                # reverse the reaction and set positive bounds,
                # as if both bounds are negative, the glp solver produces an error
                l_b, u_b = get_bounds(r)
                reverse_reaction(r)
                set_bounds(r, -u_b, -l_b)
                allowed_exchange_r_id2rev[r.id] = not rev
            r_l, r_u = get_bounds(r)
            set_bounds(r, max(min_flux, r_l), max(min_flux, r_u))
            continue
        rs, ps = set(get_reactants(r)), set(get_products(r))

        boundary_s_ids = {s_id for s_id in rs if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not rs:
            r_l, r_u = get_bounds(r)
            # if it's not only cofactors, constrain it
            set_bounds(r, min(r_l, 0), 0 if (ps - cofactors) else r_u)
            continue
        boundary_s_ids = {s_id for s_id in ps if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not ps:
            r_l, r_u = get_bounds(r)
            # if it's not only cofactors, constrain it
            set_bounds(r, 0 if (rs - cofactors) else r_l, max(r_u, 0))

    sbml = os.path.join(path, '%s.constrained.xml' % os.path.splitext(os.path.basename(sbml))[0])
    libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)

    return sbml
