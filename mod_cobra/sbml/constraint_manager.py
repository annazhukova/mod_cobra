import logging

from mod_sbml.sbml.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import reverse_reaction, get_reactants, get_products

__author__ = 'anna'


def constraint_exchange_reactions(model, forsed_r_id2rev, prohibited_r_id2rev=None, cofactors=None, min_flux=0.01):
    logging.info("Constraining input reactions...")

    for r in model.getListOfReactions():
        r_id = r.getId()
        r_l, r_u = get_bounds(r)

        if forsed_r_id2rev and r_id in forsed_r_id2rev:
            rev = forsed_r_id2rev[r_id]
            if rev:
                # reverse the reaction and set positive bounds,
                # as if both bounds are negative, the glp solver produces an error
                reverse_reaction(r)
                set_bounds(r, max(min_flux, -r_u), max(min_flux, -r_l))
                forsed_r_id2rev[r_id] = not rev
            else:
                set_bounds(r, max(min_flux, r_l), max(min_flux, r_u))
            r.setReversible(False)
            continue

        if prohibited_r_id2rev and r_id in prohibited_r_id2rev:
            rev = prohibited_r_id2rev[r_id]
            if not rev:
                reverse_reaction(r)
                set_bounds(r, 0, max(0, -r_l))
                prohibited_r_id2rev[r_id] = not rev
            else:
                set_bounds(r, 0, max(0, r_u))
            r.setReversible(False)
            continue

        if not cofactors:
            continue

        rs, ps = set(get_reactants(r)), set(get_products(r))

        boundary_s_ids = {s_id for s_id in rs if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not rs:
            if ps - cofactors:
                reverse_reaction(r)
                set_bounds(r, 0, max(-r_l, 0))
                r.setReversible(False)
            continue
        boundary_s_ids = {s_id for s_id in ps if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not ps:
            if rs - cofactors:
                set_bounds(r, 0, max(r_u, 0))
                r.setReversible(False)


def get_exchange_reactions(model):
    result = []
    for r in model.getListOfReactions():
        rs, ps = set(get_reactants(r)), set(get_products(r))

        boundary_s_ids = {s_id for s_id in rs if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not rs:
            result.append(r.getId())
            continue
        boundary_s_ids = {s_id for s_id in ps if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not ps:
            result.append(r.getId())
    return result
