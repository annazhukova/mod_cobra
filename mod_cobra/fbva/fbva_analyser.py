import logging
from mod_cobra.fbva import MAXIMIZE

from mod_cobra.fbva.serialization import format_r_id, format_r_name

from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

from mod_cobra import round_value

__author__ = 'anna'


def analyse_by_fba(cobra_model, bm_r_id, objective_sense=MAXIMIZE, threshold=0):
    cobra_bm_r_id = format_r_id(bm_r_id)
    opt_value = optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    if opt_value:
        opt_value = round_value(opt_value)
        r_id2val = get_fluxes_larger_than_threshold(cobra_model, threshold=threshold)

        return r_id2val, opt_value
    return {}, None


def analyse_by_fva(cobra_model, bm_r_id, objective_sense=MAXIMIZE, threshold=0):
    cobra_bm_r_id = format_r_id(bm_r_id)
    opt_value = optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    if opt_value:
        opt_value = round_value(opt_value)
        r_id2bounds = get_r_id2fva_bounds(cobra_model, threshold=threshold)
        return r_id2bounds, opt_value
    return {}, None


def optimise_biomass(model, bm_id, objective_sense=MAXIMIZE, level=logging.INFO):
    reaction = model.reactions.get_by_id(bm_id)
    model.change_objective([reaction])
    try:
        optimize_minimal_flux(model, objective_sense=objective_sense)
    except Exception as e:
        logging.error("Could not optimise the model %s: %s" % (model.id, e.message))
        return None
    if "optimal" != model.solution.status:
        logging.error("A problem occurred while calculating the solution for %s: %s"
                      % (format_r_name(reaction), model.solution.status))
        return None
    else:
        logging.log(level, "%s: %.4g" % (format_r_name(reaction), model.solution.f))
        return model.solution.f


def get_fluxes_larger_than_threshold(model, threshold=0):
    r_id2val = {}
    for r_id, value in model.solution.x_dict.items():
        value = round_value(value)
        if threshold is None or abs(value) > threshold:
            r_id2val[r_id] = value
    return r_id2val


def get_r_id2fva_bounds(model, threshold=None, rs=None):
    r_id2min_max = flux_variability_analysis(model, reaction_list=rs)
    r_id2bounds = {}
    for r_id, values in r_id2min_max.items():
        min_v, max_v = round_value(values['minimum']), round_value(values['maximum'])
        if threshold is None or abs(min_v) > threshold or abs(max_v) > threshold:
            r_id2bounds[r_id] = min_v, max_v
    return r_id2bounds
