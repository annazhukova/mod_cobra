import logging

from mod_cobra.fbva.fbva_manager import optimise_biomass, get_fluxes_larger_than_threshold, round_value, \
    get_r_id2fva_bounds
from mod_cobra.fbva.serialization import format_r_id

__author__ = 'anna'


def analyse_by_fba(cobra_model, bm_r_id, objective_sense='maximize', threshold=0):
    cobra_bm_r_id = format_r_id(bm_r_id)
    opt_value = optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    if opt_value:
        opt_value = round_value(opt_value)
        r_id2val = get_fluxes_larger_than_threshold(cobra_model, threshold=threshold)

        return r_id2val, opt_value
    return {}, None


def analyse_by_fva(cobra_model, bm_r_id, objective_sense='maximize', threshold=0):
    cobra_bm_r_id = format_r_id(bm_r_id)
    opt_value = optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    if opt_value:
        opt_value = round_value(opt_value)
        r_id2bounds = get_r_id2fva_bounds(cobra_model, threshold=threshold)
        return r_id2bounds, opt_value
    return {}, None