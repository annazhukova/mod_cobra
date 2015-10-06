from collections import defaultdict
import logging

from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

from mod_cobra import round_value
from mod_cobra.fbva.serialization import format_r_name

__author__ = 'anna'

NUM_STEPS = 20.0


def set_bounds(reaction, value, max_up_b, max_l_b):
    value = abs(value)
    reaction.lower_bound = max(0 - value, 0 - abs(max_l_b))
    reaction.upper_bound = min(value, abs(max_up_b))


def optimise_biomass(model, bm_id, objective_sense='maximize', level=logging.INFO):
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


def get_biomass_dependency(model, r_id, r_ids, bm_id, max_bound=None, remove_zeros=False, constrained=False,
                           minimized=False):
    reaction = model.reactions.get_by_id(r_id)
    r_up_b, r_l_b = reaction.upper_bound, reaction.lower_bound
    max_bound = abs(max_bound) if max_bound is not None else max(abs(r_up_b), abs(r_l_b))
    r_id2ys = defaultdict(list)
    step = max_bound / NUM_STEPS
    xs = []

    old_r_id2bounds = None
    if constrained:
        r_id2bounds = get_possible_reaction_bounds(model, bm_id, r_ids - {r_id, bm_id}, reaction, xs)
        old_r_id2bounds = constraint_reactions(model, r_id2bounds)

    model.change_objective([model.reactions.get_by_id(bm_id)])

    for x in [i * step for i in range(0, int(NUM_STEPS))]:
        # reaction.lower_bound = -x
        # reaction.upper_bound = -x
        set_bounds(reaction, x, abs(r_l_b), abs(r_up_b))
        model.optimize()
        if model.solution.f:
            xs.append(x)

            if minimized:
                r_id2values = flux_variability_analysis(model, the_reactions=list(set(r_ids) - {r_id, bm_id}))
                r_id_by_flux = [rn_id for (rn_id, val) in r_id2values.iteritems() if
                                rn_id not in {r_id, bm_id} and round_value(val['minimum']) != round_value(
                                    val['maximum'])]
                r_id2bounds = get_reaction_bounds(model, r_id_by_flux)
                while r_id_by_flux:
                    max_flux_r_id = max(r_id_by_flux, key=lambda rn_id: -abs(model.solution.x_dict[rn_id]))
                    min_max = r_id2values[max_flux_r_id]
                    min_v, max_v = min_max['minimum'], min_max['maximum']
                    max_flux_r = model.reactions.get_by_id(max_flux_r_id)
                    flux_value = 0 if min_v * max_v < 0 else min_v if abs(min_v) < abs(max_v) else max_v
                    max_flux_r.lower_bound = flux_value
                    max_flux_r.upper_bound = flux_value
                    model.optimize()
                    r_id2values = flux_variability_analysis(model, the_reactions=r_id_by_flux)
                    r_id_by_flux = [rn_id for (rn_id, val) in r_id2values.iteritems() if
                                    round_value(val['minimum']) != round_value(val['maximum'])]
                constraint_reactions(model, r_id2bounds)

            for r_id in r_ids:
                r_id2ys[r_id].append(round_value(model.solution.x_dict[r_id]))
    reaction.upper_bound = r_up_b
    reaction.lower_bound = r_l_b

    if constrained and old_r_id2bounds:
        constraint_reactions(model, old_r_id2bounds)

    if remove_zeros:
        return xs, {r_id: ys for (r_id, ys) in r_id2ys.iteritems() if ys and (len(set(ys)) > 1 or ys[0] != 0)}
    else:
        return xs, r_id2ys


def get_possible_reaction_bounds(model, bm_id, r_ids, r_to_vary, steps):
    model.change_objective([model.reactions.get_by_id(bm_id)])
    model.optimize()
    r_id2min_max = flux_variability_analysis(model, allow_loops=False)
    r_id2bounds = {r_id: (r_id2min_max[r_id]['minimum'], r_id2min_max[r_id]['maximum']) for r_id in r_ids}
    up_b, l_b = abs(r_to_vary.upper_bound), abs(r_to_vary.lower_bound)
    for x in steps:
        set_bounds(r_to_vary, x, up_b, l_b)
        model.optimize()
        variable_r_ids = list(r_id2bounds.keys())
        for r_id in variable_r_ids:
            v_min, v_max = r_id2bounds[r_id]
            r_id2min_max = flux_variability_analysis(model, allow_loops=False, the_reactions=variable_r_ids)[r_id]
            new_v_min, new_v_max = max(r_id2min_max['minimum'], v_min), min(r_id2min_max['maximum'], v_max)
            if new_v_min > new_v_max:
                del r_id2bounds[r_id]
            else:
                r_id2bounds[r_id] = new_v_min, new_v_max
    return r_id2bounds


def constraint_reactions(model, r_id2bounds):
    old_r_id2bounds = {}
    for r_id, (v_min, v_max) in r_id2bounds.iteritems():
        rn = model.reactions.get_by_id(r_id)
        old_r_id2bounds[r_id] = (rn.lower_bound, rn.upper_bound)
        rn.lower_bound = v_min
        rn.upper_bound = v_max
    return old_r_id2bounds


def get_reaction_bounds(model, r_ids):
    r_id2bounds = {}
    for r_id in r_ids:
        r = model.reactions.get_by_id(r_id)
        r_id2bounds[r_id] = r.lower_bound, r.upper_bound
    return r_id2bounds


def get_fluxes_larger_than_threshold(model, threshold=0):
    r_id2val = {}
    for r_id, value in model.solution.x_dict.iteritems():
        value = round_value(value)
        if threshold is None or abs(value) > threshold:
            r_id2val[r_id] = value
    return r_id2val


def get_r_id2fva_bounds(model, threshold=None, rs=None):
    r_id2min_max = flux_variability_analysis(model, reaction_list=rs)
    r_id2bounds = {}
    for r_id, values in r_id2min_max.iteritems():
        min_v, max_v = round_value(values['minimum']), round_value(values['maximum'])
        if threshold is None or abs(min_v) > threshold or abs(max_v) > threshold:
            r_id2bounds[r_id] = min_v, max_v
    return r_id2bounds


# The method find_blocked_reactions in COBRA fails as it tries to give unexpected parameters to the solver.
# This method is an analogue.
def get_blocked_reactions(cobra_model, the_reactions=None, tolerance_optimality=1e-9, open_exchanges=False, **kwargs):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for cobra_model, using flux variability
    analysis.
    """
    cobra_model = cobra_model.copy()
    if not the_reactions:
        the_reactions = cobra_model.reactions
    if open_exchanges:
        exchange_reactions = [x for x in cobra_model.reactions
                              if x.startswith('EX')]
        for the_reaction in exchange_reactions:
            if the_reaction.lower_bound >= 0:
                the_reaction.lower_bound = -1000
            if the_reaction.upper_bound >= 0:
                the_reaction.upper_bound = 1000
    flux_span_dict = flux_variability_analysis(cobra_model, the_reactions, **kwargs)
    return [k for k, v in flux_span_dict.items() if max(map(abs, v.values())) <= tolerance_optimality]


def constraint_reaction_of_interest(model, r_id, bound):
    reaction = model.reactions.get_by_id(r_id)
    up_b, l_b = abs(reaction.upper_bound), abs(reaction.lower_bound)
    set_bounds(reaction, bound, up_b, l_b)

