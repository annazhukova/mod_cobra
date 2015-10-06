from collections import defaultdict
import logging

from cobra.flux_analysis import flux_variability_analysis

from mod_cobra.fbva.cobra_model_manager import get_boundary_reactions, get_transport_reactions
from mod_cobra import round_value
from mod_cobra.fbva.serialization import format_values, format_r_id, get_cobra_r_formula

__author__ = 'anna'


def print_fluxes_larger_than_threshold_by_pathway(model, pw2r_ids, o_r_ids, pw2name, threshold=0):
    for pw in sorted(pw2r_ids.iterkeys()):
        r_ids = pw2r_ids[pw]
        logging.info("\n---------------")
        logging.info("%s (%d out of %d)" % (
            pw2name(pw), len([r_id for r_id in r_ids if r_id in model.solution.x_dict and round_value(model.solution.x_dict[r_id])]), len(r_ids)))
        logging.info("---------------\n")
        print_fluxes_larger_than_threshold(model, r_ids=r_ids, threshold=threshold)
    if o_r_ids:
        i_r_ids = get_boundary_reactions(model) & o_r_ids
        t_r_ids = get_transport_reactions(model) & o_r_ids
        ti_r_ids = i_r_ids | t_r_ids
        logging.info("\n---------------")
        logging.info("No pathway assigned: input and transport (%d out of %d)" %
              (len([r_id for r_id in ti_r_ids if r_id in model.solution.x_dict and round_value(model.solution.x_dict[r_id])]), len(ti_r_ids)))
        logging.info("---------------\n")
        print_fluxes_larger_than_threshold(model, r_ids=ti_r_ids, threshold=threshold)
        o_r_ids = o_r_ids - ti_r_ids
        logging.info("\n---------------")
        logging.info("No pathway assigned: other (%d out of %d)" %
              (len([r_id for r_id in o_r_ids if r_id in model.solution.x_dict and round_value(model.solution.x_dict[r_id])]), len(o_r_ids)))
        logging.info("---------------\n")
        print_fluxes_larger_than_threshold(model, r_ids=o_r_ids, threshold=threshold)


def print_fluxes_larger_than_threshold(model, threshold=0, r_ids=None, highlighted_r_ids=None):
    r_id2val = {}
    value2rn = defaultdict(list)
    for r_id, value in model.solution.x_dict.iteritems():
        value = round_value(value)
        if threshold is None or abs(value) > threshold:
            rn = model.reactions.get_by_id(r_id)
            # value2rn[value].append(format_r_name(rn) + ": " + rn.build_reaction_string(True))
            if not r_ids or r_id in r_ids:
                value2rn[value].append(rn.id + ": " + get_cobra_r_formula(rn, comp=False)
                                       + (" *" if highlighted_r_ids and rn.id in highlighted_r_ids else ''))
            r_id2val[r_id] = value

    prev_value = None
    for value in sorted(value2rn.iterkeys(), key=lambda v: abs(v)):
        for rn in sorted(value2rn[value]):
            if prev_value is not None and abs(prev_value) != abs(value):
                logging.info('')
            prev_value = value
            logging.info("%.4g\t%s" % (value, rn))

    return r_id2val


def serialize_fluxes(model, r_id2val, path, objective_sense, out_r_id):
    out_r_id = format_r_id(out_r_id)
    title = '%s reaction %s (%s): %.4g\n==================================\n'\
            % (objective_sense, out_r_id, get_cobra_r_formula(model.reactions.get_by_id(out_r_id), comp=True),
               r_id2val[out_r_id])
    value2rn = defaultdict(list)
    for r_id, value in r_id2val.iteritems():
        rn = model.reactions.get_by_id(r_id)
        value2rn[value].append(r_id + ": " + get_cobra_r_formula(rn, comp=True))

    prev_value = None
    with open(path, 'w+') as f:
        if title:
            f.write(title + '\n')
        for value in sorted(value2rn.iterkeys(), key=lambda v: abs(v)):
            for rn in sorted(value2rn[value]):
                if prev_value is not None and abs(prev_value) != abs(value):
                    f.write('\n')
                prev_value = value
                f.write("%.4g\t%s\n" % (value, rn))


def print_fva(model, rs=None, threshold=None, r_ids=None):
    r_id2min_max = flux_variability_analysis(model, reaction_list=rs)
    values2r_ids = defaultdict(set)
    r_id2bounds = {}
    for r_id, values in r_id2min_max.iteritems():
        min_v, max_v = round_value(values['minimum']), round_value(values['maximum'])
        if threshold is None or abs(min_v) > threshold or abs(max_v) > threshold:
            values2r_ids[(min_v, max_v)].add(r_id)
            r_id2bounds[r_id] = min_v, max_v
    keys = sorted(values2r_ids.iterkeys(), key=lambda (min_v, max_v): (-abs(max_v - min_v), min_v, max_v))
    ess_count, var_count = 0, 0
    logging.info("==============================\nESSENTIAL FLUXES\n==============================")
    for (min_v, max_v) in keys:
        if min_v * max_v > 0:
            v_r_ids = values2r_ids[(min_v, max_v)]
            if r_ids:
                v_r_ids = set(v_r_ids) & r_ids
            value = format_values(min_v, max_v)
            ess_count += len(v_r_ids)
            for r_id in sorted(v_r_ids):
                logging.info("%s\t%s: %s" % (value, r_id, model.reactions.get_by_id(r_id).build_reaction_string(True)))
            logging.info('')
    logging.info("==============================\nOTHER FLUXES\n==============================")
    for (min_v, max_v) in keys:
        if min_v * max_v == 0:
            v_r_ids = values2r_ids[(min_v, max_v)]
            if r_ids:
                v_r_ids = set(v_r_ids) & r_ids
            value = format_values(min_v, max_v)
            var_count += len(v_r_ids)
            for r_id in sorted(v_r_ids):
                logging.info("%s\t%s: %s" % (value, r_id, model.reactions.get_by_id(r_id).build_reaction_string(True)))
            logging.info('')
    logging.info("==============================")
    logging.info("In total, found %d essential reactions and %d variable reactions" % (ess_count, var_count))
    return r_id2bounds


def serialize_fva(model, r_id2bounds, path, objective_sense, out_r_id):
    out_r_id = format_r_id(out_r_id)
    title = '%s reaction %s (%s): %.4g\n==================================\n'\
              % (objective_sense, out_r_id, get_cobra_r_formula(model.reactions.get_by_id(out_r_id),
                                                                comp=True),
                 r_id2bounds[out_r_id][0])

    values2r_ids = defaultdict(set)
    for r_id, (min_v, max_v) in r_id2bounds.iteritems():
        values2r_ids[(min_v, max_v)].add(r_id)
    keys = sorted(values2r_ids.iterkeys(), key=lambda (min_v, max_v): (-abs(max_v - min_v), min_v, max_v))
    ess_count, var_count = 0, 0
    with open(path, 'w+') as f:
        if title:
            f.write(title + '\n')
        f.write("==============================\nESSENTIAL FLUXES\n==============================\n")
        for (min_v, max_v) in keys:
            if min_v * max_v > 0:
                v_r_ids = values2r_ids[(min_v, max_v)]
                value = format_values(min_v, max_v)
                ess_count += len(v_r_ids)
                for r_id in sorted(v_r_ids):
                    f.write("%s\t%s: %s\n" % (value, r_id,
                                              get_cobra_r_formula(model.reactions.get_by_id(r_id), comp=True)))
                f.write('\n')
        f.write("==============================\nOTHER FLUXES\n==============================\n")
        for (min_v, max_v) in keys:
            if min_v * max_v <= 0:
                v_r_ids = values2r_ids[(min_v, max_v)]
                value = format_values(min_v, max_v)
                var_count += len(v_r_ids)
                for r_id in sorted(v_r_ids):
                    f.write("%s\t%s: %s\n" % (value, r_id,
                                              get_cobra_r_formula(model.reactions.get_by_id(r_id), comp=True)))
                f.write('\n')
        f.write("==============================\n")
        f.write("In total, found %d essential reactions and %d variable reactions" % (ess_count, var_count))
    return ess_count, var_count
