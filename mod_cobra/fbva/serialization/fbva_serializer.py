from collections import defaultdict
from mod_cobra.fbva.serialization import format_values, format_r_id, get_cobra_r_formula

__author__ = 'anna'


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
