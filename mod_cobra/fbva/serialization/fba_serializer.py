from collections import defaultdict
import os
from mod_cobra.html import describe

from mod_cobra.fbva.serialization import format_r_id, get_cobra_r_formula

__author__ = 'anna'


def serialize_fba(model, r_id2val, path, objective_sense, out_r_id):
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


def serialize(cobra_model, opt_val, r_id2val, objective_sense, out_r_id, path, get_f_path):
    fba_file = None
    if opt_val:
        fba_file = os.path.join(path, 'fba.txt')
        serialize_fba(cobra_model, r_id2val, path=fba_file, objective_sense=objective_sense, out_r_id=out_r_id)
    return describe('fba.html', optimal_value=opt_val, r_num=len(r_id2val), description_filepath=get_f_path(fba_file))
