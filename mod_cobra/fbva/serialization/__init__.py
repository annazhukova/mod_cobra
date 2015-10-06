__author__ = 'anna'


def format_values(m, M):
    return "%g:\t%g" % (m, M) if m != M else "\tCONST\t%g" % m


def format_r_name(r):
    return "%s: %s" % (r.id, r.name) if r.name.find(r.id) == -1 else r.name


def get_cobra_r_formula(r, comp=True):
    if comp:
        return "{0} <=> {1}".format(
            " + ".join(
                ["%s%s[%s](%s)" % (('%g ' % -r.get_coefficient(m)) if r.get_coefficient(m) != -1 else '',
                                   m.name, m.compartment, m.id) for m in r.reactants]),
            " + ".join(["%s%s[%s](%s)" % (('%g ' % r.get_coefficient(m)) if r.get_coefficient(m) != 1 else '',
                                          m.name, m.compartment, m.id) for m in r.products]))
    else:
        return "{0} <=> {1}".format(
            " + ".join(
                ["%s%s(%s)" % (('%g ' % -r.get_coefficient(m)) if r.get_coefficient(m) != -1 else '',
                               m.name, m.id) for m in r.reactants]),
            " + ".join(["%s%s(%s)" % (('%g ' % r.get_coefficient(m)) if r.get_coefficient(m) != 1 else '',
                                      m.name, m.id) for m in r.products]))


def format_r_id(r_id, remove_prefix=True):
    if remove_prefix and r_id.startswith('R_'):
        r_id = r_id[2:]
    if not remove_prefix and not r_id.startswith('R_'):
        return 'R_' + r_id
    return r_id