from mod_sbml.sbml.reaction_boundary_manager import set_bounds
from mod_cobra.sbml import r_ids2sbml
from mod_cobra.fbva.serialization import format_r_id

__author__ = 'anna'


def create_fva_model(sbml, r_id2bounds, new_sbml):
    def r_updater(r):
        l_b, u_b = r_id2bounds[format_r_id(r.id)]
        set_bounds(r, l_b, max(u_b, 0))
        if l_b * u_b >= 0:
            r.setReversible(False)

    r_ids = set(r_id2bounds.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2bounds.iterkeys()}
    r_ids2sbml(r_ids, sbml, new_sbml, 'FVA', r_updater)
    return new_sbml