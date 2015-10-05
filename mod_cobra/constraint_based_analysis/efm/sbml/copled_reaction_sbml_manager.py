import libsbml

from mod_sbml.sbml.sbml_manager import create_reaction
from mod_sbml.sbml.submodel_manager import remove_unused_species

__author__ = 'anna'


def create_folded_sbml(S_coupled, S_no_duplicates, in_sbml, out_sbml):
    doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = doc.getModel()

    r_id2new_r_id = {}
    r_ids_to_remove = set()
    for mr_id, cl_id2c in S_no_duplicates.gr_id2r_id2c.iteritems():
        name = 'lumped reaction groups of type %s' % mr_id
        id_ = mr_id
        r_id2st, p_id2st = S_no_duplicates.st_matrix.get_inputs_outputs(mr_id)
        new_r_id = create_reaction(model, r_id2st, p_id2st, reversible=False, id_=id_, name=name).getId()
        for cl_id in cl_id2c.iterkeys():
            for r_id in S_coupled.gr_id2r_id2c[cl_id].iterkeys() if cl_id in S_coupled.gr_id2r_id2c else [cl_id]:
                r_id2new_r_id[r_id] = new_r_id
                r_ids_to_remove.add(r_id)

    for cl_id, r_id2c in S_coupled.gr_id2r_id2c.iteritems():
        # If already processed in the previous block, skip.
        if cl_id in S_no_duplicates.r_id2gr_id:
            continue

        name = 'lumped reaction group %s' % cl_id
        id_ = cl_id
        r_id2st, p_id2st = S_coupled.st_matrix.get_inputs_outputs(cl_id)
        new_r_id = create_reaction(model, r_id2st, p_id2st, reversible=False, id_=id_, name=name).getId()
        for r_id in r_id2c.iterkeys():
            r_id2new_r_id[r_id] = new_r_id
            r_ids_to_remove.add(r_id)

    for r_id in r_ids_to_remove:
        model.removeReaction(r_id)

    remove_unused_species(model)
    model.setId('%s_folded' % model.getId())
    model.setName('%s_folded' % model.getName())
    libsbml.SBMLWriter().writeSBMLToFile(doc, out_sbml)
    return r_id2new_r_id