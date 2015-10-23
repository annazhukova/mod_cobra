import libsbml

from mod_sbml.sbml.sbml_manager import create_reaction
from mod_sbml.sbml.submodel_manager import remove_unused_species

__author__ = 'anna'


def create_folded_sbml(S, in_sbml, out_sbml):
    doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = doc.getModel()

    r_id2new_r_id = {}
    initial_r_ids = {r_id for r_id in S.r_id2i.iterkeys() if r_id not in S.gr_id2r_id2c}
    r_ids_to_remove = [r.getId() for r in model.getListOfReactions() if r.getId() not in S.r_ids]

    def map2new_id(r_id, new_id):
        if r_id not in S.gr_id2r_id2c:
            r_id2new_r_id[r_id] = new_id
        else:
            for sub_r_id in S.gr_id2r_id2c[r_id].iterkeys():
                map2new_id(sub_r_id, new_id)

    for mr_id in S.r_ids - initial_r_ids:
        cl_id2c = S.gr_id2r_id2c[mr_id]
        name = 'lumped reaction %s' % mr_id
        id_ = mr_id
        r_id2st, p_id2st = S.st_matrix.get_inputs_outputs(mr_id)
        new_r_id = create_reaction(model, r_id2st, p_id2st, reversible=True, id_=id_, name=name).getId()
        for r_id in cl_id2c:
            map2new_id(r_id, new_r_id)

    for r_id in r_ids_to_remove:
        model.removeReaction(r_id)

    for r in model.getListOfReactions():
        for i in xrange(0, r.getNumModifiers()):
            r.removeModifier(0)

    remove_unused_species(model)
    model.setId('%s_folded' % model.getId())
    model.setName('%s_folded' % model.getName())
    libsbml.SBMLWriter().writeSBMLToFile(doc, out_sbml)
    return r_id2new_r_id
