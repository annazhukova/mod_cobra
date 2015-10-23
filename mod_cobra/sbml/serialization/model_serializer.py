import logging
import os
import libsbml
from mod_cobra.html import describe
from mod_sbml.serialization.csv_manager import serialize_model_info
from mod_sbml.serialization import get_sbml_r_formula

__author__ = 'anna'


def serialize(sbml, out_r_id, out_rev, in_r_id2rev, path, get_f_path):
    logging.info("Serializing model info...")

    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    info_prefix = os.path.join(path, 'model_info_')
    c_csv, m_csv, r_csv = serialize_model_info(model, info_prefix)

    r_string = lambda r, rev: '<b>%s</b>%s: %s' % (r.getId(), ' (reversed)' if rev else '',
                                                   get_sbml_r_formula(model, r, show_metabolite_ids=False))
    r = model.getReaction(out_r_id)
    return describe('input_data.html', model_name=(model.getName() if model.getName() else model.getId()),
                    sbml_filepath=get_f_path(sbml), obj_rn=r_string(r, out_rev),
                    c_csv=get_f_path(c_csv), m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv),
                    in_rn_len=len(in_r_id2rev) if in_r_id2rev else 0,
                    in_rns='; '.join(r_string(model.getReaction(r_id), rev)
                                     for (r_id, rev) in in_r_id2rev.iteritems()) if in_r_id2rev else '')
