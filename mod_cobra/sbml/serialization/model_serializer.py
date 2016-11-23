import logging
import os

from mod_cobra.html import describe
from mod_sbml.serialization import get_sbml_r_formula
from mod_sbml.serialization.csv_manager import serialize_model_info

__author__ = 'anna'


def serialize(sbml, model, m_name, r_id2rev, path, get_f_path):
    logging.info("Serializing model info...")

    info_prefix = os.path.join(path, 'model_info_')
    c_csv, m_csv, r_csv = serialize_model_info(model, info_prefix)

    r_string = lambda r, rev: '<b>%s</b>%s: %s' % (r.getId(), ' (reversed)' if rev else '',
                                                   get_sbml_r_formula(model, r, show_metabolite_ids=False))
    return describe('input_data.html', model_name=m_name, sbml_filepath=get_f_path(sbml),
                    c_csv=get_f_path(c_csv), m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv),
                    in_rn_len=len(r_id2rev),
                    in_rns=[r_string(model.getReaction(r_id), rev)
                            for (r_id, rev) in r_id2rev.items()] if r_id2rev else [])
