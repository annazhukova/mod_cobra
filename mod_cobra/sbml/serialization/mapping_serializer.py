import logging
import os
from mod_cobra.html import describe
from mod_sbml.serialization.csv_manager import serialize_common_elements_to_csv

__author__ = 'anna'


def serialize(model_id2dfs, model_id2c_id_groups, model_id2m_id_groups, model_id2r_id_groups, path, get_f_path):
    comp_csv, m_csv, r_csv = serialize_common_elements_to_csv(model_id2dfs, model_id2c_id_groups,
                                                              model_id2m_id_groups, model_id2r_id_groups,
                                                              os.path.join(path, 'Model_comparison_')) \
        if model_id2c_id_groups else (None, None, None)
    logging.info('Serialized the mappings.')

    return describe('model_comparison.html', 
                    c_num=len(model_id2c_id_groups), m_num=len(model_id2m_id_groups), r_num=len(model_id2r_id_groups), 
                    c_csv=get_f_path(comp_csv), m_csv=get_f_path(m_csv), r_csv=get_f_path(r_csv))
