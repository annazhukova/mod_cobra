import logging
import os

from mod_cobra.analysis_pipeline import multimodel_pipeline
from models import SR_MODEL_SER, SR_GLN_EXCHANGE, SR_SER_PRODUCTION, \
    RECON_MODEL_BOUNDARY, RECON_GLN_EXCHANGE, RECON_SER_PRODUCTION, \
    UHLEN_MODEL, UH_GLN_IN, UH_SER_PRODUCTION, \
    DATA_DIR, RECON_GLC_EXCHANGE, SR_GLC_EXCHANGE, UH_GLC_IN

sbml2parameter_gln_ser = {SR_MODEL_SER: ({SR_SER_PRODUCTION: False, SR_GLN_EXCHANGE: False}, {}),
                          UHLEN_MODEL: ({UH_SER_PRODUCTION: False, UH_GLN_IN: False}, {}),
                          RECON_MODEL_BOUNDARY: ({RECON_SER_PRODUCTION: False, RECON_GLN_EXCHANGE: True}, {})}

sbml2parameter_no_glc_no_o2 = {SR_MODEL_SER: ({SR_SER_PRODUCTION: False}, {SR_GLC_EXCHANGE: False, 'Boundary3': False}),
                               UHLEN_MODEL: ({UH_SER_PRODUCTION: False}, {UH_GLC_IN: False, 'R_HMR_9048': False,
                                                                          'R_HMR_9383': False}),
                               RECON_MODEL_BOUNDARY: ({RECON_SER_PRODUCTION: False},
                                                      {RECON_GLC_EXCHANGE: True, 'R_EX_o2_LPAREN_e_RPAREN_': True})}

TREEEFM_PATH = "/home/azhukova/Applications/TreeEFM/TreeEFMseq"

__author__ = 'anna'


if "__main__" == __name__:
    # the following block is optional, it tells the system to log the analysis process to the console
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S")
    # end of the block

    multimodel_pipeline(sbml2parameter_gln_ser, os.path.join(DATA_DIR, 'glu_ser'), treeefm_path=TREEEFM_PATH,
                        rewrite=False, max_efm_number=100)

