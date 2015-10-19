import logging
import os
import sys

from models import SR_MODEL_SER, SR_SER_B, SR_GLN_B, SR_GLN_EXCHANGE, SR_SER_PRODUCTION, \
    RECON_MODEL_BOUNDARY, RECON_SER_B, RECON_GLN_B, RECON_GLN_EXCHANGE, RECON_SER_PRODUCTION, \
    UHLEN_MODEL, UH_SER_B, UH_GLN_B, UH_GLN_IN, UH_SER_PRODUCTION, \
    DATA_DIR
from mod_cobra.pipeline import multimodel_pipeline

__author__ = 'anna'


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S")
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S")

    multimodel_pipeline({
        SR_MODEL_SER: (SR_SER_PRODUCTION, False, {SR_GLN_EXCHANGE: False}, SR_GLN_B, SR_SER_B, 'SR'),
        # RR_MODEL: (RR_SER_PRODUCTION, False, {RR_GLN_EXCHANGE: False}, RR_GLN_B, RR_SER_B, 'RR'),
        # UHLEN_MODEL: (UH_SER_PRODUCTION, False, {UH_GLN_IN: False}, UH_GLN_B, UH_SER_B, 'Uhlen'),
        RECON_MODEL_BOUNDARY: (RECON_SER_PRODUCTION, False, {RECON_GLN_EXCHANGE: True},
                               RECON_GLN_B, RECON_SER_B, 'Recon2')
    }, os.path.join(DATA_DIR, 'Human_100'), do_fba=False, do_fva=False, do_efm=True, max_efm_number=100)


if "__main__" == __name__:
    sys.exit(main())
