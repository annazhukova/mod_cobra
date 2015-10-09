import logging
import os
import sys

import mod_cobra
from mod_cobra.pipeline import multimodel_pipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(mod_cobra.__file__)), '..', 'examples', 'data')

# -----------------Uhlen----------------------------------------------
UHLEN_MODEL = os.path.join(DATA_DIR, "Uhlen_heart_MODEL1411240029.xml")

UH_SER_PRODUCTION = 'R_HMR_3843'

UH_GLN_IN = 'R_HMR_9063'
UH_SER_IN = 'R_HMR_9069'

UH_SER_B = 'M_m02896x'
UH_GLN_B = 'M_m01975x'


# -----------------Smith-Robinson--------------------
SR_MODEL_SER = os.path.join(DATA_DIR, "SR_serine.xml")

# L-Glutamine [Boundary] <=> L-Glutamine [Cytosol]
SR_GLN_EXCHANGE = 'Boundary18'

# alpha-D-Glucose [Boundary] <=> alpha-D-Glucose [Cytosol]
SR_GLC_EXCHANGE = 'Boundary24'

# L-Serine [Boundary] <=> L-Serine [Cytosol]
SR_SER_EXCHANGE = 'Boundary72'

# Oxaloacetate(C00036Cyto) + gtp(C00044Cytosol) <=>
# CO2(C00011Cyto) + Phosphoenolpyruvate(C00074Cyto) + GDP(C00035Cyto)
SR_SER_PRODUCTION = 'R_PSP_Lc'

SR_SER_B = 'C00065_b'
SR_GLN_B = 'C00064_b'

# --------------------RR--------------------
RR_MODEL = os.path.join(DATA_DIR, 'RR12.xml')

RR_GLC_EXCHANGE = 'GLUCc_exchange'
RR_GLN_EXCHANGE = 'GLNc_exchange'

RR_SER_PRODUCTION = 'AA2'

RR_SER_B = 'SERc_b'
RR_GLN_B = 'GLNc_b'


# -----------------------Recon-2---------------------------------
RECON_MODEL = os.path.join(DATA_DIR, "recon2model.v02.bounds.xml")

# L-glutamine(M_gln_L_e) <=>
RECON_GLN_EXCHANGE = 'R_EX_gln_L_LPAREN_e_RPAREN_'

# D-glucose(M_glc_D_e) <=>
RECON_GLC_EXCHANGE = 'R_EX_glc_LPAREN_e_RPAREN_'

# L-serine(M_ser_L_e) <=>
RECON_SER_L_EXCHANGE = 'R_EX_ser_L_LPAREN_e_RPAREN_'

# D-serine(M_ser_D_e) <=>
RECON_SER_D_EXCHANGE = 'R_EX_ser_D_LPAREN_e_RPAREN_'

# H2O(M_h2o_c) + proton(M_h_c) + O-Phospho-L-serine(M_pser_L_c) <=> hydrogenphosphate(M_pi_c) + L-serine(M_ser_L_c)
RECON_SER_PRODUCTION = 'R_PSP_L'

RECON_GLN_B = 'M_gln_L_e_b'
RECON_SER_B = 'M_ser_L_e_b'


__author__ = 'anna'


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S")
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S")

    multimodel_pipeline({
        SR_MODEL_SER: (SR_SER_PRODUCTION, False, {SR_GLN_EXCHANGE: False}, SR_GLN_B, SR_SER_B, 'SR'),
        RR_MODEL: (RR_SER_PRODUCTION, False, {RR_GLN_EXCHANGE: False}, RR_GLN_B, RR_SER_B, 'RR'),
        UHLEN_MODEL: (UH_SER_PRODUCTION, False, {UH_GLN_IN: False}, UH_GLN_B, UH_SER_B, 'Uhlen'),
        RECON_MODEL: (RECON_SER_PRODUCTION, False, {RECON_GLN_EXCHANGE: True}, RECON_GLN_B, RECON_SER_B, 'Recon2')
    }, os.path.join(DATA_DIR, 'Human_500'), do_fba=True, do_fva=True, do_efm=True, max_efm_number=500)

if "__main__" == __name__:
    sys.exit(main())
