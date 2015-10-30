import os
import mod_cobra

__author__ = 'anna'

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(mod_cobra.__file__)), '..', 'examples', 'data')

# -----------------Uhlen----------------------------------------------
UHLEN_MODEL = os.path.join(DATA_DIR, "MODEL1411240029.xml")

UH_SER_PRODUCTION = 'R_HMR_3843'

UH_GLN_IN = 'R_HMR_9063'
UH_GLC_IN = 'R_HMR_9034'
UH_SER_IN = 'R_HMR_9069'

UH_SER_B = 'M_m02896x'
UH_GLN_B = 'M_m01975x'


# -----------------Smith-Robinson--------------------
SR_MODEL = os.path.join(DATA_DIR, "MODEL1106160000.xml")
SR_MODEL_ANNOTATED = os.path.join(DATA_DIR, "SR.annotated.xml")
SR_MODEL_SER = os.path.join(DATA_DIR, "SR.annotated.serine.xml")

# L-Glutamine [Boundary] <=> L-Glutamine [Cytosol]
SR_GLN_EXCHANGE = 'Boundary18'

# alpha-D-Glucose [Boundary] <=> alpha-D-Glucose [Cytosol]
SR_GLC_EXCHANGE = 'Boundary24'

# L-Serine [Boundary] <=> L-Serine [Cytosol]
SR_SER_EXCHANGE = 'Boundary72'

# Oxaloacetate(C00036Cyto) + gtp(C00044Cytosol) <=>
# CO2(C00011Cyto) + Phosphoenolpyruvate(C00074Cyto) + GDP(C00035Cyto)
SR_SER_PRODUCTION = 'R_PSP_L'

SR_SER_B = 'C00065_b'
SR_GLN_B = 'C00064_b'

# --------------------RR--------------------
RR_MODEL = os.path.join(DATA_DIR, 'RR13-2.xml')

RR_GLC_EXCHANGE = 'GLUCc_exchange'
RR_GLN_EXCHANGE = 'GLNCUP'

RR_SER_PRODUCTION = 'AA2'

RR_SER_B = 'SER'
RR_GLN_B = 'GLN'


# -----------------------Recon-2---------------------------------
RECON_MODEL = os.path.join(DATA_DIR, "MODEL1109130000.xml")
RECON_MODEL_BOUNDARY = os.path.join(DATA_DIR, "MODEL1109130000.boundary.xml")

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