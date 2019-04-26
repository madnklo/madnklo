###########################################################
#
# colorful_ subtraction scheme
#
###########################################################

import os
import sys
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0,root_path)

from commons.currents_exporter import GenericCurrentsExporter

# From common resources
import commons.beam_factorization_BF as BF
import commons.beam_factorization_BS as BS

# Specific to the colorful scheme
import NLO.local_currents as NLO_local_currents
import NLO.integrated_currents as NLO_integrated_currents

__authors__ = ["valentin.hirschi@gmail.com", "simone.lionetti@gmail.com", "nicoldeu@phys.ethz.ch"]

# Mandatory variables of each subtraction module

# General properties
# ==================

requires_soft_beam_factorization = False
# Colorful does not use sectors
sector_generator = None

# Note: specifying below which resources are needed is optional
exporter = GenericCurrentsExporter(relative_resource_paths=[
    'subtraction_schemes/colorful',
    'subtraction_schemes/colorful_pp'    
])

all_subtraction_current_classes = []

# Add NLO beam factorization counterterms (BF terms)
# ==================================================
all_subtraction_current_classes.extend([
    BF.QCD_beam_factorization_F0,
    BF.QCD_beam_factorization_single_collinear,
    BF.QCD_beam_factorization_single_softcollinear
])

# Add local NLO counterterms
# ==========================
all_subtraction_current_classes.extend([
    # final-final collinears
    NLO_local_currents.QCD_final_collinear_0_qqx,
    NLO_local_currents.QCD_final_collinear_0_gq,
    NLO_local_currents.QCD_final_collinear_0_gg,
    # initial-final collinears
    #       This scheme does *not* support ISR and the DefaultCurrent
    #       implementation will be used for them with an appropriate warning.
    # soft and soft-collinears
    NLO_local_currents.QCD_soft_0_g,
    NLO_local_currents.QCD_final_softcollinear_0_gX,    
])

# Add NLO integrated counterterms
# ===============================
all_subtraction_current_classes.extend([
    # final-final collinears
    NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_qqx,
    NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gq,
    NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gg,
    # soft and soft-collinear
    NLO_integrated_currents.integrated_NLO_QCD_soft_gluon,
    NLO_integrated_currents.integrated_NLO_FF_QCD_softcollinear_gq
])
