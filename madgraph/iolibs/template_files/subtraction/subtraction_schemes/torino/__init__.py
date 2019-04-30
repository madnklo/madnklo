###########################################################
#
# torino subtraction scheme
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

# Import the sectors defining function
from sectors import SectorGenerator

# Import currents from colorful for now
import subtraction_schemes.colorful.NLO.local_currents as colorful_NLO_local_currents
import subtraction_schemes.colorful.NLO.integrated_currents as colorful_NLO_integrated_currents

__authors__ = ["valentin.hirschi@gmail.com", "simone.lionetti@gmail.com", "nicoldeu@phys.ethz.ch"]

# Mandatory variables of each subtraction module

# General properties
# ==================

# In the torino implementation, softs do not recoil symmetrically against initial states
requires_soft_beam_factorization = False 
# Colorful does not use sectors
sector_generator = SectorGenerator()
# This scheme belongs to a family of scheme using leg number information to instantiates its currents
are_current_instances_for_specific_leg_numbers = True

# Note: specifying below which resources are needed is optional
exporter = GenericCurrentsExporter(relative_resource_paths=[
    'subtraction_schemes/colorful',
    'subtraction_schemes/torino'
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
    colorful_NLO_local_currents.QCD_final_collinear_0_qqx,
    colorful_NLO_local_currents.QCD_final_collinear_0_gq,
    colorful_NLO_local_currents.QCD_final_collinear_0_gg,
    # initial-final collinears
    #       This scheme does *not* support ISR and the DefaultCurrent
    #       implementation will be used for them with an appropriate warning.
    # soft and soft-collinears
    colorful_NLO_local_currents.QCD_soft_0_g,
    colorful_NLO_local_currents.QCD_final_softcollinear_0_gX,
])

# Add NLO integrated counterterms
# ===============================
all_subtraction_current_classes.extend([
    # final-final collinears
    colorful_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_qqx,
    colorful_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gq,
    colorful_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gg,
    # soft and soft-collinear
    colorful_NLO_integrated_currents.integrated_NLO_QCD_soft_gluon,
    colorful_NLO_integrated_currents.integrated_NLO_FF_QCD_softcollinear_gq
])
