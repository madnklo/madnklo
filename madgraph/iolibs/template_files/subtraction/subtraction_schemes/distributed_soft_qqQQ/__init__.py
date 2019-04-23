###########################################################
#
# colorful_pp subtraction scheme
#
###########################################################

from commons.currents_exporter import GenericCurrentsExporter

# From common resources
import commons.beam_factorization_BF as BF
import commons.beam_factorization_BS as BS

# Imports from the colorful_pp scheme
import subtraction_schemes.colorful_pp.NLO.local_currents as colorful_pp_NLO_local_currents
import subtraction_schemes.colorful_pp.NLO.integrated_currents as colorful_pp_NLO_integrated_currents

# Specific to the colorful scheme
import NLO.local_currents as NLO_local_currents
import NLO.integrated_currents as NLO_integrated_currents

__authors__ = ["valentin.hirschi@gmail.com", "simone.lionetti@gmail.com", "nicoldeu@phys.ethz.ch"]

raise NotImplementedError
# Below everything is copy pasted from colorful, this is not yet the scheme we want


# Mandatory variables of each subtraction module

# General properties
# ==================

requires_soft_beam_factorization = False 

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

# Add NLO beam factorization counterterms of soft origin
# recoiling symmetrically against the initial state (BS)
# ======================================================
all_subtraction_current_classes.extend([
    BS.QCD_beam_factorization_single_soft,
])

# Add local NLO counterterms
# ==========================
all_subtraction_current_classes.extend([
    # final-final collinears
    colorful_pp_NLO_local_currents.QCD_final_collinear_0_qqx,
    colorful_pp_NLO_local_currents.QCD_final_collinear_0_gq,
    colorful_pp_NLO_local_currents.QCD_final_collinear_0_gg,
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
    colorful_pp_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_qqx,
    colorful_pp_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gq,
    colorful_pp_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gg,
    # soft and soft-collinear
    NLO_integrated_currents.integrated_NLO_QCD_soft_gluon,
    NLO_integrated_currents.integrated_NLO_FF_QCD_softcollinear_gq
])
