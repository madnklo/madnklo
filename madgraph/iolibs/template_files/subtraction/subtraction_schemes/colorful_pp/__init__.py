###########################################################
#
# colorful_pp subtraction scheme
#
###########################################################

__authors__ = ["valentin.hirschi@gmail.com", "simone.lionetti@gmail.com", "nicoldeu@phys.ethz.ch"]

# Mandatory variables of each subtraction module

# General properties
# ==================

# In the colorful_pp implementation softs do actually recoil against initial states
requires_soft_beam_factorization = True

# This scheme belongs to a family of scheme using leg number information to instantiates its currents
are_current_instances_for_specific_leg_numbers = True

# loaded attributes
loaded_attributes = {}

def load():

    import os
    import sys
    root_path = os.path.dirname(os.path.realpath( __file__ ))
    sys.path.insert(0,root_path)
    sys.path.insert(0,os.path.abspath(os.path.join(root_path, os.path.pardir, os.path.pardir)))

    from commons.currents_exporter import GenericCurrentsExporter
    import commons.factors_and_cuts as factors_and_cuts

    # From common resources
    import commons.beam_factorization_BF as BF
    import commons.beam_factorization_BS as BS

    # Specific to the colorful_pp scheme
    import colorful_pp_config
    import NLO.local_currents as NLO_local_currents
    import NLO.integrated_currents as NLO_integrated_currents
    import NNLO.local_currents as NNLO_local_currents

    # Imports from the colorful scheme
    import subtraction_schemes.colorful.NLO.local_currents as colorful_NLO_local_currents
    import subtraction_schemes.colorful.NLO.integrated_currents as colorful_NLO_integrated_currents

    # Colorful_pp does not use sectors
    loaded_attributes['sector_generator'] = None

    # Note: specifying below which resources are needed is optional
    loaded_attributes['exporter'] = GenericCurrentsExporter(relative_resource_paths=[
        'subtraction_schemes/colorful',
        'subtraction_schemes/colorful_pp'
    ])

    all_subtraction_current_classes = []

    ###########
    # NLO
    ###########

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
        # initial-final collinears
        NLO_local_currents.QCD_initial_collinear_0_qg,
        NLO_local_currents.QCD_initial_collinear_0_gq,
        NLO_local_currents.QCD_initial_collinear_0_qq,
        NLO_local_currents.QCD_initial_collinear_0_gg,
        # soft and soft-collinears
        NLO_local_currents.QCD_soft_0_g,
        NLO_local_currents.QCD_final_softcollinear_0_gX,
        NLO_local_currents.QCD_initial_softcollinear_0_Xg,
    ])

    NLO_final_collinears = [
        NLO_local_currents.QCD_final_collinear_0_qqx,
        NLO_local_currents.QCD_final_collinear_0_gq,
        NLO_local_currents.QCD_final_collinear_0_gg,
    ]
    # final-final collinears
    all_subtraction_current_classes.extend(NLO_final_collinears)

    # Add NLO integrated counterterms
    # ===============================

    integrated_final_collinears_from_colorful = [
        colorful_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_qqx,
        colorful_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gq,
        colorful_NLO_integrated_currents.integrated_NLO_FF_QCD_collinear_gg,
    ]
    for integrated_final_collinear in integrated_final_collinears_from_colorful:
        # We must make sure that the alpha_0 virtuality bound is now computed dynamically so as to achieve
        # mapping independence, hence making our LorentzMapping for final final collinear consistent with these integrated CTs.
        integrated_final_collinear.get_alpha_virtuality_upper_bound = \
                                                staticmethod(NLO_integrated_currents.dynamic_alpha_virtuality_upper_bound)

    # final-final collinears
    all_subtraction_current_classes.extend(integrated_final_collinears_from_colorful)
    all_subtraction_current_classes.extend([
        # soft and soft-collinear
        NLO_integrated_currents.integrated_NLO_QCD_soft_gluon,
        NLO_integrated_currents.integrated_NLO_FF_QCD_softcollinear_gq
    ])

    ###########
    # NNLO
    ###########

    # For now we are only trying an elementary IFF q > q q' q' collinear
    all_subtraction_current_classes.extend([
        # IF
        NNLO_local_currents.QCD_initial_collinear_0_qqpqp,
        # FF
        NNLO_local_currents.QCD_final_collinear_0_qqpqp,
        # S(FF)
        NNLO_local_currents.QCD_soft_0_qqp,
        # IS(FF)
        NNLO_local_currents.QCD_initial_soft_collinear_0_qqpqp,
        # C(I,C(FF))
        NNLO_local_currents.QCD_C_FqFqx_C_IqpFqFqx,
        # S(C(FF))
        NNLO_local_currents.QCD_S_FqFqx_C_FqFqx,
        # C(S(C(FF)),I)
        NNLO_local_currents.QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx
    ])

    # Finally register the subtraction current classes loaded
    loaded_attributes['all_subtraction_current_classes'] = all_subtraction_current_classes
