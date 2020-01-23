###########################################################
#
# colorful_pp subtraction scheme
#
###########################################################

import os

__authors__ = ["valentin.hirschi@gmail.com", "nicoldeu@phys.ethz.ch"]

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
    subtraction_scheme_name = os.path.split(os.path.dirname(os.path.realpath(__file__)))[-1]

    from commons.currents_exporter import GenericCurrentsExporter

    class ColorfulPPExporter(GenericCurrentsExporter):

        def does_require_correlated_beam_convolution(self, singular_structure):
            """ The subtraction scheme module must also specify which type of singular structure yield integrated
            counterterm that require a correlated convolution with both beams.
            For the colorful_pp scheme this is necessary whenever a soft limit or a collinear limit of only final state
            partons is involved.
            """

            requires_correlated_beam_convolution = any(any(
                    (c.name() == 'S' or all(l.state == l.FINAL for l in c.legs))
                for c in sub_ss.decompose() ) for sub_ss in singular_structure.substructures)

            return requires_correlated_beam_convolution

    import factors_and_cuts as factors_and_cuts

    # Specific to the colorful_pp scheme
    import colorful_pp_config
    import NLO.local_currents as NLO_local_currents
    import NLO.integrated_currents as NLO_integrated_currents
    import NNLO.local_currents as NNLO_local_currents

    # Colorful_pp does not use sectors
    loaded_attributes['sector_generator'] = None

    # Note: specifying below which resources are needed is optional
    loaded_attributes['exporter'] = ColorfulPPExporter(relative_resource_paths=[
        'subtraction_schemes/%s'%subtraction_scheme_name
    ])

    all_subtraction_current_classes = []

    ###########
    # NLO
    ###########

    # Add NLO beam factorization counterterms (BF terms)
    # ==================================================
    all_subtraction_current_classes.extend([
        # PDF counterterms
        NLO_integrated_currents.QCD_F0,
        # collinear IF counterterms
        NLO_integrated_currents.QCD_integrated_C_IF,
        # collinear FF counterterms
        NLO_integrated_currents.QCD_integrated_C_FqFq,
        NLO_integrated_currents.QCD_integrated_C_FqFg,
        NLO_integrated_currents.QCD_integrated_C_FgFg,
        # soft F counterterms
        NLO_integrated_currents.QCD_integrated_S_Fg,
        # soft-collinear counterterms
        NLO_integrated_currents.QCD_integrated_CS_IF,
        NLO_integrated_currents.QCD_integrated_CS_FF,
    ])

    # Add NLO beam factorization counterterms of soft origin
    # recoiling symmetrically against the initial state (BS)
    # ======================================================
    all_subtraction_current_classes.extend([ ])

    # Add local NLO counterterms
    # ==========================
    all_subtraction_current_classes.extend([
        # final-final collinears
        NLO_local_currents.QCD_C_FqFqx,
        NLO_local_currents.QCD_C_FgFg,
        NLO_local_currents.QCD_C_FqFg,
        # initial-final collinears
        NLO_local_currents.QCD_C_IgFq,
        NLO_local_currents.QCD_C_IqFg,
        NLO_local_currents.QCD_C_IqFq,
        NLO_local_currents.QCD_C_IgFg,
        # soft and soft-collinears
        NLO_local_currents.QCD_S_g,
        NLO_local_currents.QCD_CS_FgFg,
        NLO_local_currents.QCD_CS_FgFq,
        NLO_local_currents.QCD_CS_IgFg,
        NLO_local_currents.QCD_CS_IqFg,
    ])

    # Add NLO integrated counterterms
    # ===============================

    all_subtraction_current_classes.extend([
        # There are None in this scheme since *all* local counterterms
        # recoil against the initial states.
    ])

    ###########
    # NNLO
    ###########

    # For now we are only trying an elementary IFF q > q q' q' collinear
    all_subtraction_current_classes.extend([
        # IF
        # FF
        # S(FF)
        # IS(FF)
        # C(I,C(FF))
        NNLO_local_currents.QCD_C_FqFqx_C_IqpFqFqx,
        # S(C(FF))
        NNLO_local_currents.QCD_S_FqFqx_C_FqFqx,
        # C(S(C(FF)),I)
        NNLO_local_currents.QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx
    ])

    # Finally register the subtraction current classes loaded
    loaded_attributes['all_subtraction_current_classes'] = all_subtraction_current_classes
