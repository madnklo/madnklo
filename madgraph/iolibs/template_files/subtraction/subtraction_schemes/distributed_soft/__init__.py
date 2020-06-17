###########################################################
#
# distributed softs subtraction scheme
#
###########################################################

import os
import madgraph.various.misc as misc

__authors__ = ["valentin.hirschi@gmail.com", "hainzc@student.ethz.ch"]

# Mandatory variables of each subtraction module

# General properties
# ==================

# Here there is no recoil against the initial states
requires_soft_beam_factorization = False

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
    from madgraph.core.subtraction import SubtractionLeg

    class DistributedSoftExporter(GenericCurrentsExporter):

        def does_require_correlated_beam_convolution(self, singular_structure):
            """ The subtraction scheme module must also specify which type of singular structure yield integrated
            counterterm that require a correlated convolution with both beams.
            For the colorful_pp scheme this is necessary whenever a soft limit or a collinear limit of only final state
            partons is involved.
            """

            return False

        def get_n_non_factorisable_double_sided_convolution(self, counterterm):
            """ Test if this counterterm involves a non-factorisable double-sided convolution which appears when one has
            a combination of a correlated convolution with single-sided ones (like for IF integrated collinear CT or PDF
            counterterms). For now it can only return 2 convolutions if non-factorisable; it may need to be extended for N^3LO."""

            return False

    # Specific to the distibuted softs scheme
    # NLO local
    import NLO.local_currents as NLO_local_currents
    # NLO integrated
    import NLO.integrated_currents as NLO_integrated_currents

    # Colorful_pp does not use sectors
    loaded_attributes['sector_generator'] = None

    # Note: specifying below which resources are needed is optional
    loaded_attributes['exporter'] = DistributedSoftExporter(relative_resource_paths=[
        'subtraction_schemes/%s'%subtraction_scheme_name
    ])

    all_subtraction_current_classes = []

    ###########
    # NLO
    ###########

    # Add NLO beam factorization counterterms (BF terms)
    # ==================================================
    all_subtraction_current_classes.extend([
    ])

    # Add NLO beam factorization counterterms of soft origin
    # recoiling symmetrically against the initial state (BS)
    # ======================================================
    all_subtraction_current_classes.extend([ ])

    # Add local NLO counterterms
    # ==========================
    all_subtraction_current_classes.extend([
        # final-final collinears
        NLO_local_currents.QCD_C_FqFg,
        NLO_local_currents.QCD_C_FqFqx,
        NLO_local_currents.QCD_C_FgFg,
        # initial-final collinears
        # soft and soft-collinears
        NLO_local_currents.QCD_S_g,
        NLO_local_currents.QCD_CS_FgFg,
        NLO_local_currents.QCD_CS_FgFq,
        #NLO_local_currents.QCD_CS_IgFg, # Are these inital state currents?
        #NLO_local_currents.QCD_CS_IqFg,
    ])

    # Add NLO integrated counterterms
    # ===============================

    all_subtraction_current_classes.extend([
        # TODO Place here integrated counterterms
        NLO_integrated_currents.QCD_integrated_C_FqFq,
        NLO_integrated_currents.QCD_integrated_C_FqFg,
        NLO_integrated_currents.QCD_integrated_C_FgFg,
        NLO_integrated_currents.QCD_integrated_CS_FF,
        NLO_integrated_currents.QCD_integrated_S_Fg,#Need soft term
    ])

    ###########
    # NNLO
    ###########

    ###########
    # NNNLO
    ###########


    # Finally register the subtraction current classes loaded
    loaded_attributes['all_subtraction_current_classes'] = all_subtraction_current_classes
