###########################################################
#
# Torino subtraction scheme
#
###########################################################


__authors__ = ["valentin.hirschi@gmail.com", "nicoldeu@phys.ethz.ch"]

# Mandatory variables of each subtraction module

# General properties
# ==================

# In the Torino implementation softs do not recoil against initial states
requires_soft_beam_factorization = True

#gl
torino_sub_BS = True

# This scheme belongs to a family of scheme using leg number information to instantiates its currents
are_current_instances_for_specific_leg_numbers = True

# loaded attributes
loaded_attributes = {}

def load():

    import os
    import sys
#Ez
    import logging

    logger1 = logging.getLogger('madgraph')

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
            #print('__initi__.py - requires_correlated_beam_convolution : ' + str(requires_correlated_beam_convolution))    

            return requires_correlated_beam_convolution

#gl
        def does_require_torino_sub_BS(self, singular_structure):

            requires_torino_sub_BS = torino_sub_BS
            #print('__initi__.py - requires_torino_sub_BS : ' + str(requires_torino_sub_BS))
            return requires_torino_sub_BS


    import factors_and_cuts as factors_and_cuts

    # From common resources
    ##import commons.beam_factorization_BF as BF
    ##import commons.beam_factorization_BS as BS

    # Specific to the Torino scheme
    import torino_config
    import NLO.local_currents as NLO_local_currents
    import NLO.integrated_currents as NLO_integrated_currents
    ##import NNLO.local_currents as NNLO_local_currents

    import sectors as sectors
    loaded_attributes['sector_generator'] = sectors.SectorGenerator()
#Ez Doesn't work
#    logger1.info("Sectors: "+str(list(loaded_attributes['sector_generator'])))

    # Note: specifying below which resources are needed is optional
    loaded_attributes['exporter'] = ColorfulPPExporter(relative_resource_paths=[
        'subtraction_schemes/%s'%subtraction_scheme_name
    ])

    all_subtraction_current_classes = []

    ###########
    # NLO
    ###########
    # Soft included in final part. Includes initial case.
    NLO_final_collinears = [
        NLO_local_currents.QCD_TRN_C_FgFq,
        NLO_local_currents.QCD_TRN_C_FqFqx,
        NLO_local_currents.QCD_TRN_C_FgFg,
        NLO_local_currents.QCD_TRN_S_g,
        NLO_local_currents.QCD_TRN_CS_FgFq,
        NLO_local_currents.QCD_TRN_CS_FgFg,
    ]
    # final-final collinears
    all_subtraction_current_classes.extend(NLO_final_collinears)

    NLO_initial_collinears = [
        NLO_local_currents.QCD_TRN_C_IgFq,
        NLO_local_currents.QCD_TRN_C_IqFg,
        NLO_local_currents.QCD_TRN_C_IqFq,
        NLO_local_currents.QCD_TRN_C_IgFg,
        NLO_local_currents.QCD_TRN_CS_IqFg,
        NLO_local_currents.QCD_TRN_CS_IgFg,
    ]
    # initial-final collinears
    all_subtraction_current_classes.extend(NLO_initial_collinears)

    # Add NLO integrated counterterms
    # ===============================
    all_subtraction_current_classes.extend([
                NLO_integrated_currents.QCD_integrated_TRN_C_FqFg,
                NLO_integrated_currents.QCD_integrated_TRN_C_FqFqx,
                NLO_integrated_currents.QCD_integrated_TRN_C_FgFg,
                NLO_integrated_currents.QCD_integrated_TRN_S_g,
                NLO_integrated_currents.QCD_integrated_TRN_CS_FqFg,
                NLO_integrated_currents.QCD_integrated_TRN_CS_FgFg
    ])

    NLO_integrated_initial_collinears = [
        NLO_integrated_currents.QCD_F0,
        NLO_integrated_currents.QCD_F0_lepton,
        NLO_integrated_currents.QCD_integrated_TRN_C_IgFq,
        NLO_integrated_currents.QCD_integrated_TRN_C_IqFg,
        NLO_integrated_currents.QCD_integrated_TRN_C_IqFq,
        NLO_integrated_currents.QCD_integrated_TRN_C_IgFg,
        NLO_integrated_currents.QCD_integrated_TRN_CS_IqFg,
        NLO_integrated_currents.QCD_integrated_TRN_CS_IgFg,
    ]
    # initial-final collinears
    all_subtraction_current_classes.extend(NLO_integrated_initial_collinears)

    # Finally register the subtraction current classes loaded
    loaded_attributes['all_subtraction_current_classes'] = all_subtraction_current_classes