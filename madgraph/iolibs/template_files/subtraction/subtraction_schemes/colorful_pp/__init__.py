###########################################################
#
# colorful_pp subtraction scheme
#
###########################################################

import os
import madgraph.various.misc as misc

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
    from madgraph.core.subtraction import SubtractionLeg

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

            #misc.sprint([[str(c) for c in sub_ss.decompose()] for sub_ss in singular_structure.substructures])
            #misc.sprint(requires_correlated_beam_convolution)
            return requires_correlated_beam_convolution

        def get_n_non_factorisable_double_sided_convolution(self, counterterm):
            """ Test if this counterterm involves a non-factorisable double-sided convolution which appears when one has
            a combination of a correlated convolution with single-sided ones (like for IF integrated collinear CT or PDF
            counterterms). For now it can only return 2 convolutions if non-factorisable; it may need to be extended for N^3LO."""

            #misc.sprint("doing %s"%str(counterterm))
            #misc.sprint(str(counterterm.get_integrated_current()))
            #misc.sprint(type(counterterm.get_integrated_current()),counterterm.get_integrated_current().__class__.__name__)
            #misc.sprint(str(counterterm.get_integrated_current().get('singular_structure')))
            #misc.sprint(self.does_require_correlated_beam_convolution(counterterm.get_integrated_current().get('singular_structure')))

            #misc.sprint([leg.n for node in counterterm.nodes for leg in node.current['singular_structure'].get_all_legs()])
            all_legs_states = [leg.state for node in counterterm.nodes for leg in node.current['singular_structure'].get_all_legs()]
            #misc.sprint(all_legs_states)

            n_non_factorisable_double_sided_convolution = (    2 if (
                    # DO NOT substitute the line below by 'counterterm.does_require_correlated_beam_convolution()'
                    # since this would then be incorrect as it would force to False for IntegratedBeamCurrents
                    self.does_require_correlated_beam_convolution(
                                    counterterm.get_integrated_current().get('singular_structure')) and
                    all_legs_states.count(SubtractionLeg.INITIAL)>0 and
                    # The CT ((C(2,4),S(5),),) requires a non-factorisable beam convolution.
                    # but the soft-collinear CT ((C(S((5, 21, f)),(2, 21, i)),),) requires a one-sided one.
                    # We therefore need the condition below
                    #( (counterterm.count_unresolved()-all_legs_states.count(SubtractionLeg.FINAL))>0 or
                    #  any(
                    #      any(len(ss.substructures)>1 for ss in node.current['singular_structure'].substructures)
                    #      for node in counterterm.nodes
                    #     )
                    #  )
                    #)
                    len(counterterm.nodes)>1 or any(
                        len(node.current['singular_structure'].substructures)>1 for node in counterterm.nodes)
            ) else 0 )
            #misc.sprint(counterterm.does_require_correlated_beam_convolution())
            #misc.sprint(all_legs_states.count(SubtractionLeg.INITIAL))
            #misc.sprint(counterterm.count_unresolved())
            #misc.sprint(all_legs_states.count(SubtractionLeg.FINAL))
            #misc.sprint([[len(ss.substructures) for ss in node.current['singular_structure'].substructures] for node in counterterm.nodes])
            #misc.sprint(str(counterterm),n_non_factorisable_double_sided_convolution)
            return n_non_factorisable_double_sided_convolution

    import factors_and_cuts as factors_and_cuts

    # Specific to the colorful_pp scheme
    import colorful_pp_config
    # NLO local
    import NLO.local_currents as NLO_local_currents
    # NLO integrated
    import NLO.integrated_currents as NLO_integrated_currents
    # NNLO local
    import NNLO.local_currents as NNLO_local_currents
    import NNLO.local_one_loop_currents as NNLO_local_one_loop_currents
    # NNLO integrated
    import NNLO.integrated_currents as NNLO_integrated_currents
    import NNLO.integrated_one_loop_currents as NNLO_integrated_one_loop_currents
    # N3LO local
    import NNNLO.local_currents as NNNLO_local_currents

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
        # Add zero PDF counterterms for supporting lepton collisions with this scheme too
        NLO_integrated_currents.QCD_F0_lepton,
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

    # Add beam currents which are also integrated counterterpart of all local counterterms recoiling against the
    # intiial states
    # First the integrated counterterms featuring one-loop currents
    all_subtraction_current_classes.extend([
        # PDF counterterms
        NNLO_integrated_one_loop_currents.TODO_QCD_1_F0,
        # Add zero PDF counterterms for supporting lepton collisions with this scheme too
        NNLO_integrated_one_loop_currents.QCD_1_F0_lepton,
        # collinear IF counterterms
        NNLO_integrated_one_loop_currents.TODO_QCD_1_integrated_C_IF,
        # collinear FF counterterms
        NNLO_integrated_one_loop_currents.TODO_QCD_1_integrated_C_FqFq,
        NNLO_integrated_one_loop_currents.TODO_QCD_1_integrated_C_FqFg,
        NNLO_integrated_one_loop_currents.TODO_QCD_1_integrated_C_FgFg,
        # soft F counterterms
        NNLO_integrated_one_loop_currents.TODO_QCD_1_integrated_S_Fg,
        # soft-collinear counterterms
        NNLO_integrated_one_loop_currents.QCD_1_integrated_CS_IF,
        NNLO_integrated_one_loop_currents.QCD_1_integrated_CS_FF,
    ])

    # Add compound kernels
    all_subtraction_current_classes.extend([
        # Disable counterterms of the form ([((F,X,)]) (type: (IntegratedBeamCurrent) ) since they are already accounted
        # for in colorful_pp by counterterms of the form: ([((X,))], :(F,):)
        NNLO_integrated_currents.QCD_disable_single_block_compound_PDF_counterterms,
        # Combined endpoint convolution of the NLO PDF counterterms of both beams
        NNLO_integrated_currents.TOTEST_QCD_integrated_Fx_Fx,
        # Single soft integrated together with a PDF counterterm
        NNLO_integrated_currents.TODO_QCD_integrated_S_Fg_F,
        # Single soft-collinear integrated together with a PDF counterterm, they are all zero since already accounted
        # for in the current above
        NNLO_integrated_currents.QCD_integrated_CS_X_F
    ])

    # Then the integrated counterterms of doubly unresolved limits
    all_subtraction_current_classes.extend([
        # Triple collinear limits
        NNLO_integrated_currents.TODO_QCD_integrated_C_IqFqpFqpx,
        NNLO_integrated_currents.TODO_QCD_integrated_C_IqFqFqx,
        # Nested collinear limits
        NNLO_integrated_currents.TODO_QCD_integrated_C_FqFqx_C_IqpFqFqx,
        NNLO_integrated_currents.TODO_QCD_integrated_C_FqFqx_C_IqFqFqx,
        # Soft limit (containing also all soft-collinear overlaps)
        NNLO_integrated_currents.TODO_QCD_integrated_S_FqFqx,
        # Soft-collinear integrated kernels that are zero since already all contained in the integrated soft
        NNLO_integrated_currents.QCD_integrated_S_FqFqx_C_IgFqFqx,
        NNLO_integrated_currents.QCD_integrated_S_FqFqx_C_IqFqFqx,
        NNLO_integrated_currents.QCD_integrated_S_FqFqx_C_FqFqx,
        NNLO_integrated_currents.QCD_integrated_S_FqFqx_C_FqFqx_C_IqpFqFqx,
        NNLO_integrated_currents.QCD_integrated_S_FqFqx_C_FqFqx_C_IqFqFqx,
    ])

    # Double-unresolved local counterterms
    # For now we are only trying an elementary IFF q > q q' q' collinear
    all_subtraction_current_classes.extend([
        # S(FF)
        NNLO_local_currents.QCD_S_FqFqx,
        NNLO_local_currents.QCD_S_FgFg,
        # C(IFF)
        NNLO_local_currents.QCD_C_IqFqpFqpx,
        NNLO_local_currents.QCD_C_IqFqFqx,
        NNLO_local_currents.QCD_C_IqFgFg,
        # C(FFF)
        NNLO_local_currents.QCD_C_FqFqpFqpx,
        NNLO_local_currents.QCD_C_FqFqFqx,
        NNLO_local_currents.QCD_C_FqFgFg,
        # C(I,S(FF))
        NNLO_local_currents.QCD_S_FqFqx_C_IgFqFqx,
        NNLO_local_currents.QCD_S_FqFqx_C_IqFqFqx,
        # S(C(FF))
        NNLO_local_currents.QCD_S_FqFqx_C_FqFqx,
        # C(I,C(FF))
        NNLO_local_currents.QCD_C_FqFqx_C_IqpFqFqx,
        NNLO_local_currents.QCD_C_FqFqx_C_IqFqFqx,
        # C(C(IF),F)
        NNLO_local_currents.QCD_C_IqFqx_C_IqFqFqx,
        # C(S(C(FF)),I)
        NNLO_local_currents.QCD_S_FqFqx_C_FqFqx_C_IqpFqFqx,
        NNLO_local_currents.QCD_S_FqFqx_C_FqFqx_C_IqFqFqx,
    ])

    # One-loop singly unresolved local counterterms
    all_subtraction_current_classes.extend([
        # final-final collinears
        NNLO_local_one_loop_currents.TODO_QCD_1_C_FqFqx,
        NNLO_local_one_loop_currents.TODO_QCD_1_C_FgFg,
        NNLO_local_one_loop_currents.TODO_QCD_1_C_FqFg,
        # initial-final collinears
        NNLO_local_one_loop_currents.TODO_QCD_1_C_IgFq,
        NNLO_local_one_loop_currents.TODO_QCD_1_C_IqFg,
        NNLO_local_one_loop_currents.TODO_QCD_1_C_IqFq,
        NNLO_local_one_loop_currents.TODO_QCD_1_C_IgFg,
        # soft and soft-collinears
        NNLO_local_one_loop_currents.TODO_QCD_1_S_g,
        NNLO_local_one_loop_currents.TODO_QCD_1_CS_FgFg,
        NNLO_local_one_loop_currents.TODO_QCD_1_CS_FgFq,
        NNLO_local_one_loop_currents.TODO_QCD_1_CS_IgFg,
        NNLO_local_one_loop_currents.TODO_QCD_1_CS_IqFg,
    ])

    # Add NNLO integrated counterterms
    all_subtraction_current_classes.extend([
        # There are None in this scheme since *all* local counterterms
        # recoil against the initial states.
    ])

    ###########
    # NNNLO
    ###########

    all_subtraction_current_classes.extend([
        # S(FFF)
        NNNLO_local_currents.QCD_S_FgFgFg,
        # C(FFFF)
        NNNLO_local_currents.QCD_C_FqFqxFqpFqpx,
        NNNLO_local_currents.QCD_C_FgFgFgFg,
    ])


    # Finally register the subtraction current classes loaded
    loaded_attributes['all_subtraction_current_classes'] = all_subtraction_current_classes