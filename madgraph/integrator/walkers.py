##########################################################################################
#
# Copyright (c) 2017 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
##########################################################################################

import logging

import madgraph.various.misc as misc
import madgraph.integrator.mappings as mappings
import madgraph.core.subtraction as sub
from madgraph import MadGraph5Error

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')

#=========================================================================================
# Low-level approach limit function
#=========================================================================================

def low_level_approach_limit(
    PS_point, mappings_sequence, scaling_parameter, momenta_dict,
    verbose=False):
    """Parametrically approach a limit
    starting from a given resolved kinematic configuration.

    :param PS_point: starting phase-space point

    :param mappings_sequence: a list of tuples (mapping, structure, power)
    that specifies a sequence of mappings to be applied with their input structure,
    and the weight of the rescaling to be applied at each step

    :param scaling_parameter: the value in (0,1) of the dimensionless
    scaling parameter which regulate the distance from the limit

    :param momenta_dict: the momentum dictionary used by all mappings
    to determine the leg numbers of parent and children

    :param verbose: False suppresses any debugging output
    """

    tmp_PS_point = PS_point.get_copy()
    kin_variables_list = list()
    if verbose:
        misc.sprint(PS_point.__str__(n_initial=1))
        # misc.sprint(PS_point)
        misc.sprint(momenta_dict)
        misc.sprint("*** Walking down ***")
    for mapping, structure, power in mappings_sequence:
        if verbose:
            misc.sprint(mapping.__class__.__name__, str(structure))
        kin_variables = {}
        sub.IRSubtraction.update_momenta_dict(momenta_dict, structure)
        tmp_PS_point, _ = mapping.map_to_lower_multiplicity(
            tmp_PS_point, structure, momenta_dict, None, kin_variables)
        if verbose:
            misc.sprint(tmp_PS_point.__str__(n_initial=1))
            # misc.sprint(tmp_PS_point)
            misc.sprint(momenta_dict)
            misc.sprint(kin_variables)
        mapping.rescale_kinematic_variables(
            structure, momenta_dict, kin_variables, scaling_parameter ** power)
        if verbose:
            misc.sprint(power)
            misc.sprint(kin_variables)
        kin_variables_list.append(kin_variables)
    if verbose:
        misc.sprint("*** Walking up ***")
    for mapping, structure, power in reversed(mappings_sequence):
        tmp_PS_point, _ = mapping.map_to_higher_multiplicity(
            tmp_PS_point, structure, momenta_dict, kin_variables_list[-1])
        kin_variables_list = kin_variables_list[:-1]
        if verbose:
            misc.sprint(mapping.__name__, str(structure))
            misc.sprint(tmp_PS_point.__str__(n_initial=1))
            # misc.sprint(tmp_PS_point)
    return tmp_PS_point

#=========================================================================================
# Generic walker parent class
#=========================================================================================

class VirtualWalker(object):
    """Base class for walker implementations."""
    
    def __new__(cls, walker=None, **opts):
        """Factory class to make plugin easy."""

        if cls is VirtualWalker:
            if walker is None:
                raise MadGraph5Error(
                    "VirtualWalker called without a walker name.")
            if not walker_classes_map.has_key(walker):
                raise MadGraph5Error(
                    "Unknown mapping walker of type '%s'." % walker )
            target_class = walker_classes_map[walker]
            return super(VirtualWalker, cls).__new__(target_class, **opts)
        else:
            return super(VirtualWalker, cls).__new__(cls, **opts)

    def __init__(self, model=None, **opts):
        """General initialization of a walker.
        The model is an optional specification,
        which can be useful to know properties of the leg mapped.
        """

        self.model = model

    @classmethod
    def determine_mapping(cls, structure):
        """Determine which elementary mapping to use for a given singular structure."""

        raise NotImplemented

    @classmethod
    def get_recoilers(cls, counterterm, exclude=None):
        """Select particles to be used as recoilers for a given counterterm."""

        raise NotImplemented

    @classmethod
    def decompose_counterterm(cls, counterterm, counterterms):
        """Determine the sequence of elementary mappings that must be applied
        in order to approach the limit.
        """

        complete_ss = counterterm.reconstruct_complete_singular_structure()
        decomposed_sss = sorted(
            complete_ss.decompose(), key=lambda x: x.__str__(True, True, True) )
        decomposed_cts = []
        for ss in decomposed_sss:
            found = False
            for ct in counterterms:
                if len(ct.nodes) != 1: continue
                if ct.nodes[0].nodes: continue
                ss2 = ct.nodes[0].current['singular_structure']
                if ss != ss2: continue
                decomposed_cts.append(ct)
                found = True
                break
            if not found: raise MadGraph5Error('Counterterm not found')
        return decomposed_cts

    def approach_limit(
        self, PS_point, structure, scaling_parameter, process ):
        """Produce a higher multiplicity phase-space point from PS_point,
        according to kinematic_variables that approach the limit of structure
        parametrically with scaling_parameter.
        """

        # Decompose the counterterm
        decomposed = structure.decompose()        
        # The rescaling of the convolution variables is done independently of the mapping
        # and should therefore not be considered
        decomposed = [step for step in decomposed if step.name() != "F"]
        
        # Always approach the limit at the same speed
        base = scaling_parameter ** (1. / max(len(decomposed),1))
        #base = scaling_parameter
        # Prepare a momentum dictionary for each mapping
        mom_dict = sub.bidict()
        for leg in process['legs']:
            mom_dict[leg['number']] = frozenset([leg['number'], ])
        parent_index = max(leg['number'] for leg in process['legs']) + 1
        fake_ct = sub.Counterterm(process=process, momenta_dict=mom_dict)
        closer_PS_point = PS_point.get_copy()
        # Walk the hike up and down
        for step in decomposed:
            mapping = self.determine_mapping(step)
            all_children = frozenset([leg.n for leg in step.get_all_legs()])
            recoilers = self.get_recoilers(fake_ct, exclude=all_children)
            # Below is a hack to recoil against the Higgs for C(1,3),C(2,4) of g g > d d~ h.
            #recoilers = [sub.SubtractionLeg(5,25,sub.SubtractionLeg.FINAL),]
            new_ss = sub.SingularStructure(substructures=[step, ], legs=recoilers)
            if step.name() == "C":
                mom_dict[parent_index] = all_children
            elif step.name() == "S":
                pass
            else:
                raise MadGraph5Error("Unrecognized structure of type " + step.name())
            kin_variables = {}
            #misc.sprint('Now doing step: %s'%str(step))
            #misc.sprint('Starting PS point:\n',str(closer_PS_point))
            low_PS_point, _ = mapping.map_to_lower_multiplicity(
                closer_PS_point, new_ss, mom_dict, None, kin_variables )
            #misc.sprint('Mapped down PS point:\n',str(low_PS_point))
            #misc.sprint('kin_variables=',kin_variables)
            mapping.rescale_kinematic_variables(
                new_ss, mom_dict, kin_variables, base)
            #misc.sprint('rescaled kin_variables=',base,kin_variables)
            closer_PS_point, _ = mapping.map_to_higher_multiplicity(
                low_PS_point, new_ss, mom_dict, kin_variables )
            #misc.sprint('Mapped up PS point:\n',str(closer_PS_point))
            #misc.sprint('kin_variables=',kin_variables)
            if parent_index in mom_dict.keys():
                del mom_dict[parent_index]
        return closer_PS_point

#=========================================================================================
# Generic walker for one-level counterterms
#=========================================================================================

class OneNodeWalker(VirtualWalker):
    """Implement a generic determine_hike for one-node counterterms."""

    @classmethod
    def cannot_handle_msg(cls, obj):

        name = obj.__class__.__name__
        return cls.__name__ + " cannot handle the " + name + ": " + str(obj)

    @classmethod
    def good_recoiler(cls, model, leg):
        """Indicate if a particle is apt to be a recoiler according to this walker."""

        raise NotImplemented

    @classmethod
    def not_excluded(cls, leg, exclude):

        return not(exclude and leg['number'] in exclude)

    @classmethod
    def get_recoilers(cls, counterterm, exclude=None):
        """Select particles in the reduced process to be used as recoilers."""

        legs = counterterm.process['legs']
        model = counterterm.process['model']
        recoilers = [
            sub.SubtractionLeg(leg)
            for leg in legs
            if cls.good_recoiler(model, leg) and cls.not_excluded(leg, exclude) ]
        return recoilers

# FinalCollinearOneWalker
#=========================================================================================

class FinalCollinearOneWalker(OneNodeWalker):

    collinear_map = None
    only_colored_recoilers = None

    @classmethod
    def good_recoiler(cls, model, leg):

        return (
            leg['state'] == leg.FINAL and not
            (cls.only_colored_recoilers and model.get_particle(leg['id'])['color'] == 1) )

    @classmethod
    def determine_mapping(cls, structure):

        if structure.name() == 'S' or structure.substructures:
            raise MadGraph5Error(cls.cannot_handle_msg(structure))
        else:
            return cls.collinear_map

class FinalRescalingOneWalker(FinalCollinearOneWalker):

    collinear_map = mappings.FinalRescalingOneMapping()
    only_colored_recoilers = True

class FinalLorentzOneWalker(FinalCollinearOneWalker):

    collinear_map = mappings.FinalLorentzOneMapping()
    only_colored_recoilers = False

# FinalNLOWalker
#=========================================================================================

class FinalNLOWalker(OneNodeWalker):

    collinear_map = None
    soft_map = None
    soft_collinear_map = None
    only_colored_recoilers = None

    @classmethod
    def good_recoiler(cls, model, leg):

        return (
            leg['state'] == leg.FINAL and not
            (cls.only_colored_recoilers and model.get_particle(leg['id'])['color'] == 1) )

    @classmethod
    def determine_mapping(cls, structure):

        if structure.name() == 'S':
            return cls.soft_map
        elif structure.name() == 'C' and not structure.substructures:
            return cls.collinear_map
        elif (structure.name() == 'C'
            and len(structure.substructures) == 1
            and structure.substructures[0].name() == 'S'
            and not structure.substructures[0].substructures
            and len(structure.legs) == 1 ):
            return cls.soft_collinear_map
        else:
            raise MadGraph5Error(cls.cannot_handle_msg(structure))

class FinalRescalingNLOWalker(FinalNLOWalker):

    collinear_map = mappings.FinalRescalingOneMapping()
    soft_map = mappings.SoftVsFinalPureRescalingMapping()
    soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, collinear_map)
    only_colored_recoilers = True

class FinalLorentzNLOWalker(FinalNLOWalker):

    collinear_map = mappings.FinalLorentzOneMapping()
    soft_map = mappings.SoftVsFinalPureRescalingMapping()
    soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, collinear_map) 
    # Beware that integrated counterterms are only correct
    # when recoiling against *all* final states
    only_colored_recoilers = True

# General NLO walker
#=========================================================================================

class NLOWalker(OneNodeWalker):

    f_collinear_map = None
    i_collinear_map = None
    soft_map = None
    f_soft_collinear_map = None
    i_soft_collinear_map = None
    only_colored_recoilers = None

    @classmethod
    def good_recoiler(cls, model, leg):

        return (
            leg['state'] == leg.FINAL and not
            (cls.only_colored_recoilers and model.get_particle(leg['id'])['color'] == 1) )

    @classmethod
    def determine_mapping(cls, structure):

        if structure.name() == 'S':
            return cls.soft_map
        elif structure.name() == 'C' and not structure.substructures:
            if structure.legs.has_initial_state_leg():
                return cls.i_collinear_map
            else:
                return cls.f_collinear_map
        elif (structure.name() == 'C'
              and len(structure.substructures) == 1
              and structure.substructures[0].name() == 'S'
              and not structure.substructures[0].substructures
              and len(structure.legs) == 1):
            if structure.legs.has_initial_state_leg():
                return cls.i_soft_collinear_map
            else:
                return cls.f_soft_collinear_map
        else:
            logger.critical("Error while processing %s" % structure)
            raise MadGraph5Error(cls.cannot_handle_msg("SingularStructure"))

class LorentzNLOWalker(NLOWalker):

    f_collinear_map = mappings.FinalLorentzOneMapping()
    i_collinear_map = mappings.InitialLorentzOneMapping()
    soft_map = mappings.SoftVsFinalPureRescalingMapping()
    f_soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, f_collinear_map)
    i_soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, i_collinear_map)
    # The integrated counterterms are only correct when recoiling against *all* final states
    # Take care that the soft mapping only works for massless particles instead
    only_colored_recoilers = False

class SoftBeamsRecoilNLOWalker(NLOWalker):
    """ Set of mappings designed to work for the NLO topology pp > X(color-singlet) + at most
    one jet. The collinear mapping is left untouched compared to LorentzNLOWalker, but the 
    soft one is not the original Colorful mapping where recoilers are individually rescaled but
    instead it is a mapping that rescales both initial states without boosting the c.o.m frame.
    This is well-suited for integrating 'p p > X (+j)' where X is a color-singlet. """
    
    f_collinear_map = mappings.FinalLorentzOneMapping()
    i_collinear_map = mappings.InitialLorentzOneMapping()
    
    # The two lines below yield the difference w.r.t LorentzNLOWalker
    soft_map = mappings.SoftVsInitialMapping()
    only_colored_recoilers = False
    
    f_soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, f_collinear_map)
    i_soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, i_collinear_map)

# Dictionary of available walkers
#=========================================================================================

# Mapping classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Mapping.
# Note that this must be placed after all Mapping daughter classes in this module
# have been declared.
walker_classes_map = {
    'FinalRescalingOne': FinalRescalingOneWalker,
    'FinalLorentzOne': FinalLorentzOneWalker,
    'FinalRescalingNLO': FinalRescalingNLOWalker,
    'FinalLorentzNLO': FinalLorentzNLOWalker,
    'LorentzNLO': LorentzNLOWalker,
    'SoftBeamsRecoilNLO': SoftBeamsRecoilNLOWalker
}
