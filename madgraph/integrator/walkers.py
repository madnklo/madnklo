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

import madgraph
import madgraph.various.misc as misc
import madgraph.integrator.mappings as mappings
import madgraph.core.subtraction as sub
from madgraph import InvalidCmd, MadGraph5Error

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')

#=========================================================================================
# Stroll
#=========================================================================================

class Stroll(object):
    """Container for a mapping call."""

    def __init__(self, mapping, structure, squared_masses, currents, variables=None):

        super(Stroll, self).__init__()
        assert isinstance(mapping, mappings.VirtualMapping)
        self.mapping = mapping
        assert isinstance(structure, sub.SingularStructure)
        self.structure = structure
        self.squared_masses = squared_masses
        for current in currents:
            assert isinstance(current, sub.Current)
        self.currents = currents
        self.variables = variables

    def __str__(self):

        foo = self.mapping.__class__.__name__
        foo += " with structure " + str(self.structure)
        if self.squared_masses:
            foo += str(self.squared_masses)
        foo += " for currents: "
        foo += ", ".join(str(current) for current in self.currents)
        if self.variables is not None:
            foo += " (with extra variables %s)" % str(self.variables)
        return foo

#=========================================================================================
# Hike
#=========================================================================================

class Hike(list):
    """Container for a sequence of mapping calls."""

    def __new__(cls, *args, **opts):

        for arg in args:
            assert isinstance(arg, Stroll)
        return super(Hike, cls).__new__(cls, *args, **opts)

    def __str__(self):

        foo = "--- Hike start ---\n"
        foo += "\n".join(
            ["Stroll "+str(i+1)+": "+str(stroll) for (i, stroll) in enumerate(self)] )
        foo += "\n---  Hike end  ---"
        return foo

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
    def determine_hike(cls, counterterm):
        """Determine the sequence of elementary mappings that must be applied
        in order to map to lower multiplicity.
        """

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

    def get_hike(self, counterterm):
        """Try to recycle the hike from a cached previous determination
        using the counterterm hash.
        """

        # TODO Implement
        # Remember to make the cache instance-dependent (self.###)

        return self.determine_hike(counterterm)

    def walk_to_lower_multiplicity(
        self, PS_point, counterterm,
        compute_kinematic_variables=False, compute_jacobian=False, verbose=False ):
        """Starting from the highest-multiplicity phase-space point,
        generate all lower-multiplicity phase-space points
        that are necessary for the evaluation of the given counterterm.

        :param PS_point: highest-multiplicity starting phase-space point,
        as a dictionary that associates integers to Lorentz vectors.

        :param counterterm: Counterterm object that specifies
        clusters of particle and recoilers recursively.
        Momenta_dict will obtained directly from it.

        :param compute_kinematic_variables: flag that specifies whether to compute
        all kinematic variables needed to recover the starting phase-space point
        from the lowest multiplicity one.

        :param compute_jacobian: flag that specifies whether to compute
        the phase-space jacobian.

        :return: a dictionary with the following entries:
        'currents', a list of stroll output dictionaries
            {stroll_currents, higher_PS_point, lower_PS_point, stroll_vars}
            for each stroll/mapping applied;
        'matrix_element', the reduced matrix element,
            paired with its corresponding reduced phase-space point;
        'kinematic_variables', a dictionary of all variables needed to recover
            the starting phase-space point from the lowest multiplicity one,
            or None if such variables were not requested.
        """

        # Identify the starting phase-space point
        point = PS_point
        if verbose: logger.debug(str(point))
        # Initialize return variables
        currents_4_eval = []
        kinematic_variables = dict() if compute_kinematic_variables else None
        # Determine the hike
        hike = self.get_hike(counterterm)
        # Walk along the hike
        for stroll in hike:
            # Map to lower multiplicity
            new_point, vars = stroll.mapping.map_to_lower_multiplicity(
                point, stroll.structure, counterterm.momenta_dict, stroll.squared_masses,
                kinematic_variables=kinematic_variables,
                compute_jacobian=compute_jacobian )
            if verbose: logger.debug(str(new_point)+"\n"+str(vars))
            # Append the current and the momenta
            stroll_output = {
                'stroll_currents': stroll.currents,
                'higher_PS_point': point,
                'lower_PS_point': new_point,
                'stroll_vars': vars }
            currents_4_eval.append(stroll_output)
            point = new_point
        # Identify reduced matrix element,
        # computed in the point which has received all mappings
        ME_PS_pair = [counterterm.process, point]
        # Return
        return {
            'currents': currents_4_eval,
            'matrix_element': ME_PS_pair,
            'kinematic_variables': kinematic_variables }

    def walk_to_higher_multiplicity(
        self, PS_point, counterterm, kinematic_variables,
        compute_jacobian=False, verbose=False ):
        """Starting from the lowest-multiplicity phase-space point,
        generate all higher-multiplicity phase-space points
        that are necessary for the evaluation of the given counterterm.

        :param PS_point: lowest-multiplicity starting phase-space point,
        as a dictionary that associates integers to Lorentz vectors.

        :param counterterm: Counterterm object that specifies
        clusters of particle and recoilers recursively.
        Momenta_dict will obtained directly from it.

        :param kinematic_variables: dictionary of all variables needed to recover
        the highest-multiplicity phase-space point from the starting one.

        :param compute_jacobian: flag that specifies whether to compute
        the phase-space jacobian.

        :return: a dictionary with the following entries:
        'currents', a list of stroll output dictionaries
            {stroll_currents, higher_PS_point, lower_PS_point, stroll_vars}
            for each stroll/mapping applied;
        'matrix_element', the reduced matrix element,
            paired with its corresponding reduced phase-space point.
        """

        # Identify reduced matrix element,
        # computed in the lowest multiplicity point
        ME_PS_pair = [counterterm.process, PS_point]
        # Identify the starting phase-space point
        point = PS_point
        if verbose: logger.debug(str(point))
        # Initialize return variables
        currents_4_eval = []
        # Determine the hike
        hike = self.get_hike(counterterm)
        for stroll in reversed(hike):
            # Compute jacobian and map to higher multiplicity
            new_point, vars = stroll.mapping.map_to_higher_multiplicity(
                point, stroll.structure, counterterm.momenta_dict, kinematic_variables,
                compute_jacobian=compute_jacobian )
            if verbose: logger.debug(str(new_point)+"\n"+str(vars))
            # Update jacobian and mapping variables
            stroll_output = {
                'stroll_currents': stroll.currents,
                'higher_PS_point': new_point,
                'lower_PS_point': point,
                'stroll_vars': vars }
            currents_4_eval.insert(0, stroll_output)
            point = new_point
        # Return
        return {
            'currents': currents_4_eval,
            'matrix_element': ME_PS_pair, }

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

    @classmethod
    def determine_hike(cls, counterterm):

        # Initialize return variables
        hike = Hike()
        all_nodes = [ node for node in counterterm.nodes if not isinstance(node.current, 
                    (sub.BeamCurrent, sub.IntegratedCurrent, sub.IntegratedBeamCurrent) ) ]
        
        # If the counterterm is not trivial
        if len(all_nodes)>=1:
            # Check it is only one
            if len(counterterm.nodes) != 1:
                raise MadGraph5Error(cls.cannot_handle_msg(counterterm))
            # Alias node and singular structure for convenience
            node = counterterm.nodes[0]
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.nodes:
                raise MadGraph5Error(cls.cannot_handle_msg(counterterm))
            # Get parent and children numbers
            parent, children, _ = mappings.get_structure_numbers(
                ss, counterterm.momenta_dict )
            # Pick recoilers
            recoilers = cls.get_recoilers(counterterm, (parent, ))
            structure = sub.SingularStructure(legs=recoilers, substructures=(ss,))
            # Choose mapping to use
            mapping = cls.determine_mapping(ss)
            # Determine the parent mass
            if parent is not None:
                try:
                    parent_pdg = counterterm.find_leg(parent)['id']
                except KeyError:
                    raise MadGraph5Error(
                        "Impossible to find parent " + str(parent) +
                        " within counterterm " + str(counterterm) )
                parent_particle = counterterm.process['model'].get_particle(parent_pdg)
                if parent_particle['mass'].lower() != 'zero':
                    raise MadGraph5Error("DEVELOPER: retrieve parent mass!")
            squared_masses = None
            # Compute jacobian and map to lower multiplicity
            hike.append(Stroll(mapping, structure, squared_masses, (node.current, )))
        # Return
        return hike

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
    soft_map = mappings.SoftVsFinalMapping()
    soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, collinear_map)
    only_colored_recoilers = True

class FinalLorentzNLOWalker(FinalNLOWalker):

    collinear_map = mappings.FinalLorentzOneMapping()
    soft_map = mappings.SoftVsFinalMapping()
    soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, collinear_map)
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
    soft_map = mappings.SoftVsFinalMapping()
    f_soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, f_collinear_map)
    i_soft_collinear_map = mappings.SoftCollinearVsFinalMapping(soft_map, i_collinear_map)
    # The integrated counterterms are only correct when recoiling against *all* final states.
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

#=========================================================================================
# Walker for disjoint counterterms
#=========================================================================================

class DisjointWalker(OneNodeWalker):
    """Implement a generic determine_hike for disjoint counterterms."""

    @classmethod
    def determine_hike(cls, counterterm):

        # Initialize return variables
        hike = Hike()
        # If the counterterm is not trivial
        if not counterterm.nodes:
            return hike
        parents = []
        currents = []
        substructures = []
        soft_particles = []
        for node in counterterm.nodes:
            # Do not accept nested currents
            if node.nodes:
                raise MadGraph5Error(cls.cannot_handle_msg(counterterm))
            # Alias singular structure for convenience
            ss = node.current['singular_structure']
            # Get parent and children numbers
            parent, _, _ = mappings.get_structure_numbers(ss, counterterm.momenta_dict )
            # Collinear structures
            if parent is not None:
                # Exclude parents from recoilers
                parents.append(parent)
                currents.append(node.current)
                substructures.append(ss)
                # Determine the parent mass
                # TODO: For now, check that it is zero
                try:
                    parent_pdg = counterterm.find_leg(parent)['id']
                except KeyError:
                    raise MadGraph5Error(
                        "Impossible to find parent " + str(parent) +
                        " within counterterm " + str(counterterm) )
                parent_particle = counterterm.process['model'].get_particle(parent_pdg)
                if parent_particle['mass'].lower() != 'zero':
                    raise MadGraph5Error("DEVELOPER: retrieve parent mass!")
            # Soft structures
            else:
                # Group all soft particles in one structure
                soft_particles += ss.get_all_legs()
                currents.append(node.current)
        squared_masses = None
        if soft_particles:
            substructures.append(sub.SoftStructure(legs=soft_particles))
        # Pick recoilers as everything in the final state of the reduced process
        recoilers = cls.get_recoilers(counterterm, parents)
        structure = sub.SingularStructure(legs=recoilers, substructures=substructures)
        # Choose mapping to use
        mapping = cls.determine_mapping(structure)
        # Compute jacobian and map to lower multiplicity
        hike.append(Stroll(mapping, structure, squared_masses, currents))
        # Returns
        return hike

# Walker for disjoint final-collinear counterterms
#=========================================================================================

class FinalCollinearDisjointWalker(DisjointWalker):
    """Implement a generic determine_hike for disjoint final-collinear counterterms."""

    collinear_map = None

    @classmethod
    def good_recoiler(cls, model, leg):

        return leg['state'] == leg.FINAL

    @classmethod
    def determine_mapping(cls, structure):

        for substructure in structure.substructures:
            if substructure.name() != 'C':
                raise MadGraph5Error(cls.cannot_handle_msg(structure))
            if substructure.substructures:
                raise MadGraph5Error(cls.cannot_handle_msg(structure))
            if structure.legs.has_initial_state_leg():
                raise MadGraph5Error(cls.cannot_handle_msg(structure))
        return cls.collinear_map

class FinalLorentzDisjointWalker(FinalCollinearDisjointWalker):

    collinear_map = mappings.FinalLorentzMapping()

class FinalGroupingDisjointWalker(FinalCollinearDisjointWalker):

    collinear_map = mappings.FinalGroupingMapping()

# Walker for disjoint final-collinear counterterms
#=========================================================================================

class SoftDisjointWalker(DisjointWalker):
    """Implement a generic determine_hike for disjoint soft counterterms."""

    soft_map = None

    @classmethod
    def good_recoiler(cls, model, leg):

        return leg['state'] == leg.FINAL

    @classmethod
    def determine_mapping(cls, structure):

        for substructure in structure.substructures:
            if substructure.name() != 'S':
                raise MadGraph5Error(cls.cannot_handle_msg(structure))
            if substructure.substructures:
                raise MadGraph5Error(cls.cannot_handle_msg(structure))
            if structure.legs.has_initial_state_leg():
                raise MadGraph5Error(cls.cannot_handle_msg(structure))
        return cls.soft_map

class SoftVsFinalDisjointWalker(SoftDisjointWalker):

    soft_map = mappings.SoftVsFinalMapping()

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
    'FinalLorentzDisjoint': FinalLorentzDisjointWalker,
    'FinalGroupingDisjoint': FinalGroupingDisjointWalker,
    'SoftVsFinalDisjoint': SoftVsFinalDisjointWalker,
    'SoftBeamsRecoilNLO': SoftBeamsRecoilNLOWalker
}
