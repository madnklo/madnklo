#!/usr/bin/env python
################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
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
################################################################################
"""Definition of all the classes and features relevant to the handling of 
higher order IR subtraction.
"""

import copy
import itertools
import collections
import logging
import types
import math
import sys
import os
if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(
        os.path.realpath(__file__)),
        os.path.pardir,
        os.path.pardir
    ))
from madgraph import MadGraph5Error, MG5DIR, InvalidCmd
import madgraph.core.base_objects as base_objects
import madgraph.fks.fks_common as fks_common
import madgraph.various.misc as misc
from bidict import bidict

logger = logging.getLogger('madgraph')
pjoin = os.path.join

#===============================================================================
# Multinomial function
#===============================================================================

def multinomial(i_s):
    """Compute the multinomial (i_1 + ... + i_n)!/(i_1! ... i_n!)."""

    itot = 0
    denom = 1
    for i in i_s:
        itot += i
        denom *= math.factorial(i)
    return math.factorial(itot)/denom

#===============================================================================
# SubtractionLeg
#===============================================================================
class SubtractionLeg(tuple):
    """Leg object specialized for subtraction."""

    FINAL = base_objects.Leg.FINAL
    INITIAL = base_objects.Leg.INITIAL

    def __new__(cls, *args, **opts):
        """Initialize a SubtractionLeg object from various representations."""

        target = (0, SubtractionLeg.FINAL)
        # One argument passed
        if len(args) == 1:
            leg = args[0]
            # Initialization from Leg object
            if isinstance(leg, base_objects.Leg):
                target = (leg.get('number'), leg.get('id'), leg.get('state'))
            # Initialization from dictionary
            elif isinstance(leg, dict):
                target = (leg['number'], leg['id'], leg['state'])
            # Initialization from iterable
            elif (
                hasattr(leg, '__len__') and len(leg) == 3 and
                isinstance(leg[0], int) and isinstance(leg[1], int) and
                isinstance(leg[2], bool)
            ):
                target = leg
            # Invalid argument
            else:
                raise MadGraph5Error(
                        "SubtractionLeg cannot be initialized with "
                        "%s of type %s." % (leg, str(type(leg)))
                )
        # Initialize with 3 arguments
        elif (
            (len(args) == 3) and
            isinstance(args[0], int) and isinstance(args[1], int) and
            isinstance(args[2], bool)
        ):
            target = (args)
        # Invalid number of arguments
        else:
            raise MadGraph5Error(
                    "SubtractionLeg cannot be initialized with argument %s." % str(args)
            )

        return super(SubtractionLeg, cls).__new__(cls, target)

    @property
    def n(self):

        return self[0]

    @property
    def pdg(self):

        return self[1]

    @property
    def state(self):

        return self[2]

    @staticmethod
    def state_str(state):

        if state == SubtractionLeg.FINAL:
            return 'f'
        return 'i'

    def __str__(self, print_n = True, print_pdg = False, print_state = False):
        """Return a string representation of this subtraction leg."""

        tmp = [str(self.n), str(self.pdg), SubtractionLeg.state_str(self.state)]
        if not print_state:
            del tmp[2]
        if not print_pdg:
            del tmp[1]
        if not print_n:
            del tmp[0]
        if tmp:
            if len(tmp) == 1:
                return tmp[0]
            else:
                return "(" + ", ".join(tmp) + ")"
        return ""

#===============================================================================
# SubtractionLegSet
#===============================================================================
class SubtractionLegSet(frozenset):
    """Set of SubtractionLeg objects."""

    def __new__(cls, *args, **opts):
        """Initialize set, trying to convert arguments into SubtractionLeg's."""    
        
        if not args:
            return super(SubtractionLegSet, cls).__new__(cls)

        if (
            isinstance(args[0], collections.Iterable) and
            not isinstance(args[0],(dict, SubtractionLeg))
        ):
            legs = args[0]
        else:
            legs = args

        return super(SubtractionLegSet, cls).__new__(
                cls,
                [SubtractionLeg(leg) for leg in legs]
        )


    def __init__(self, *args, **opts):
        """Initialize set, trying to convert arguments into SubtractionLeg's."""

        if not args:
            super(SubtractionLegSet, self).__init__()
            return

        if (
            isinstance(args[0], collections.Iterable) and
            not isinstance(args[0],(dict, SubtractionLeg))
        ):
            legs = args[0]
        else:
            legs = args

        super(SubtractionLegSet, self).__init__(
                [SubtractionLeg(leg) for leg in legs]
        )
        return

    def has_initial_state_leg(self):
        """Returns True if this leg set has at least one initial state leg,
        False otherwise.
        """

        for leg in self:
            if leg.state == SubtractionLeg.INITIAL:
                return True
        return False

#===============================================================================
# SingularStructure
#===============================================================================
class SingularStructure(object):
    """Object that represents a hierarchical structure of IR singularities."""

    def __init__(self, *args, **opts):
        """Initialize a hierarchical singular structure."""

        components = None
        if args:
            if isinstance(args[0], types.GeneratorType):
                components = tuple(a for a in args[0])
            elif (
                isinstance(args[0], collections.Iterable) and
                not isinstance(args[0], (dict, SubtractionLeg) )
            ):
                components = args[0]
        if not components:
            components = args

        # If this Structure is annihilated, any operator acting on it would
        # not apply and this structure evaluated to False.
        self.is_annihilated = False

        # Type check
        for a in components:
            if not isinstance(a, (SingularStructure, SubtractionLeg)):
                raise MadGraph5Error(
                        "SingularStructure initialized with invalid argument "
                        "%s of type %s." % (a, str(type(a)))
                )

        # List of substructures that this SubtractionOperators acts on
        self.substructures = [
            a for a in components if isinstance(a, SingularStructure)
        ]
        # Set of simple legs this SubtractionOperators acts on
        self.legs = SubtractionLegSet(
            # HACK without tuple and expansion fails for empty components
            *(a for a in components if isinstance(a, SubtractionLeg))
        )

    def __str__(self, print_n = True, print_pdg = False, print_state = False):
        """Return a string representation of the singular structure."""

        if self.is_annihilated:
            return 'Non-existing structure'

        tmp_str = self.name() + "("
        tmp_str += ",".join(sorted(
            sub.__str__(print_n, print_pdg, print_state)
            for sub in self.substructures
        ))
        if self.substructures:
            tmp_str += ","
        tmp_str += ",".join(sorted(
            leg.__str__(print_n, print_pdg, print_state)
            for leg in self.legs
        ))
        tmp_str += ")"
        return tmp_str

    def is_void(self):

        return self.is_annihilated

    def annihilate(self):
        """When an operator cannot act on this structure,
        remove this structure altogether
        by flagging it as unapplicable to any other operator.
        """

        # Really clear structures and legs
        del self.substructures[:]
        self.legs = SubtractionLegSet()
        self.is_annihilated = True

    def get_all_legs(self):
        """Return all legs involved in this singular structure."""

        all_legs = set(sub.get_all_legs() for sub in self.substructures)
        all_legs.add(self.legs)
        return SubtractionLegSet().union(*all_legs)

    def discard_leg_numbers(self):
        """Set all leg numbers to zero, discarding this information forever."""

        self.legs = SubtractionLegSet(
                SubtractionLeg(0, leg.pdg, leg.state)
                for leg in self.legs
        )
        for substructure in self.substructures:
            substructure.discard_leg_numbers()
        return

    def get_canonical_representation(self, track_leg_numbers = True):
        """Creates a canonical hashable representation of self."""
        
        canonical = {}

        canonical['is_annihilated'] = self.is_annihilated
        if track_leg_numbers:
            canonical['legs'] = self.legs
        else:
            canonical['legs'] = SubtractionLegSet(
                    SubtractionLeg(0,l.pdg,l.state)
                    for l in self.legs
            )
        canonical['name'] = self.name()
        canonical['substructures'] = tuple(
                structure.get_canonical_representation()
                for structure in self.substructures
        )

        return tuple(sorted(canonical.items()))

    def prefactor(self):
        """Determine the prefactor of the counterterm
        associated with this singular structure.
        """

        pref = -1
        if type(self) == SingularStructure:
            pref = 1
        # Account for simplification of soft operators
        # TODO Check this works for collinears inside softs
        softs = []
        for sub in self.substructures:
            pref *= sub.prefactor()
            if sub.name() == "S":
                softs += [len(sub.substructures) + len(sub.legs), ]
        return pref * multinomial(softs)

    def __eq__(self, other):
        """Check if two singular structures are the same."""

        assert isinstance(other, SingularStructure)
        self_can = self.get_canonical_representation()
        other_can = other.get_canonical_representation()
        return self_can == other_can

    def name(self):

        return ""

class SoftStructure(SingularStructure):

    def name(self):

        return "S"

class CollStructure(SingularStructure):

    def name(self):

        return "C"

#===============================================================================
# SingularOperator
#===============================================================================
class SingularOperator(SubtractionLegSet):
    """Virtual base class for elementary singular operators."""

    def __init__(self, *args, **opts):

        super(SingularOperator, self).__init__(*args, **opts)

    def __str__(self):
        """Return a simple string representation of the singular operator."""

        return self.name() + str(sorted((leg.n for leg in self)))

    def name(self):
        """Symbol used to represent this operator within output."""
        raise MadGraph5Error(
            "name called in SingularOperator of unspecified type."
        )

    def get_structure(self):
        """Singular structure corresponding to this operator
        acting on hard particles.
        """

        raise MadGraph5Error(
            "structure called in SingularOperator of unspecified type."
        )

    def act_here_needed(self, structure):
        """Return True if the counterterm obtained acting with this operator
        on the given structure at the first level is needed for subtraction.
        """

        # WARNING Rules may not extend to non-QCD cases
        raise MadGraph5Error(
            "act_here_needed called in SingularOperator of unspecified type."
        )

    def non_overlapping_with(self, structure):
        """Determines if this operator does not affect any leg
        that appears within the given structure, at any level.
        """

        assert isinstance(structure, SingularStructure)

        if not self.isdisjoint(structure.legs):
            return False

        for substructure in structure.substructures:
            if not self.non_overlapping_with(substructure):
                return False

        return True

    def act_on(self, structure):
        """Act with an elementary operator on a non-elementary structure."""

        assert isinstance(structure, SingularStructure)

        # If the limit does not overlap with the existing structure at all,
        # just append it at the end
        if self.non_overlapping_with(structure):
            structure.substructures.append(self.get_structure())
            return

        # If the limit acts at least partly at the current level
        if not self.isdisjoint(structure.legs):
            # If the limit acts at different levels,
            # it is not needed for subtraction
            if not self.issubset(structure.legs):
                structure.annihilate()
                return
            # If the limit acts completely within this level,
            # it may be needed or not
            else:
                if self.act_here_needed(structure):
                    structure.substructures.append(self.get_structure())
                    structure.legs = structure.legs.difference(self)
                    return
                else:
                    structure.annihilate()
                    return

        # The limit acts deeper
        for substructure in structure.substructures:
            # Find the first substructure hit by the limit, even partially
            if not self.non_overlapping_with(substructure):

                self.act_on(substructure)
                # If the action on the substructure annihilated it,
                # then annihilate this structure as well
                if substructure.is_void():
                    structure.annihilate()
                return

class SoftOperator(SingularOperator):
    """Object that represents a soft elementary singular operator."""

    def __init__(self, *args, **opts):

        super(SoftOperator, self).__init__(*args, **opts)
        return

    def name(self):

        return "S"

    def get_structure(self):

        return SoftStructure(self)

    def act_here_needed(self, structure):

        assert isinstance(structure, SingularStructure)

        # Soft limits are not needed within soft structures
        if structure.name() == "S":
            return False

        # Soft limits are needed within collinear structures
        # if there is at least one hard particle left
        if structure.name() == "C":
            # If there is a hard particle at this level, the limit is needed
            if len(self) < len(structure.legs):
                return True
            # If there is a hard collinear substructure, the limit is needed
            else:
                for substructure in structure.substructures:
                    if not substructure.name() == "S":
                        return True
            # If did not return yet, there is no hard particle left
            return False

        # By default, not skipping overlaps ensures correct subtraction
        return True

class CollOperator(SingularOperator):
    """Object that represents a collinear elementary singular operator."""

    def name(self):

        return "C"

    def get_structure(self):

        return CollStructure(self)

    def __init__(self, *args, **opts):

        super(CollOperator, self).__init__(*args, **opts)
        return

    def act_here_needed(self, structure):

        assert isinstance(structure, SingularStructure)

        # Collinear limits are needed within soft _and_ collinear structures
        # WARNING: this statement has to be checked
        # By default, not skipping overlaps ensures correct subtraction
        return True

#===============================================================================
# SingularOperatorList
#===============================================================================
class SingularOperatorList(list):
    """List of subtraction operators."""

    def __str__(self):
        """Print this singular operator list."""

        return "[" + ", ".join([str(leg) for leg in self]) + "]"
    
    def act_on(self, structure):
        """Act with a list of operators on a substructure."""

        # Empty list of operators to apply or invalid structure: done
        if structure.is_void() or not self:
            return structure
        # Apply last operator and advance recursion
        all_but_last_operators = SingularOperatorList(self[:-1])
        # Act using the last operator
        self[-1].act_on(structure)
        # And now act using the rest of the operators, if needed
        if not structure.is_void():
            return all_but_last_operators.act_on(structure)

    def simplify(self):
        """Simplify a list of operators,
        returning the real structure of the corresponding singularities.
        """

        structure = SingularStructure()
        self.act_on(structure)
        return structure

#===============================================================================
# Current
#===============================================================================
class Current(base_objects.Process):

    def default_setup(self):
        """Default values for all properties specific to current,
        additional to Process.
        """

        super(Current, self).default_setup()
        self['n_loops'] = None
        self['singular_structure'] = SingularStructure()
        return

    def discard_leg_numbers(self):
        """Discard all leg numbers in the singular_structure
        and in the parent_subtraction_leg,
        effectively getting the type of current
        without information about momenta for which it should be evaluated.
        """

        self['singular_structure'].discard_leg_numbers()
        return

    def __str__(self, print_n = True, print_pdg = True, print_state = True):
        """Convert this current to a nice readable string."""

        readable_string = self['singular_structure'].__str__(
            print_n, print_pdg, print_state
        )
        n_loops = self.get('n_loops')
        if not n_loops is None:
            readable_string += " @ " + str(n_loops) + " loop"
            if n_loops != 1:
                readable_string += "s"
        return readable_string

    def get_key(self):
        """Return the ProcessKey associated to this current."""

        from contributions import ProcessKey
        return ProcessKey(
                self,
                allowed_attributes = [
                    'singular_structure',
                    'n_loops',
                    'squared_orders'
                ]
            )

    def __eq__(self, other):
        """Compare two currents using their ProcessKey's."""

        return self.get_key().key_dict == other.get_key().key_dict

#===============================================================================
# CountertermNode
#===============================================================================
class CountertermNode(object):
    """Class representing a current attaching to legs and other currents."""

    def __init__(
        self,
        current = None,
        subcurrents = None
    ):

        if current:
            assert isinstance(current, base_objects.Process)
            self.current = current
        else:
            self.current = Current()
        if subcurrents:
            assert isinstance(subcurrents, list)
            self.subcurrents = subcurrents
        else:
            self.subcurrents = []

    def __str__(self, level = 0):

        tmp_str = "    " * level + str(self.current) + "\n"
        for subcurrent in self.subcurrents:
            tmp_str += subcurrent.__str__(level + 1)
        return tmp_str

    def n_loops(self):
        """If number of loops are assigned,
        return the total number of loops in this node and its children.
        """

        result = self.current.get('n_loops')
        if result is not None:
            for subcurrent in self.subcurrents:
                sub_n_loops = subcurrent.n_loops()
                if sub_n_loops is None:
                    result = None
                else:
                    result += sub_n_loops
        return result

    def get_copy(self, copied_attributes = ()):
        """Make sure that attributes args of the object returned
        and all its children can be modified without changing the original.
        """

        subcurrents_copy = [
            subcurrent.get_copy(copied_attributes)
            for subcurrent in self.subcurrents
        ]
        return CountertermNode(
            self.current.get_copy(copied_attributes),
            subcurrents_copy
        )

    def get_all_currents(self):
        """Return a list of all currents involved in this elementary node
        and its children.
        """

        currents = []
        for subcurrent in self.subcurrents:
            currents += subcurrent.get_all_currents()
        currents += [self.current, ]
        return currents

#===============================================================================
# Counterterm
#===============================================================================
class Counterterm(CountertermNode):
    """Class representing a tree of currents multiplying a matrix element."""

    def __init__(
        self,
        process = None,
        subcurrents = None,
        momenta_dict = None
    ):

        super(Counterterm, self).__init__(process, subcurrents)
        if momenta_dict:
            assert isinstance(momenta_dict, bidict)
            self.momenta_dict = momenta_dict
        else:
            self.momenta_dict = bidict()

    @property
    def process(self):
        return self.current

    def __str__(self, level = 0):

        tmp_str  = "    " * level + self.process.nice_string(0, True, False)
        tmp_str += " ("
        tmp_str += " ".join(
            str(leg['number'])
            for leg in self.process['legs']
            if leg['state'] == SubtractionLeg.INITIAL
        )
        tmp_str += " > "
        tmp_str += " ".join(
            str(leg['number'])
            for leg in self.process['legs']
            if leg['state'] == SubtractionLeg.FINAL
        )
        tmp_str += ")"
        process_n_loops = self.current.get('n_loops')
        if process_n_loops is None:
            pass
        else:
            tmp_str += " @ " + str(process_n_loops) + " loop"
            if process_n_loops != 1:
                tmp_str += "s"
        tmp_str += "\n"
        for subcurrent in self.subcurrents:
            tmp_str += subcurrent.__str__(level + 1)
        tmp_str += "    " * level + "Pseudoparticles: {"
        tmp_str += "; ".join(
            str(key) + ": (" +
            ",".join(str(n) for n in self.momenta_dict[key]) + ")"
            for key in self.momenta_dict
        )
        tmp_str += "}"
        return tmp_str

    def get_copy(self, copied_attributes = (), copy_momenta_dict = False):
        """Make sure that attributes args of the object returned
        and all its children can be modified without changing the original.
        """

        node = super(Counterterm, self).get_copy(copied_attributes)
        momenta_dict = self.momenta_dict
        if copy_momenta_dict:
            momenta_dict = bidict(momenta_dict)
        return Counterterm(node.current, node.subcurrents, momenta_dict)

    def get_all_currents(self):
        """Return a list of all currents involved in this counterterm."""

        currents = []
        for subcurrent in self.subcurrents:
            currents += subcurrent.get_all_currents()
        return currents

#===============================================================================
# order_2_string
#===============================================================================

def order_2_string(order):
    """Convert powers of the squared coupling to (N...)LO string."""

    assert isinstance(order, int)
    return order * 'N' + 'LO'

#===============================================================================
# IRSubtraction
#===============================================================================
class IRSubtraction(object):

    def __init__(self, model, orders = {'QCD': 1, 'QED': 0}):
        """Initialize a IR subtractions for a given model,
        correction order and type.
        """
        
        self.model = model
        self.orders = orders
        # Map perturbed coupling orders to the corresponding relevant interactions and particles.
        # The values of the dictionary are 'interactions', 'pert_particles' and 'soft_particles'
        self.IR_quantities_for_corrections_types = dict(
            (
                coupling,
                fks_common.find_pert_particles_interactions(
                    self.model,
                    pert_order = coupling
                )
            )
            for coupling in orders.keys()
        )


    def global_order(self):
        """Return the total perturbative order in all couplings."""

        return sum(pow for pow in self.orders.values())

    def can_be_IR_unresolved(self, PDG):
        """Check whether a particle given by its PDG can become unresolved
        and lead to singular behavior.
        """

        return any(
            (PDG in self.IR_quantities_for_corrections_types[coupling]['pert_particles'])
            for coupling in self.orders.keys()
        )

    def parent_PDGs(self, legs):
        """List all possible parent PDGs to a given set of legs."""

        # BALDY: parent_PDGs hardcoded to SM QCD (with some sanity checks)
        # A generic way to build histories of possible combination of external legs is using the 
        # model ref_dict_to1
        # misc.sprint(self.model.get('ref_dict_to1')[(-2,2)])
        # misc.sprint(self.model.get_particle(-2).get('spin'),
        #             self.model.get_particle(-2).get_color(),
        #             self.model.get_particle(-2).get('mass')=='zero'
        #             )

        if not self.model.get('name').startswith('sm'):
            raise InvalidCmd(
                "The function parent_PDGs is implemented for the SM only."
            )
        if any(order != 'QCD' for order in self.orders.keys()):
            raise InvalidCmd(
                "The function parent_PDGs is implemented for QCD only."
            )

        # Get parton flavors, eliminating gluons
        flavored_legs = [leg for leg in legs if leg.pdg != 21]

        # If all daughters were gluons, the only parent is a gluon
        if not flavored_legs:
            return [21]

        # Consider last particle
        last_leg = flavored_legs.pop()
        ll_state = last_leg.state
        ll_id    = last_leg.pdg
        # Look for a corresponding anti-particle
        for leg in range(len(flavored_legs)):
            cur_state = flavored_legs[leg].state
            cur_id    = flavored_legs[leg].pdg
            if (
                        (cur_state == ll_state and cur_id == -ll_id)
                        or
                        (cur_state != ll_state and cur_id == ll_id)
            ):
                # Eliminate it and start over
                flavored_legs.pop(leg)
                return self.parent_PDGs(flavored_legs)

        # If there was no anti-particle,
        # check if all other legs have been emitted by the last one
        if self.parent_PDGs(flavored_legs) == [21]:
            # Return the 'initial-state PDG' of the particle
            if ll_state == SubtractionLeg.INITIAL:
                return [ll_id]
            else:
                return [-ll_id]

        # At this point, there is no valid parent: return empty list
        return []
        
    def can_become_soft(self, legs):
        """Check whether a bunch of legs going simultaneously soft 
        lead to singular behavior.
        """

        for pdg in self.parent_PDGs(legs):
            particle = self.model.get_particle(pdg)
            if (particle.get('spin') == 3 and particle.get('mass') == 'zero'):
                return True
        return False
    
    def can_become_collinear(self, legs):
        """Check whether a bunch of legs going collinear to each other 
        lead to singular behavior.
        """
        
        for pdg in self.parent_PDGs(legs):
            if self.can_be_IR_unresolved(pdg):
                return True
        return False        
    
    def get_all_elementary_operators(self, process):
        """Generate all 'building blocks' operators relevant for the process
        passed in argument.
        """
        
        elementary_operator_list = SingularOperatorList([])
        
        # Eliminate particles that do not have a role in the subtraction
        legs = SubtractionLegSet(
            SubtractionLeg(leg) for leg in process.get('legs')
            if self.can_be_IR_unresolved(leg.get('id'))
        )
        fs_legs = SubtractionLegSet(
                leg for leg in legs if leg.state == SubtractionLeg.FINAL
        )
        is_legs = SubtractionLegSet(legs.difference(fs_legs))

        # Loop over number of unresolved particles
        for unresolved in range(1, self.global_order()+1):
            # Get iterators at the start of the final-state list
            it = iter(fs_legs)
            soft_it, coll_final_it, coll_initial_it = itertools.tee(it, 3)
            # Final state particle sets going soft
            for soft_set in itertools.combinations(soft_it, unresolved):
                if self.can_become_soft(soft_set):
                    elementary_operator_list.append(SoftOperator(soft_set))
            # Final state particle sets going collinear
            for coll_final_set in itertools.combinations(coll_final_it, unresolved + 1):
                if self.can_become_collinear(coll_final_set):
                    elementary_operator_list.append(CollOperator(coll_final_set))
            # Initial-final collinear
            # For any final-state set with one less particle
            for coll_initial_set in itertools.combinations(coll_initial_it, unresolved):
                # Combine with all initial legs
                for coll_initial_leg in is_legs:
                    coll_set = coll_initial_set + (coll_initial_leg, )
                    if self.can_become_collinear(coll_set):
                        elementary_operator_list.append(CollOperator(coll_set))

        return SingularOperatorList(elementary_operator_list)

    def get_all_raw_combinations(self, elementary_operators):
        """Determine all combinations of elementary operators,
        without applying any simplification.
        """

        combos = []
        for nop in range(len(elementary_operators) + 1):
            combos += [
                SingularOperatorList(combo)
                for combo in itertools.combinations(
                        iter(elementary_operators),
                        nop
                )
            ]
        return combos

    def get_all_combinations(self, elementary_operators):
        """Determine all combinations of elementary operators,
        applying simplification and discarding the ones that vanish.
        """

        return [
            simple_combo
            for simple_combo in [
                combo.simplify()
                for combo in self.get_all_raw_combinations(elementary_operators)
                ]
            if not simple_combo.is_void()
        ]

    def count_unresolved(self, structure):
        """Count the number of unresolved particles in some structure."""

        assert(type(structure) is SingularStructure)
        number_of_unresolved_particles = 0
        for sub in structure.substructures:
            if sub.name() == "S":
                number_of_unresolved_particles += len(sub.get_all_legs())
            elif sub.name() == "C":
                number_of_unresolved_particles += len(sub.get_all_legs()) - 1
            else:
                raise MadGraph5Error(
                        "Unrecognised structure of type %s"
                        "in IRSubtraction.count_unresolved" %
                        str(type(structure))
                )
        return number_of_unresolved_particles

    def filter_combinations(self, combinations):
        """Filter out combinations with too many unresolved particles
        for the order that has been set.
        """

        max_unresolved = self.global_order()
        return [
            combo
            for combo in combinations
            if self.count_unresolved(combo) <= max_unresolved
        ]

    def get_counterterm(
        self,
        structure,
        process,
        momenta_dict_so_far = None
    ):
        """Build the product of a set of currents and a matrix element
        that approximates the matrix element for process
        in the singular limit specified by structure.
        """

        # 1. Initialize variables

        assert isinstance(structure, SingularStructure)
        assert isinstance(process, base_objects.Process)

        reduced_process = process
        reduced_process['n_loops'] = None

        # If no momenta dictionary was passed
        if not momenta_dict_so_far:
            # Initialize it with process legs
            momenta_dict_so_far = bidict()
            for leg in process['legs']:
                # Check that legs are numbered progressively
                # from 1 to len(process['legs']),
                # else a more elaborate treatment of indices is needed
                assert leg['number'] == len(momenta_dict_so_far) + 1
                momenta_dict_so_far[leg['number']] = frozenset((leg['number'],))
            reduced_process = reduced_process.get_copy(['legs', 'n_loops'])

        subcurrents = []

        # 2. Recursively look into substructures

        current_args = set(structure.legs)
        for substructure in structure.substructures:
            node = self.get_counterterm(
                substructure,
                reduced_process,
                momenta_dict_so_far
            )
            current_structure = node.current['singular_structure']
            current_legs = current_structure.get_all_legs()
            current_leg_ns = frozenset(
                leg.n for leg in current_legs
            )
            # Replace collinear substructures with their parent
            if current_structure.name() == "C":
                # The parent has already been generated by recursion,
                # retrieve its number
                parent_index = momenta_dict_so_far.inv[current_leg_ns]
                # Build the SubtractionLeg that appears in the current arguments
                parent_PDGs = self.parent_PDGs(current_legs)
                assert len(parent_PDGs) == 1
                parent_PDG = parent_PDGs[0]
                parent_state = SubtractionLeg.FINAL
                if current_legs.has_initial_state_leg():
                    parent_state = SubtractionLeg.INITIAL
                current_args.add(
                    SubtractionLeg(parent_index, parent_PDG, parent_state)
                )
                # Eliminate soft sub-nodes without losing their children
                for subcurrent in node.subcurrents:
                    if subcurrent.current['singular_structure'].name() == "S":
                        node.subcurrents += subcurrent.subcurrents
                        node.subcurrents.remove(subcurrent)
            # Replace soft structures with their flattened versions
            elif current_structure.name() == "S":
                current_args.add(current_structure)
            # Other structures need to be implemented
            else:
                raise MadGraph5Error(
                    "Unrecognized current of type %s" %
                    str(type(current_structure))
                )
            # Add this node
            subcurrents.append(node)

        # If this is the outermost level,
        # the recursion was all that needed to be done
        if type(structure) is SingularStructure:
            return Counterterm(
                reduced_process,
                subcurrents,
                momenta_dict_so_far
            )

        # 3. Else build the current and update
        #    the reduced process as well as the dictionary
        current_type = type(structure)(current_args)
        current = Current({
            'singular_structure': current_type
        })
        structure_legs = current_type.get_all_legs()
        structure_leg_ns = frozenset(leg.n for leg in structure_legs)
        parent = None
        if structure.name() == "C":
            # Add entry to dictionary
            parent_index = len(momenta_dict_so_far)+1
            momenta_dict_so_far[parent_index] = structure_leg_ns
            # Work out the complete SubtractionLeg for the parent
            parent_PDGs = self.parent_PDGs(structure_legs)
            assert len(parent_PDGs) == 1
            parent_PDG = parent_PDGs[0]
            parent_state = SubtractionLeg.FINAL
            if structure_legs.has_initial_state_leg():
                parent_state = SubtractionLeg.INITIAL
            parent = SubtractionLeg(parent_index, parent_PDG, parent_state)
        elif structure.name() == "S":
            # No parent
            pass
        else:
            raise MadGraph5Error(
                "Building unrecognized current of type %s" %
                str(type(structure))
            )
        # Remove legs of this structure
        legs_to_remove = []
        for leg in reduced_process['legs']:
            if leg['number'] in structure_leg_ns:
                legs_to_remove.append(leg)
        for leg in legs_to_remove:
            reduced_process['legs'].remove(leg)
        # Add parent of this structure
        if parent:
            reduced_process['legs'].append(
                base_objects.Leg({
                    'number': parent.n,
                    'id': parent.pdg,
                    'state': parent.state
                })
            )
        # Sort preserving the initial state order
        rp_legs = reduced_process['legs']
        rp_legs.sort(key = lambda x: (x['state'], x['number']))
        if rp_legs[0]['number'] == 2:
            rp_legs[0], rp_legs[1] = rp_legs[1], rp_legs[0]

        # Finally return the counterterm node
        return CountertermNode(
            current,
            subcurrents
        )

    def split_loops(self, counterterm, n_loops):
        """Split a counterterm in several ones according to the individual
        loop orders of its currents.
        """

        assert isinstance(counterterm, CountertermNode)
        assert isinstance(n_loops, int)

        # Generate all combinations of subcurrents with total number of loops
        # less or equal to n_loops
        subcurrent_combinations = []
        subcurrents = counterterm.subcurrents
        # If this Counterterm(Node) has no subcurrents, nothing to do
        if not subcurrents:
            subcurrent_combinations = [[], ]
        # Else construct all combinations recursively
        else:
            # Initialize subcurrent_combinations with the first subcurrent
            # at any loop number between 0 and n_loops
            first_without_loops = subcurrents[0]
            for loop_n in range(0, n_loops + 1):
                for first_with_loops in self.split_loops(
                    first_without_loops, loop_n
                ):
                    subcurrent_combinations += [[first_with_loops, ], ]
            # Add the other subcurrents one by one recursively,
            # taking care of never exceeding n_loops
            n_subcurrents = len(subcurrents)
            i_current = 1
            while i_current < n_subcurrents:
                new_subcurrent_combinations = []
                for combination in subcurrent_combinations:
                    combination_loops = sum(
                        cur.n_loops()
                        for cur in combination
                    )
                    for new_loop_n in range(0, n_loops + 1 - combination_loops):
                        # Distribute the order of this subcurrent
                        # in all possible ways within the current itself
                        # (recursion step)
                        for ith_subcurrent in self.split_loops(
                            subcurrents[i_current], new_loop_n
                        ):
                            new_subcurrent_combinations.append(
                                combination + [ith_subcurrent, ]
                            )
                subcurrent_combinations = new_subcurrent_combinations
                i_current += 1

        # This current should be evaluated with the number of loops
        # that is missing to reach exactly n_loops
        result = []
        for combination in subcurrent_combinations:
            combination_loops = sum(
                cur.n_loops()
                for cur in combination
            )
            if type(counterterm) == Counterterm:
                result.append(
                    Counterterm(
                        counterterm.current.get_copy(('n_loops')),
                        combination,
                        counterterm.momenta_dict
                    )
                )
            else:
                result.append(
                    CountertermNode(
                        counterterm.current.get_copy(('n_loops')),
                        combination
                    )
                )
            result[-1].current['n_loops'] = n_loops - combination_loops

        return result

#===============================================================================
# Standalone main for debugging / standalone trials
#===============================================================================
if __name__ == '__main__':
    misc.sprint("Put your standalone subtraction code here.")
