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
import math

try:
    import madgraph
except ImportError:
    MADEVENT = True
    import internal.misc as misc
    import internal.integrands as integrands
    import internal.base_objects as base_objects
    import internal.subtraction as subtraction
    from internal import InvalidCmd, MadGraph5Error

else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.integrator.integrands as integrands
    import madgraph.core.base_objects as base_objects
    import madgraph.core.subtraction as subtraction
    from madgraph import InvalidCmd, MadGraph5Error

import copy

from madgraph.integrator.vectors import Vector, LorentzVector
from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')

#=========================================================================================
# Kinematic functions
#=========================================================================================

def Kaellen(*args):

    l = len(args)
    foo = 0.
    for i in range(l):
        foo += args[i]**2
        for j in range(i+1, l):
            foo -= 2*args[i]*args[j]
    return foo

#=========================================================================================
# Compute parent & children numbers
#=========================================================================================

def get_structure_numbers(structure, momenta_dict):
    """Return the number of the parent and children of a given structure
    according to some momenta dictionary.
    """

    legs = structure.get_all_legs()
    children = frozenset((leg.n for leg in legs))
    if structure.name() == "S":
        return None, children, None
    else:
        parent = momenta_dict.inv[children]
        is_legs = tuple(
            leg.n for leg in legs
            if leg.state == subtraction.SubtractionLeg.INITIAL )
        if not is_legs:
            return parent, children, None
        is_leg = is_legs[0]
        fs_children = frozenset((child for child in children if child != is_leg))
        return parent, fs_children, is_leg

#=========================================================================================
# Final-collinear variables
#=========================================================================================

class FinalCollinearVariables(object):

    precision_loss_message = "Precision loss detected when computing collinear variables."

    @staticmethod
    def names(parent, children):
        """Get the names of variables describing particles going unresolved."""

        names = ['s' + str(parent), ]
        for child in children:
            names += ['z' + str(child), 'p' + str(child) + '2', 'kt' + str(child), ]
        return names

    @staticmethod
    def collinear_and_reference(p):
        """Given a momentum, return normalized vectors on the light-cone."""

        n = Vector(p.space())
        n.normalize()
        return LorentzVector([1, ] + list(n)), LorentzVector([1, 0, 0, 1])

    @staticmethod
    def get(
        PS_point, children, na, nb, kinematic_variables,
        precision=1e-6 ):
        """Given unmapped momenta and reference vectors, compute the kinematic variables
        that describe the internal structure of particles going unresolved.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        # Compute the sum of momenta
        p = LorentzVector()
        for i in children:
            p += PS_point[i]
        # Pre-compute scalar products
        nap  = na.dot(p)
        nbp  = nb.dot(p)
        nanb = na.dot(nb)
        pt = p - (nbp*na + nap * nb) / nanb
        # Initialize variables for sum rules check
        zsum = 0
        ktsum = LorentzVector()
        ktabssum = LorentzVector()
        # Compute all kinematic variables
        for i in children:
            pi = PS_point[i]
            napi = na.dot(pi)
            nbpi = nb.dot(pi)
            zi = nbpi / nbp
            kti = pi - (nbpi*na+napi*nb) / nanb - zi*pt
            kinematic_variables['z'  + str(i)] = zi
            kinematic_variables['kt' + str(i)] = kti
            kinematic_variables['p'  + str(i) + '2'] = pi.square()
            zsum += zi
            ktsum += kti
            for j in range(len(kti)):
                ktabssum[j] += abs(kti[j])
        # Check numerical accuracy
        # TODO Ideally switch to quadruple precision if the check fails
        if len(children) < 2: return
        ktsum_abs = abs(ktsum.view(Vector))
        ktabssum_abs = abs(ktabssum.view(Vector))
        ktsum_ratio = ktsum_abs / ktabssum_abs
        if (abs(zsum - 1) > precision) or (ktsum_ratio > precision):
            logger.critical(FinalCollinearVariables.precision_loss_message)
            logger.critical("The sum of z's is %.16e" % zsum)
            logger.critical("The sum of kt's is %s" % str(ktsum))
            logger.critical("abs(sum(kt's)) / sum(abs(kt's)) =  %s" % ktsum_ratio)
            logger.critical("Inputs for CollinearVariables.get():")
            logger.critical("na = %s, nb = %s" % (str(na), str(nb)))
            for i in children:
                logger.critical("child %d: %s" % (i, str(PS_point[i])))
            logger.critical("Output of CollinearVariables.get():")
            logger.critical(str(kinematic_variables))
        return

    @staticmethod
    def set(
        PS_point, children, total_momentum, na, nb, kinematic_variables,
        precision=1e-6 ):
        """Given the lower multiplicity parent momentum,
        the total momentum of children and collinear variables,
        compute and set the children momenta.
        Parent and children are indices that already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        # Rename the sum of momenta
        p = total_momentum
        # Pre-compute scalar products
        nap  = na.dot(p)
        nbp  = nb.dot(p)
        nanb = na.dot(nb)
        pt = p - (nbp*na + nap*nb) / nanb
        # Variables for sums
        p_sum = LorentzVector()
        # Set momenta for all children
        for i in children:
            zi  = kinematic_variables['z'  + str(i)]
            kti = kinematic_variables['kt' + str(i)]
            pi2 = kinematic_variables['p'  + str(i) + '2']
            nbpi = zi*nbp
            pti = kti + zi*pt
            napi = (pi2-pti.square())*nanb/(2*nbpi)
            PS_point[i] = (nbpi*na+napi*nb) / nanb + pti
            p_sum += PS_point[i]
        # Check how well the parent's momentum is reproduced
        # TODO Ideally switch to quadruple precision if the check fails
        if len(children) < 2: return
        deviation = abs((p - p_sum).view(Vector))
        benchmark = abs(p.view(Vector))
        if deviation / benchmark > precision:
            logger.critical(FinalCollinearVariables.precision_loss_message)
            logger.critical("The sum of children momenta is %s" % str(p_sum))
            logger.critical("Inputs for CollinearVariables.set():")
            logger.critical("total momentum = %s" % str(p))
            logger.critical("na = %s, nb = %s" % (str(na), str(nb)))
            logger.critical("kinematic variables:")
            logger.critical(str(kinematic_variables))
            logger.critical("Output of CollinearVariables.set():")
            for i in children:
                logger.critical("child %d: %s" % (i, str(PS_point[i])))
        return

#=========================================================================================
# Initial-collinear variables
#=========================================================================================

class InitialCollinearVariables(object):

    @staticmethod
    def names(parent, fs_children, is_child):
        """Get the names of variables describing particles going unresolved."""

        names = ['z' + str(parent), ]
        for child in fs_children:
            names += ['z' + str(child), 'p' + str(child) + '2', 'kt' + str(child), ]
        return names

    @staticmethod
    def collinear_and_reference(p):
        """Given a momentum, return normalized vectors on the light-cone."""

        # In this case taking the anti-collinear direction as a reference poses no risks,
        # because the direction of the incoming parton is fixed

        n = Vector(p.space())
        n.normalize()
        # For good phase space points, p[0] >= 0, but test PS points might not be valid
        if p[0] >= 0:
            return LorentzVector([1, ] + list(+n)), LorentzVector([1, ] + list(-n))
        else:
            return LorentzVector([1, ] + list(-n)), LorentzVector([1, ] + list(+n))

    @staticmethod
    def get(
        PS_point, fs_children, is_child, na, nb, kinematic_variables,
        precision=1e-6 ):
        """Given unmapped momenta and reference vectors, compute the kinematic variables
        that describe the internal structure of particles going unresolved.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        pa = PS_point[is_child]
        # Compute the sum of momenta
        pA = LorentzVector(pa)
        for i in fs_children:
            pA -= PS_point[i]
        # Pre-compute variables
        napA = na.dot(pA)
        nbpA = nb.dot(pA)
        nanb = na.dot(nb)
        nbpa = nb.dot(pa)
        zA = nbpA / nbpa
        ktA = pA - (nbpA * na + napA * nb) / nanb
        # ptA = pA - (nbpA * na + napA * nb) / nanb
        # ktA = ptA / zA
        # Initialize variables for sum rules check
        zsum = 0
        ktsum = LorentzVector()
        ktabssum = LorentzVector()
        # Fill in A data, using child number improperly
        kinematic_variables['z'  + str(is_child)] = zA
        kinematic_variables['kt' + str(is_child)] = ktA
        kinematic_variables['s'  + str(is_child)] = pA.square()
        zsum += zA
        ktsum += ktA
        for j in range(len(ktA)):
            ktabssum[j] += abs(ktA[j])
        # Compute all kinematic variables
        for i in fs_children:
            pi = PS_point[i]
            napi = na.dot(pi)
            nbpi = nb.dot(pi)
            zi = nbpi / nbpa
            kti = pi - (nbpi*na+napi*nb) / nanb
            # pti = pi - (nbpi*na+napi*nb) / nanb
            # kti = pti - zi * ktA
            kinematic_variables['z'  + str(i)] = zi
            kinematic_variables['kt' + str(i)] = kti
            kinematic_variables['p'  + str(i) + '2'] = pi.square()
            zsum += zi
            ktsum += kti
            for j in range(len(kti)):
                ktabssum[j] += abs(kti[j])
        # Check numerical accuracy
        # TODO Ideally switch to quadruple precision if the check fails
        if not fs_children: return
        ktsum_abs = abs(ktsum.view(Vector))
        ktabssum_abs = abs(ktabssum.view(Vector))
        ktsum_ratio = ktsum_abs / ktabssum_abs
        if (abs(zsum - 1) > precision) or (ktsum_ratio > precision):
            logger.critical(FinalCollinearVariables.precision_loss_message)
            logger.critical("The sum of z's is %.16e" % zsum)
            logger.critical("The sum of kt's is %s" % str(ktsum))
            logger.critical("abs(sum(kt's)) / sum(abs(kt's)) =  %s" % ktsum_ratio)
            logger.critical("Inputs for InitialCollinearVariables.get():")
            logger.critical("na, nb = %s, %s" % (str(na), str(nb)))
            for i in fs_children:
                logger.critical("fs_child %d: %s" % (i, str(PS_point[i])))
            logger.critical("is_child %d: %s" % (is_child, str(PS_point[is_child])))
            logger.critical("Output of InitialCollinearVariables.get():")
            logger.critical(str(kinematic_variables))
        return

    @staticmethod
    def set(
        PS_point, fs_children, is_child, pa, na, nb, kinematic_variables,
        precision=1e-6 ):
        """Given the lower multiplicity momentum of the incoming parton
        and collinear variables compute and set the children momenta.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        PS_point[is_child] = pa
        # Get A data
        zA  = kinematic_variables['z'  + str(is_child)]
        ktA = kinematic_variables['kt' + str(is_child)]
        pA2 = kinematic_variables['s'  + str(is_child)]
        # Pre-compute variables
        nanb = na.dot(nb)
        nbpa = nb.dot(pa)
        ptA = ktA
        # ptA = zA * ktA
        nbpA = zA * nbpa
        napA = (pA2 - ptA.square()) * nanb / (2*nbpA)
        pA = (nbpA*na + napA*nb) / nanb + ptA
        # Variables for sums
        p_sum = LorentzVector(pa)
        # Set momenta for all children
        for i in fs_children:
            zi  = kinematic_variables['z'  + str(i)]
            kti = kinematic_variables['kt' + str(i)]
            pi2 = kinematic_variables['p'  + str(i) + '2']
            pti = kti
            # pti = kti + zi * ktA
            nbpi = zi * nbpa
            napi = (pi2 - pti.square()) * nanb / (2 * nbpi)
            PS_point[i] = (nbpi*na + napi*nb) / nanb + pti
            p_sum -= PS_point[i]
        # Check how well the parent's momentum is reproduced
        # TODO Ideally switch to quadruple precision if the check fails
        if not fs_children: return
        deviation = abs((pA - p_sum).view(Vector))
        benchmark = abs(pA.view(Vector))
        if deviation / benchmark > precision:
            logger.critical(FinalCollinearVariables.precision_loss_message)
            logger.critical("The sum of children momenta is %s" % str(p_sum))
            logger.critical("vs the total: %s" % str(pA))
            logger.critical("Inputs for InitialCollinearVariables.set():")
            logger.critical("pa = %s" % str(pa))
            logger.critical("na = %s, nb = %s" % (str(na), str(nb)))
            logger.critical("kinematic variables:")
            logger.critical(str(kinematic_variables))
            logger.critical("Output of InitialCollinearVariables.set():")
            for i in fs_children:
                logger.critical("fs_child %d: %s" % (i, str(PS_point[i])))
        return

#=========================================================================================
# Soft variables
#=========================================================================================

class SoftVariables(object):

    @staticmethod
    def names(children):
        """Get the names of variables describing particles going unresolved."""

        return ['p' + str(child) for child in children]

    @staticmethod
    def get(PS_point, children, kinematic_variables):
        """Given unmapped and mapped momenta, compute the kinematic variables
        that describe the internal structure of particles going unresolved.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        """

        # For soft particles, just pass the whole momentum
        for i in children:
            kinematic_variables['p' + str(i)] = PS_point[i].get_copy()
        return

    @staticmethod
    def set(PS_point, children, kinematic_variables):
        """Given a dictionary of variables that describe the unresolved partons,
        compute and set the children momenta.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        """

        # Set momenta for all children
        for i in children:
            PS_point[i] = kinematic_variables['p' + str(i)]
        return

#=========================================================================================
# VirtualMapping
#=========================================================================================

class VirtualMapping(object):
    """Base class for elementary mapping implementations."""

    # TODO Add a method is_valid_recoiler?

    @classmethod
    def is_valid_structure(cls, singular_structure):
        """Return true if it makes sense to apply this mapping
        to the given singular structure.
        """

        raise NotImplemented

    @classmethod
    def get_kinematic_variables_names(cls, singular_structure, momenta_dict):
        """For a valid singular structure,
        returns a list of variable names necessary to apply this mapping.
        """

        raise NotImplemented

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None ):
        """Map a given phase-space point to lower multiplicity,
        by clustering the substructures and recoiling against the legs
        specified in singular_structure.

        :param PS_point: higher-multiplicity phase-space point,
        as a dictionary that associates integers to Lorentz vectors,
        which will be modified to the lower-multiplicity one

        :param singular_structure: SingularStructure object that specifies
        sets of unresolved particles and recoilers recursively

        :param momenta_dict: two-way dictionary that associates a unique label
        to each set of one or more unresolved particles identified by their number

        :param kinematic_variables: if a non-empty dictionary is passed,
        the kinematic variables that are necessary to reproduce the higher-multiplicity
        phase-space point from the lower-multiplicity one will be set

        :param compute_jacobian: if False, will not compute the jacobian for the mapping

        :param masses: masses of the parents of the collinear sets,
        as a dictionary with (number of parent, mass)

        :return: dictionary containing the jacobian weight due to the mapping
                 and eventually other characteristic variables of the mapping
        """
        
        raise NotImplemented

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):
        """Map a given phase-space point to higher multiplicity,
        by splitting the (pseudo)particles and recoiling against the legs
        specified in singular_structure.

        :param PS_point: lower-multiplicity phase-space point,
        as a dictionary that associates integers to Lorentz vectors,
        which will be modified to the higher-multiplicity one

        :param singular_structure: SingularStructure object that specifies
        sets of unresolved particles and recoilers recursively

        :param momenta_dict: two-way dictionary that associates a unique label
        to each set of one or more unresolved particles identified by their number

        :param kinematic_variables: variables describing the splitting,
        as a dictionary that associates variable names to values

        :param compute_jacobian: if False, will not compute the jacobian for the mapping

        :return: dictionary containing the jacobian weight due to the mapping
                 and eventually other characteristic variables of the mapping
        """
        
        raise NotImplemented


    @classmethod
    def rescale_kinematic_variables(
        cls, singular_structure, momenta_dict, kinematic_variables, scaling_parameter ):
        """Rescale in-place the given kinematic variables so as to approach the limit.

        :param singular_structure: SingularStructure object that specifies
        sets of unresolved particles and recoilers recursively

        :param momenta_dict: two-way dictionary that associates a unique label
        to each set of one or more unresolved particles identified by their number

        :param kinematic_variables: variables describing the splitting,
        as a dictionary that associates variable names to values

        :param scaling_parameter: a parameter from 0 to 1 that indicates the distance
        from the singular limit

        :return: the modified kinematic variables
        """

        raise NotImplemented

#=========================================================================================
# Final-collinear mappings
#=========================================================================================

# Generic final-collinear mapping parent class
#=========================================================================================

class FinalCollinearMapping(VirtualMapping):
    """Common functions for final-collinear elementary mappings."""

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, subtraction.SingularStructure)
        # Valid only for sets of final-state particles going collinear,
        # with no recursive substructure
        for substructure in singular_structure.substructures:
            if not substructure.name() == "C":
                return False
            if substructure.substructures:
                return False
            if substructure.get_all_legs().has_initial_state_leg():
                return False
        return True

    @classmethod
    def get_kinematic_variables_names(cls, singular_structure, momenta_dict):
        """Get the names of variables describing particles going unresolved."""

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        names = []
        for substructure in singular_structure.substructures:
            parent, children, _ = get_structure_numbers(substructure, momenta_dict)
            FinalCollinearVariables.names(parent, children)
        return names


    @classmethod
    def rescale_kinematic_variables(
        cls, singular_structure, momenta_dict, kinematic_variables, scaling_parameter):

        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict))
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        # Determine the correct scaling for the divergence to go like 1/parameter
        base = scaling_parameter ** (0.5 * (len(children) - 1))
        kinematic_variables['s' + str(parent)] *= base ** 2
        for child in children:
            kinematic_variables['kt' + str(child)] *= base
        return kinematic_variables

# Final-collinear rescaling mapping, one set
#=========================================================================================

class FinalRescalingMappingOne(FinalCollinearMapping):
    """Implementation of the rescaling mapping
    for one bunch of collinear particles (with massless parent)
    and an arbitrary number of massless recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one bunch of final-state particles going collinear,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(FinalRescalingMappingOne, cls).is_valid_structure(
                singular_structure )
        return False

    @staticmethod
    def alpha(pC, Q):
        """Return the parameter alpha of the rescaling mapping."""

        Q2  = Q.square()
        pC2 = pC.square()
        QpC = pC.dot(Q)
        return (QpC - math.sqrt(QpC**2 - Q2*pC2)) / Q2

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        # Build collective momenta
        pC = LorentzVector()
        for j in children:
            pC += PS_point[j]
        pR = LorentzVector()
        for leg in singular_structure.legs:
            pR += PS_point[leg.n]
        Q = pC + pR
        # Compute the parameter alpha
        alpha = FinalRescalingMappingOne.alpha(pC, Q)
        # Map all recoilers' momenta
        for recoiler in recoilers:
            PS_point[recoiler] /= (1-alpha)
        # Map the set's momentum
        PS_point[parent] = (pC - alpha * Q) / (1-alpha)
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = FinalCollinearVariables.collinear_and_reference(PS_point[parent])
            FinalCollinearVariables.get(PS_point, children, na, nb, kinematic_variables)
            kinematic_variables['s' + str(parent)] = pC.square()
        # Eliminate children momenta from the mapped phase-space point
        for j in children:
            if j != parent: # Bypass degenerate case of 1->1 splitting
                del PS_point[j]
        # Compute the jacobian for this mapping
        jacobian = (1-alpha)**(2*(len(recoilers)-1))
        # Return characteristic variables
        return {
            'jacobian':           jacobian,
            'alpha'+str(parent):  alpha,
            'pC'+str(parent):     pC,
            'Q':                  Q, }

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        # Build collective momenta
        qC = PS_point[parent]
        qR = LorentzVector()
        for leg in singular_structure.legs:
            qR += PS_point[leg.n]
        Q = qR + qC
        # Compute scalar products
        assert abs(qC.square()) < math.sqrt(qC.eps())
        qR2 = qR.square()
        QqC = Q.dot(qC)
        sC  = kinematic_variables['s' + str(parent)]
        # Obtain parameter alpha
        alpha = 0
        if (qR2*sC)/(QqC**2) < 1e-6:
            alpha = 0.5 * sC / QqC
        else:
            alpha = (math.sqrt(QqC**2 + sC*qR2) - QqC) / qR2
        # Compute reverse-mapped momentum
        pC = (1-alpha) * qC + alpha * Q
        # Map recoil momenta
        for recoiler in recoilers:
            PS_point[recoiler] *= 1-alpha
        # Set children momenta
        na, nb = FinalCollinearVariables.collinear_and_reference(qC)
        FinalCollinearVariables.set(PS_point, children, pC, na, nb, kinematic_variables)
        # Remove parent's momentum
        if parent not in children: # Bypass degenerate case of 1->1 splitting
            del PS_point[parent]
        # Return the jacobian for this mapping
        jacobian = (1-alpha)**(2*(len(recoilers)-1))
        return {
            'jacobian':           jacobian,
            'alpha'+str(parent):  alpha,
            'pC'+str(parent):     pC,
            'Q':                  Q, }

# Final-collinear Lorentz mapping, one set
#=========================================================================================

class FinalLorentzMappingOne(FinalCollinearMapping):
    """Implementation of the Lorentz transformation mapping
    for one bunch of collinear particles (with massless parent)
    and an arbitrary number of (eventually massive) recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one bunch of final-state particles going collinear,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(FinalLorentzMappingOne, cls).is_valid_structure(
                singular_structure )
        return False

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        # Build collective momenta
        pC = LorentzVector()
        for j in children:
            pC += PS_point[j]
        pR = LorentzVector()
        for leg in singular_structure.legs:
            pR += PS_point[leg.n]
        Q = pC + pR
        # Compute parameters
        Q2  = Q.square()
        pC2 = pC.square()
        pR2 = pR.square()
        pC_perp = pC - ((Q2+pC2-pR2)/(2*Q2)) * Q
        alpha = (Q2-pR2) / math.sqrt(Kaellen(Q2, pR2, pC2))
        # Map the set's momentum
        PS_point[parent] = alpha*pC_perp + ((Q2-pR2)/(2*Q2))*Q
        # Map all recoilers' momenta
        qR = Q - PS_point[parent]
        for recoiler in recoilers:
            # TODO Move this try/except to higher level
            try:
                PS_point[recoiler].rotoboost(pR, qR)
            except:
                logger.critical("Problem encountered for %s" % str(singular_structure))
                logger.critical("The full phase space point was\n%s" % str(PS_point))
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = FinalCollinearVariables.collinear_and_reference(PS_point[parent])
            kinematic_variables['s' + str(parent)] = pC2
            FinalCollinearVariables.get(PS_point, children, na, nb, kinematic_variables)
        # Eliminate children momenta from the mapped phase-space point
        for j in children:
            if j != parent: # Bypass degenerate case of 1->1 splitting
                del PS_point[j]
        jacobian = alpha**(-3) * (Q2-pR2)/(Q2+pC2-pR2)
        return {
            'jacobian':           jacobian,
            'pC'+str(parent):     pC,
            'Q':                  Q, }

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        # Build collective momenta
        qC = PS_point[parent]
        qR = LorentzVector()
        for recoiler in recoilers:
            qR += PS_point[recoiler]
        Q = qR + qC
        # Compute scalar products
        pC2 = kinematic_variables['s' + str(parent)]
        # Compute parameters
        assert abs(qC.square()) < math.sqrt(qC.eps())
        Q2  = Q.square()
        qR2 = qR.square()
        qC_perp = qC - ((Q2-qR2)/(2*Q2)) * Q
        alpham1 = math.sqrt(Kaellen(Q2, qR2, pC2)) / (Q2-qR2)
        # Compute reverse-mapped momentum
        pC = alpham1*qC_perp + ((Q2+pC2-qR2)/(2*Q2))*Q
        pR = Q - pC
        # Map recoil momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n].rotoboost(qR, pR)
        # Set children momenta
        na, nb = FinalCollinearVariables.collinear_and_reference(qC)
        FinalCollinearVariables.set(PS_point, children, pC, na, nb, kinematic_variables)
        # Remove parent's momentum
        if parent not in children: # Bypass degenerate case of 1->1 splitting
            del PS_point[parent]
        jacobian = alpham1**3 * (Q2-qR2)/(Q2+pC2-qR2)
        return {
            'jacobian':           jacobian,
            'pC'+str(parent):     pC,
            'Q':                  Q, }

#=========================================================================================
# Initial-collinear mappings
#=========================================================================================

# Generic initial-collinear mapping parent class
#=========================================================================================

class InitialCollinearMapping(VirtualMapping):
    """Common functions for initial-collinear elementary mappings."""

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, subtraction.SingularStructure)
        # Valid only for sets of initial-state particles going collinear,
        # with no recursive substructure
        for substructure in singular_structure.substructures:
            if not substructure.name() == "C":
                return False
            if substructure.substructures:
                return False
            if not substructure.get_all_legs().has_initial_state_leg():
                return False
        # There cannot be more than two collinear sets with initial-state particles
        if len(singular_structure.substructures) > 2:
            return False
        return True

    @classmethod
    def get_kinematic_variables_names(cls, singular_structure, momenta_dict):
        """Get the names of variables describing particles going unresolved."""

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        names = []
        for substructure in singular_structure.substructures:
            parent, fs_children, is_child = get_structure_numbers(
                substructure, momenta_dict )
            InitialCollinearVariables.names(parent, fs_children, is_child)
        return names

    @classmethod
    def rescale_kinematic_variables(
        cls, singular_structure, momenta_dict, kinematic_variables, scaling_parameter):

        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict))
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        _, fs_children, is_child = get_structure_numbers(substructure, momenta_dict)
        # Determine the correct scaling for the divergence to go like 1/parameter
        base = scaling_parameter ** (0.5 * (len(fs_children) - 1))
        kinematic_variables['s' + str(is_child)] *= base ** 2
        kinematic_variables['kt' + str(is_child)] *= base
        for child in fs_children:
            kinematic_variables['kt' + str(child)] *= base
        return kinematic_variables

# Initial-collinear Lorentz mapping, one set
#=========================================================================================

class InitialLorentzMappingOne(InitialCollinearMapping):
    """Implementation of the Lorentz transformation mapping
    for one set of collinear particles (with massless parent)
    and an arbitrary number of (eventually massive) recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one set of particles going collinear to an initial-state parton,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(InitialLorentzMappingOne, cls).is_valid_structure(
                singular_structure )
        return False

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, fs_children, is_child = get_structure_numbers(substructure, momenta_dict)
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        # Build collective momenta
        pCa = LorentzVector()
        for j in fs_children:
            pCa += PS_point[j]
        pR = LorentzVector()
        for recoiler in recoilers:
            pR += PS_point[recoiler]
        pa = PS_point[is_child]
        pA = pa - pCa
        pAmpR = pA - pR
        # Compute parameters
        xia = (pAmpR.square() - pR.square())/(2*pa.dot(pAmpR))
        # Map the set's momentum
        qA = xia * pa
        PS_point[parent] = qA
        # Map all recoilers' momenta
        qR = qA - pAmpR
        for recoiler in singular_structure.legs:
            # TODO Move this try/except to higher level
            try:
                PS_point[recoiler.n].rotoboost(pR, qR)
            except:
                logger.critical("Problem encountered for %s" % str(singular_structure))
                logger.critical("The full phase space point was\n%s" % str(PS_point))
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = InitialCollinearVariables.collinear_and_reference(qA)
            InitialCollinearVariables.get(
                PS_point, fs_children, is_child, na, nb, kinematic_variables )
        # Eliminate children momenta from the mapped phase-space point
        for j in fs_children:
            del PS_point[j]
        if is_child != parent:  # Bypass degenerate case of 1->1 splitting
            del is_child
        # TODO Compute the jacobian for this mapping
        jacobian = 1.0
        return {
            'jacobian':           jacobian,
            'pC'+str(parent):     pA,
            'Q':                  pAmpR, }

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False):

        # print "Mapping up the PS point:\n", PS_point
        # print "With variables:\n", kinematic_variables
        # print "and structure", singular_structure

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict))
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, fs_children, is_child = get_structure_numbers(substructure, momenta_dict)
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        # Build collective momenta
        qA = PS_point[parent]
        na, nb = InitialCollinearVariables.collinear_and_reference(qA)
        nanb = na.dot(nb)
        qR = LorentzVector()
        for recoiler in recoilers:
            qR += PS_point[recoiler]
        qRmqA = qR - qA
        zA = kinematic_variables['z' + str(is_child)]
        ktA = kinematic_variables['kt' + str(is_child)]
        pA2 = kinematic_variables['s' + str(is_child)]
        ptA = ktA
        # ptA = zA * ktA
        quad_a = na.dot(qRmqA)
        quad_b = 2*ptA.dot(qRmqA) + pA2 + qRmqA.square() - qR.square()
        quad_c = (pA2 - ptA.square()) * 2*nb.dot(qRmqA) / nanb
        sqrt_delta = math.sqrt(quad_b**2 - 4*quad_a*quad_c)
        # WARNING
        # The plus sign in front of sqrt_delta is the correct one
        # for > 90% of physical phase-space points, but not all the times.
        # It seems like the minus sign has to bo chosen _mostly_
        # when there is a single recoiler and no unchanged partons.
        # To reproduce the issue, within test_mappings.py, go to test_invertible,
        # set random.seed(7) and skip the first 6 PS points.
        pAp = (-quad_b+sqrt_delta)/(2*quad_a)
        xia = 2*zA*nb.dot(qA)/(pAp*nanb)
        # xia = kinematic_variables['xi' + str(parent)]
        # Compute parameters
        pa = qA / xia
        # Set children momenta
        InitialCollinearVariables.set(
            PS_point, fs_children, is_child, pa, na, nb, kinematic_variables )
        # Build collective momenta
        pCa = LorentzVector()
        for j in fs_children:
            pCa += PS_point[j]
        pA = pa - pCa
        # Map recoil momenta
        pR = qR + pA - qA
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n].rotoboost(qR, pR)
        # Remove parent's momentum
        if parent != is_child:  # Bypass degenerate case of 1->1 splitting
            del PS_point[parent]
        # TODO Compute the jacobian for this mapping
        jacobian = 1.0
        return {
            'jacobian':           jacobian,
            'pC'+str(parent):     pA,
            'Q':                  pA - pR, }

#=========================================================================================
# Soft mappings
#=========================================================================================

# Generic soft mapping parent class
#=========================================================================================

class ElementaryMappingSoft(VirtualMapping):
    """Common functions for soft elementary mappings."""

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, subtraction.SingularStructure)
        # Valid only for bunches of final-state particles going soft,
        # with no recursive substructure
        for substructure in singular_structure.substructures:
            if not substructure.name() == "S":
                return False
            if substructure.substructures:
                return False
        return True

    @classmethod
    def get_kinematic_variables_names(cls, singular_structure, momenta_dict):
        """Get the names of variables describing particles going unresolved."""

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        # For every soft particle, just return its momentum
        names = []
        for substructure in singular_structure.substructures:
            _, children, _ = get_structure_numbers(substructure, momenta_dict)
            for legn in sorted(children):
                names += ['p' + str(legn), ]
        return names

    @classmethod
    def rescale_kinematic_variables(
        cls, singular_structure, momenta_dict, kinematic_variables, scaling_parameter):

        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict))
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        _, children, _ = get_structure_numbers(substructure, momenta_dict)
        # Determine the correct scaling for the divergence to go like 1/parameter
        base = scaling_parameter ** (0.5 * len(children))
        for child in children:
            kinematic_variables['p' + str(child)] *= base
        return kinematic_variables

# Soft mapping, massless final-state recoilers
#=========================================================================================

class MappingSomogyietalSoft(ElementaryMappingSoft):
    """Implementation of the mapping used by Somogyi et al.
    in arXiv:hep-ph/0609042 for soft particles.
    It applies when there are at least two recoilers in the final state
    and all recoilers are massless.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only if there are at least two recoilers (assumed in the final state)
        if len(singular_structure.legs) < 2: return False
        return super(MappingSomogyietalSoft, cls).is_valid_structure(singular_structure)

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        # Build the total soft momentum,
        # save the soft momenta in variables and eliminate them from PS_point
        pS = LorentzVector()
        for substructure in singular_structure.substructures:
            children = tuple(leg.n for leg in substructure.legs)
            if kinematic_variables is not None:
                SoftVariables.get(PS_point, children, kinematic_variables)
            for child in children:
                pS += PS_point.pop(child)
        # Build the total momentum of recoilers
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        pR = LorentzVector()
        for recoiler in recoilers:
            pR += PS_point[recoiler]
        # Build the total momentum Q
        Q = pS + pR
        # Compute the parameter la
        pR2_Q2 = pR.square() / Q.square()
        y = 1. - pR2_Q2
        la = math.sqrt(pR2_Q2)
        P = pR / la
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] /= la
            PS_point[recoiler.n].rotoboost(P, Q)
        jacobian = (pR2_Q2)**(len(recoilers)-2)
        return {
            'jacobian': jacobian,
            'y':        y,
            'pS':       pS,
            'Q':        Q, }

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict ) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Build the total soft momentum,
        # get the soft momenta from variables and save them in PS_point
        pS = LorentzVector()
        for substructure in singular_structure.substructures:
            children = tuple(leg.n for leg in substructure.legs)
            SoftVariables.set(PS_point, children, kinematic_variables)
            for child in children:
                pS += PS_point[child]
        # Build the total momentum, which is equal to the mapped recoilers'
        Q = LorentzVector()
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        for recoiler in recoilers:
            Q += PS_point[recoiler]
        # Build the recoilers' momentum
        pR = Q - pS
        # Compute the parameter la
        pR2_Q2 = pR.square() / Q.square()
        y = 1. - pR2_Q2
        la = math.sqrt(pR2_Q2)
        P = pR / la
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] *= la
            PS_point[recoiler.n].rotoboost(Q, P)
        jacobian = (pR2_Q2)**(len(recoilers)-2)
        return {
            'jacobian': jacobian,
            'y':        y,
            'pS':       pS,
            'Q':        Q, }

#=========================================================================================
# Soft-collinear mappings
#=========================================================================================

# Generic soft-collinear mapping
#=========================================================================================

class SoftCollinearMapping(VirtualMapping):
    """Proxy to dispatch soft-collinear mappings."""

    def __init__(self, soft_mapping, collinear_mapping):

        self.soft_mapping = soft_mapping
        self.coll_mapping = collinear_mapping
        super(SoftCollinearMapping, self).__init__()

    @staticmethod
    def coll_structure(structure):

        all_collinear_structures = [
            subtraction.CollStructure(legs=substructure.legs)
            for substructure in structure.substructures ]
        return subtraction.SingularStructure(
            substructures=all_collinear_structures,
            legs=structure.legs )

    @staticmethod
    def soft_structure(structure):

        all_soft_structures = []
        all_soft_recoils    = copy.copy(structure.legs)
        for substructure in structure.substructures:
            all_soft_structures += substructure.substructures
            all_soft_recoils    += substructure.legs
        return subtraction.SingularStructure(
            substructures=all_soft_structures,
            legs=all_soft_recoils )

    @staticmethod
    def coll_momenta_dict(structure, momenta_dict):

        coll_momenta_dict = subtraction.bidict()
        for substructure in structure.substructures:
            coll_ns = frozenset(leg.n for leg in substructure.legs)
            parent, _, _ = get_structure_numbers(substructure, momenta_dict)
            coll_momenta_dict[parent] = coll_ns
        return coll_momenta_dict

    def is_valid_structure(
        self, singular_structure,
        i_soft_structure=None, i_coll_structure=None ):

        assert isinstance(singular_structure, subtraction.SingularStructure)
        # Check collinear structure
        if i_coll_structure: coll_structure = i_coll_structure
        else: coll_structure = self.coll_structure(singular_structure)
        if not self.coll_mapping.is_valid_structure(coll_structure):
            return False
        # Check soft structure
        if i_soft_structure: soft_structure = i_soft_structure
        else: soft_structure = self.soft_structure(singular_structure)
        return self.soft_mapping.is_valid_structure(soft_structure)

    def get_kinematic_variables_names(
        self, singular_structure, momenta_dict,
        i_soft_structure=None, i_coll_structure=None ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        if i_coll_structure: coll_structure = i_coll_structure
        else: coll_structure = self.coll_structure(singular_structure)
        if i_soft_structure: soft_structure = i_soft_structure
        else: soft_structure = self.soft_structure(singular_structure)
        assert self.is_valid_structure(
            singular_structure,
            i_soft_structure=soft_structure, i_coll_structure=coll_structure )

        # For every soft particle, just return its momentum
        soft_names = self.soft_mapping.get_kinematic_variables_names(soft_structure)
        coll_names = self.coll_mapping.get_kinematic_variables_names(coll_structure)
        return soft_names + coll_names

    def map_to_lower_multiplicity(
        self, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        coll_structure = self.coll_structure(singular_structure)
        soft_structure = self.soft_structure(singular_structure)
        assert self.is_valid_structure(
            singular_structure,
            i_soft_structure=soft_structure, i_coll_structure=coll_structure )

        soft_result = self.soft_mapping.map_to_lower_multiplicity(
            PS_point, soft_structure, momenta_dict,
            kinematic_variables, compute_jacobian, masses )
        coll_momenta_dict = self.coll_momenta_dict(singular_structure, momenta_dict)
        coll_result = self.coll_mapping.map_to_lower_multiplicity(
            PS_point, coll_structure, coll_momenta_dict,
            kinematic_variables, compute_jacobian, masses )
        soft_result['jacobian'] *= coll_result['jacobian']
        return soft_result

    def map_to_higher_multiplicity(
        self, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        coll_structure = self.coll_structure(singular_structure)
        soft_structure = self.soft_structure(singular_structure)
        assert self.is_valid_structure(
            singular_structure,
            i_soft_structure=soft_structure, i_coll_structure=coll_structure )

        coll_momenta_dict = self.coll_momenta_dict(singular_structure, momenta_dict)
        coll_result = self.coll_mapping.map_to_higher_multiplicity(
            PS_point, coll_structure, coll_momenta_dict,
            kinematic_variables, compute_jacobian )
        soft_result = self.soft_mapping.map_to_higher_multiplicity(
            PS_point, soft_structure, momenta_dict,
            kinematic_variables, compute_jacobian )
        soft_result['jacobian'] *= coll_result['jacobian']
        return soft_result

#=========================================================================================
# Mapping walkers
#=========================================================================================

# Stroll
#=========================================================================================

class Stroll(object):
    """Container for a mapping call."""

    def __init__(self, mapping, structure, current):

        super(Stroll, self).__init__()
        assert isinstance(mapping, VirtualMapping)
        self.mapping = mapping
        assert isinstance(structure, subtraction.SingularStructure)
        self.structure = structure
        assert isinstance(current, subtraction.Current)
        self.current = current

    def __str__(self):

        foo = self.mapping.__class__.__name__
        foo += " with structure " + str(self.structure)
        foo += " for current " + str(self.current)
        return foo

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

# Generic walker parent class
#=========================================================================================

class VirtualWalker(object):
    """Base class for walker implementations."""
    
    def __new__(cls, **opts):
        """Factory class to make plugin easy."""

        if cls is VirtualWalker:
            map_type = opts.pop('map_type') if 'map_type' in opts else 'Unknown'
            if (map_type not in mapping_walker_classes_map
                or not mapping_walker_classes_map[map_type] ):
                raise MadGraph5Error(
                    "Unknown mapping walker of type '%s'." % str(map_type) )
            target_class = mapping_walker_classes_map[map_type]
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
            str1 = ss.__str__(True, True, True)
            for ct in counterterms:
                if len(ct.nodes) != 1: continue
                if ct.nodes[0].nodes: continue
                str2 = ct.nodes[0].current['singular_structure'].__str__(True, True, True)
                if str1 != str2: continue
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
        self, PS_point, counterterm, compute_kinematic_variables=False, verbose=False ):
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

        :return: a dictionary with the following entries:
        'currents', a list of all currents that need to be evaluated,
            paired with the momenta needed for their evaluation;
        'matrix_element', the reduced matrix element,
            paired with its corresponding reduced phase-space point;
        'mapping_variables', the cumulative variables of all mappings that were applied;
        'kinematic_variables', a dictionary of all variables needed to recover
            the starting phase-space point from the lowest multiplicity one,
            or None if such variables were not requested;
        'resulting_PS_point', the final lowest multiplicity PS point.
        """

        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        if verbose: print point
        # Initialize return variables
        current_PS_pairs = []
        mapping_variables = {'jacobian': 1.}
        kinematic_variables = dict() if compute_kinematic_variables else None
        # Determine the hike
        hike = self.get_hike(counterterm)
        # Walk along the hike
        for stroll in hike:
            # Deep copy the momenta
            momenta = point.get_copy()
            # Map to lower multiplicity
            result = stroll.mapping.map_to_lower_multiplicity(
                point, stroll.structure, counterterm.momenta_dict, kinematic_variables )
            # Transcribe parents' momenta
            for key in point.keys():
                if not momenta.has_key(key):
                    momenta[key] = LorentzVector(point[key])
            # Update jacobian and mapping variables
            mapping_variables['jacobian'] *= result.pop('jacobian')
            mapping_variables.update(result)
            # Append the current and the momenta
            current_PS_pairs.append((stroll.current, momenta))
        # Identify reduced matrix element,
        # computed in the point which has received all mappings
        ME_PS_pair = [counterterm.process, point]
        if verbose: print point
        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'mapping_variables': mapping_variables,
            'kinematic_variables': kinematic_variables,
            'resulting_PS_point': point}

    def walk_to_higher_multiplicity(
        self, PS_point, counterterm, kinematic_variables, verbose=False ):
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

        :return: a dictionary with the following entries:
        'currents', a list of all currents that need to be evaluated,
            paired with the momenta needed for their evaluation;
        'matrix_element', the reduced matrix element,
            paired with its corresponding reduced phase-space point;
        'mapping_variables', the cumulative variables of all mappings that were applied;
        'resulting_PS_point', the final lowest multiplicity PS point.
        """

        # Identify reduced matrix element,
        # computed in the lowest multiplicity point
        ME_PS_pair = [counterterm.process, PS_point]
        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        if verbose: print point
        # Initialize return variables
        current_PS_pairs = []
        mapping_variables = {'jacobian': 1}
        # Determine the hike
        hike = self.get_hike(counterterm)
        # Walk the hike backwards
        prev_point = PS_point
        prev_parents = []
        for stroll in reversed(hike):
            # Compute jacobian and map to lower multiplicity
            result = stroll.mapping.map_to_higher_multiplicity(
                point, stroll.structure, counterterm.momenta_dict, kinematic_variables )
            if verbose: print point
            # Update jacobian and mapping variables
            mapping_variables['jacobian'] *= result.pop('jacobian')
            mapping_variables.update(result)
            # Prepend pair of this current and the momenta,
            # deep copy wanted
            momenta = point.get_copy()
            # Transcribe parents' momenta
            parents = []
            for key in prev_point.keys():
                if not momenta.has_key(key) and key not in prev_parents:
                    parents.append(key)
                    momenta[key] = LorentzVector(prev_point[key])
            current_PS_pairs.insert(0, (stroll.current, momenta))
            prev_point = momenta
            prev_parents = parents
        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'mapping_variables': mapping_variables,
            'resulting_PS_point': point}

    def approach_limit(
        self, PS_point, counterterm, scaling_parameter, all_counterterms ):
        """Produce a higher multiplicity phase-space point from PS_point,
        according to kinematic_variables that approach the limit of counterterm
        parametrically with scaling_parameter.
        """

        # Decompose the counterterm
        decomposed = self.decompose_counterterm(counterterm, all_counterterms)
        # Always approach the limit at the same speed
        base = scaling_parameter ** (1. / len(decomposed))
        # Walk the hike up and down
        for ct in decomposed:
            hike = self.determine_hike(ct)
            assert len(hike) == 1
            stroll = hike[0]
            mapping = stroll.mapping
            kin_variables = {}
            mapping.map_to_lower_multiplicity(
                PS_point, stroll.structure, counterterm.momenta_dict, kin_variables )
            mapping.rescale_kinematic_variables(
                stroll.structure, counterterm.momenta_dict, kin_variables, base)
            mapping.map_to_higher_multiplicity(
                PS_point, stroll.structure, counterterm.momenta_dict, kin_variables )
        return

# Flat collinear walker
#=========================================================================================

class FlatCollinearWalker(VirtualWalker):

    cannot_handle = """The Flat Collinear walker found a singular structure
    it is not capable to handle.
    """
    collinear_map = FinalRescalingMappingOne()

    @staticmethod
    def get_recoilers(counterterm, parents=None):

        legs = counterterm.process['legs']
        model = counterterm.process['model']

        FINAL = base_objects.Leg.FINAL

        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in legs
            if leg['state'] == FINAL ]
        if parents:
            for recoiler in recoilers:
                if recoiler.n in parents:
                    recoilers.remove(recoiler)
                    break
        return recoilers

    @classmethod
    def determine_hike(cls, counterterm):

        # Initialize return variables
        hike = Hike()
        # If the counterterm is not trivial
        if counterterm.nodes:
            # Check it is only one
            if len(counterterm.nodes) != 1:
                raise MadGraph5Error(cls.cannot_handle)
            # Alias node and singular structure for convenience
            node = counterterm.nodes[0]
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.nodes:
                raise MadGraph5Error(cls.cannot_handle)
            # Get parent and children numbers
            parent, children, _ = get_structure_numbers(
                ss, counterterm.momenta_dict )
            # Pick recoilers
            recoilers = cls.get_recoilers(counterterm, (parent, ))
            # Compute jacobian and map to lower multiplicity
            if ss.name() == 'S' or node.nodes:
                raise MadGraph5Error(cls.cannot_handle)
            hike.append(Stroll(
                cls.collinear_map,
                subtraction.SingularStructure(legs=recoilers, substructures=(ss, )),
                node.current))
        # Return
        return hike

# Final-state only NLO walker
#=========================================================================================

class FinalNLOWalker(VirtualWalker):

    cannot_handle = """FinalNLOWalker found a singular structure
    it is not capable to handle.
    """

    # collinear_map = FinalLorentzMappingOne()
    collinear_map = FinalRescalingMappingOne()
    soft_map = MappingSomogyietalSoft()
    soft_collinear_map = SoftCollinearMapping(soft_map, collinear_map)

    @staticmethod
    def get_recoilers(counterterm, parents=None):

        legs = counterterm.process['legs']
        model = counterterm.process['model']

        FINAL = base_objects.Leg.FINAL

        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in legs
            if leg['state'] == FINAL and model.get_particle(leg['id'])['color'] != 1 ]
        if parents:
            for recoiler in recoilers:
                if recoiler.n in parents:
                    recoilers.remove(recoiler)
                    break
        return recoilers

    @classmethod
    def determine_hike(cls, counterterm):

        # Initialize return variables
        hike = Hike()
        # If the counterterm is not trivial
        if counterterm.nodes:
            # Check it is only one
            if len(counterterm.nodes) != 1:
                raise MadGraph5Error(cls.cannot_handle)
            # Alias node and singular structure for convenience
            node = counterterm.nodes[0]
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.nodes:
                raise MadGraph5Error(cls.cannot_handle)
            # Get parent and children numbers
            parent, children, _ = get_structure_numbers(
                ss, counterterm.momenta_dict )
            # Pick recoilers
            recoilers = cls.get_recoilers(counterterm, (parent, ))
            structure = subtraction.SingularStructure(legs=recoilers, substructures=(ss,))
            # Compute jacobian and map to lower multiplicity
            if ss.name() == 'S':
                hike.append(Stroll(cls.soft_map, structure, node.current))
            elif ss.name() == 'C' and not ss.substructures:
                hike.append(Stroll(cls.collinear_map, structure, node.current))
            elif (ss.name() == 'C'
                and len(ss.substructures) == 1
                and ss.substructures[0].name() == 'S'
                and not ss.substructures[0].substructures
                and len(ss.legs) == 1 ):
                hike.append(Stroll(cls.soft_collinear_map, structure, node.current))
            else:
                raise MadGraph5Error(cls.cannot_handle)
        # Return
        return hike

# General NLO walker
#=========================================================================================

class NLOWalker(VirtualWalker):

    cannot_handle = "NLOWalker found a singular structure it cannot handle."

    # NOTE: If rescaling mappings are used, only_colored_recoilers should be set to True.
    #       This might fail for some processes,
    #       e.g. gluon fusion Higgs production and Drell--Yan.

    f_collinear_map = FinalLorentzMappingOne()
    i_collinear_map = InitialLorentzMappingOne()
    soft_map = MappingSomogyietalSoft()
    f_soft_collinear_map = SoftCollinearMapping(soft_map, f_collinear_map)
    i_soft_collinear_map = SoftCollinearMapping(soft_map, i_collinear_map)
    only_colored_recoilers = True

    @classmethod
    def good_recoiler(cls, model, leg):

        return (
            leg['state'] == base_objects.Leg.FINAL and not
            (cls.only_colored_recoilers and model.get_particle(leg['id'])['color'] == 1) )

    @classmethod
    def get_recoilers(cls, counterterm, parents=None):

        legs = counterterm.process['legs']
        model = counterterm.process['model']

        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in legs if cls.good_recoiler(model, leg) ]
        if parents:
            for recoiler in recoilers:
                if recoiler.n in parents:
                    recoilers.remove(recoiler)
                    break
        return recoilers

    @classmethod
    def determine_hike(cls, counterterm):

        # Initialize return variables
        hike = Hike()
        # If the counterterm is not trivial
        if counterterm.nodes:
            # Check it is only one
            if len(counterterm.nodes) != 1:
                raise MadGraph5Error(cls.cannot_handle)
            # Alias node and singular structure for convenience
            node = counterterm.nodes[0]
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.nodes:
                raise MadGraph5Error(cls.cannot_handle)
            # Get parent number
            parent, _, _ = get_structure_numbers(ss, counterterm.momenta_dict)
            # Pick recoilers
            recoilers = cls.get_recoilers(counterterm, (parent, ))
            structure = subtraction.SingularStructure(legs=recoilers, substructures=(ss,))
            # Map to lower multiplicity
            if ss.name() == 'S':
                hike.append(Stroll(cls.soft_map,structure, node.current))
            elif ss.name() == 'C' and not ss.substructures:
                if ss.legs.has_initial_state_leg():
                    hike.append(Stroll(cls.i_collinear_map, structure, node.current))
                else:
                    hike.append(Stroll(cls.f_collinear_map, structure, node.current))
            elif (ss.name() == 'C'
                  and len(ss.substructures) == 1
                  and ss.substructures[0].name() == 'S'
                  and not ss.substructures[0].substructures
                  and len(ss.legs) == 1):
                if ss.legs.has_initial_state_leg():
                    hike.append(Stroll(cls.i_soft_collinear_map, structure, node.current))
                else:
                    hike.append(Stroll(cls.f_soft_collinear_map, structure, node.current))
            else:
                raise MadGraph5Error(cls.cannot_handle)
        # Return
        return hike

# Dictionary of available walkers
#=========================================================================================

# Mapping classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Mapping.
# Note that this must be placed after all Mapping daughter classes in this module
# have been declared.
mapping_walker_classes_map = {
    'FlatCollinear': FlatCollinearWalker,
    'FinalNLO': FinalNLOWalker,
    'NLO': NLOWalker,
    'Unknown': None
}
