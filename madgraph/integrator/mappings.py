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

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')

from madgraph.integrator.vectors import Vector, LorentzVector
from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList

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

    children = frozenset((leg.n for leg in structure.get_all_legs()))
    if structure.name() == "S":
        return None, children
    else:
        return momenta_dict.inv[children], children

#=========================================================================================
# VirtualMapping
#=========================================================================================

class VirtualMapping(object):
    """Base class for elementary mapping implementations."""

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
        clusters of particle and recoilers recursively

        :param momenta_dict: two-way dictionary that associates a unique label
        to each cluster of one or more particles identified by their number

        :param kinematic_variables: if a non-empty dictionary is passed,
        the kinematic variables that are necessary to reproduce the higher-multiplicity
        phase-space point from the lower-multiplicity one will be set

        :param compute_jacobian: if False, will not compute the jacobian for the mapping

        :param masses: masses of the parents of the collinear sets,
        as a dictionary with (number of parent, mass)

        :return: the jacobian weight due to the mapping
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
        clusters of particle and recoilers recursively

        :param momenta_dict: two-way dictionary that associates a unique label
        to each cluster of one or more particles identified by their number

        :param kinematic_variables: variables describing the splitting,
        as a dictionary that associates variable names to values

        :param compute_jacobian: if False, will not compute the jacobian for the mapping

        :return: the jacobian weight due to the mapping
        """
        
        raise NotImplemented

#===============================================================================
# Final-final collinear variables
#===============================================================================

class CollinearVariables(object):

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
        nap = na.dot(p)
        nbp = nb.dot(p)
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
            kinematic_variables['z' + str(i)] = zi
            kinematic_variables['kt' + str(i)] = kti
            kinematic_variables['p' + str(i) + '2'] = pi.square()
            zsum += zi
            ktsum += kti
            for j in range(len(kti)):
                ktabssum[j] += abs(kti[j])
        # Check numerical accuracy
        # TODO Ideally switch to quadruple precision if the check fails
        ktsum_abs = abs(ktsum.view(Vector))
        ktabssum_abs = abs(ktabssum.view(Vector))
        ktsum_ratio = ktsum_abs / ktabssum_abs
        if (abs(zsum - 1) > precision) or (ktsum_ratio > precision):
            logger.critical(CollinearVariables.precision_loss_message)
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
        nap = na.dot(p)
        nbp = nb.dot(p)
        nanb = na.dot(nb)
        pt = p - (nbp*na + nap * nb) / nanb
        # Variables for sums
        p_sum = LorentzVector()
        # Set momenta for all children
        for i in children:
            zi = kinematic_variables['z' + str(i)]
            kti = kinematic_variables['kt' + str(i)]
            pi2 = kinematic_variables['p' + str(i) + '2']
            nbpi = zi*nbp
            pti = kti + zi*pt
            napi = (pi2-pti.square())*nanb/(2*nbpi)
            PS_point[i] = (nbpi*na+napi*nb) / nanb + pti
            p_sum += PS_point[i]
        # Check how well the parent's momentum is reproduced
        # TODO Ideally switch to quadruple precision if the check fails
        deviation = abs((p - p_sum).view(Vector))
        benchmark = abs(p.view(Vector))
        if deviation / benchmark > precision:
            logger.critical(CollinearVariables.precision_loss_message)
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

    @staticmethod
    def azimuth_reference_vector(n):
        """Given a unit vector n,
        return a reference vector in the plane orthogonal to n.
        """

        n_ref = Vector([1., 0., 0.])
        return Vector(n_ref - n.dot(n_ref) * n).normalize()

    @staticmethod
    def to_random(mapped_momentum, children, variables):
        """Get [0,1] random variables
        that represent the internal kinematic variables.
        """

        # WARNING This function is a stub

        # Empty dictionary for randoms
        randoms = dict()

        # Compute reference vectors for this collinear subset
        nvec = Vector([mapped_momentum[j + 1] for j in range(3)]).normalize()
        # Get a reference phi=0 direction
        n_phi = CollinearVariables.azimuth_reference_vector(nvec)
        n_triple = nvec.cross(n_phi)

        # TODO Include virtuality
        # Compute all kinematic variables
        for i in children[:-1]:
            # Actual kinematic variables
            zi = variables['z'+str(i)]
            kti = variables['kt'+str(i)]
            kti_norm = math.sqrt(-kti.square())
            nti = Vector([kti[1:3]]).normalize()
            phi = math.acos(nti.dot(n_phi))
            if nti.dot(n_triple) < 0:
                phi *= -1
            # TODO zi and kti_norm need to be mapped to [0,1]
            randoms['z'+str(j)] = zi
            randoms['kt'+str(i)] = kti_norm
            randoms['phi'+str(i)] = phi

        return randoms

    @staticmethod
    def from_random(mapped_momentum, children, randoms):
        """Compute collinear variables from randoms in the range [0,1].
        """

        # TODO implement this
        variables = dict()
        return variables

#===============================================================================
# Final-final collinear mappings
#===============================================================================

class FFCollinearMapping(VirtualMapping):
    """Common functions for final-final collinear elementary mappings."""

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
            parent, children = get_structure_numbers(substructure, momenta_dict)
            CollinearVariables.names(parent, children)
        return names

class FFRescalingMappingOne(FFCollinearMapping):
    """Implementation of the rescaling mapping
    for one bunch of collinear particles (with massless parent)
    and an arbitrary number of massless recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one bunch of final-state particles going collinear,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(FFRescalingMappingOne, cls).is_valid_structure(singular_structure)
        return False

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        # Precompute sets and numbers
        cluster = singular_structure.substructures[0]
        parent, children = get_structure_numbers(cluster, momenta_dict)
        # Build collective momenta
        pC = LorentzVector()
        for j in children:
            pC += PS_point[j]
        pR = LorentzVector()
        for leg in singular_structure.legs:
            pR += PS_point[leg.n]
        Q = pC + pR
        # Compute the parameter alpha
        Q2  = Q.square()
        pC2 = pC.square()
        QpC = pC.dot(Q)
        alpha = (QpC - math.sqrt(QpC**2 - Q2*pC2)) / Q2
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] /= (1-alpha)
        # Map the cluster's momentum
        PS_point[parent] = (pC - alpha * Q) / (1-alpha)
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = CollinearVariables.collinear_and_reference(PS_point[parent])
            kinematic_variables['s'+str(parent)] = pC2
            CollinearVariables.get(PS_point, children, na, nb, kinematic_variables)
        # Eliminate children momenta from the mapped phase-space point
        for j in children:
            if j != parent: # Bypass degenerate case of 1->1 splitting
                del PS_point[j]

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

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
        cluster = singular_structure.substructures[0]
        parent, children = get_structure_numbers(cluster, momenta_dict)
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
        sC  = kinematic_variables['s'+str(parent)]
        # Obtain parameter alpha
        alpha = 0
        if (qR2*sC)/(QqC**2) < 1e-6:
            alpha = 0.5 * sC / QqC
        else:
            alpha = (math.sqrt(QqC**2 + sC*qR2) - QqC) / qR2
        # Compute reverse-mapped momentum
        pC = (1-alpha) * qC + alpha * Q
        # Map recoil momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] *= 1-alpha
        # Set children momenta
        na, nb = CollinearVariables.collinear_and_reference(qC)
        CollinearVariables.set(PS_point, children, pC, na, nb, kinematic_variables)
        # Remove parent's momentum
        if parent not in children: # Bypass degenerate case of 1->1 splitting
            del PS_point[parent]

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

class FFLorentzMappingOne(FFCollinearMapping):
    """Implementation of the Lorentz transformation mapping
    for one bunch of collinear particles (with massless parent)
    and an arbitrary number of (eventually massive) recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one bunch of final-state particles going collinear,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(FFLorentzMappingOne, cls).is_valid_structure(singular_structure)
        return False

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None, compute_jacobian=False, masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert cls.is_valid_structure(singular_structure)

        # Precompute sets and numbers
        cluster = singular_structure.substructures[0]
        parent, children = get_structure_numbers(cluster, momenta_dict)
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
        # Map the cluster's momentum
        PS_point[parent] = alpha*pC_perp + ((Q2-pR2)/(2*Q2))*Q
        # Map all recoilers' momenta
        qR = Q - PS_point[parent]
        for recoiler in singular_structure.legs:
            # TODO Move this try/except to higher level
            try:
                PS_point[recoiler.n].rotoboost(pR, qR)
            except:
                logger.critical("Problem encountered for %s" % str(singular_structure))
                logger.critical("The full phase space point was\n%s" % str(PS_point))
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = CollinearVariables.collinear_and_reference(PS_point[parent])
            kinematic_variables['s'+str(parent)] = pC2
            CollinearVariables.get(PS_point, children, na, nb, kinematic_variables)
        # Eliminate children momenta from the mapped phase-space point
        for j in children:
            if j != parent: # Bypass degenerate case of 1->1 splitting
                del PS_point[j]

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

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
        cluster = singular_structure.substructures[0]
        parent, children = get_structure_numbers(cluster, momenta_dict)
        # Build collective momenta
        qC = PS_point[parent]
        qR = LorentzVector()
        for leg in singular_structure.legs:
            qR += PS_point[leg.n]
        Q = qR + qC
        # Compute scalar products
        pC2 = kinematic_variables['s'+str(parent)]
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
        na, nb = CollinearVariables.collinear_and_reference(qC)
        CollinearVariables.set(PS_point, children, pC, na, nb, kinematic_variables)
        # Remove parent's momentum
        if parent not in children: # Bypass degenerate case of 1->1 splitting
            del PS_point[parent]

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

#===============================================================================
# Final-final soft variables
#===============================================================================

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
            kinematic_variables['p' + str(i)] = PS_point[i]
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

    @staticmethod
    def to_randoms(mapped_momentum, children, variables):
        """Get random variables in the range [0,1]
        that correspond to the internal kinematic variables.
        """

        raise NotImplemented

    @staticmethod
    def from_randoms(mapped_momentum, children, randoms):
        """Compute internal kinematic variables
        from random numbers in the range [0,1].
        """

        raise NotImplemented

#=========================================================================================
# Soft mappings
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
            parent, children = get_structure_numbers(substructure, momenta_dict)
            for legn in sorted(children):
                names += ['p'+str(legn), ]
        return names

class MappingSomogyietalSoft(ElementaryMappingSoft):
    """Implementation of the mapping used by Somogyi et al.
    in arXiv:hep-ph/0609042 for soft particles.
    It applies when there is at least one recoiler in the final state
    and all recoilers are massless.
    """

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
            children = [leg.n for leg in substructure.legs]
            if kinematic_variables is not None:
                SoftVariables.get(PS_point, children, kinematic_variables)
            for child in children:
                pS += PS_point.pop(child)
        # Build the total momentum of recoilers
        pR = LorentzVector()
        for leg in singular_structure.legs:
            pR += PS_point[leg.n]
        # Build the total momentum Q
        Q = pS + pR
        # Compute the parameter la
        la = math.sqrt(pR.square() / Q.square())
        P = pR / la
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] /= la
            PS_point[recoiler.n].rotoboost(P, Q)

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

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
            children = [leg.n for leg in substructure.legs]
            SoftVariables.set(PS_point, children, kinematic_variables)
            for child in children:
                pS += PS_point[child]
        # Build the total momentum, which is equal to the mapped recoilers'
        Q = LorentzVector()
        for leg in singular_structure.legs:
            Q += PS_point[leg.n]
        # Build the total momentum Q
        pR = Q - pS
        # Compute the parameter la
        la = math.sqrt(pR.square() / Q.square())
        P = pR / la
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] *= la
            PS_point[recoiler.n].rotoboost(Q, P)

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

#=========================================================================================
# Mapping walkers
#=========================================================================================

class VirtualWalker(object):
    """Base class for walker implementations."""
    
    def __new__(cls, **opts):
        """Factory class to make plugin easy."""

        if cls is VirtualWalker:
            map_type = opts.pop('map_type') if 'map_type' in opts else 'Unknown'
            if map_type not in mapping_walker_classes_map or not mapping_walker_classes_map[map_type]:
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

    def walk_to_lower_multiplicity(
        self, PS_point, counterterm, kinematic_variables=False
    ):
        """Starting from the highest-multiplicity phase-space point,
        generate all lower-multiplicity phase-space points
        that are necessary for the evaluation of the given counterterm.

        :param PS_point: highest-multiplicity starting phase-space point,
        as a dictionary that associates integers to Lorentz vectors.

        :param counterterm: Counterterm object that specifies
        clusters of particle and recoilers recursively.
        Momenta_dict will obtained directly from it.

        :param kinematic_variables: flag that specifies whether to compute
        all kinematic variables needed to recover the starting phase-space point
        from the lowest multiplicity one.

        :return: a dictionary with the following entries:
        'currents', a list of all currents that need to be evaluated,
            paired with the momenta needed for their evaluation;
        'matrix_element', the reduced matrix element,
            paired with its corresponding reduced phase-space point;
        'jacobian', the cumulative jacobian of all mappings that were applied;
        'kinematic_variables', a dictionary of all variables needed to recover
            the starting phase-space point from the lowest multiplicity one,
            or None if such variables were not requested;
        'resulting_PS_point', the final lowest multiplicity PS point.
        """

        raise NotImplementedError

    def walk_to_higher_multiplicity(
        self, PS_point, counterterm, kinematic_variables
    ):
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
        'jacobian', the cumulative jacobian of all mappings that were applied;
        'resulting_PS_point', the final lowest multiplicity PS point.
        """
         
        raise NotImplementedError
        
    def rescale_kinematic_variables(
        self, counterterm, kinematic_variables, scaling_parameter
    ):
        """Rescale kinematic_variables with the scaling_parameter provided,
        so as to progressively approach the IR limit specified by counterterm.
        """

        raise NotImplementedError
        
    def approach_limit(
        self, PS_point, counterterm, kinematic_variables, scaling_parameter
    ):
        """Produce a higher multiplicity phase-space point from PS_point,
        according to kinematic_variables that approach the limit of counterterm
        parametrically with scaling_parameter.
        """

        new_kin_variables = self.rescale_kinematic_variables(
            counterterm, kinematic_variables, scaling_parameter )
        return self.walk_to_higher_multiplicity(
            PS_point, counterterm, new_kin_variables )

class FlatCollinearWalker(VirtualWalker):

    cannot_handle = """The Flat Collinear walker found a singular structure
    it is not capable to handle.
    """
    collinear_map = FFRescalingMappingOne()

    def walk_to_lower_multiplicity(
        self, PS_point, counterterm, compute_kinematic_variables=False, verbose=False
    ):

        point = PS_point.get_copy()
        if verbose:
            print point
        # Initialize return variables
        current_PS_pairs = []
        jacobian = 1.
        kinematic_variables = dict() if compute_kinematic_variables else None
        # Recoil against all final-state particles that are not going singular
        # TODO Recoilers and numbers for a given counterterm should be cached
        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in counterterm.process['legs']
            if leg['state'] == base_objects.Leg.FINAL
        ]
        for node in counterterm.nodes:
            parent, _ = get_structure_numbers(
                node.current['singular_structure'],
                counterterm.momenta_dict
            )
            for recoiler in recoilers:
                if recoiler.n == parent:
                    recoilers.remove(recoiler)
        # Loop over the counterterm's first-level currents
        for node in counterterm.nodes:
            # Alias singular structure for convenience
            ss = node.current['singular_structure']
            # Do not accept soft currents or nested ones
            if ss.name == 'S' or node.nodes:
                raise MadGraph5Error(self.cannot_handle)
            # Get parent and children numbers
            # TODO Read these from cache
            parent, children = get_structure_numbers(
                ss, counterterm.momenta_dict
            )
            # Deep copy the momenta
            momenta = point.get_copy()
            # Compute the jacobian and map to lower multiplicity
            jacobian *= self.collinear_map.map_to_lower_multiplicity(
                point,
                subtraction.SingularStructure(legs=recoilers, substructures=(ss, )),
                counterterm.momenta_dict,
                kinematic_variables
            )
            if verbose:
                print point
            # Append the current and the momenta,
            # deep copy wanted
            momenta[parent] = LorentzVector(point[parent])
            current_PS_pairs.append((node.current, momenta))
        # Identify reduced matrix element,
        # computed in the point which has received all mappings
        ME_PS_pair = [counterterm.process, point]
        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'jacobian': jacobian,
            'kinematic_variables': kinematic_variables,
            'resulting_PS_point': point
        }

    def walk_to_higher_multiplicity(
        self, PS_point, counterterm, kinematic_variables, verbose=False
    ):

        # Identify reduced matrix element,
        # computed in the lowest multiplicity point
        ME_PS_pair = [counterterm.process, PS_point]
        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        if verbose:
            print point
        # Initialize return variables
        current_PS_pairs = []
        jacobian = 1.
        # Recoil against all final-state particles that are not going singular
        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in counterterm.process['legs']
            if leg['state'] == base_objects.Leg.FINAL
        ]
        for node in counterterm.nodes:
            parent, _ = get_structure_numbers(
                node.current['singular_structure'],
                counterterm.momenta_dict
            )
            for recoiler in recoilers:
                if recoiler.n == parent:
                    recoilers.remove(recoiler)
        # Loop over the counterterm's first-level currents
        for node in reversed(counterterm.nodes):
            # Alias singular structure for convenience
            ss = node.current['singular_structure']
            # Do not accept soft currents or nested ones
            if ss.name == 'S' or node.nodes:
                raise MadGraph5Error(self.cannot_handle)
            # Get parent and children numbers
            # TODO Read these from cache
            parent, children = get_structure_numbers(
                ss, counterterm.momenta_dict
            )
            # Deep copy the parent momentum
            momenta = LorentzVectorDict({parent: LorentzVector(point[parent])})
            # Compute the jacobian and map to higher multiplicity
            jacobian *= self.collinear_map.map_to_higher_multiplicity(
                point,
                subtraction.SingularStructure(legs=recoilers, substructures=(ss, )),
                counterterm.momenta_dict,
                kinematic_variables )
            if verbose:
                print point
            # Prepend pair of this current and the momenta,
            # deep copy wanted
            momenta.update(point.get_copy())
            current_PS_pairs.insert(0, (node.current, momenta))

        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'jacobian': jacobian,
            'resulting_PS_point': point }

    def rescale_kinematic_variables(
        self, counterterm, kinematic_variables, scaling_parameter ):

        # For collinear clusters, just rescale the virtuality
        new_kinematic_variables = {}
        for var, value in kinematic_variables.items():
            if var.startswith('s'):
                new_kinematic_variables[var] = value*scaling_parameter**2
            elif var.startswith('k'):
                new_kinematic_variables[var] = value*scaling_parameter
            else:
                new_kinematic_variables[var] = value
        return new_kinematic_variables

class FFNLOWalker(VirtualWalker):

    cannot_handle = """FFNLOWalker found a singular structure
    it is not capable to handle.
    """

    collinear_map = FFLorentzMappingOne()
#    collinear_map = FFRescalingMappingOne()
    soft_map = MappingSomogyietalSoft()

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

    def walk_to_lower_multiplicity(
        self, PS_point, counterterm, compute_kinematic_variables=False, verbose=False ):

        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        if verbose:
            print point
        # Initialize return variables
        current_PS_pairs = []
        jacobian = 1.
        kinematic_variables = dict() if compute_kinematic_variables else None
        # If the counterterm is not trivial
        if counterterm.nodes:
            # Check it is only one
            if len(counterterm.nodes) != 1:
                raise MadGraph5Error(self.cannot_handle)
            # Alias node and singular structure for convenience
            node = counterterm.nodes[0]
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.nodes:
                raise MadGraph5Error(self.cannot_handle)
            # Get parent and children numbers
            parent, children = get_structure_numbers(
                ss, counterterm.momenta_dict )
            # Deep copy the momenta
            momenta = point.get_copy()
            # Pick recoilers
            recoilers = self.get_recoilers(counterterm, (parent,))
            # Compute jacobian and map to lower multiplicity
            if ss.name() == 'S':
                jacobian *= self.soft_map.map_to_lower_multiplicity(
                    point,
                    subtraction.SingularStructure(legs=recoilers, substructures=(ss,)),
                    counterterm.momenta_dict,
                    kinematic_variables )
            elif ss.name() == 'C' and not ss.substructures:
                jacobian *= self.collinear_map.map_to_lower_multiplicity(
                    point,
                    subtraction.SingularStructure(legs=recoilers, substructures=(ss,) ),
                    counterterm.momenta_dict,
                    kinematic_variables )
                momenta[parent] = LorentzVector(point[parent])
            elif (ss.name() == 'C'
                and len(ss.substructures) == 1
                and ss.substructures[0].name() == 'S'
                and not ss.substructures[0].substructures
                and len(ss.legs) == 1 ):
                # Make a fake soft structure for soft-collinears
                # in order to use the mapping for soft substructures
                new_re = recoilers + list(ss.legs)
                new_ss = subtraction.SoftStructure(legs=ss.substructures[0].legs)
                # kinematic_variables['s'+str(parent)] = point[parent].square()
                jacobian *= self.soft_map.map_to_lower_multiplicity(
                    point,
                    subtraction.SingularStructure(legs=new_re, substructures=(new_ss,) ),
                    counterterm.momenta_dict,
                    kinematic_variables )
                point[parent] = point[ss.legs[0].n]
                del point[ss.legs[0].n]
                momenta[parent] = LorentzVector(point[parent])
            else:
                raise MadGraph5Error(self.cannot_handle)
            # Append the current and the momenta
            current_PS_pairs.append((node.current, momenta))
        # Identify reduced matrix element,
        # computed in the point which has received all mappings
        ME_PS_pair = [counterterm.process, point]
        if verbose:
            print point
        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'jacobian': jacobian,
            'kinematic_variables': kinematic_variables,
            'resulting_PS_point': point }

    def walk_to_higher_multiplicity(
        self, PS_point, counterterm, kinematic_variables, verbose=False ):

        # Identify reduced matrix element,
        # computed in the lowest multiplicity point
        ME_PS_pair = [counterterm.process, PS_point]
        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        if verbose:
            print point
        # Initialize return variables
        current_PS_pairs = []
        jacobian = 1.
        # If the counterterm is not trivial
        if counterterm.nodes:
            # Check it is only one
            if len(counterterm.nodes) != 1:
                raise MadGraph5Error(self.cannot_handle)
            # Alias node and singular structure for convenience
            node = counterterm.nodes[0]
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.nodes:
                raise MadGraph5Error(self.cannot_handle)
            # Get parent and children numbers
            parent, children = get_structure_numbers(
                ss, counterterm.momenta_dict )
            # Deep copy the momenta
            momenta = LorentzVectorDict()
            # Pick recoilers
            recoilers = self.get_recoilers(counterterm, (parent,))
            # Compute jacobian and map to lower multiplicity
            if ss.name() == 'S':
                jacobian *= self.soft_map.map_to_higher_multiplicity(
                    point,
                    subtraction.SingularStructure(legs=recoilers, substructures=(ss,)),
                    counterterm.momenta_dict,
                    kinematic_variables )
            elif ss.name() == 'C' and not ss.substructures:
                momenta[parent] = LorentzVector(point[parent])
                jacobian *= self.collinear_map.map_to_higher_multiplicity(
                    point,
                    subtraction.SingularStructure(legs=recoilers, substructures=(ss,) ),
                    counterterm.momenta_dict,
                    kinematic_variables )
            elif (ss.name() == 'C'
                and len(ss.substructures) == 1
                and ss.substructures[0].name() == 'S'
                and not ss.substructures[0].substructures
                and len(ss.legs) == 1 ):
                # Make a fake soft structure for soft-collinears
                # in order to use the mapping for soft substructures
                new_re = recoilers + list(ss.legs)
                new_ss = subtraction.SoftStructure(legs=ss.substructures[0].legs)
                momenta[parent] = LorentzVector(point[parent])
                point[ss.legs[0].n] = point[parent]
                del point[parent]
                jacobian *= self.soft_map.map_to_higher_multiplicity(
                    point,
                    subtraction.SingularStructure(legs=new_re, substructures=(new_ss,) ),
                    counterterm.momenta_dict,
                    kinematic_variables )
            else:
                raise MadGraph5Error(self.cannot_handle)
            # Prepend pair of this current and the momenta,
            # deep copy wanted
            momenta.update(point.get_copy())
            current_PS_pairs.insert(0, (node.current, momenta))
            if verbose:
                print point
        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'jacobian': jacobian,
            'resulting_PS_point': point }

    def rescale_kinematic_variables(
        self, counterterm, kinematic_variables, scaling_parameter ):

        new_variables = {
            key: type(val)(val) for (key, val) in kinematic_variables.items() }
        
        # Determine 'equivalent power' of rescaling
        # This structure has at most two levels - checked later
        total_limits = 0
        for node in counterterm.nodes:
            total_limits += 1 + len(node.current['singular_structure'].substructures)
        scaling_parameter **= 1. / total_limits
        for node in counterterm.nodes:
            ss = node.current['singular_structure']
            if len(ss.substructures) != 0:
                raise MadGraph5Error(self.cannot_handle)
            if ss.name() == 'C':
                # For collinear clusters, rescale virtuality and transverse momenta
                for var in new_variables.keys():
                    if var.startswith('s'):
                        new_variables[var] *= scaling_parameter**2
                    elif var.startswith('k'):
                        new_variables[var] *= scaling_parameter
            elif ss.name() == 'S':
                # For soft clusters, rescale the whole momenta
                for leg in ss.legs:
                    new_variables['p' + str(leg.n)] *= scaling_parameter
            else:
                raise MadGraph5Error(self.cannot_handle)
        return new_variables

# Mapping classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Mapping.
# Note that this must be placed after all the Mapping daughter classes in this module have been declared.
mapping_walker_classes_map = {
    'FlatCollinear': FlatCollinearWalker,
    'FFNLO': FFNLOWalker,
    'Unknown': None
}
