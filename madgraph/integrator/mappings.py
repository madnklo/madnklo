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
    import internal.subtraction as sub
    from internal import InvalidCmd, MadGraph5Error
else:
    MADEVENT= False
    import madgraph.various.misc as misc
    import madgraph.core.subtraction as sub
    from madgraph import InvalidCmd, MadGraph5Error

import copy

from madgraph.integrator.vectors import Vector, LorentzVector
from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')

class MappingError(Exception):
    """Exception raised if an exception is triggered in implementation of the currents.""" 
    pass

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
            if leg.state == leg.INITIAL )
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
            names += ['z' + str(child), 'm2' + str(child), 'kt' + str(child), ]
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

        if len(children) < 2: return
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
            kinematic_variables['m2' + str(i)] = pi.square()
            zsum += zi
            ktsum += kti
            for j in range(len(kti)):
                ktabssum[j] += abs(kti[j])
        # Check numerical accuracy
        # TODO Ideally switch to quadruple precision if the check fails
        ktsum_abs = abs(ktsum)
        ktabssum_abs = abs(ktabssum)
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
        PS_point, parent, children, na, nb, kinematic_variables,
        precision=1e-6 ):
        """Given a phase-space point with an off-shell parent and collinear variables,
        compute and set the children momenta.
        The parent's momentum is removed.
        Parent and children are indices that already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        if len(children) < 2:
            for child in children:
                if child != parent:
                    PS_point[child] = PS_point[parent]
                    del PS_point[parent]
            return
        # Rename the sum of momenta
        p = PS_point[parent]
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
            pi2 = kinematic_variables['m2' + str(i)]
            nbpi = zi*nbp
            pti = kti + zi*pt
            napi = (pi2-pti.square())*nanb/(2*nbpi)
            PS_point[i] = (nbpi*na+napi*nb) / nanb + pti
            p_sum += PS_point[i]
        # Check how well the parent's momentum is reproduced
        # TODO Ideally switch to quadruple precision if the check fails
        deviation = abs(p - p_sum)
        benchmark = abs(p)
        if deviation / benchmark > precision:
            logger.critical(FinalCollinearVariables.precision_loss_message)
            logger.critical("The sum of children momenta is %s" % str(p_sum))
            logger.critical("Inputs for FinalCollinearVariables.set():")
            logger.critical("total momentum = %s" % str(p))
            logger.critical("na = %s, nb = %s" % (str(na), str(nb)))
            logger.critical("kinematic variables:")
            logger.critical(str(kinematic_variables))
            logger.critical("Output of FinalCollinearVariables.set():")
            for i in children:
                logger.critical("child %d: %s" % (i, str(PS_point[i])))
        del PS_point[parent]
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
            names += ['z' + str(child), 'm2' + str(child), 'kt' + str(child), ]
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
            kinematic_variables['m2' + str(i)] = pi.square()
            zsum += zi
            ktsum += kti
            for j in range(len(kti)):
                ktabssum[j] += abs(kti[j])
        # Check numerical accuracy
        # TODO Ideally switch to quadruple precision if the check fails
        if not fs_children: return
        ktsum_abs = abs(ktsum)
        ktabssum_abs = abs(ktabssum)
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
        PS_point, is_child, fs_children, na, nb, kinematic_variables,
        precision=1e-6 ):
        """Given the lower multiplicity momentum of the incoming parton
        as PS_point[is_child] and collinear variables,
        compute and set the final-state children momenta.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        if not fs_children: return
        pa = PS_point[is_child]
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
            pi2 = kinematic_variables['m2' + str(i)]
            pti = kti
            # pti = kti + zi * ktA
            nbpi = zi * nbpa
            napi = (pi2 - pti.square()) * nanb / (2 * nbpi)
            PS_point[i] = (nbpi*na + napi*nb) / nanb + pti
            p_sum -= PS_point[i]
        # Check how well the parent's momentum is reproduced
        # TODO Ideally switch to quadruple precision if the check fails
        deviation = abs(pA - p_sum)
        benchmark = abs(pA)
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
        """Check if the mapping can be applied to a given singular structure."""

        raise NotImplemented

    @classmethod
    def get_kinematic_variables_names(cls, singular_structure, momenta_dict):
        """Get the names of the variables that describe unresolved particles."""

        raise NotImplemented

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):
        """Map a phase-space point to lower multiplicity,
        by clustering the substructures and recoiling against the legs
        specified in singular_structure.

        :param PS_point: higher-multiplicity phase-space point,
        as a dictionary that associates integers to Lorentz vectors;
        this will not be modified
        :type PS_point: LorentzVectorDict

        :param singular_structure: SingularStructure object that specifies
        sets of unresolved particles and recoilers recursively
        :type PS_point: SingularStructure

        :param momenta_dict: two-way dictionary that associates a unique label
        to each set of one or more unresolved particles identified by their number
        :type momenta_dict: sub.bidict

        :param squared_masses: squared masses of parents of particle sets,
        as a dictionary {'m2i': $m_i^2$} where i is the parent number

        :param kinematic_variables: if a non-empty dictionary is passed,
        the kinematic variables that are necessary to reproduce the higher-multiplicity
        phase-space point from the lower-multiplicity one will be set

        :param compute_jacobian: if False, the jacobian of the mapping will be set to 1
        :type compute_jacobian: bool

        :return: lower-multiplicity phase-space point and mapping variables,
        including the jacobian weight of the mapping and the total momentum involved
        """
        
        raise NotImplemented

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):
        """Map a phase-space point to higher multiplicity,
        by splitting the (pseudo)particles and recoiling against the legs
        specified in singular_structure.

        :param PS_point: lower-multiplicity phase-space point
        which will not be modified
        :type PS_point: LorentzVectorDict

        :param singular_structure: SingularStructure object that specifies
        sets of unresolved particles and recoilers recursively
        :type PS_point: SingularStructure

        :param momenta_dict: two-way dictionary that associates a unique label
        to each set of one or more unresolved particles identified by their number
        :type momenta_dict: sub.bidict

        :param kinematic_variables: variables describing the splitting,
        as a dictionary that associates variable names to values

        :param compute_jacobian: if False, the jacobian of the mapping will be set to 1
        :type compute_jacobian: bool

        :return: higher-multiplicity phase-space point and jacobian of the mapping
        """
        
        raise NotImplemented

    @classmethod
    def can_map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables ):
        """Map a given phase-space point to higher multiplicity,
        by splitting the (pseudo)particles and recoiling against the legs
        specified in singular_structure.

        :param PS_point: lower-multiplicity phase-space point
        which will be modified to the higher-multiplicity one
        :type PS_point: LorentzVectorDict

        :param singular_structure: SingularStructure object that specifies
        sets of unresolved particles and recoilers recursively
        :type PS_point: SingularStructure

        :param momenta_dict: two-way dictionary that associates a unique label
        to each set of one or more unresolved particles identified by their number
        :type momenta_dict: sub.bidict

        :param kinematic_variables: variables describing the splitting,
        as a dictionary that associates variable names to values

        :return: boolean that specifies if the lower-multiplicity phase-space point
        has a corresponding higher-multiplicity point with the given kinematic variables
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


class FailedMapping(MadGraph5Error):
    """Exception raised when a mapping cannot be applied."""

    pass

#=========================================================================================
# Final mappings changing invariant masses
#=========================================================================================

# Final mapping to/from zero masses
#=========================================================================================

class FinalZeroMassesMapping(VirtualMapping):
    """Mapping that sends massive particles into massless particles."""

    plot = False

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, sub.SingularStructure)
        # Valid only for final-state particles with no recursive substructure
        for substructure in singular_structure.substructures:
            if substructure.substructures:
                return False
            if len(substructure.legs) != 1:
                return False
            if substructure.get_all_legs().has_initial_state_leg():
                return False
        # At least two particles
        assert len(singular_structure.substructures) > 1
        # No recoilers
        if singular_structure.legs:
            return False
        return True

    @classmethod
    def get_kinematic_variables_names(cls, singular_structure, momenta_dict):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        names = []
        for substructure in singular_structure.substructures:
            names.append('s' + str(tuple(substructure.legs)[0].n))
        return names

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Determine leg numbers and check that target masses are zero
        js = []
        for substructure in singular_structure.substructures:
            j = tuple(substructure.legs)[0].n
            js.append(j)
            assert (squared_masses is None) or (squared_masses['m2' + str(j)] == 0)
        # Build total momentum
        Q = LorentzVector()
        for j in js:
            Q += PS_point[j]
        Q2 = Q.square()
        # Compute the parameters gamma, alpha, beta
        gamma = {j: (Q.dot(PS_point[j]) / Q2) for j in js}
        # assert abs(sum(gamma.values()) - 1.) < 1.e-6
        mu2 = {j: (PS_point[j].square() / Q2) for j in js}
        alpha = 0.
        for j in js:
            alpha += (gamma[j] ** 2 - mu2[j]) ** 0.5
        beta = {j: ((gamma[j] ** 2 - mu2[j]) ** 0.5 / alpha) for j in js}
        # assert abs(sum(beta.values())  - 1.) < 1.e-6
        # Map all momenta
        new_PS_point = PS_point.get_copy()
        for j in js:
            new_PS_point[j] /= alpha
            new_PS_point[j] += (beta[j] - gamma[j]/alpha) * Q
            # assert abs(PS_point[j].square() / Q2) < 1.e-6
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            for j in js:
                kinematic_variables['s' + str(j)] = mu2[j] * Q2
        # Compute the jacobian for this mapping
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        sum_beta2_gamma = 0.
        prod_beta_gamma = 1.
        for j in js:
            sum_beta2_gamma += beta[j] ** 2 / gamma[j]
            prod_beta_gamma *= beta[j] / gamma[j]
        jacobian = alpha**(3*len(js)-5) * prod_beta_gamma / sum_beta2_gamma
        mapping_variables = {'jacobian': jacobian, 'Q': Q}
        # Return characteristic variables
        return new_PS_point, mapping_variables

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Determine leg numbers
        js = []
        for substructure in singular_structure.substructures:
            js.append(tuple(substructure.legs)[0].n)
        # Build total momentum
        Q = LorentzVector()
        for j in js:
            Q += PS_point[j]
        Q2 = Q.square()
        # Compute the parameters beta
        beta = {j: (Q.dot(PS_point[j]) / Q2) for j in js}
        # assert abs(sum(beta.values())  - 1.) < 1.e-6
        mu2 = {j: (kinematic_variables['s' + str(j)] / Q2) for j in js}
        if len(js) == 2:
            alpha = (Kaellen(1., mu2[js[0]], mu2[js[1]])) ** 0.5
        else:
            # Solve the equation for alpha numerically
            from scipy.optimize import newton
            func = lambda a: sum(((a*beta[j])**2 + mu2[j]) ** 0.5 for j in js) - 1.
            alpha = newton(func, 1.)
            # misc.sprint(alpha)
            # Optional plot to monitor the numerical solution
            if cls.plot:
                import matplotlib.pyplot as plt
                import numpy as np
                zero = lambda a: 0.*a
                tau = np.linspace(0., 1., 101)
                plt.plot(tau, func(tau), "-", tau, zero(tau), "-")
                plt.grid()
                plt.show()
        # Compute the parameters gamma
        gamma = {j: ((alpha*beta[j])**2 + mu2[j]) ** 0.5 for j in js}
        # assert abs(sum(gamma.values()) - 1.) < 1.e-6
        # Map all momenta
        new_PS_point = PS_point.get_copy()
        for j in js:
            new_PS_point[j] *= alpha
            new_PS_point[j] += (gamma[j] - alpha*beta[j]) * Q
            # assert abs(PS_point[j].square() / Q2 - mu2[j]) < 1.e-6
        # Compute the jacobian for this mapping
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        sum_beta2_gamma = 0.
        prod_beta_gamma = 1.
        for j in js:
            sum_beta2_gamma += beta[j] ** 2 / gamma[j]
            prod_beta_gamma *= beta[j] / gamma[j]
        jacobian = alpha**(3*len(js)-5) * prod_beta_gamma / sum_beta2_gamma
        mapping_variables = {'jacobian': jacobian, 'Q': Q}
        # Return characteristic variables
        return new_PS_point, mapping_variables

# Final masses mapping
#=========================================================================================

class FinalMassesMapping(FinalZeroMassesMapping):
    """Mapping that changes momenta invariant masses."""

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        zero_PS_point, zero_vars = FinalZeroMassesMapping.map_to_lower_multiplicity(
            PS_point, singular_structure, momenta_dict, None,
            kinematic_variables, compute_jacobian )
        fake_kinematic_variables = None
        if squared_masses is not None:
            fake_kinematic_variables = {
                's'+key[2:]: value for (key, value) in squared_masses.items()}
        mass_PS_point, mass_vars = FinalZeroMassesMapping.map_to_higher_multiplicity(
            zero_PS_point, singular_structure, momenta_dict, fake_kinematic_variables,
            compute_jacobian=compute_jacobian )
        if compute_jacobian:
            zero_vars['jacobian'] /= mass_vars['jacobian']
        # Return characteristic variables
        return mass_PS_point, zero_vars

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        fake_squared_masses = {
            'm2' + key[1:]: value
            for (key, value) in kinematic_variables.items()
            if key.startswith('s') }
        mapped_PS_point, vars = cls.map_to_lower_multiplicity(
            PS_point, singular_structure, momenta_dict, fake_squared_masses,
            None, compute_jacobian )
        if compute_jacobian:
            vars['jacobian'] **= -1
        return mapped_PS_point, vars

#=========================================================================================
# Final-collinear mappings
#=========================================================================================

# Generic final-collinear mapping parent class
#=========================================================================================

class FinalCollinearMapping(VirtualMapping):
    """Common functions for final-collinear elementary mappings."""

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, sub.SingularStructure)
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
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        names = []
        for substructure in singular_structure.substructures:
            parent, children, _ = get_structure_numbers(substructure, momenta_dict)
            FinalCollinearVariables.names(parent, children)
        return names


    @classmethod
    def rescale_kinematic_variables(
        cls, singular_structure, momenta_dict, kinematic_variables, scaling_parameter):

        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict))
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        # Determine the correct scaling of variables for the phase-space volume
        # to be proportional to scaling_parameter
        base = scaling_parameter ** (0.5 / (len(children) - 1))
        kinematic_variables['s' + str(parent)] *= base ** 2
        for child in children:
            kinematic_variables['kt' + str(child)] *= base
        return kinematic_variables

# Final-collinear rescaling mapping, one set
#=========================================================================================

class FinalRescalingOneMapping(FinalCollinearMapping):
    """Implementation of the rescaling mapping
    for one bunch of collinear particles (with massless parent)
    and an arbitrary number of massless recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one bunch of final-state particles going collinear,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(FinalRescalingOneMapping, cls).is_valid_structure(
                singular_structure )
        return False

    @staticmethod
    def abc(pC, Q):
        """Return the parameters alpha, beta, gamma and mu2 of the rescaling mapping."""

        Q2  = Q.square()
        mu2 = pC.square() / Q2
        gamma = pC.dot(Q) / Q2
        alpha = gamma - (gamma ** 2 - mu2) ** 0.5
        beta = (gamma - alpha) / (1 - alpha)
        return alpha, beta, gamma, mu2

    @staticmethod
    def alpha(pC, Q):
        """Return the parameter alpha of the rescaling mapping."""

        return FinalRescalingOneMapping.abc(pC, Q)[0]

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        assert (squared_masses is None) or (squared_masses['m2' + str(parent)] == 0)
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
        alpha, beta, gamma, mu2 = FinalRescalingOneMapping.abc(pC, Q)
        # Create new PS point
        new_PS_point = PS_point.get_copy()
        # Map all recoilers' momenta
        for recoiler in recoilers:
            new_PS_point[recoiler] /= (1-alpha)
        # Map the set's momentum
        new_PS_point[parent] = (pC - alpha * Q) / (1-alpha)
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = FinalCollinearVariables.collinear_and_reference(new_PS_point[parent])
            FinalCollinearVariables.get(PS_point, children, na, nb, kinematic_variables)
            kinematic_variables['s' + str(parent)] = pC.square()
        # Eliminate children momenta from the mapped phase-space point
        for j in children:
            if j != parent: # Bypass degenerate case of 1->1 splitting
                del new_PS_point[j]
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        # Compute the jacobian for this mapping
        jacobian = (1-alpha)**(2*len(recoilers)) * beta / (gamma - mu2)
        # Return characteristic variables
        return new_PS_point, {'jacobian': jacobian, 'Q': Q}

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
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
        assert qC.square_almost_zero()
        Q2 = Q.square()
        beta = Q.dot(qC) / Q2
        mu2 = kinematic_variables['s' + str(parent)] / Q2
        # Obtain parameter alpha
        if (abs(beta - 0.5) < 1.e-6):
            alpha = mu2
        else:
            alpha = ((beta ** 2 + (1-2*beta)*mu2) ** 0.5 - beta) / (1-2*beta)
        gamma = alpha + (1-alpha) * beta
        # Compute reverse-mapped momentum
        pC = (1-alpha) * qC + alpha * Q
        # Create new PS point
        new_PS_point = PS_point.get_copy()
        new_PS_point[parent] = pC
        # Map recoil momenta
        for recoiler in recoilers:
            new_PS_point[recoiler] *= 1-alpha
        # Set children momenta
        na, nb = FinalCollinearVariables.collinear_and_reference(qC)
        FinalCollinearVariables.set(
            new_PS_point, parent, children, na, nb, kinematic_variables)
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        # Return the jacobian for this mapping
        jacobian = (1-alpha)**(2*len(recoilers)) * beta / (gamma - mu2)
        # Return characteristic variables
        return new_PS_point, {'jacobian': jacobian, 'Q': Q}

# Final-collinear Lorentz mapping, one set
#=========================================================================================

class FinalLorentzOneMapping(FinalCollinearMapping):
    """Implementation of the Lorentz transformation mapping
    for one bunch of collinear particles (with massless parent)
    and an arbitrary number of (eventually massive) recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one bunch of final-state particles going collinear,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(FinalLorentzOneMapping, cls).is_valid_structure(
                singular_structure )
        return False

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, children, _ = get_structure_numbers(substructure, momenta_dict)
        assert (squared_masses is None) or (squared_masses['m2' + str(parent)] == 0)
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
        alpha = (Q2-pR2) / math.sqrt(Kaellen(Q2, pR2, pC2))
        # Create new PS point
        new_PS_point = PS_point.get_copy()
        # Map the set's momentum
        pC_perp = pC - ((Q2+pC2-pR2)/(2*Q2)) * Q
        new_PS_point[parent] = alpha*pC_perp + ((Q2-pR2)/(2*Q2))*Q
        # Map all recoilers' momenta
        qR = Q - new_PS_point[parent]
        for recoiler in recoilers:
            # TODO Move this try/except to higher level
            try:
                new_PS_point[recoiler].rotoboost(pR, qR)
                # new_PS_point[recoiler].boost_from_to(pR, qR)
            except:
                logger.critical("Problem encountered for %s" % str(singular_structure))
                logger.critical("The full phase space point was\n%s" % str(new_PS_point))
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = FinalCollinearVariables.collinear_and_reference(new_PS_point[parent])
            kinematic_variables['s' + str(parent)] = pC2
            FinalCollinearVariables.get(PS_point, children, na, nb, kinematic_variables)
        # Eliminate children momenta from the mapped phase-space point
        for j in children:
            if j != parent: # Bypass degenerate case of 1->1 splitting
                del new_PS_point[j]
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        jacobian = 1. / alpha
        # Return characteristic variables
        return new_PS_point, {'jacobian': jacobian, 'Q': Q}

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
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
        assert qC.square_almost_zero()
        Q2  = Q.square()
        qR2 = qR.square()
        qC_perp = qC - ((Q2-qR2)/(2*Q2)) * Q
        kaellen = Kaellen(Q2, qR2, pC2)
        if kaellen < 0:
            raise FailedMapping
        alpham1 = kaellen ** 0.5 / (Q2-qR2)
        # Compute reverse-mapped momentum
        pC = alpham1*qC_perp + ((Q2+pC2-qR2)/(2*Q2))*Q
        pR = Q - pC
        # Create new PS point
        new_PS_point = PS_point.get_copy()
        new_PS_point[parent] = pC
        # Map recoil momenta
        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n].rotoboost(qR, pR)
            # new_PS_point[recoiler.n].boost_from_to(qR, pR)
        # Set children momenta
        na, nb = FinalCollinearVariables.collinear_and_reference(qC)
        FinalCollinearVariables.set(
            new_PS_point, parent, children, na, nb, kinematic_variables)
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        # Return characteristic variables
        return new_PS_point, {'jacobian': alpham1, 'Q': Q}

    @classmethod
    def can_map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables, ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
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
        pC2 = kinematic_variables['s'+str(parent)]
        # Compute parameters
        assert qC.square_almost_zero()
        Q2  = Q.square()
        qR2 = qR.square()
        kaellen = Kaellen(Q2, qR2, pC2)
        return kaellen >= 0

# Final-collinear grouping mapping
#=========================================================================================

class FinalGroupingMapping(FinalCollinearMapping):
    """Implementation of the mapping that groups sets of collinear particles
    and assigns them a new mass.
    """

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Build reduced PS point, singular structure & masses to call FinalMassesMapping
        reduced_PS_point = PS_point.get_copy()
        reduced_singular_structure = sub.SingularStructure()
        reduced_squared_masses = {}
        parents = []
        for substructure in singular_structure.substructures:
            parent, children, _ = get_structure_numbers(substructure, momenta_dict)
            for child in children: del reduced_PS_point[child]
            reduced_PS_point[parent] = sum(PS_point[child] for child in children)
            parent_leg = sub.SubtractionLeg(parent, 0, sub.SubtractionLeg.FINAL)
            reduced_singular_structure.substructures.append(sub.CollStructure(parent_leg))
            m2i = 'm2' + str(parent)
            if squared_masses is None:
                reduced_squared_masses[m2i] = 0
            else:
                reduced_squared_masses[m2i] = squared_masses[m2i]
            parents.append(parent)
            if kinematic_variables is not None:
                kinematic_variables['s' + str(parent)] = reduced_PS_point[parent].square()
        for leg in singular_structure.legs:
            # reduced_PS_point[leg.n] = LorentzVector(PS_point[leg.n])
            reduced_singular_structure.substructures.append(sub.CollStructure(leg))
            reduced_squared_masses['m2'+str(leg.n)] = PS_point[leg.n].square()
        # Perform mapping of parent momenta to target masses
        new_PS_point, vars = FinalMassesMapping.map_to_lower_multiplicity(
            reduced_PS_point, reduced_singular_structure,
            momenta_dict, reduced_squared_masses,
            kinematic_variables, compute_jacobian )
        # If need be, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            for parent in parents:
                na, nb = FinalCollinearVariables.collinear_and_reference(
                    new_PS_point[parent] )
                FinalCollinearVariables.get(
                    PS_point, momenta_dict[parent], na, nb, kinematic_variables)
        return new_PS_point, vars

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Build reduced PS point, singular structure & masses to call FinalMassesMapping
        reduced_singular_structure = sub.SingularStructure()
        reduced_kinematic_variables = dict()
        parents = []
        for substructure in singular_structure.substructures:
            parent, _, _ = get_structure_numbers(substructure, momenta_dict)
            parent_leg = sub.SubtractionLeg(parent, 0, sub.SubtractionLeg.FINAL)
            reduced_singular_structure.substructures.append(sub.CollStructure(parent_leg))
            s_i = 's' + str(parent)
            reduced_kinematic_variables[s_i] = kinematic_variables[s_i]
            parents.append(parent)
        for leg in singular_structure.legs:
            reduced_singular_structure.substructures.append(sub.CollStructure(leg))
            reduced_kinematic_variables['s' + str(leg.n)] = PS_point[leg.n].square()
        # Perform mapping of parent momenta to target masses
        new_PS_point, vars = FinalMassesMapping.map_to_higher_multiplicity(
            PS_point, reduced_singular_structure, momenta_dict,
            reduced_kinematic_variables, compute_jacobian )
        # Set the kinematic_variables
        for parent in parents:
            children = momenta_dict[parent]
            na, nb = FinalCollinearVariables.collinear_and_reference(PS_point[parent])
            FinalCollinearVariables.set(
                new_PS_point, parent, children, na, nb, kinematic_variables)
        return new_PS_point, vars

# Final-collinear Lorentz mapping
#=========================================================================================

class FinalLorentzMapping(FinalCollinearMapping):
    """Implementation of the Lorentz mapping for multiple collinear sets."""

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Build reduced PS point, singular structure & masses to call FinalMassesMapping
        reduced_PS_point = PS_point.get_copy()
        reduced_singular_structure = sub.SingularStructure()
        reduced_squared_masses = {}
        parents = []
        Q = LorentzVector()
        # First build pseudo-particles for the collinear sets
        for substructure in singular_structure.substructures:
            parent, children, _ = get_structure_numbers(substructure, momenta_dict)
            for child in children: del reduced_PS_point[child]
            reduced_PS_point[parent] = sum(PS_point[child] for child in children)
            Q += reduced_PS_point[parent]
            parent_leg = sub.SubtractionLeg(parent, 0, sub.SubtractionLeg.FINAL)
            reduced_singular_structure.substructures.append(sub.CollStructure(parent_leg))
            m2i = 'm2' + str(parent)
            if squared_masses is None:
                reduced_squared_masses[m2i] = 0
            else:
                reduced_squared_masses[m2i] = squared_masses[m2i]
            parents.append(parent)
            if kinematic_variables is not None:
                kinematic_variables['s' + str(parent)] = reduced_PS_point[parent].square()
        mass_sum = sum(s ** 0.5 for s in reduced_squared_masses.values())
        # Then treat recoilers collectively if any
        recoilers = [leg.n for leg in singular_structure.legs]
        R = None
        s_R = 0.
        if recoilers:
            for recoiler in recoilers:
                del reduced_PS_point[recoiler]
            if len(recoilers) == 1:
                R = recoilers[0]
            else:
                R = max(momenta_dict.keys()) + 1
            reduced_PS_point[R] = LorentzVector()
            for recoiler in recoilers:
                reduced_PS_point[R] += PS_point[recoiler]
            recoiler_leg = sub.SubtractionLeg(R, 0, sub.SubtractionLeg.FINAL)
            reduced_singular_structure.substructures.append(
                sub.CollStructure(recoiler_leg))
            s_R = reduced_PS_point[R].square()
            reduced_squared_masses['m2'+str(R)] = s_R
            Q += reduced_PS_point[R]
        # If the mapping is not valid, raise
        Q2 = Q.square()
        if (abs(s_R) ** 0.5) > (Q2 ** 0.5 - mass_sum):
            raise FailedMapping
        # Perform mapping of parent momenta to target masses
        new_PS_point, vars = FinalMassesMapping.map_to_lower_multiplicity(
            reduced_PS_point, reduced_singular_structure,
            momenta_dict, reduced_squared_masses,
            kinematic_variables, compute_jacobian )
        # If need be, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            for parent in parents:
                na, nb = FinalCollinearVariables.collinear_and_reference(
                    new_PS_point[parent] )
                FinalCollinearVariables.get(
                    PS_point, momenta_dict[parent], na, nb, kinematic_variables)
        # Boost recoilers
        if len(recoilers) > 1:
            pR = reduced_PS_point[R]
            qR = new_PS_point[R]
            for recoiler in recoilers:
                new_PS_point[recoiler] = LorentzVector(PS_point[recoiler])
                # TODO Move this try/except to higher level
                try:
                    new_PS_point[recoiler].rotoboost(pR, qR)
                    # new_PS_point[recoiler].boost_from_to(pR, qR)
                except:
                    logger.critical(
                        "Problem encountered for " + str(singular_structure))
                    logger.critical("The full phase space point was\n" + str(PS_point))
            del new_PS_point[R]
        return new_PS_point, vars

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Build reduced PS point, singular structure & masses to call FinalMassesMapping
        reduced_PS_point = PS_point.get_copy()
        reduced_singular_structure = sub.SingularStructure()
        reduced_kinematic_variables = dict()
        parents = []
        Q = LorentzVector()
        # First build pseudo-particles for the collinear sets
        for substructure in singular_structure.substructures:
            parent, _, _ = get_structure_numbers(substructure, momenta_dict)
            parent_leg = sub.SubtractionLeg(parent, 0, sub.SubtractionLeg.FINAL)
            reduced_singular_structure.substructures.append(sub.CollStructure(parent_leg))
            s_i = 's' + str(parent)
            reduced_kinematic_variables[s_i] = kinematic_variables[s_i]
            parents.append(parent)
            Q += PS_point[parent]
        # Then treat recoilers collectively if any
        recoilers = [leg.n for leg in singular_structure.legs]
        mass_sum = sum(s ** 0.5 for s in reduced_kinematic_variables.values())
        R = None
        qR = None
        s_R = 0.
        if recoilers:
            if len(recoilers) == 1:
                R = recoilers[0]
            else:
                R = max(momenta_dict.keys()) + 1
                qR = sum(PS_point[r] for r in recoilers)
                reduced_PS_point[R] = LorentzVector(qR)
            recoiler_leg = sub.SubtractionLeg(R, 0, sub.SubtractionLeg.FINAL)
            reduced_singular_structure.substructures.append(
                sub.CollStructure(recoiler_leg) )
            Q += reduced_PS_point[R]
            s_R = reduced_PS_point[R].square()
            reduced_kinematic_variables['s'+str(R)] = s_R
        # If the mapping is not valid, undo all changes and raise
        Q2 = Q.square()
        if (s_R ** 0.5) > (Q2 ** 0.5 - mass_sum):
            raise FailedMapping
        # Perform mapping of parent momenta to target masses
        new_PS_point, vars = FinalMassesMapping.map_to_higher_multiplicity(
            reduced_PS_point, reduced_singular_structure, momenta_dict,
            reduced_kinematic_variables, compute_jacobian )
        # Set the kinematic_variables
        for parent in parents:
            children = momenta_dict[parent]
            na, nb = FinalCollinearVariables.collinear_and_reference(PS_point[parent])
            FinalCollinearVariables.set(
                new_PS_point, parent, children, na, nb, kinematic_variables)
        # Boost recoilers
        if len(recoilers) > 1:
            pR = new_PS_point[R]
            for recoiler in recoilers:
                # TODO Move this try/except to higher level
                try:
                    new_PS_point[recoiler].rotoboost(qR, pR)
                    # new_PS_point[recoiler].boost_from_to(qR, pR)
                except:
                    logger.critical(
                        "Problem encountered for %s" % str(singular_structure))
                    logger.critical("The full phase space point was\n%s" % str(PS_point))
            del new_PS_point[R]
        return new_PS_point, vars

    @classmethod
    def can_map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables, ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict) )
        assert needed_variables.issubset(kinematic_variables.keys())

        Q = LorentzVector()
        mass_sum = 0.
        for substructure in singular_structure.substructures:
            parent, _, _ = get_structure_numbers(substructure, momenta_dict)
            Q += PS_point[parent]
            sqr = kinematic_variables['s'+str(parent)]
            mass_sum += sqr ** 0.5
        # Then treat recoilers collectively if any
        recoilers = [leg.n for leg in singular_structure.legs]
        qR = LorentzVector()
        if recoilers:
            qR = sum(PS_point[r] for r in recoilers)
        Q += qR
        qR2 = qR.square()
        Q2 = Q.square()
        return (qR2 ** 0.5) <= (Q2 ** 0.5 - mass_sum)

    @classmethod
    def can_map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        Q = LorentzVector()
        mass_sum = 0.
        for substructure in singular_structure.substructures:
            parent, _, _ = get_structure_numbers(substructure, momenta_dict)
            Q += PS_point[parent]
            if squared_masses is None:
                sqr = 0
            else:
                sqr = squared_masses['m2' + str(parent)]
            mass_sum += sqr ** 0.5
        # Then treat recoilers collectively if any
        recoilers = [leg.n for leg in singular_structure.legs]
        qR = LorentzVector()
        if recoilers:
            qR = sum(PS_point[r] for r in recoilers)
        Q += qR
        qR2 = qR.square()
        Q2 = Q.square()
        return (qR2 ** 0.5) <= (Q2 ** 0.5 - mass_sum)

#=========================================================================================
# Initial-collinear mappings
#=========================================================================================

# Generic initial-collinear mapping parent class
#=========================================================================================

class InitialCollinearMapping(VirtualMapping):
    """Common functions for initial-collinear elementary mappings."""

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, sub.SingularStructure)
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
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        names = []
        for substructure in singular_structure.substructures:
            parent, fs_children, is_child = get_structure_numbers(
                substructure, momenta_dict )
            InitialCollinearVariables.names(parent, fs_children, is_child)
        return names

    @classmethod
    def rescale_kinematic_variables(
        cls, singular_structure, momenta_dict, kinematic_variables, scaling_parameter):

        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict))
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        _, fs_children, is_child = get_structure_numbers(substructure, momenta_dict)
        # Determine the correct scaling for the divergence to go like 1/parameter
        base = scaling_parameter ** (0.5 / len(fs_children))
        
        kinematic_variables['s' + str(is_child)] *= base ** 2
        kinematic_variables['kt' + str(is_child)] *= base
        for child in fs_children:
            kinematic_variables['kt' + str(child)] *= base
        return kinematic_variables

# Initial-collinear Lorentz mapping, one set
#=========================================================================================

class InitialLorentzOneMapping(InitialCollinearMapping):
    """Implementation of the Lorentz transformation mapping
    for one set of collinear particles (with massless parent)
    and an arbitrary number of (eventually massive) recoilers.
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only for one set of particles going collinear to an initial-state parton,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(InitialLorentzOneMapping, cls).is_valid_structure(
                singular_structure )
        return False

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Precompute sets and numbers
        substructure = singular_structure.substructures[0]
        parent, fs_children, is_child = get_structure_numbers(substructure, momenta_dict)
        assert (squared_masses is None) or (squared_masses['m2' + str(parent)] == 0)
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
        # Create new PS point
        new_PS_point = PS_point.get_copy()
        # Map the set's momentum
        qA = xia * pa
        new_PS_point[parent] = qA
        # Map all recoilers' momenta
        qR = qA - pAmpR
        for recoiler in singular_structure.legs:
            # TODO Move this try/except to higher level
            try:
                new_PS_point[recoiler.n].rotoboost(pR, qR)
                # new_PS_point[recoiler.n].boost_from_to(pR, qR)
            except:
                logger.critical("Problem encountered for " + str(singular_structure))
                logger.critical("The full phase space point was\n" + str(PS_point))
        # Eliminate children momenta from the mapped phase-space point
        for j in fs_children:
            del new_PS_point[j]
        if is_child != parent:  # Bypass degenerate case of 1->1 splitting
            del new_PS_point[is_child]
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            na, nb = InitialCollinearVariables.collinear_and_reference(qA)
            InitialCollinearVariables.get(
                PS_point, fs_children, is_child, na, nb, kinematic_variables )
        # TODO Check if the jacobian for this mapping is really 1
        jacobian = 1.
        

        # WARNING: for ISR counterterms, the quantity Q below may *not* match with 
        # PS_point[1] + PS_point[2]. So it should be used only for what it actually is!
        Q = pAmpR
        
        mapping_variables = {'jacobian': jacobian, 'Q': Q}
        # Return characteristic variables
        return new_PS_point, mapping_variables

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False):

        #misc.sprint("Mapping up the PS point:\n", str(PS_point))
        #misc.sprint("with momenta dict: %s",momenta_dict)
        #misc.sprint("with variables:\n", kinematic_variables)
        #misc.sprint("and structure", singular_structure)

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
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
        # Create new PS point
        new_PS_point = PS_point.get_copy()
        del new_PS_point[parent]
        new_PS_point[is_child] = pa
        # Set children momenta
        InitialCollinearVariables.set(
            new_PS_point, is_child, fs_children, na, nb, kinematic_variables )
        # Build collective momenta
        pCa = LorentzVector()
        for j in fs_children:
            pCa += new_PS_point[j]
        pA = pa - pCa
        # Map recoil momenta
        pR = qR + pA - qA
        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n].rotoboost(qR, pR)
            # new_PS_point[recoiler.n].boost_from_to(qR, pR)
        # TODO Check if the jacobian for this mapping is really 1
        jacobian = 1.

        # VH :: START :: TEMPORARY FIX OF MODIFYING Q DEFINITION
        #Q = -qRmqA
        Q = new_PS_point[1] + new_PS_point[2]
        # VH :: END :: TEMPORARY FIX OF MODIFYING Q DEFINITION
        
        mapping_variables = {'jacobian': jacobian, 'Q': -qRmqA}
        # Return characteristic variables
        return new_PS_point, mapping_variables

#=========================================================================================
# Soft mappings
#=========================================================================================

# Generic soft mapping parent class
#=========================================================================================

class ElementaryMappingSoft(VirtualMapping):
    """Common functions for soft elementary mappings."""

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, sub.SingularStructure)
        # Valid only for bunches of final-state particles going soft,
        # with no recursive substructure
        for substructure in singular_structure.substructures:
            if not substructure.name() == "S":
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
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

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

        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict))
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers
        base = scaling_parameter ** (1. / len(singular_structure.substructures))
        for substructure in singular_structure.substructures:
            _, children, _ = get_structure_numbers(substructure, momenta_dict)
            # Determine the correct scaling for the divergence to go like 1/parameter
            this_base = base ** (0.5 / len(children))
            for child in children:
                kinematic_variables['p' + str(child)] *= this_base
        return kinematic_variables

# Final mapping to/from zero masses with soft
#=========================================================================================

class FinalAssociativeSoftMappingZero(ElementaryMappingSoft):
    """Mapping that removes momenta from a set of massless ones,
    producing momenta that are still massless.
    """

    plot = False

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only if there are at least two recoilers in the final state
        if len(singular_structure.legs) < 2:
            return False
        if singular_structure.legs.has_initial_state_leg():
            return False
        return super(FinalAssociativeSoftMappingZero, cls).is_valid_structure(
            singular_structure)

    @classmethod
    def get_kinematic_variables_names(cls, singular_structure, momenta_dict):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        names = []
        for substructure in singular_structure.substructures:
            _, children, _ = get_structure_numbers(substructure, momenta_dict)
            names += SoftVariables.names(children)
        return names

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Determine leg numbers and check that target masses are zero
        pS = LorentzVector()
        recoilers = [leg.n for leg in singular_structure.legs]
        all_soft_children = []
        new_PS_point = PS_point.get_copy()
        for substructure in singular_structure.substructures:
            _, children, _ = get_structure_numbers(substructure, momenta_dict)
            all_soft_children.append(children)
            pS += sum(new_PS_point.pop(child) for child in children)
        # Build total momentum
        Q = LorentzVector(pS)
        for r in recoilers:
            assert PS_point[r].square_almost_zero()
            Q += PS_point[r]
        Q2 = Q.square()
        # Compute the parameters ys, y, ybar, mubar2, ytilde, alpha
        ps = pS / len(recoilers)
        ys = Q.dot(ps) / Q2
        y = {r: (Q.dot(new_PS_point[r]) / Q2) for r in recoilers}
        ybar = {r: (y[r] + ys) for r in recoilers}
        assert abs(sum(ybar.values()) - 1.) < 1.e-6
        mubar2 = {r: ((new_PS_point[r] + ps).square() / Q2) for r in recoilers}
        ytilde = {r: (ybar[r] ** 2 - mubar2[r]) ** 0.5 for r in recoilers}
        alpha = sum(ytilde.values())
        # misc.sprint(alpha)
        for r in recoilers:
            ytilde[r] /= alpha
        assert abs(sum(ytilde.values())  - 1.) < 1.e-6
        # Map all momenta
        for r in recoilers:
            new_PS_point[r] = (new_PS_point[r] + ps - ybar[r] * Q) / alpha + ytilde[r] * Q
            assert abs(new_PS_point[r].square() / Q2) < 1.e-6
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            for soft_children in all_soft_children:
                SoftVariables.get(PS_point, soft_children, kinematic_variables)
        # Compute the jacobian for this mapping
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        sum_deriv = 0.
        prod_ytilde_y = 1.
        for r in recoilers:
            space_scalar_prod = ytilde[r] * ys - ps.dot(new_PS_point[r]) / Q2
            sum_deriv += (alpha * ytilde[r]**2 - space_scalar_prod) / y[r]
            prod_ytilde_y *= ytilde[r] / y[r]
        jacobian = alpha**(3*len(recoilers)-4) * prod_ytilde_y / sum_deriv
        mapping_variables = {'jacobian': jacobian, 'Q': Q}
        # Return characteristic variables
        return new_PS_point, mapping_variables

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Determine leg numbers
        pS = LorentzVector()
        recoilers = [leg.n for leg in singular_structure.legs]
        all_soft_children = []
        new_PS_point = PS_point.get_copy()
        for substructure in singular_structure.substructures:
            _, children, _ = get_structure_numbers(substructure, momenta_dict)
            all_soft_children.append(children)
            SoftVariables.set(new_PS_point, children, kinematic_variables)
            pS += sum(new_PS_point[child] for child in children)
        # Build total momentum
        Q = LorentzVector()
        for r in recoilers:
            assert PS_point[r].square_almost_zero()
            Q += PS_point[r]
        Q2 = Q.square()
        # Compute the parameters
        ps = pS / len(recoilers)
        # misc.sprint(ps)
        yS = pS.dot(Q) / Q2
        ys = yS / len(recoilers)
        mus2 = ps.square() / Q2
        ytilde = {}
        a = {}
        b = {}
        c = {}
        for r in recoilers:
            ytilde[r] = Q.dot(PS_point[r]) / Q2
            a[r] = ytilde[r] ** 2
            b[r] = -2. * (ytilde[r] * ys - ps.dot(PS_point[r]) / Q2)
            c[r] = ys ** 2 - mus2
        assert abs(sum(ytilde.values())  - 1.) < 1.e-6
        # if len(rs) == 2:
        #     alpha = (Kaellen(1., mu2[js[0]], mu2[js[1]])) ** 0.5
        # else:
        # Solve the equation for alpha numerically
        y_fun = lambda r, x: (a[r]*(x**2) + b[r]*x + c[r])**0.5
        from scipy.optimize import newton
        fun = lambda x: sum(y_fun(r, x) for r in recoilers) - (1. - yS)
        # Optional plot to monitor the numerical solution
        if cls.plot:
            import numpy as np
            import matplotlib.pyplot as plt
            zero = lambda a: 0.*a
            tau = np.linspace(0., 1., 501)
            plt.plot(tau, fun(tau), "-", tau, zero(tau), "-")
            plt.grid()
            plt.show()
        alpha = newton(fun, 1.)
        # misc.sprint(alpha)
        # Compute the parameters gamma
        y = {r: y_fun(r, alpha) for r in recoilers}
        assert abs(sum(y.values()) - (1. - yS)) < 1.e-6
        # Map all momenta
        for r in recoilers:
            new_PS_point[r] *= alpha
            new_PS_point[r] += (y[r] + ys - alpha*ytilde[r]) * Q - ps
            assert abs(new_PS_point[r].square() / Q2) < 1.e-6
        # Compute the jacobian for this mapping
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        sum_deriv = 0.
        prod_ytilde_y = 1.
        for r in recoilers:
            space_scalar_prod = ytilde[r] * ys - ps.dot(PS_point[r]) / Q2
            sum_deriv += (alpha * ytilde[r]**2 - space_scalar_prod) / y[r]
            prod_ytilde_y *= ytilde[r] / y[r]
        jacobian = alpha**(3*len(recoilers)-4) * prod_ytilde_y / sum_deriv
        mapping_variables = {'jacobian': jacobian, 'Q': Q}
        # Return characteristic variables
        return new_PS_point, mapping_variables

class FinalAssociativeSoftMapping(FinalAssociativeSoftMappingZero):
    """Associative soft mapping with final-state recoilers."""

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        recoilers = [leg.n for leg in singular_structure.legs]
        target_masses = {'s' + str(r): PS_point[r].square() for r in recoilers}
        if kinematic_variables is not None:
            target_masses.update(kinematic_variables)
        singular_structure_to_zero = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[leg]) for leg in singular_structure.get_all_legs() ])
        singular_structure_from_zero = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[leg]) for leg in singular_structure.legs ])
        zero_PS_point, zero_vars = FinalZeroMassesMapping.map_to_lower_multiplicity(
            PS_point, singular_structure_to_zero, momenta_dict,
            kinematic_variables=kinematic_variables, compute_jacobian=compute_jacobian )
        soft_PS_point, soft_vars = FinalAssociativeSoftMappingZero.map_to_lower_multiplicity(
            zero_PS_point, singular_structure, momenta_dict,
            kinematic_variables=kinematic_variables, compute_jacobian=compute_jacobian )
        mass_PS_point, mass_vars = FinalZeroMassesMapping.map_to_higher_multiplicity(
            soft_PS_point, singular_structure_from_zero, momenta_dict,
            kinematic_variables=target_masses, compute_jacobian=compute_jacobian )
        if compute_jacobian:
            zero_vars['jacobian'] *= soft_vars['jacobian'] / mass_vars['jacobian']
        # Return characteristic variables
        return mass_PS_point, zero_vars

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # recoilers = [leg.n for leg in singular_structure.legs]
        # target_masses = {'s' + str(r): PS_point[r].square() for r in recoilers}
        singular_structure_from_zero = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[leg]) for leg in singular_structure.get_all_legs() ])
        singular_structure_to_zero = sub.SingularStructure(substructures=[
            sub.CollStructure(legs=[leg]) for leg in singular_structure.legs ])
        zero_PS_point, zero_vars = FinalZeroMassesMapping.map_to_lower_multiplicity(
            PS_point, singular_structure_to_zero, momenta_dict,
            compute_jacobian=compute_jacobian )
        soft_PS_point, soft_vars = FinalAssociativeSoftMappingZero.map_to_higher_multiplicity(
            zero_PS_point, singular_structure, momenta_dict,
            kinematic_variables=kinematic_variables, compute_jacobian=compute_jacobian )
        mass_PS_point, mass_vars = FinalZeroMassesMapping.map_to_higher_multiplicity(
            soft_PS_point, singular_structure_from_zero, momenta_dict,
            kinematic_variables=kinematic_variables, compute_jacobian=compute_jacobian )
        if compute_jacobian:
            mass_vars['jacobian'] *= soft_vars['jacobian'] / zero_vars['jacobian']
        # Return characteristic variables
        return mass_PS_point, mass_vars

# Soft mapping, massless final-state recoilers
#=========================================================================================
class SoftVsFinalPureRescalingMapping(ElementaryMappingSoft):
    """Implementation of the mapping used by Somogyi et al.
    in arXiv:hep-ph/0609042 for soft particles.
    It applies when there are at least two recoilers in the final state
    and all recoilers are massless.
    """

    @staticmethod
    def y(pS, Q):
        """Return the parameter y of the SoftVsFinal mapping."""
        pR = Q - pS
        pR2_Q2 = pR.square() / Q.square()
        return 1. - pR2_Q2

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only if there are at least two recoilers (assumed in the final state)
        if len(singular_structure.legs) < 2: return False
        return super(SoftVsFinalPureRescalingMapping, cls).is_valid_structure(singular_structure)

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        # Build the total soft momentum,
        # save the soft momenta in kinematic_variables
        # and eliminate them from new_PS_point
        new_PS_point = PS_point.get_copy()
        pS = LorentzVector()
        for substructure in singular_structure.substructures:
            children = tuple(leg.n for leg in substructure.legs)
            if kinematic_variables is not None:
                SoftVariables.get(PS_point, children, kinematic_variables)
            for child in children:
                pS += new_PS_point.pop(child)
        # Build the total momentum of recoilers
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        pR = LorentzVector()
        for recoiler in recoilers:
            pR += PS_point[recoiler]
        # Build the total momentum Q
        Q = pS + pR
        # Compute the parameter la
        pR2_Q2 = pR.square() / Q.square()
        la = math.sqrt(pR2_Q2)
        P = pR / la
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n] /= la
            new_PS_point[recoiler.n].rotoboost(P, Q)
            # new_PS_point[recoiler.n].boost_from_to(P, Q)
        mapping_variables = {'Q': Q}
        if compute_jacobian:
            mapping_variables['jacobian'] = (pR2_Q2)**(len(recoilers)-2)
        # Return characteristic variables
        return new_PS_point, mapping_variables

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict ) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Build the total soft momentum,
        # get the soft momenta from kinematic_variables
        # and save them in new_PS_point
        new_PS_point = PS_point.get_copy()
        pS = LorentzVector()
        for substructure in singular_structure.substructures:
            children = tuple(leg.n for leg in substructure.legs)
            SoftVariables.set(new_PS_point, children, kinematic_variables)
            for child in children:
                pS += new_PS_point[child]
        # Build the total momentum, which is equal to the mapped recoilers'
        Q = LorentzVector()
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        for recoiler in recoilers:
            Q += PS_point[recoiler]
        # Build the recoilers' momentum
        pR = Q - pS
        # Compute the parameter la
        pR2_Q2 = pR.square() / Q.square()
        la = math.sqrt(pR2_Q2)
        P = pR / la
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n] *= la
            new_PS_point[recoiler.n].rotoboost(Q, P)
            # new_PS_point[recoiler.n].boost_from_to(Q, P)
        mapping_variables = {'Q': Q}
        if compute_jacobian:
            mapping_variables['jacobian'] = (pR2_Q2)**(len(recoilers)-2)
        # Return characteristic variables
        return new_PS_point, mapping_variables

# Soft mapping with final-state recoilers, based on rescaling + boost
#=========================================================================================
class SoftVsFinalRescaleThenBoostMapping(ElementaryMappingSoft):
    """Generalisation of the soft mapping used by Somogyi et al. in arXiv:hep-ph/0609042.
    It combines a generalised rescaling and a Lorentz transformation,
    and it applies when there are at least two recoilers in the final state.
    """

    plot = False

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only if there are at least two recoilers (assumed in the final state)
        if len(singular_structure.legs) < 2: return False
        mysuper = super(SoftVsFinalRescaleThenBoostMapping, cls) 
        return mysuper.is_valid_structure(singular_structure)

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError(
                "Singular structure '%s' is not supported by mapping '%s'"
                % (str(singular_structure), cls.__name__) )

        # Build the total soft momentum,
        # save the soft momenta in kinematic_variables
        # and eliminate them from new_PS_point
        new_PS_point = PS_point.get_copy()
        pS = LorentzVector()
        for substructure in singular_structure.substructures:
            children = tuple(leg.n for leg in substructure.legs)
            if kinematic_variables is not None:
                SoftVariables.get(PS_point, children, kinematic_variables)
            for child in children:
                pS += new_PS_point.pop(child)
        # Build the total momentum of recoilers
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        pR = LorentzVector()
        for recoiler in recoilers:
            pR += PS_point[recoiler]
        # Build the total momentum Q
        Q = pS + pR
        Q2 = Q.square()
        # Compute the parameters
        ER = Q.dot(pR)/Q2
        pRe = pR - ER * Q
        pRe2 = -pRe.square() / Q2
        Ers = {r: (Q.dot(PS_point[r]) / Q2) for r in recoilers}
        pres = {r: PS_point[r] - Ers[r] * Q for r in recoilers}
        pre2s = {r: -pres[r].square() / Q2 for r in recoilers}
        mr2s = {r: PS_point[r].square() / Q2 for r in recoilers}
        # Solve the equation for alpha numerically
        from scipy.optimize import newton
        func = lambda a2: 1. + pRe2 * a2 - sum(
            (pre2s[r] * a2 + mr2s[r]) ** 0.5 for r in recoilers) ** 2
        alpha2 = newton(func, 1.)
        alpha = alpha2 ** 0.5
        # misc.sprint(alpha2)
        # Optional plot to monitor the numerical solution
        if cls.plot:
            import matplotlib.pyplot as plt
            import numpy as np
            zero = lambda a: 0.*a
            tau = np.linspace(1., 100., 101)
            plt.plot(tau, func(tau), "-", tau, zero(tau), "-")
            plt.grid()
            plt.show()
        assert alpha2 >= 1.
        # Compute the rescaled momenta
        prs_new = {
            r: pres[r] * alpha + ((pre2s[r] * alpha2 + mr2s[r]) ** 0.5) * Q
            for r in recoilers}
        pR_new = sum(prs_new.values())
        pR2_new = pR_new.square()
        assert abs(pR2_new - Q2) < pR_new.eps
        # Map all momenta
        for r in recoilers:
            new_PS_point[r] = prs_new[r].rotoboost(pR_new, Q)
            # new_PS_point[r] = prs_new[r].boost_from_to(pR_new, Q)
            assert abs(new_PS_point[r].square() / Q2 - mr2s[r]) < 1.e-6
        # Compute the jacobian for this mapping
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        raise NotImplementedError

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        raise NotImplementedError

# Soft mapping with final-state recoilers, based on boost + rescaling
#=========================================================================================
class SoftVsFinalBoostThenRescaleMapping(ElementaryMappingSoft):
    """Generalisation of the soft mapping used by Somogyi et al. in arXiv:hep-ph/0609042.
    It combines a generalized a Lorentz transformation and a rescaling,
    and it applies when there are at least two recoilers in the final state.
    """

    plot = False

    @classmethod
    def is_valid_structure(cls, singular_structure):

        # Valid only if there are at least two recoilers (assumed in the final state)
        if len(singular_structure.legs) < 2: return False
        mysuper = super(SoftVsFinalBoostThenRescaleMapping, cls)
        return mysuper.is_valid_structure(singular_structure)

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError(
                "Singular structure '%s' is not supported by mapping '%s'"
                % (str(singular_structure), cls.__name__) )

        # Build the total soft momentum,
        # save the soft momenta in variables and eliminate them from PS_point
        new_PS_point = PS_point.get_copy()
        pS = LorentzVector()
        for substructure in singular_structure.substructures:
            children = tuple(leg.n for leg in substructure.legs)
            if kinematic_variables is not None:
                SoftVariables.get(PS_point, children, kinematic_variables)
            for child in children:
                pS += new_PS_point.pop(child)
        # Build the total momentum of recoilers
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        pR = LorentzVector()
        for recoiler in recoilers:
            pR += PS_point[recoiler]
        pR2 = pR.square()
        # Build the total momentum Q
        Q = pS + pR
        Q2 = Q.square()
        # Boost all momenta
        pRhat = ((pR2 / Q2) ** 0.5) * Q
        for r in recoilers:
            new_PS_point[r] = new_PS_point[r].rotoboost(pR, pRhat)
            # new_PS_point[r] = new_PS_point[r].boost_from_to(pR, pRhat)
        # Rescale all three momenta
        Ers = {r: (Q.dot(new_PS_point[r]) / Q2) for r in recoilers}
        pres = {r: new_PS_point[r] - Ers[r] * Q for r in recoilers}
        pre2s = {r: -pres[r].square() / Q2 for r in recoilers}
        mr2s = {r: PS_point[r].square() / Q2 for r in recoilers}
        # Solve the equation for alpha numerically
        from scipy.optimize import newton
        func = lambda a2: 1. - sum((pre2s[r] * a2 + mr2s[r]) ** 0.5 for r in recoilers)
        alpha2 = newton(func, 1.)
        alpha = alpha2 ** 0.5
        misc.sprint(alpha2)
        # Optional plot to monitor the numerical solution
        if cls.plot:
            import matplotlib.pyplot as plt
            import numpy as np
            zero = lambda a: 0.*a
            tau = np.linspace(1., 100., 101)
            plt.plot(tau, func(tau), "-", tau, zero(tau), "-")
            plt.grid()
            plt.show()
        assert alpha2 >= 1.
        # Compute the rescaled momenta
        for r in recoilers:
            new_PS_point[r] = pres[r] * alpha + ((pre2s[r] * alpha2 + mr2s[r]) ** 0.5) * Q
        # Compute the jacobian for this mapping
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        raise NotImplementedError

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError(
                "Singular structure '%s' is not supported by mapping '%s'"
                % (str(singular_structure), cls.__name__) )
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict ) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Build the total soft momentum,
        # get the soft momenta from kinematic_variables
        # and save them in new_PS_point
        new_PS_point = PS_point.get_copy()
        pS = LorentzVector()
        for substructure in singular_structure.substructures:
            children = tuple(leg.n for leg in substructure.legs)
            SoftVariables.set(new_PS_point, children, kinematic_variables)
            for child in children:
                pS += new_PS_point[child]
        # Build the total momentum, which is equal to the mapped recoilers'
        Q = LorentzVector()
        recoilers = tuple(leg.n for leg in singular_structure.legs)
        for recoiler in recoilers:
            Q += PS_point[recoiler]
        Q2 = Q.square()
        # Build the recoilers' momentum
        pR = Q - pS
        pR2 = pR.square()
        # Compute the parameters
        pR2_Q2 = pR2 / Q2
        Erts = {r: (Q.dot(PS_point[r]) / Q2) for r in recoilers}
        prtes = {r: PS_point[r] - Erts[r] * Q for r in recoilers}
        prte2s = {r: -prtes[r].square() / Q2 for r in recoilers}
        mr2s = {r: PS_point[r].square() / Q2 for r in recoilers}
        # Solve the equation for alpha numerically
        from scipy.optimize import newton
        func = lambda a2: pR2_Q2 ** 0.5 - sum(
            (prte2s[r] / a2 + mr2s[r]) ** 0.5 for r in recoilers)
        alpha2 = newton(func, 1.)
        alpha = alpha2 ** 0.5
        misc.sprint(alpha2)
        # Optional plot to monitor the numerical solution
        if cls.plot:
            import matplotlib.pyplot as plt
            import numpy as np
            zero = lambda a: 0.*a
            tau = np.linspace(1., 100., 101)
            plt.plot(tau, func(tau), "-", tau, zero(tau), "-")
            plt.grid()
            plt.show()
        assert alpha2 >= 1.
        # Compute the rescaled momenta
        for r in recoilers:
            new_PS_point[r] = prtes[r] / alpha + ((prte2s[r] / alpha2 + mr2s[r]) ** 0.5) * Q
        # Boost back the rescaled momenta
        Qresc = (pR2_Q2 ** 0.5) * Q
        for r in recoilers:
            new_PS_point[r] = new_PS_point[r].rotoboost(Qresc, pR)
            # new_PS_point[r] = new_PS_point[r].boost_from_to(Qresc, pR)
        # Compute the jacobian for this mapping
        if not compute_jacobian:
            return new_PS_point, {'Q': Q}
        raise NotImplementedError

# Final collinear mapping recoiling against initial states, suited for p p > X collisions.
#=========================================================================================
class FinalCollinearVsInitialMapping(FinalCollinearMapping):
    """Implementation of the mapping described in Eqs. 2.13-2.15 of Gabor's notes (in the FF section).
    It is applicable for hadronic collisions where the initial states can be rescaled.
    - is_valid_structure is entirely inherited from the parent.
    This mapping works for any hadron-hadron collider.
    The initial state should be checked at initialization
    """

    @classmethod
    def is_valid_structure(cls, singular_structure):

        assert isinstance(singular_structure, sub.SingularStructure)
        # Valid only for sets of final-state particles going collinear,
        # with no recursive substructure
        # Also for now we support only on substructure
        if len(singular_structure.substructures)>1:
            return False
        for substructure in singular_structure.substructures:
            if not substructure.name() == "C":
                return False
            if substructure.substructures:
                return False
            if substructure.get_all_legs().has_initial_state_leg():
                return False
        return True

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None, kinematic_variables=None, compute_jacobian=False ):
        """Map the resolved momenta of a real emission process to on-shell kinematics by recoling against the initial state partons.

        :param PS_point: dictionary of Real kinematics
        :param singular_structure: SingularStructure describing the limit we want to factorize with this mapping
        :param momenta_dict: Routing of momenta label
        :param squared_masses: dictionary
        :param kinematic_variables: dictionary to store the variables of the unresolved PS
        :param compute_jacobian: Boolean, specify that the jacobian should be computed if it is necessary to divide by it in the local CT.
        :return:
            * new_PS_point: a dictionary with Born kinematics
            * mapping_variables: a dictionary with variables, possible used in the current evaluation.
        """

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Recoilers must be the two initial state legs:
        if len(singular_structure.legs)!=2 and not all(leg.state==leg.INITIAL for leg in singular_structure.legs):
            raise MappingError("Singular structure '%s' legs should be two initial states when using the mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Construct a new phase space point
        new_PS_point = PS_point.get_copy()

        # Build the total collinear momentum and remove all collinear legs from the mapped PS point
        # Below must be upgraded for nested collinears
        pC = LorentzVector()
        collinear_children = set([])
        for substructure in singular_structure.substructures: # S(3,4).legs returns [3,4]
            children = tuple(leg.n for leg in substructure.legs)
            for child in children:
                pC += new_PS_point.pop(child)
            collinear_children |= set(children)

        # Build the total momentum Q from the initial states. This is of course only applicable for 2>N kinematics.
        # Fetch the initial states numbers using the recoiler ids.
        Q = LorentzVector()
        Q += PS_point[singular_structure.legs[0].n]
        Q += PS_point[singular_structure.legs[1].n]
        Q_square = Q.square()
        pa_dot_pb = Q_square/2.

        # Compute the alpha parameter
        alpha = 0.5*( pC.dot(Q)/pa_dot_pb - math.sqrt( (pC.dot(Q)/pa_dot_pb)**2 - 4.* (pC.square()/Q_square) ) )

        # Recoil against initial states
        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n] = PS_point[recoiler.n]*(1.-alpha)

        # Find the leg mother label
        # Note: in the case of nested collinear limits, the fetching below must be upgraded
        mother_label = momenta_dict.inv[frozenset(collinear_children)]
        new_PS_point[mother_label] = pC - alpha*Q

        mapping_variables = {'Q': Q, 'alpha': alpha}

        # The jacobian may be one, but not sure; haven't computed it yet.
        if compute_jacobian:
            raise NotImplementedError

        return new_PS_point, mapping_variables

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables, compute_jacobian=False ):
        """ Inverse function of map_to_lower_multiplicity using the kinematic_variables supplied as generating variables."""

        raise NotImplementedError

# Soft mapping recoiling against initial states, suited for p p > X collisions.
#=========================================================================================
class SoftVsInitialMapping(ElementaryMappingSoft):
    """Implementation of the mapping described in sect. 5.3.2. of Vittorio's notes.
    It is applicable for hadronic collisions where the initial states can be rescaled.
    - is_valid_structure is entirely inherited from the parent.
    This mapping works for any hadron-hadron collider.
    The initial state should be checked at initialization
    - map_to_lower_multiplicity does TODO write description
    - map_to_higher_multiplicity does TODO write description
    """
    #TODO implement a check for Had-Had initial state when calling the schemes that call SoftVsInitialMapping
    #TODO validate map_to_lower | checked with unit_test, need to check test_IR_limits
    #TODO validate map_to_higher | checked with unit_test, need to check test_IR_limits
    #TODO once validated remove NotImplementedError from the methods

    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False ):
        """Map the resolved momenta of a real emission process to on-shell kinematics by recoling against the initial state partons.
        This function applies the transformation described in equation 5.44 of Vittorio's notes. 
        This is also described in Simone Lionetti's thesis at 5.3.3 but the notation here is aligned with Vittorio's notes.

        :param PS_point: dictionary of Real kinematics
        :param singular_structure: SingularStructure describing the limit we want to factorize with this mapping
        :param momenta_dict: ?
        :param squared_masses: dictionary
        :param kinematic_variables: dictionary to store the variables of the unresolved PS
        :param compute_jacobian: Boolean, should the jacobian be computed for each PS point?
        :return:
            * new_PS_point: a dictionary with Born kinematics
            * mapping_variables: a dictionary with variables needed to compute splitting functions
        """
        #TODO review documentation

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Construct a new phase space point
        new_PS_point = PS_point.get_copy()

        # Build the total soft momentum and remove all soft legs from the mapped PS point
        # Save the soft legs in kinematic_variables
        pS = LorentzVector()
        for substructure in singular_structure.substructures: # S(3,4).legs returns [3,4]
            children = tuple(leg.n for leg in substructure.legs)
            if kinematic_variables is not None:
                SoftVariables.get(PS_point, children, kinematic_variables)
            for child in children:
                pS += new_PS_point.pop(child)

        # Build the total momentum Q from the intial state
        Q = LorentzVector()
        Q += PS_point[1]
        Q += PS_point[2]

        # Compute the parameter la (lambda_s1...sn in Vittorio's 5.44)
        la = math.sqrt((Q - pS).square() / Q.square())

        # All recoilers are boosted by a Lorentz transform K->Ktilde
        Ktilde = Q * la
        K = Q - pS

        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n].rotoboost(K, Ktilde)  # this is Lambda(Ktilde,K) p_n
            # new_PS_point[recoiler.n].boost_from_to(K, Ktilde)  # this is Lambda(Ktilde,K) p_n
        # The initial state is rescaled by lambda
        new_PS_point[1]*=la
        new_PS_point[2]*=la

        mapping_variables = {'Q': Q}

        # There is no jacobian for this transformation as it is based on a pure Lorentz transform
        if compute_jacobian:
            mapping_variables['jacobian'] = 1
        return new_PS_point, mapping_variables

    @classmethod
    def map_to_higher_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):
        """Create a Real PS point from a Born PS point and unresolved PS variables by recoling against the initial state partons.
        This function applies the inverse transformation of that described in equation 5.44 of Vittorio's notes. This is also described in Simone Lionetti's thesis at 5.3.3 but the notation here is aligned with Vittorio's notes.

        :param PS_point: dictionary of Born kinematics + unresolved
        :param singular_structure: SingularStructure describing the limit we want to factorize with this mapping
        :param momenta_dict: ?
        :param kinematic_variables: dictionary to store the variables of the unresolved PS
        :param compute_jacobian: Boolean, should the jacobian be computed for each PS point?
        :return:
                * new_PS_point: a dictionary with Real kinematics
                * mapping_variables: a dictionary with variables needed to compute splitting functions
        """

        # Consistency checks
        # Type check
        assert isinstance(momenta_dict, sub.bidict)
            # Are we using the right kind of mapping for the limit?
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))
        # Do we have all the input we need?
        needed_variables = set(
            cls.get_kinematic_variables_names(singular_structure, momenta_dict ) )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Construct a new phase space point
        new_PS_point = PS_point.get_copy()

        # Build the total soft momentum,
        pS = LorentzVector()
        for substructure in singular_structure.substructures: # singular_structure = (S(...),): get the SoftStructure
            children = tuple(leg.n for leg in substructure.legs) # Get all unresolved legs: S(1,2).legs = (1,2)
            # the kinematic variables describing the unresolved legs are copied from kinematic_variables
            SoftVariables.set(new_PS_point, children, kinematic_variables)
            for child in children:
                pS += new_PS_point[child]

        # Build the total momentum
        # The sum of the momenta of initial state particles yields lambda*Q as they have been mapped
        laQ = LorentzVector()
        laQ += PS_point[1]
        laQ += PS_point[2]

        # Compute lambda (la)
        # In the direct mapping: lambda = sqrt(1 - 2*pS.Q/Q^2 + pS^2/Q^2) = sqrt(1 - y + muS)
        # Here the total momentum is laQ with Q = laQ / lambda.
        # Defining y' = pS.laQ / laQ^2 = y/lambda and muS' = pS^2/Qla^2 = muS / lambda^2
        # => lambda = sqrt(1 - lambda*y' + lambda^2*muS')
        # This has a unique positive solution
        yprime = 2*laQ.dot(pS)/laQ.square()
        muSprime = pS.square()/laQ.square()
        la = .5*(math.sqrt(4.*(1.-muSprime)+yprime**2)-yprime)/(1.-muSprime) # Sanity analytic check: limit of la in pS->0 is 1

        # Reconstruct the momenta of the inverse Lorentz transform
        Q = laQ.get_copy()/la
        K = Q-pS
        Ktilde = laQ

        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n].rotoboost(Ktilde,K)
            # new_PS_point[recoiler.n].boost_from_to(Ktilde, K)
        # Map the initial state
        new_PS_point[1] *= 1. / la
        new_PS_point[2] *= 1. / la

        # Save the variables needed for the splitting functions
        mapping_variables = {'Q': Q}
        if compute_jacobian:
            mapping_variables['jacobian'] = 1.

        return new_PS_point, mapping_variables

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

    @classmethod
    def good_soft_recoiler(cls, recoiler_leg):

        raise NotImplemented

    @classmethod
    def coll_structure(cls, structure):

        all_collinear_structures = [
            sub.CollStructure(legs=substructure.legs)
            for substructure in structure.substructures ]
        return sub.SingularStructure(
            substructures=all_collinear_structures,
            legs=structure.legs )

    @classmethod
    def soft_structure(cls, structure):

        all_soft_structures = []
        all_soft_recoils    = copy.copy(structure.legs)
        for substructure in structure.substructures:
            all_soft_structures += substructure.substructures
            all_soft_recoils    += tuple(
                leg for leg in substructure.legs
                if cls.good_soft_recoiler(leg) )
        return sub.SingularStructure(
            substructures=all_soft_structures,
            legs=all_soft_recoils )

    @classmethod
    def coll_momenta_dict(cls, structure, momenta_dict):

        coll_momenta_dict = sub.bidict()
        for substructure in structure.substructures:
            coll_ns = frozenset(leg.n for leg in substructure.legs)
            parent, _, _ = get_structure_numbers(substructure, momenta_dict)
            coll_momenta_dict[parent] = coll_ns
        return coll_momenta_dict

    def is_valid_structure(
        self, singular_structure,
        i_soft_structure=None, i_coll_structure=None ):

        assert isinstance(singular_structure, sub.SingularStructure)
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
        assert isinstance(momenta_dict, sub.bidict)
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
        self, PS_point, singular_structure, momenta_dict, squared_masses=None,
        kinematic_variables=None, compute_jacobian=False, masses=None ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        coll_structure = self.coll_structure(singular_structure)
        soft_structure = self.soft_structure(singular_structure)
        assert self.is_valid_structure(
            singular_structure,
            i_soft_structure=soft_structure, i_coll_structure=coll_structure )

        soft_PS_point, soft_vars = self.soft_mapping.map_to_lower_multiplicity(
            PS_point, soft_structure, momenta_dict, squared_masses,
            kinematic_variables, compute_jacobian )
        coll_momenta_dict = self.coll_momenta_dict(singular_structure, momenta_dict)
        coll_PS_point, coll_vars = self.coll_mapping.map_to_lower_multiplicity(
            soft_PS_point, coll_structure, coll_momenta_dict, squared_masses,
            kinematic_variables, compute_jacobian )
        # We adopt the definition of Q from the soft mapping
        coll_vars['Q'] = soft_vars['Q']
        if compute_jacobian:
            coll_vars['jacobian'] *= soft_vars['jacobian']
        return coll_PS_point, coll_vars

    def map_to_higher_multiplicity(
        self, PS_point, singular_structure, momenta_dict, kinematic_variables,
        compute_jacobian=False ):

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        coll_structure = self.coll_structure(singular_structure)
        soft_structure = self.soft_structure(singular_structure)
        assert self.is_valid_structure(
            singular_structure,
            i_soft_structure=soft_structure, i_coll_structure=coll_structure )

        coll_momenta_dict = self.coll_momenta_dict(singular_structure, momenta_dict)
        coll_PS_point, coll_vars = self.coll_mapping.map_to_higher_multiplicity(
            PS_point, coll_structure, coll_momenta_dict,
            kinematic_variables, compute_jacobian )
        soft_PS_point, soft_vars = self.soft_mapping.map_to_higher_multiplicity(
            coll_PS_point, soft_structure, momenta_dict,
            kinematic_variables, compute_jacobian )
        # We adopt the definition of Q from the soft mapping
        coll_vars['Q'] = soft_vars['Q']
        if compute_jacobian:
            coll_vars['jacobian'] *= soft_vars['jacobian']
        return soft_PS_point, coll_vars

class SoftCollinearVsFinalMapping(SoftCollinearMapping):
    """Soft-collinear mapping that selects only final-state recoilers
    for the soft mapping.
    """

    @classmethod
    def good_soft_recoiler(cls, recoiler_leg):

        return recoiler_leg.state == recoiler_leg.FINAL

class RelabellingMapping(VirtualMapping):
    """ Mapping that does nothing but relabel momenta."""

    @classmethod
    def is_valid_structure(cls, singular_structure):
        """Check if the mapping can be applied to a given singular structure."""

        raise True

    @classmethod
    def needs_relabelling(cls, momenta_dict):
        """ Verfies if specified momenta_dict contains relabelling instructions."""
        return any(len(v)==1 for v in momenta_dict.values())

    @classmethod
    def map_to_lower_multiplicity(
            cls, PS_point, singular_structure, momenta_dict, **opts):
        """ Apply relabelling specified from momenta_dict with values containing a single entry. """

        out_PS_point = PS_point.get_copy()
        for mapped_to, mapped_from in momenta_dict.items():
            if len(mapped_from)==1:
                mapped_from = list(mapped_from)[0]
                out_PS_point[mapped_to] = out_PS_point[mapped_from]
                del out_PS_point[mapped_from]
        return out_PS_point

    @classmethod
    def map_to_higher_multiplicity(
            cls, PS_point, singular_structure, momenta_dict, **opts):
        """ Apply relabelling specified from momenta_dict with values containing a single entry. """

        out_PS_point = PS_point.get_copy()
        for mapped_to, mapped_from in momenta_dict.items():
            if len(mapped_from)==1:
                mapped_from = list(mapped_from)[0]
                out_PS_point[mapped_from] = out_PS_point[mapped_to]
                del out_PS_point[mapped_to]
        return out_PS_point