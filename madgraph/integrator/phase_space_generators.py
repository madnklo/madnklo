################################################################################
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
################################################################################

import sys
import os

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))

import logging
import math
import numpy
import copy

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

class PhaseSpaceGeneratorError(Exception):
    """Exception raised if an exception is triggered in integrators.""" 

#===============================================================================
# Vectors
#===============================================================================

class Vector(numpy.ndarray):

    def __new__(cls, *args, **opts):

        if args and isinstance(args[0], Vector):
            foo = args[0].get_copy()
        else:
            foo = numpy.asanyarray(*args, **opts).view(cls)
        return foo

    def eps(self):

        if numpy.issubdtype(self.dtype, numpy.inexact):
            return numpy.finfo(self.dtype).eps
        else:
            return 0

    def huge(self):

        if numpy.issubdtype(self.dtype, numpy.inexact):
            return numpy.finfo(self.dtype).max
        elif numpy.issubdtype(self.dtype, numpy.integer):
            return numpy.iinfo(self.dtype).max
        else:
            raise ValueError

    def __eq__(self, other):

        eps = max(self.eps(), other.eps())
        # numpy.allclose uses abs(self-other) which would compute the norm
        # for LorentzVector, thus view self and other as numpy.ndarray's
        return numpy.allclose(
            self.view(type=numpy.ndarray), other.view(type=numpy.ndarray),
            math.sqrt(eps), 0.
        )

    def __ne__(self, other):

        return not self.__eq__(other)

    def get_copy(self):

        # The vector instantiated by get_copy() should be modified
        # without changing the previous instance, irrespectively of the
        # (presumably few) layers that compose entries of the vector
        return copy.deepcopy(self)

    def dot(self, v):

        assert len(self) == len(v)
        return sum(self[i] * v[i] for i in range(len(self)))

    def square(self):

        return self.dot(self)

    def __abs__(self):

        return math.sqrt(self.square())

    def normalize(self):

        self /= abs(self)
        return

    def project_onto(self, v):

        return (self.dot(v) / v.square()) * v

    def component_orthogonal_to(self, v):

        return self - self.project_onto(v)

    # Specific to 3D vectors
    def cross(self, v):

        assert len(self) == 3
        assert len(v) == 3
        return Vector([
            self[1] * v[2] - self[2] * v[1],
            self[2] * v[0] - self[0] * v[2],
            self[0] * v[1] - self[1] * v[0]
        ])

class LorentzVector(Vector):

    def __new__(cls, *args, **opts):

        if len(args)==0:
            return super(LorentzVector, cls).__new__(cls, [0.,0.,0.,0.], **opts)
        return super(LorentzVector, cls).__new__(cls, *args, **opts)

    def __str__(self):
        """Nice representation of a Lorentz vector."""

        return "(E = %s, p_x = %s, p_y = %s, p_z = %s, M^2 = %s)"%(
            self[0],self[1],self[2],self[3],self.square()
        )

    def space(self):
        """Return the spatial part of this LorentzVector."""

        return self[1:].view(type=Vector)

    def dot(self, v):
        """Compute the Lorentz scalar product."""

        pos = self[0]*v[0]
        neg = self.space().dot(v.space())
        if pos+neg != 0 and abs(2*(pos-neg)/(pos+neg)) < 100.*self.eps():
            return 0
        return pos - neg

    def rho2(self):
        """Compute the radius squared."""

        return self.space().square()

    def rho(self):
        """Compute the radius."""

        return abs(self.space())

    def set_square(self, square, negative=False):
        """Change the time component of this LorentzVector
        in such a way that self.square() = square.
        If negative is True, set the time component to be negative,
        else assume it is positive.
        """

        # Note: square = self[0]**2 - self.rho2(),
        # so if (self.rho2() + square) is negative, self[0] is imaginary.
        # Letting math.sqrt fail if data is not complex on purpose in this case.
        self[0] = math.sqrt(self.rho2() + square)
        if negative:
            self[0] *= -1
        return self

    def rotoboost(self, p, q):
        """Apply the Lorentz transformation that sends p in q to this vector."""

        # Check that the two invariants are close,
        # else the transformation is invalid
        p2 = p.square()
        q2 = q.square()
        assert abs(p2-q2) < 1e-9*(abs(p2)+abs(q2))
        # Compute scalar products
        pq = p + q
        pq2 = pq.square()
        p_s = self.dot(p)
        pq_s = self.dot(pq)
        # Assemble vector
        self += 2 * ((p_s/q2) * q - (pq_s/pq2) * pq)
        return self

    def pt(self, axis=3):
        """Compute transverse momentum."""

        return math.sqrt(
            sum(self[i]**2 for i in range(1, len(self)) if i != axis)
        )

    def pseudoRap(self):
        """Compute pseudorapidity."""

        pt = self.pt()
        if pt < self.eps() and abs(self[3]) < self.eps():
            return self.huge()*(self[3]/abs(self[3]))
        th = math.atan2(pt, self[3])
        return -math.log(math.tan(th/2.))

    def getdelphi(self,p2):
        """Compute the phi-angle separation with p2."""

        pt1 = self.pt()
        pt2 = p2.pt()
        if pt1 == 0. or pt2 == 0.:
            return self.huge()
        tmp = self[1]*p2[1] + self[2]*p2[2]
        tmp /= (pt1*pt2)
        if abs(tmp) > (1.0+self.eps()):
            raise PhaseSpaceGeneratorError(
                "Cosine larger than 1. in phase-space cuts."
            )
        if abs(tmp) > 1.0:
            return math.acos(tmp/abs(tmp))
        return math.acos(tmp)

    def deltaR(self, p2):
        """Compute the deltaR separation with momentum p2"""

        delta_eta = self.pseudoRap() - p2.pseudoRap()
        delta_phi = self.getdelphi(p2)

        return math.sqrt(delta_eta**2 + delta_phi**2)

    def boostVector(self):

        if self[0] == 0.:
            if self.rho2() == 0.:
                return (0.,0.,0.)
            else:
                raise PhaseSpaceGeneratorError(
                    "Attempting to compute a boost"
                    "from a reference vector with zero energy."
                )
        if abs(self) < 0.:
            raise PhaseSpaceGeneratorError(
                    "Attempting to compute a boost"
                    "from a reference vector with negative mass."
            )

        return tuple(_*(1./self[0]) for _ in self[1:4])

    def cosTheta(self):

        ptot = self.rho()
        assert (ptot > 0.)
        return self[3] / ptot

    def phi(self):

        return math.atan2(self[2], self[1])
    
    def boost(self, boost_vector, gamma=-1.):
        """Transport self into the rest frame of the boost_vector in argument.
        This means that the following command, for any vector p=(E, px, py, pz)
            p.boost(p.boostVector())
        transforms p to (M,0,0,0).
        """

        # TODO Rewrite this using the nice functions of Vector

        bx, by, bz = boost_vector[0], boost_vector[1], boost_vector[2]
        b2 = bx**2 + by**2 + bz**2
        if gamma < 0.:
            gamma = 1.0 / math.sqrt(1.0 - b2)

        bp = bx*self[1] + by*self[2] + bz*self[3]
        gamma2 = (gamma-1.0) / b2 if b2 > 0 else 0.
        self[1] = self[1] + gamma2*bp*bx + gamma*bx*self[0]
        self[2] = self[2] + gamma2*bp*by + gamma*by*self[0]
        self[3] = self[3] + gamma2*bp*bz + gamma*bz*self[0]
        self[0] = gamma*(self[0] + bp)

class LorentzVectorDict(dict):
    """A simple class wrapping dictionaries that store Lorentz vectors."""

    def __str__(self, n_initial=2):
        """Nice printout of the momenta."""

        # Use padding for minus signs
        def special_float_format(float):
            return '%s%.16e' % ('' if float < 0.0 else ' ', float)

        cols_widths = [4, 25, 25, 25, 25, 25]
        template = ' '.join(
            '%%-%ds' % col_width for col_width in cols_widths
        )
        line = '-' * (sum(cols_widths) + len(cols_widths) - 1)

        out_lines = [template % ('#', ' E', ' p_x', ' p_y', ' p_z', ' M',)]
        out_lines.append(line)
        running_sum = LorentzVector()
        for i in sorted(self.keys()):
            # BALDY why are you copying mom here? (pun came out for free)
            mom = LorentzVector(self[i])
            if i <= n_initial:
                running_sum += mom
            else:
                running_sum -= mom
            out_lines.append(template % tuple(
                ['%d' % i] + [
                    special_float_format(el) for el in (list(mom) + [abs(mom)])
                ]
            ))
        out_lines.append(line)
        out_lines.append(template % tuple(
            ['Sum'] + [special_float_format(el) for el in running_sum] + ['']
        ))

        return '\n'.join(out_lines)

    def get_copy(self):
        """Return a copy that can be freely modified
        without changing the current instance.
        """

        return type(self)({i: LorentzVector(self[i]) for i in self.keys()})

class LorentzVectorList(list):
    """A simple class wrapping lists that store Lorentz vectors."""

    def __str__(self, n_initial=2):
        """Nice printout of the momenta."""

        return LorentzVectorDict(
            (i + 1, v) for i, v in enumerate(self)
        ).__str__(n_initial=n_initial)

    def get_copy(self):
        """Return a copy that can be freely modified
        without changing the current instance.
        """

        return type(self)([LorentzVector(p) for p in self])

#===============================================================================
# Kinematic functions
#===============================================================================

def Kaellen(*args):

    l = len(args)
    foo = 0.
    for i in range(l):
        foo += args[i]**2
        for j in range(i+1, l):
            foo -= 2*args[i]*args[j]
    return foo

#===============================================================================
# Phase space generation
#===============================================================================

class VirtualPhaseSpaceGenerator(object):

    def __init__(self, initial_masses, final_masses,
                 beam_Es,
                 beam_types=(1,1)
                ):
        
        self.initial_masses  = initial_masses
        self.masses          = final_masses
        self.n_initial       = len(initial_masses)
        self.n_final         = len(final_masses)
        self.beam_Es         = beam_Es
        self.collider_energy = sum(beam_Es)
        self.beam_types      = beam_types
        self.dimensions      = self.get_dimensions()
        self.dim_ordered_names = [d.name for d in self.dimensions]
        # BALDY shouldn't you use bidict here, given we need it anyway?
        self.dim_name_to_position = dict((d.name,i) for i, d in enumerate(self.dimensions))
        self.position_to_dim_name = dict((v,k) for (k,v) in self.dim_name_to_position.items())
        
    def generateKinematics(self, E_cm, random_variables):
        """ Generates a phase-space point with fixed center of mass energy."""
        raise NotImplementedError
    
    def get_PS_point(self, random_variables):
        """ Generates a complete PS point, including Bjorken x's, dictating a specific choice
        of incoming particle's momenta,"""
        raise NotImplementedError

    def boost_to_COM_frame(self, PS_point, xb_1, xb_2):
        """ Boost a phase-space point to the rest-frame, given Bjorken x's."""
        
        if self.n_initial != 2 and (xb_1!=1. or xb_2!=1.):
            ref_lab = PS_point[0]*xb_1 + PS_point[1]*xb_2
            if ref_lab.rho2() != 0.:
                for p in PS_point:
                    p.boost(ref_lab.boostVector())

    def nDimPhaseSpace(self):
        """ Return the number of random numbers required to produce a given
        multiplicity final state. """

        if self.n_final == 1:
            return 0
        return 3*self.n_final -4

    def get_dimensions(self):
        """ Generates a list of dimensions for this integrand."""
        
        dims = integrands.DimensionList()

        # Add the PDF dimensions if necessary
        if self.beam_types[0]==self.beam_types[1]==1:
            dims.append(integrands.ContinuousDimension('ycms',lower_bound=0.0, upper_bound=1.0))
            # The 2>1 topology requires a special treatment
            if not (self.n_initial==2 and self.n_final==1):
                dims.append(integrands.ContinuousDimension('tau',lower_bound=0.0, upper_bound=1.0)) 

        # Add the phase-space dimensions
        dims.extend([ integrands.ContinuousDimension('x_%d'%i,lower_bound=0.0, upper_bound=1.0) 
                                     for i in range(1, self.nDimPhaseSpace()+1) ])
        
        return dims

class MultiChannelPhasespace(VirtualPhaseSpaceGenerator):
    """Implementation of a phase-space generator that lines up integration variables
    with s- and t-channels specifying a paticular provided topology.
    """
    
    def __init__(self, *args, **opts):
        
        if 'model' not in opts:
            raise PhaseSpaceGeneratorError("A model must be specified with the option 'model' when"+
                " instantiating the class %s."%self.__class__.__name__)
        self.model = opts.pop('model')
        
        if 'topology' not in opts:
            raise PhaseSpaceGeneratorError("A list of s- and t-channels must be specified with the "+
                "option 'topology' when instantiating the class %s."%self.__class__.__name__)
        self.topology = opts.pop('topology')
        
        super(MultiChannelPhasespace, self).__init__(*args, **opts)
        
        # One can do additional business here upon instantiating this PS generator, like renaming
        # the random variables to names describing the "propagators generated" with them. 

    def get_PS_point(self, random_variables, **opts):
        """ Generates a complete PS point, including Bjorken x's, dictating a specific choice
        of incoming particle's momenta,"""

        #
        # IMPLEMENTATION TO DO
        # 
        # For now just return a random PS point from flat generation
        return FlatInvertiblePhasespace(self.initial_masses,self.masses,self.beam_Es,
                                        beam_types  = self.beam_types).get_PS_point(random_variables)
        raise NotImplementedError

class FlatInvertiblePhasespace(VirtualPhaseSpaceGenerator):
    """ Implementation following S. Platzer: arxiv:1308.2922 """

    # This parameter defines a thin layer around the boundary of the unit hypercube of the 
    # random variables generating the phase-space, so as to avoid extrema which are an issue in most
    # PS generators.
    epsilon_border = 1e-10

    # The lowest value that the center of mass energy can take.
    # We take here 1 GeV, as anyway below this non-perturbative effects dominate and factorization does not
    # make sense anymore
    absolute_Ecm_min = 1.

    # For reference here we put the flat weights that Simon uses in his
    # Herwig implementation. I will remove them once I will have understood
    # why they don't match the physical PS volume.
    # So these are not used for now, and get_flatWeights() is used instead.
    flatWeights =  { 2 :  0.039788735772973833942,
                     3 :  0.00012598255637968550463,
                     4 :  1.3296564302788840628e-7,
                     5 :  7.0167897579949011130e-11,
                     6 :  2.2217170114046130768e-14 
                   }

    def __init__(self, *args, **opts):
        super(FlatInvertiblePhasespace, self).__init__(*args, **opts)
        
        if self.n_initial == 1:
            raise InvalidCmd("This basic generator does not support decay topologies.")

    def get_dimensions(self):
        """ Make sure the collider setup is supported."""

        # Add the PDF dimensions if necessary
        if (not abs(self.beam_types[0])==abs(self.beam_types[1])==1) and \
           (not self.beam_types[0]==self.beam_types[1]==0):
            dims.append(integrands.ContinuousDimension('ycms',lower_bound=0.0, upper_bound=1.0))
            raise InvalidCmd(
                "This basic generator does not support the collider configuration: (lpp1=%d, lpp2=%d)"%
                             (self.run_card['lpp1'], self.run_card['lpp2']))
        
        if self.beam_Es[0]!=self.beam_Es[1]:
            raise InvalidCmd(
                "This basic generator only supports colliders with incoming beams equally energetic.")

        return super(FlatInvertiblePhasespace,self).get_dimensions()

    @classmethod
    def get_flatWeights(cls, E_cm, n, mass=None):
        """ Return the phase-space volume for a n massless final states.
          Vol(E_cm, n) = (pi/2)^(n-2) *  (E_cm^2)^(n-2) / ((n-1)!*(n-2)!)
        """
        if n==1: 
            # The jacobian from \delta(s_hat - m_final**2) present in 2->1 convolution
            # must typically be accounted for in the MC integration framework since we
            # don't have access to shat here, so we just return 1.
            return 1.

        return math.pow((math.pi/2.0),n-1)*\
            (math.pow((E_cm**2),n-2)/(math.factorial(n-1)*math.factorial(n-2)))

    @classmethod
    def bisect(cls, v, n, target = -16., maxLevel = 80):
        """ Solve v = (n+2) * u^(n+1) - (n+1) * u^(n+2) for u. """
        
        if (v == 0.0 or v == 1.0):
            return v

        level = 0
        left  = 0.
        right = 1.
            
        checkV = -1.
        u = -1.
        while ( level < maxLevel ):
            u = (left + right)*math.pow(0.5,level+1.)
            checkV = math.pow(u,n+1)*(n+2.-(n+1.)*u)
            error = abs(1.-checkV/v) 
            if ( error==0. or math.log10(error) <= target ):
                break

            left *= 2.
            right *= 2.
            if (v <= checkV ):
                right -= 1.
            else:
                left += 1.

            level += 1

        return u
    
    @staticmethod
    def rho(M, N, m):
        """ Returns sqrt((sqr(M)-sqr(N+m))*(sqr(M)-sqr(N-m)))/(8.*sqr(M)) """
        Msqr = math.pow(M,2)
        return math.sqrt( (Msqr-math.pow((N+m),2))*(Msqr-math.pow((N-m),2)) ) / (8.*Msqr)

    def setInitialStateMomenta(self, output_momenta, E_cm):
        """ Generates the initial state momenta."""

        if self.n_initial not in [1,2]:
            raise InvalidCmd(
               "This PS generator only supports 1 or 2 initial states")

        if self.n_initial == 1:
            if self.initial_masses[0]==0.:
                raise PhaseSpaceGeneratorError(
                    "Cannot generate the decay phase-space of a massless particle.")
            if self.E_cm != self.initial_masses[0]:
                raise PhaseSpaceGeneratorError(
                    "Can only generate the decay phase-space of a particle at rest.")

        if self.n_initial == 2:
            if self.initial_masses[0] != 0. or self.initial_masses[1] != 0.:
                raise InvalidCmd(
                    "Phase-space of massive 2-particle collision not implemented yet (trivial though).")

        if self.n_initial == 1:
            output_momenta[0] = LorentzVector([self.initial_masses[0] , 0., 0., 0.])
            return

        output_momenta[0] = LorentzVector([E_cm/2.0 , 0., 0., E_cm/2.0])
        output_momenta[1] = LorentzVector([E_cm/2.0 , 0., 0., -E_cm/2.0])
        return

    def get_PS_point(self, random_variables):
        """ Generates a complete PS point, including Bjorken x's, dictating a specific choice
        of incoming particle's momenta,"""

        # if random_variables are not defined, than just throw a completely random point
        if random_variables is None:
            random_variables = self.dimensions.random_sample()
        
        # Phase-space point weight to return
        wgt = 1.0
        
        # Avoid extrema since the phase-space generation algorithm doesn't like it
        random_variables = [min(max(rv,self.epsilon_border),1.-self.epsilon_border) for rv in random_variables]

        # Assign variables to their meaning.
        if 'ycms' in self.dim_name_to_position:
            PDF_ycm = random_variables[self.dim_name_to_position['ycms']]
        else:
            PDF_ycm = None
        if 'tau' in self.dim_name_to_position:
            PDF_tau = random_variables[self.dim_name_to_position['tau']]
        else:
            PDF_tau = None
        PS_random_variables  = [rv for i, rv in enumerate(random_variables) if self.position_to_dim_name[i].startswith('x') ]

        # Now take care of the Phase-space generation:
        # Set some defaults for the variables to be set further
        xb_1 = 1.
        xb_2 = 1.
        E_cm = self.collider_energy
        
        # We generate the PDF from two variables \tau = x1*x2 and ycm = 1/2 * log(x1/x2), so that:
        #  x_1 = sqrt(tau) * exp(ycm)
        #  x_2 = sqrt(tau) * exp(-ycm)
        # The jacobian of this transformation is 1.
        
        if abs(self.beam_types[0])==abs(self.beam_types[1])==1:
            
            tot_final_state_masses = sum(self.masses)
            if tot_final_state_masses > self.collider_energy:
                raise PhaseSpaceGeneratorError("Collider energy is not large enough, there is no phase-space left.")
            
            # Keep a hard cut at 1 GeV, which is the default for absolute_Ecm_min
            tau_min = (max(tot_final_state_masses, self.absolute_Ecm_min)/self.collider_energy)**2
            tau_max = 1.0

            if self.n_initial == 2 and self.n_final == 1:
                # Here tau is fixed by the \delta(xb_1*xb_2*s - m_h**2) which sets tau to 
                PDF_tau = tau_min
                # Account for the \delta(xb_1*xb_2*s - m_h**2) and corresponding y_cm matching to unit volume
                wgt *= (1./self.collider_energy**2)
            else:
                # Rescale tau appropriately
                PDF_tau = tau_min+(tau_max-tau_min)*PDF_tau
                # Including the corresponding Jacobian
                wgt *= (tau_max-tau_min)

            # And we can now rescale ycm appropriately
            ycm_min = 0.5 * math.log(PDF_tau)
            ycm_max = -ycm_min
            PDF_ycm = ycm_min + (ycm_max - ycm_min)*PDF_ycm            
            # and account for the corresponding Jacobina
            wgt *= (ycm_max - ycm_min)

            xb_1 = math.sqrt(PDF_tau) * math.exp(PDF_ycm)
            xb_2 = math.sqrt(PDF_tau) * math.exp(-PDF_ycm)
            E_cm = math.sqrt(xb_1*xb_2)*self.collider_energy

        elif self.beam_types[0]==self.beam_types[1]==0:
            xb_1 = 1.
            xb_2 = 1.
            E_cm = self.collider_energy
        else:
            raise InvalidCmd("This basic PS generator does not yet support collider mode (%d,%d)."%self.beam_types)

        # Make sure the Bjorken x's are physical:
        if xb_1>1. or xb_2>1.:
            return None, 0.0, xb_1, xb_2

        # Now generate a PS point
        PS_point, PS_weight = self.generateKinematics(E_cm, PS_random_variables)
        
        # Apply the phase-space weight
        wgt *= PS_weight
        
        return LorentzVectorList(PS_point), wgt, xb_1, xb_2

    def generateKinematics(self, E_cm, random_variables):
        """ Generate a self.n_initial -> self.n_final phase-space point using 
        the random variables passed in argument."""

        # Make sure the right number of random variables are passed
        assert (len(random_variables)==self.nDimPhaseSpace())

        # Make sure that none of the random_variables is NaN.
        if any(math.isnan(rv) for rv in random_variables):
            raise PhaseSpaceGeneratorError("Some of the random variables passed "+
              "to the phase-space generator are NaN: %s"%str(random_variables))

        # The distribution weight of the generate PS point
        weight = 1.
        
        output_momenta = []

        if self.n_final == 1:
            if self.n_initial == 1:
                raise InvalidCmd("1 > 1 phase-space generation not supported.")
            if self.masses[0]/E_cm < 1.e-7 or ((E_cm-self.masses[0])/self.masses[0]) > 1.e-7:
                raise PhaseSpaceGeneratorError("1 > 2 phase-space generation needs a final state mass equal to E_c.o.m.")
            output_momenta.append(LorentzVector([self.masses[0]/2.,0.,0.,self.masses[0]/2.]))
            output_momenta.append(LorentzVector([self.masses[0]/2.,0.,0.,-self.masses[0]/2.]))
            output_momenta.append(LorentzVector([self.masses[0]   ,0.,0.,0.]))
            weight = self.get_flatWeights(E_cm,1)
            return output_momenta, weight
  
        M    = [ 0. ]*(self.n_final-1)
        M[0] = E_cm

        weight *= self.generateIntermediatesMassive(M, E_cm, random_variables)
        M.append(self.masses[-1])

        Q     = LorentzVector([M[0],0.,0.,0.])
        nextQ = LorentzVector(4*[0., ])

        for i in range(self.n_initial+self.n_final-1):
            
            if i < self.n_initial:
                output_momenta.append(LorentzVector(4*[0., ]))
                continue

            q = 4.*M[i-self.n_initial]*\
                 self.rho(M[i-self.n_initial],M[i-self.n_initial+1],self.masses[i-self.n_initial])
            cos_theta = 2.*random_variables[self.n_final-2+2*(i-self.n_initial)]-1.;
            sin_theta = math.sqrt(1.-math.pow(cos_theta,2))
            phi = 2.*math.pi*random_variables[self.n_final-1+2*(i-self.n_initial)]
            cos_phi = math.cos(phi)
            sin_phi = math.sqrt(1.-math.pow(cos_phi,2))

            if (phi > math.pi):
                sin_phi = -sin_phi; 
            
            p = LorentzVector([0.,q*cos_phi*sin_theta, 
                                 q*sin_phi*sin_theta,
                                 q*cos_theta])
            p.set_square(self.masses[i-self.n_initial]**2)
            p.boost(Q.boostVector())
            p.set_square(self.masses[i-self.n_initial]**2)
            output_momenta.append(p)

            nextQ = Q - p
            nextQ.set_square(M[i-self.n_initial+1]**2)

            Q = nextQ
       
        output_momenta.append(Q)

        self.setInitialStateMomenta(output_momenta, E_cm)

        return LorentzVectorList(output_momenta), weight

    def generateIntermediatesMassless(self, M, E_cm, random_variables):
        """ Generate intermediate masses for a massless final state."""
        
        for i in range(2,self.n_final):
            u = self.bisect(random_variables[i-2],self.n_final-1-i)
            M[i-1] = math.sqrt(u*math.pow(M[i-2],2))

        return self.get_flatWeights(E_cm,self.n_final)
   

    def generateIntermediatesMassive(self, M, E_cm, random_variables):
        """ Generate intermediate masses for a massive final state."""

        K = list(M)
        K[0] -= sum(self.masses)

        weight = self.generateIntermediatesMassless(K, E_cm, random_variables)
        del M[:]
        M.extend(K)
        
        for i in range(1,self.n_final):
            for k in range(i,self.n_final+1):
                M[i-1] += self.masses[k-1]
        
        weight *= 8.*self.rho(M[self.n_final-2],self.masses[self.n_final-1],self.masses[self.n_final-2])

        for i in range(2,self.n_final):
            weight *= (self.rho(M[i-2],M[i-1],self.masses[i-2])/self.rho(K[i-2],K[i-1],0.)) * (M[i-1]/K[i-1])

        weight *= math.pow(K[0]/M[0],2*self.n_final-4)

        return weight

    def invertKinematics(self, E_cm, momenta):
        """ Returns the random variables that yields the specified momenta configuration."""

        # Make sure the right number of momenta are passed
        assert (len(momenta) == (self.n_initial + self.n_final) )
        moms = momenta.get_copy()

        # The weight of the corresponding PS point
        weight = 1.

        if self.n_final == 1:
            if self.n_initial == 1:
                raise PhaseSpaceGeneratorError("1 > 1 phase-space generation not supported."%str(random_variables))
            return [], self.get_flatWeights(E_cm,1) 

        # The random variables that would yield this PS point.
        random_variables = [-1.0]*self.nDimPhaseSpace()
        
        M    = [ 0. ]*(self.n_final-1)
        M[0] = E_cm

        Q     = [LorentzVector(4*[0., ])]*(self.n_final-1)
        Q[0]  = LorentzVector([M[0],0.,0.,0.])

        for i in range(2,self.n_final):
            for k in range(i, self.n_final+1):
                Q[i-1] = Q[i-1] + moms[k+self.n_initial-1]
            M[i-1] = abs(Q[i-1])

        weight = self.invertIntermediatesMassive(M, E_cm, random_variables)

        for i in range(self.n_initial,self.n_final+1):
            # BALDY another copy? moms not used afterwards
            p = LorentzVector(moms[i])
            # Take the opposite boost vector
            boost_vec = tuple([-_ for _ in Q[i-self.n_initial].boostVector()])
            p.boost(boost_vec)
            random_variables[self.n_final-2+2*(i-self.n_initial)] = (p.cosTheta()+1.)/2.
            phi = p.phi()
            if (phi < 0.):
                phi += 2.*math.pi
            random_variables[self.n_final-1+2*(i-self.n_initial)] = phi / (2.*math.pi)
        
        return random_variables, weight

    def invertIntermediatesMassive(self, M, E_cm, random_variables):
        """ Invert intermediate masses for a massive final state."""

        K = list(M)
        for i in range(1, self.n_final):
            K[i-1] -= sum(self.masses[i-1:])
        
        weight = self.invertIntermediatesMassless(K, E_cm, random_variables)
        weight *= 8.*self.rho(M[self.n_final-2],
                              self.masses[self.n_final-1],
                              self.masses[self.n_final-2])
        for i in range(2, self.n_final):
            weight *= (self.rho(M[i-2],M[i-1],self.masses[i-2])/self.rho(K[i-2],K[i-1],0.)) \
                                                                      * (M[i-1]/K[i-1])
        
        weight *= math.pow(K[0]/M[0],2*self.n_final-4)

        return weight

    def invertIntermediatesMassless(self, K, E_cm, random_variables):
        """ Invert intermediate masses for a massless final state."""

        for i in range(2, self.n_final):
            u = (K[i-1]/K[i-2])**2
            random_variables[i-2] = \
                (self.n_final+1-i)*math.pow(u,self.n_final-i) - \
                (self.n_final-i)*math.pow(u, self.n_final+1-i)
        
        return self.get_flatWeights(E_cm, self.n_final)

#===============================================================================
# Mappings
#===============================================================================

def get_structure_numbers(structure, momenta_dict):
    """Return the number of the parent and children of a given structure
    according to some momenta dictionary.
    """

    children = frozenset((leg.n for leg in structure.get_all_legs()))
    if structure.name() == "S":
        return None, children
    else:
        return momenta_dict.inv[children], children

class VirtualMapping(object):
    """Base class for elementary mapping implementations."""

    def __init__(self, model = None, **opts):
        """General initialization of any mapping.
        The model is an optional specification,
        which can be useful to know properties of the leg mapped.
        """        
        self.model = model
        
    def is_valid_structure(self, singular_structure):
        """Return true if it makes sense to apply this mapping
        to the given singular structure.
        """

        raise NotImplemented

    def get_kinematic_variables_names(self, singular_structure, momenta_dict):
        """For a valid singular structure,
        returns a list of variable names necessary to apply this mapping.
        """

        raise NotImplemented

    def map_to_lower_multiplicity(
        self, PS_point, singular_structure, momenta_dict,
        kinematic_variables = None
    ):
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

        :return: the jacobian weight due to the mapping
        """
        
        raise NotImplemented

    def map_to_higher_multiplicity(
        self, PS_point, singular_structure, momenta_dict, kinematic_variables
    ):
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

        :return: the jacobian weight due to the mapping
        """
        
        raise NotImplemented

    def approach_limit(
        self, PSpoint, singular_structure, kinematic_variables,
        scaling_parameter
    ):
        """Map to higher multiplicity, with the distance to the limit expressed
        by scaling_parameter.
        """

        # TODO This is a stub

        raise NotImplemented

#===============================================================================
# Final-final collinear mappings
#===============================================================================

class ElementaryMappingCollinearFinal(VirtualMapping):
    """Common functions for final-final collinear elementary mappings."""

    def __init__(self, *args, **opts):
        """Additional options for a final-final collinear elementary mapping."""

        super(ElementaryMappingCollinearFinal, self).__init__(*args, **opts)

    def is_valid_structure(self, singular_structure):

        assert isinstance(singular_structure, subtraction.SingularStructure)
        # Valid only for bunches of final-state particles going collinear,
        # with no recursive substructure
        for substructure in singular_structure.substructures:
            if not substructure.name() == "C":
                return False
            if substructure.substructures:
                return False
            if substructure.get_all_legs().has_initial_state_leg():
                return False
        return True

    def get_kinematic_variables_names(self, singular_structure, momenta_dict):
        """Get the names of variables describing particles going unresolved."""

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert self.is_valid_structure(singular_structure)

        names = []
        # For every collinear subset of N particles,
        # a set of (N-1) za's and zb's and 2-component normalised nt's is needed
        # plus the set's virtuality
        for substructure in singular_structure.substructures:
            # Add direct legs, last one not needed because of sum rules
            parent, children = get_structure_numbers(substructure, momenta_dict)
            for legn in sorted(children)[:-1]:
                names += ['za'+str(legn), 'zb'+str(legn), 'nt'+str(legn), ]
            # Add virtuality of this subset
            names += ['s'+str(parent), ]

        return names

    precision_loss_message = "Precision loss detected when computing collinear variables."

    def get_collinear_variables(
        self, PS_point, parent, children, kinematic_variables,
        precision=1e-6
    ):
        """Given unmapped and mapped momenta, compute the kinematic variables
        that describe the internal structure of particles going unresolved.
        Parent and children are indices that already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        # WARNING This function assumes massless particles

        # Retrieve the parent's momentum
        p = PS_point[parent]
        # Compute the sum of momenta
        q = LorentzVector(4*[0., ])
        for i in children:
            q += PS_point[i]
        # Pre-compute scalar products
        q2 = q.square()
        pq = p.dot(q)
        sqrtq2 = math.sqrt(q2)
        na = (sqrtq2/pq) * p
        nb = (2./sqrtq2) * q - na
        # Initialize variables for sum rules check
        zasum = 0
        ktsum = LorentzVector()
        ktabssum = LorentzVector()
        # Compute all kinematic variables
        for i in children:
            pi = PS_point[i]
            napi = na.dot(pi)
            nbpi = nb.dot(pi)
            zai = nbpi/nb.dot(q)
            zbi = napi/na.dot(q)
            kti = pi - 0.5*(nbpi*na+napi*nb)
            nti = kti / math.sqrt(napi*nbpi)
            kinematic_variables['za' + str(i)] = zai
            kinematic_variables['zb' + str(i)] = zbi
            kinematic_variables['nt' + str(i)] = nti
            zasum += zai
            ktsum += kti
            for j in range(len(kti)):
                ktabssum[j] += abs(kti[j])
        # Check numerical accuracy
        # TODO Ideally switch to quadruple precision if the check fails
        if abs(zasum - 1) > precision:
            misc.sprint(self.precision_loss_message)
            misc.sprint("Sum of z's is %.16e" % zasum)
        ktsum_abs = abs(ktsum.view(Vector))
        ktabssum_abs = abs(ktabssum.view(Vector))
        if ktsum_abs / ktabssum_abs > precision:
            misc.sprint(self.precision_loss_message)
            misc.sprint("Sum of kt's is "+str(ktsum))
        return

    def set_collinear_variables(
        self, PS_point, parent, children, total_momentum, kinematic_variables,
        precision=1e-6
    ):
        """Given the lower multiplicity parent momentum,
        the total momentum of children and collinear variables,
        compute and set the children momenta.
        Parent and children are indices that already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        Sum rules are checked to assess numerical accuracy.
        """

        # WARNING This function assumes massless particles

        # Retrieve the parent's momentum
        p = PS_point[parent]
        # Rename the sum of momenta
        q = total_momentum
        # Pre-compute scalar products
        q2 = q.square()
        pq = p.dot(q)
        sqrtq2 = math.sqrt(q2)
        na = (sqrtq2/pq) * p
        nb = (2./sqrtq2) * q - na
        # Variables for sums
        p_sum = LorentzVector()
        # Set momenta for all children
        for i in children:
            zai = kinematic_variables['za' + str(i)]
            zbi = kinematic_variables['zb' + str(i)]
            nti = kinematic_variables['nt' + str(i)]
            nbpi = zai*nb.dot(q)
            napi = zbi*na.dot(q)
            PS_point[i] = nti*math.sqrt(napi*nbpi) + 0.5*(napi*nb+nbpi*na)
            p_sum += PS_point[i]
        # Check how well the parent's momentum is reproduced
        # TODO Ideally switch to quadruple precision if the check fails
        deviation = abs((q - p_sum).view(Vector))
        benchmark = abs(q.view(Vector))
        if deviation / benchmark > precision:
            misc.sprint(self.precision_loss_message)
            misc.sprint(
                "Sum of children differs from parent momentum by "+str(q-p_sum)
            )
        return

    def azimuth_reference_vector(self, n):
        """Given a unit vector n,
        return a reference vector in the plane orthogonal to n.
        """

        n_ref = Vector([1., 0., 0.])
        return Vector(n_ref - n.dot(n_ref) * n).normalize()

    def get_random_variables(self, mapped_momentum, children, variables):
        """Get [0,1] random variables
        that represent the internal kinematic variables.
        """

        # WARNING This function is a stub

        # Empty dictionary for randoms
        randoms = dict()

        # Compute reference vectors for this collinear subset
        nvec = Vector([mapped_momentum[j + 1] for j in range(3)]).normalize()
        # Get a reference phi=0 direction
        n_phi = self.azimuth_reference_vector(nvec)
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

    def set_random_variables(self, mapped_momentum, children, randoms):
        """Compute collinear variables from randoms in the range [0,1].
        """

        # TODO implement this

        variables = dict()
        return variables

class MappingCataniSeymourFFOne(ElementaryMappingCollinearFinal):
    """Implementation of the Catani--Seymour mapping
    for one bunch of collinear particles
    generalized to an arbitrary number of spectators
    but restricted to massless particles.
    See arXiv:hep-ph/0609042 and references therein.
    """

    def __init__(self, *args, **opts):
        """Additional options for the Catani--Seymour final-final mapping."""

        super(MappingCataniSeymourFFOne, self).__init__(*args, **opts)

    def is_valid_structure(self, singular_structure):

        # Valid only for one bunch of final-state particles going collinear,
        # with no recursive substructure
        if len(singular_structure.substructures) == 1:
            return super(
                MappingCataniSeymourFFOne, self
            ).is_valid_structure(singular_structure)
        return False

    def map_to_lower_multiplicity(
        self, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None
    ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        print singular_structure
        assert self.is_valid_structure(singular_structure)

        # Precompute sets and numbers
        cluster = singular_structure.substructures[0]
        parent, children = get_structure_numbers(cluster, momenta_dict)
        # Build the collinear momentum
        pC = LorentzVector(4*[0., ])
        for j in children:
            pC += PS_point[j]
        # Build the total momentum of recoilers
        pR = LorentzVector(4*[0., ])
        for leg in singular_structure.legs:
            pR += PS_point[leg.n]
        Q = pC + pR
        Q2 = Q.square()
        # Compute the parameter alpha for this collinear subset
        QpC = pC.dot(Q)
        sC = pC.square()
        alpha = (QpC - math.sqrt(QpC**2 - Q2*sC)) / Q2
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] /= (1-alpha)
        # Map the cluster's momentum
        PS_point[parent] = (pC - alpha * Q) / (1-alpha)
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables is not None:
            kinematic_variables['s'+str(parent)] = sC
            self.get_collinear_variables(
                PS_point, parent, sorted(children), kinematic_variables
            )
        # Eliminate children momenta from the mapped phase-space point
        for j in children:
            if j != parent: # Bypass degenerate case of 1->1 splitting
                del PS_point[j]

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

    def map_to_higher_multiplicity(
        self, PS_point, singular_structure, momenta_dict, kinematic_variables
    ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert self.is_valid_structure(singular_structure)
        needed_variables = set(
            self.get_kinematic_variables_names(
                singular_structure, momenta_dict
            )
        )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Precompute sets and numbers...
        cluster = singular_structure.substructures[0]
        parent, children = get_structure_numbers(cluster, momenta_dict)
        # Find the parent's momentum
        qC = PS_point[parent]
        # Compute the momentum of recoilers and the total
        qR = LorentzVector(4*[0., ])
        for leg in singular_structure.legs:
            qR += PS_point[leg.n]
        Q = qR + qC
        # Compute scalar products
        QqC = Q.dot(qC)
        qC2 = qC.square()
        Q2 = Q.square()
        qR2 = qR.square()
        sC = kinematic_variables['s'+str(parent)]
        # Obtain parameter alpha
        alpha = 0
        if qR2 == 0:
            alpha = 0.5 * (sC - qC2) / (QqC - qC2)
        else:
            alpha = (math.sqrt(QqC**2 - Q2*qC2 + sC*qR2) - (QqC - qC2)) / qR2
        # Compute reverse-mapped momentum
        pC = (1-alpha) * qC + alpha * Q
        # Set children momenta
        self.set_collinear_variables(
            PS_point, parent, sorted(children), pC, kinematic_variables
        )
        # Map recoil momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler.n] *= 1-alpha
        # Remove parent's momentum
        if parent not in children: # Bypass degenerate case of 1->1 splitting
            del PS_point[parent]

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

class MappingCataniSeymourFFMany(ElementaryMappingCollinearFinal):
    """Implementation of the Catani--Seymour final-final collinear mapping,
    for multiple collinear bunches and any number of spectators
    (but restricted to massless particles).
    See arXiv:hep-ph/0609042 and references therein.
    """
    
    def __init__(self, *args, **opts):
        """Additional options for the Catani--Seymour mapping."""

        super(ElementaryMappingCollinearFinal, self).__init__(*args, **opts)

    def map_to_lower_multiplicity(
        self, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None
    ):

        # WARNING Needs checking

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert self.is_valid_structure(singular_structure)

        # Build the total momentum of collinear sets and recoilers
        Q = sum(
            PS_point[leg.n]
            for leg in singular_structure.get_all_legs()
        )
        Q2 = Q.square()
        # Prepare containers
        n_clusters = len(singular_structure.substructures)
        p_C = [LorentzVector(), ] * n_clusters
        alpha_C = [0., ] * n_clusters
        # For every cluster of collinear particles
        for i in range(n_clusters):
            cluster = singular_structure.substructures[i]
            # Build the collective momentum
            leg_ns = frozenset((leg.n for leg in cluster.legs))
            p_C[i] = sum(
                PS_point[momenta_dict.inv[(collinear_particle, )]]
                for collinear_particle in leg_ns
            )
            # Compute the parameter alpha for this collinear subset
            y_CQ = p_C[i].dot(Q) / Q2
            y_C = p_C[i].square() / Q2
            alpha_C[i] = y_CQ - math.sqrt(y_CQ**2 - y_C)
        # Compute the common scale factor
        scale_factor = 1. / (1. - sum(alpha_C))
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler] *= scale_factor
        # For every cluster of collinear particles
        for i in range(n_clusters):
            cluster = singular_structure.substructures[i]
            leg_ns = frozenset((leg.n for leg in cluster.legs))
            parent_number = momenta_dict.inv[leg_ns]
            # Map the cluster's momentum
            p_tilde_C = scale_factor * (p_C[i] - alpha_C[i] * Q)
            PS_point[parent_number] = p_tilde_C
            # If needed, update the kinematic_variables dictionary
            if kinematic_variables:
                kinematic_variables['s' + str(parent_number)] = p_C[i].square()
                children_numbers = [ momenta_dict.inv[(j,)] for j in leg_ns ]
                self.get_collinear_variables(
                    PS_point, parent_number, children_numbers, kinematic_variables
                )
        # TODO Compute the jacobian for this mapping
        jacobian = 1.0
#
        return jacobian

class Mapping_NagySoper(VirtualMapping):
    """Implementation of the Nagy--Soper mapping.
    See arXiv:hep-ph/0706.0017.
    """
    # TODO, see example above
    pass

#===============================================================================
# Soft mappings
#===============================================================================

class ElementaryMappingSoft(VirtualMapping):
    """Common functions for soft elementary mappings."""

    def __init__(self, *args, **opts):
        """Additional options for a soft elementary mapping."""

        super(ElementaryMappingSoft, self).__init__(*args, **opts)

    def is_valid_structure(self, singular_structure):

        assert isinstance(singular_structure, subtraction.SingularStructure)
        # Valid only for bunches of final-state particles going soft,
        # with no recursive substructure
        for substructure in singular_structure.substructures:
            if not substructure.name() == "S":
                return False
            if substructure.substructures:
                return False
        return True

    def get_kinematic_variables_names(self, singular_structure, momenta_dict):
        """Get the names of variables describing particles going unresolved."""

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert self.is_valid_structure(singular_structure)

        # For every soft particle, just return its momentum
        names = []
        for substructure in singular_structure.substructures:
            parent, children = get_structure_numbers(substructure, momenta_dict)
            for legn in sorted(children):
                names += ['p'+str(legn), ]
        return names

    def get_soft_variables(self, PS_point, children, kinematic_variables):
        """Given unmapped and mapped momenta, compute the kinematic variables
        that describe the internal structure of particles going unresolved.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        """

        # For soft particles, just pass the whole momentum
        for i in children:
            kinematic_variables['p' + str(i)] = PS_point[i]
        return


    def set_soft_variables(self, PS_point, children, kinematic_variables):
        """Given a dictionary of variables that describe the unresolved partons,
        compute and set the children momenta.
        Children indices should already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        """

        # Set momenta for all children
        for i in children:
            PS_point[i] = kinematic_variables['p' + str(i)]
        return


    def get_random_variables(self, mapped_momentum, children, variables):
        """Get random variables in the range [0,1]
        that correspond to the internal kinematic variables.
        """

        raise NotImplemented


    def set_random_variables(self, mapped_momentum, children, randoms):
        """Compute internal kinematic variables
        from random numbers in the range [0,1].
        """

        raise NotImplemented

class MappingSomogyietalSoft(ElementaryMappingSoft):
    """Implementation of the mapping used by Somogyi et al.
    in arXiv:hep-ph/0609042 for soft particles.
    It applies when there is at least one recoiler in the final state
    and all recoilers are massless.
    """

    def __init__(self, *args, **opts):
        """Additional options for the Somogyi et al. soft mapping."""

        super(MappingSomogyietalSoft, self).__init__(*args, **opts)

    def map_to_lower_multiplicity(
        self, PS_point, singular_structure, momenta_dict,
        kinematic_variables=None
    ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert self.is_valid_structure(singular_structure)

        # Build the total soft momentum,
        # save the soft momenta in variables and eliminate them from PS_point
        pS = LorentzVector(4 * [0., ])
        for substructure in singular_structure.substructures:
            children = [leg.n for leg in substructure.legs]
            if kinematic_variables is not None:
                self.get_soft_variables(PS_point, children, kinematic_variables)
            for child in children:
                pS += PS_point.pop(child)
        # Build the total momentum of recoilers
        pR = LorentzVector(4 * [0., ])
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

    def map_to_higher_multiplicity(
        self, PS_point, singular_structure, momenta_dict, kinematic_variables
    ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert self.is_valid_structure(singular_structure)
        needed_variables = set(
            self.get_kinematic_variables_names(
                singular_structure, momenta_dict
            )
        )
        assert needed_variables.issubset(kinematic_variables.keys())

        # Build the total soft momentum,
        # get the soft momenta from variables and save them in PS_point
        pS = LorentzVector(4 * [0., ])
        for substructure in singular_structure.substructures:
            children = [leg.n for leg in substructure.legs]
            self.set_soft_variables(PS_point, children, kinematic_variables)
            for child in children:
                pS += PS_point[child]
        # Build the total momentum, which is equal to the mapped recoilers'
        Q = LorentzVector(4 * [0., ])
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

#===============================================================================
# Mapping walkers
#===============================================================================

class VirtualWalker(object):
    """Base class for walker implementations."""
    
    def __new__(cls, **opts):
        """Factory class to make plugin easy."""

        if cls is VirtualWalker:
            map_type = opts.pop('map_type') if 'map_type' in opts else 'Unknown'
            if map_type not in mapping_walker_classes_map or not mapping_walker_classes_map[map_type]:
                raise MadGraph5Error(
                    "Unknown mapping walker of type '%s'."%str(map_type)
                )
            target_class = mapping_walker_classes_map[map_type]
            return super(VirtualWalker, cls).__new__(target_class, **opts)
        else:
            return super(VirtualWalker, cls).__new__(cls, **opts)

    def __init__(self, model=None, **opts):
        """General initialization of any mapping.
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
            counterterm, kinematic_variables, scaling_parameter
        )
        
        return self.walk_to_higher_multiplicity(
            PS_point, counterterm, new_kin_variables
        )

class FlatCollinearWalker(VirtualWalker):

    cannot_handle = """The Flat Collinear walker found a singular structure
    it is not capable to handle.
    """
    collinear_map = MappingCataniSeymourFFOne()

    def walk_to_lower_multiplicity(
        self, PS_point, counterterm, compute_kinematic_variables=False
    ):

        point = PS_point.get_copy()
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
        for node in counterterm.subcurrents:
            parent, _ = get_structure_numbers(
                node.current['singular_structure'],
                counterterm.momenta_dict
            )
            for recoiler in recoilers:
                if recoiler.n == parent:
                    recoilers.remove(recoiler)
        # Loop over the counterterm's first-level currents
        for node in counterterm.subcurrents:
            # Alias singular structure for convenience
            ss = node.current['singular_structure']
            # Do not accept soft currents or nested ones
            if ss.name == 'S' or node.subcurrents:
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
                subtraction.SingularStructure(ss, *recoilers),
                counterterm.momenta_dict,
                kinematic_variables
            )
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
        self, PS_point, counterterm, kinematic_variables
    ):

        # Identify reduced matrix element,
        # computed in the lowest multiplicity point
        ME_PS_pair = [counterterm.process, PS_point]
        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        # Initialize return variables
        current_PS_pairs = []
        jacobian = 1.
        # Recoil against all final-state particles that are not going singular
        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in counterterm.process['legs']
            if leg['state'] == base_objects.Leg.FINAL
        ]
        for node in counterterm.subcurrents:
            parent, _ = get_structure_numbers(
                node.current['singular_structure'],
                counterterm.momenta_dict
            )
            for recoiler in recoilers:
                if recoiler.n == parent:
                    recoilers.remove(recoiler)
        # Loop over the counterterm's first-level currents
        for node in reversed(counterterm.subcurrents):
            # Alias singular structure for convenience
            ss = node.current['singular_structure']
            # Do not accept soft currents or nested ones
            if ss.name == 'S' or node.subcurrents:
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
                subtraction.SingularStructure(ss, *recoilers),
                counterterm.momenta_dict,
                kinematic_variables
            )
            # Prepend pair of this current and the momenta,
            # deep copy wanted
            momenta.update(point.get_copy())
            current_PS_pairs.insert(0, (node.current, momenta))

        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'jacobian': jacobian,
            'resulting_PS_point': point
        }

    def rescale_kinematic_variables(
        self, counterterm, kinematic_variables, scaling_parameter
    ):

        # For collinear clusters, just rescale the virtuality
        new_kinematic_variables = {}
        for var, value in kinematic_variables.items():
            if var.startswith('s'):
                new_kinematic_variables[var] = value*scaling_parameter**2
            else:
                new_kinematic_variables[var] = value
        return new_kinematic_variables

class SimpleNLOWalker(VirtualWalker):

    cannot_handle = """The Simple NLO Walker found a singular structure
    it is not capable to handle.
    """
    collinear_map = MappingCataniSeymourFFOne()
    soft_map = MappingSomogyietalSoft()

    def walk_to_lower_multiplicity(
        self, PS_point, counterterm, compute_kinematic_variables=False
    ):

        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        # Initialize return variables
        current_PS_pairs = []
        jacobian = 1.
        kinematic_variables = dict() if compute_kinematic_variables else None
        # Recoil against all final-state particles that are not going unresolved
        # TODO Recoilers and numbers for a given counterterm should be cached
        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in counterterm.process['legs']
            if leg['state'] == base_objects.Leg.FINAL
        ]
        for node in counterterm.subcurrents:
            parent, children = get_structure_numbers(
                node.current['singular_structure'],
                counterterm.momenta_dict
            )
            if parent is None:
                for recoiler in recoilers:
                    if recoiler.n in children:
                        recoilers.remove(recoiler)
            else:
                for recoiler in recoilers:
                    if recoiler.n == parent:
                        recoilers.remove(recoiler)
                        break
        # Loop over the counterterm's first-level currents
        # Remember soft-collinear is a separate elementary current
        for node in counterterm.subcurrents:
            # Alias singular structure for convenience
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.subcurrents:
                raise MadGraph5Error(self.cannot_handle)
            # Get parent and children numbers
            # TODO Read these from cache
            parent, children = get_structure_numbers(
                ss, counterterm.momenta_dict
            )
            # Deep copy the momenta
            momenta = point.get_copy()
            # Compute jacobian and map to lower multiplicity
            if ss.name() == 'C':
                # Make a fake collinear structure for soft-collinears
                # in order to use the same mapping
                new_ss = ss
                if new_ss.substructures:
                    new_ss = subtraction.CollStructure(ss.get_all_legs())
                jacobian *= self.collinear_map.map_to_lower_multiplicity(
                    point,
                    subtraction.SingularStructure(new_ss, *recoilers),
                    counterterm.momenta_dict,
                    kinematic_variables
                )
                momenta[parent] = LorentzVector(point[parent])
            elif ss.name() == 'S':
                jacobian *= self.soft_map.map_to_lower_multiplicity(
                    point,
                    subtraction.SingularStructure(ss, *recoilers),
                    counterterm.momenta_dict,
                    kinematic_variables
                )
            else:
                raise MadGraph5Error(self.cannot_handle)
            # Append the current and the momenta
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
        self, PS_point, counterterm, kinematic_variables
    ):

        # Identify reduced matrix element,
        # computed in the lowest multiplicity point
        ME_PS_pair = [counterterm.process, PS_point]
        # This phase-space point will be destroyed, deep copy wanted
        point = PS_point.get_copy()
        # Initialize return variables
        current_PS_pairs = []
        jacobian = 1.
        # Recoil against all final-state particles that are not going unresolved
        # TODO Recoilers and numbers for a given counterterm should be cached
        recoilers = [
            subtraction.SubtractionLeg(leg)
            for leg in counterterm.process['legs']
            if leg['state'] == base_objects.Leg.FINAL
        ]
        for node in counterterm.subcurrents:
            parent, children = get_structure_numbers(
                node.current['singular_structure'],
                counterterm.momenta_dict
            )
            if parent is None:
                for recoiler in recoilers:
                    if recoiler.n in children:
                        recoilers.remove(recoiler)
            else:
                for recoiler in recoilers:
                    if recoiler.n == parent:
                        recoilers.remove(recoiler)
                        break
        # Loop over the counterterm's first-level currents
        for node in reversed(counterterm.subcurrents):
            # Alias singular structure for convenience
            ss = node.current['singular_structure']
            # Do not accept nested currents
            if node.subcurrents:
                raise MadGraph5Error(self.cannot_handle)
            # Get parent and children numbers
            # TODO Read these from cache
            parent, children = get_structure_numbers(
                ss, counterterm.momenta_dict
            )
            # Initialize the phase-space point
            momenta = LorentzVectorDict()
            # Compute the jacobian and map to higher multiplicity
            if ss.name() == 'C':
                # Make a fake collinear structure for soft-collinears
                # in order to use the same mapping
                new_ss = ss
                if new_ss.substructures:
                    new_ss = subtraction.CollStructure(ss.get_all_legs())
                momenta[parent] = LorentzVector(point[parent])
                jacobian *= self.collinear_map.map_to_higher_multiplicity(
                    point,
                    subtraction.SingularStructure(new_ss, *recoilers),
                    counterterm.momenta_dict,
                    kinematic_variables
                )
            elif ss.name() == 'S':
                jacobian *= self.soft_map.map_to_higher_multiplicity(
                    point,
                    subtraction.SingularStructure(ss, *recoilers),
                    counterterm.momenta_dict,
                    kinematic_variables
                )
            else:
                raise MadGraph5Error(self.cannot_handle)
            # Prepend pair of this current and the momenta,
            # deep copy wanted
            momenta.update(point.get_copy())
            current_PS_pairs.insert(0, (node.current, momenta))

        # Return
        return {
            'currents': current_PS_pairs,
            'matrix_element': ME_PS_pair,
            'jacobian': jacobian,
            'resulting_PS_point': point
        }

    def rescale_kinematic_variables(
        self, counterterm, kinematic_variables, scaling_parameter
    ):

        new_variables = {
            key: type(val)(val) for (key, val) in kinematic_variables.items()
        }
        # Determine 'equivalent power' of rescaling
        # This structure has at most two levels - checked later
        total_limits = 0
        for node in counterterm.subcurrents:
            total_limits += 1 + len(node.subcurrents)
        scaling_parameter **= 1./total_limits
        for node in counterterm.subcurrents:
            ss = node.current['singular_structure']
            if ss.name() == 'C':
                # For collinear clusters, rescale the virtuality
                parent, children = get_structure_numbers(
                    ss, counterterm.momenta_dict
                )
                new_variables['s' + str(parent)] *= scaling_parameter**2
                # For soft-collinears, rescale z's on top
                for subnode in node.subcurrents:
                    subss = subnode.current['singular_structure']
                    if subss.name() == 'S':
                        # Check depth
                        assert len(subss.subcurrents) == 0
                        for leg in subss.legs:
                            new_variables['z' + str(leg.n)] *= scaling_parameter
                    else:
                        raise MadGraph5Error(self.cannot_handle)
            elif ss.name() == 'S':
                # For soft clusters, rescale the whole momenta
                # Check depth
                assert len(ss.subcurrents) == 0
                for leg in ss.legs:
                    new_variables['p' + str(leg.n)] *= scaling_parameter
            else:
                raise MadGraph5Error(self.cannot_handle)
        return new_variables

class CataniSeymourWalker(VirtualWalker):
    pass

class NagySoperWalker(VirtualWalker):
    pass

# Mapping classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Mapping.
# Notice that this must be placed after all the Mapping daughter classes in this module have been declared.
mapping_walker_classes_map = {
    'FlatCollinear': FlatCollinearWalker,
    'SimpleNLO': SimpleNLOWalker,
    'CataniSeymour': CataniSeymourWalker,
    'NagySoper': NagySoperWalker,
    'Unknown': None
}

#===============================================================================
# Standalone main for debugging / standalone trials
#===============================================================================
if __name__ == '__main__':

    import random

    E_cm  = 5000.0

    # Try to run the above for a 2->8.
    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [100. + 10.*i for i in range(8)],
                                            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
    # Try to run the above for a 2->1.    
    #    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [5000.0])
    
    random_variables = [random.random() for _ in range(my_PS_generator.nDimPhaseSpace())]

    momenta, wgt = my_PS_generator.generateKinematics(E_cm, random_variables)
   
    print "\n ========================="
    print " ||    PS generation    ||"
    print " ========================="

    print "\nRandom variables :\n",random_variables
    print "\n%s\n"%momenta.__str__(n_initial=my_PS_generator.n_initial)
    print "Phase-space weight : %.16e\n"%wgt,

    variables_reconstructed, wgt_reconstructed = \
                                         my_PS_generator.invertKinematics(E_cm, momenta)

    print "\n ========================="
    print " || Kinematic inversion ||"
    print " ========================="
    print "\nReconstructed random variables :\n",variables_reconstructed
    differences = [
        abs(variables_reconstructed[i]-random_variables[i])
        for i in range(len(variables_reconstructed))
    ]
    print "Reconstructed weight = %.16e"%wgt_reconstructed
    if differences:
        print "\nMax. relative diff. in reconstructed variables = %.3e"%\
            max(differences[i]/random_variables[i] for i in range(len(differences)))
    print "Rel. diff. in PS weight = %.3e\n"%((wgt_reconstructed-wgt)/wgt)
