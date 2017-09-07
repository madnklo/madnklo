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
import array

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
pjoin = os.path.join

class PhaseSpaceGeneratorError(Exception):
    """Exception raised if an exception is triggered in integrators.""" 

#===============================================================================
# Vector
#===============================================================================

class Vector(array.array):

    def __new__(cls, arg1, arg2 = 0):

        if isinstance(arg1, int):
            return super(Vector, cls).__new__(cls, 'd', (arg2, ) * arg1)
        assert len(arg1) > 0
        return super(Vector, cls).__new__(cls, 'd', arg1)

    def __iadd__(self, y):

        for i in range(min(len(self), len(y))):
            self[i] += y[i]
        return self

    def __add__(self, y):

        new_vector = type(self)(self)
        new_vector += y
        return new_vector

    def __isub__(self, y):

        for i in range(min(len(self), len(y))):
            self[i] -= y[i]
        return self

    def __sub__(self, y):

        new_vector = type(self)(self)
        new_vector -= y
        return new_vector

    def __imul__(self, x):

        for i in range(len(self)):
            self[i] *= x
        return self

    def __mul__(self, x):

        tmp = type(self)(self)
        tmp *= x
        return tmp

    __rmul__ = __mul__

    def __idiv__(self, x):

        self *= 1. / x
        return self

    def __div__(self, x):

        return self * (1./x)

    def __neg__(self):

        tmp = type(self)(self)
        tmp *= -1
        return tmp

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

#===============================================================================
# LorentzVector
#===============================================================================

class LorentzVector(Vector):

    def __new__(cls, *args):

        return super(LorentzVector, cls).__new__(cls, *args)

    def dot(self, v):

        assert len(self) == len(v)
        return self[0]*v[0] - sum(self[i]*v[i] for i in range(1, len(self)))

#===============================================================================
# Lorentz5Vector
#===============================================================================
class Lorentz5Vector(list):
    """ A convenient class to manipulate Lorentz 4-vectors while keeping
    track of their mass component.
    Format used: (E, p_x, p_y, p_z, mass)
    """

    _TINY = 1.0e-5
    _MAX_OUTPUT = 1.0e8

    def __init__(self, *args, **opts):
        super(Lorentz5Vector,self).__init__(*args, **opts)
        
        if len(self)==0:
            self.extend([0., 0., 0., 0., 0.])
        if len(self)==4:
            M2 = self[0] ** 2 - self.rho2()
            if M2 >= 0:
                self.append(math.sqrt(M2))
            else:
                self.append(0)
    
    def rescaleEnergy(self):
        self[0]=math.sqrt(self[1]**2+self[2]**2+self[3]**2+self[4]**2)

    def rho2(self):
        """ Radius squared."""
        return self[1]**2+self[2]**2+self[3]**2

    def rho(self):
        """ Radius """
        return math.sqrt(self[1]**2+self[2]**2+self[3]**2)

    def pt(self):
        """ Compute transverse momentum."""
        return math.sqrt(self[1]**2 + self[2]**2)

    def deltaR(self, p2):
        """ Compute the deltaR separation with momentum p2"""

        delta_eta = self.pseudoRap() - p2.pseudoRap()
        delta_phi = self.getdelphi(p2)

        return math.sqrt(delta_eta**2 + delta_phi**2)

    def pseudoRap(self):
        """ Return pseudo-rapidity."""

        pt = self.pt()
        if pt < self._TINY and abs(self[3]) < self._TINY:
            return self._MAX_OUTPUT*(self[3]/abs(self[3]))
        th = math.atan2(pt, self[3])
        return -math.log(math.tan(th/2.))

    def getdelphi(self,p2):
        """ Return the phi-angle separation with p2."""
        pt1 = self.pt()
        pt2 = p2.pt()
        if pt1 == 0. or pt2 == 0.:
            return self._MAX_OUTPUT
        tmp = self[1]*p2[1] + self[2]*p2[2]
        tmp /= (pt1*pt2)
        if abs(tmp) > (1.0+self._TINY):
            raise PhaseSpaceGeneratorError("Cosine larger than 1. in phase-space cuts.")
        if abs(tmp) > 1.0:
            return math.acos(tmp/abs(tmp))
        return math.acos(tmp)

    def boostVector(self):
        if self[0] == 0.:
            if self.rho2() == 0.:
                return (0.,0.,0.)
            else:
                raise PhaseSpaceGeneratorError(
                    "Attempting to compute a boost from a reference vector with zero energy.") 
        if self[4] < 0.:
            raise PhaseSpaceGeneratorError(
                    "Attempting to compute a boost from a reference vector with negative mass.") 

        return tuple(_*(1./self[0]) for _ in self[1:4]) 

    def calculateMass(self):
        return math.sqrt(self[0]**2-self.rho2())

    def setMass(self, mass):
        self[4] = mass

    def getMass(self):
        return self[4]

    def cosTheta(self):
        ptot = self.rho()
        assert( ptot > 0. )
        return self[3] / ptot

    def phi(self):
        return math.atan2(self[2],self[1])
    
    def boost(self, boost_vector, gamma=-1.):
        """ Transport self into the restframe of the boostvector in argument.
        This means that the following command, for any vector p=(E, px, py, pz, M)
            p.boost(p.boostVector())
        transforms p to (M,0,0,0,M).
        """
        bx, by, bz = boost_vector[0], boost_vector[1], boost_vector[2]
        b2 = bx**2 + by**2 + bz**2
        if(gamma < 0.):
            gamma = 1.0 / math.sqrt(1.0 - b2)

        bp = bx*self[1] + by*self[2] + bz*self[3]
        gamma2 = (gamma -1.0) / b2 if b2 > 0 else 0.
        self[1] = self[1] + gamma2*bp*bx + gamma*bx*self[0]
        self[2] = self[2] + gamma2*bp*by + gamma*by*self[0]
        self[3] = self[3] + gamma2*bp*bz + gamma*bz*self[0]
        self[0] = gamma*(self[0] + bp)

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
                ref_lab.setMass(ref_lab.calculateMass())
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

    @classmethod
    def nice_momenta_string(cls, momenta,recompute_mass=False, n_initial=2):
        """ Nice printout of the momenta."""

        # Use padding for minus signs
        def special_float_format(float):
            return '%s%.16e'%('' if float<0.0 else ' ',float)
        
        cols_widths = [4,25,25,25,25,25]
        template = ' '.join('%%-%ds'%col_width for col_width in cols_widths)
        line     = '-'*(sum(cols_widths)+len(cols_widths)-1)

        out_lines = [template%('#', ' E', ' p_x',' p_y', ' p_z', ' M',)]
        out_lines.append(line)
        running_sum = Lorentz5Vector()
        for i, m in enumerate(momenta):
            mom = Lorentz5Vector(m)
            if i < n_initial:
                running_sum = running_sum + mom
            else:
                running_sum = running_sum - mom
            out_lines.append(template%tuple(['%d'%(i+1)]+
               [special_float_format(el) for el in mom[:4]+
               ([math.sqrt(abs(mom[0]**2-mom.rho2()))] if recompute_mass else [mom[4]]) ]))
        out_lines.append(line)
        out_lines.append(template%tuple(['Sum']+
               [special_float_format(el) for el in running_sum[:4]]+['']))

        return '\n'.join(out_lines)


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
            output_momenta[0] = Lorentz5Vector([self.initial_masses[0] , 0., 0., 0., self.initial_masses[0]])
            return

        output_momenta[0] = Lorentz5Vector([E_cm/2.0 , 0., 0., E_cm/2.0, 0.])
        output_momenta[1] = Lorentz5Vector([E_cm/2.0 , 0., 0., -E_cm/2.0, 0.])
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
        
        return PS_point, wgt, xb_1, xb_2

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
            output_momenta.append(Lorentz5Vector([self.masses[0]/2.,0.,0.,self.masses[0]/2.,0.]))
            output_momenta.append(Lorentz5Vector([self.masses[0]/2.,0.,0.,-self.masses[0]/2.,0.]))
            output_momenta.append(Lorentz5Vector([self.masses[0]   ,0.,0.,0.,self.masses[0]]))
            weight = self.get_flatWeights(E_cm,1)
            return output_momenta, weight
  
        M    = [ 0. ]*(self.n_final-1)
        M[0] = E_cm

        weight *= self.generateIntermediatesMassive(M, E_cm, random_variables)
        M.append(self.masses[-1])

        Q     = Lorentz5Vector([M[0],0.,0.,0.,M[0]])
        nextQ = Lorentz5Vector()

        for i in range(self.n_initial+self.n_final-1):
            
            if i < self.n_initial:
                output_momenta.append(Lorentz5Vector())
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
            
            p = Lorentz5Vector([0.,
                                q*cos_phi*sin_theta, 
                                q*sin_phi*sin_theta,
                                q*cos_theta,
                                self.masses[i-self.n_initial]])
            p.rescaleEnergy()
            p.boost(Q.boostVector())
            p.rescaleEnergy()
            output_momenta.append(p)

            nextQ = Q - p
            nextQ.setMass(M[i-self.n_initial+1])
            nextQ.rescaleEnergy()

            Q = nextQ
       
        output_momenta.append(Q)

        self.setInitialStateMomenta(output_momenta, E_cm)

        return output_momenta, weight

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
        moms = [Lorentz5Vector(mom) for mom in momenta]

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

        Q     = [Lorentz5Vector()]*(self.n_final-1)
        Q[0]  = Lorentz5Vector([M[0],0.,0.,0.,M[0]])

        for i in range(2,self.n_final):
            for k in range(i, self.n_final+1):
                Q[i-1] = Q[i-1] + moms[k+self.n_initial-1]
            M[i-1] = Q[i-1].calculateMass()

        weight = self.invertIntermediatesMassive(M, E_cm, random_variables)

        for i in range(self.n_initial,self.n_final+1):
            p = Lorentz5Vector(moms[i])
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

def get_parent_number(structure, momenta_dict):
    """Return the number of the parent of a given structure
    according to some momenta dictionary.
    """

    if structure.name() == "S":
        return None
    return momenta_dict.inv[
        frozenset((leg.n for leg in structure.get_all_legs()))
    ]

class VirtualMapping(object):
    """ A virtual class from which all Mapping implementations must inherit."""
    
    def __new__(cls, **opts):

        if cls is VirtualMapping:
            map_type = opts.pop('map_type') if 'map_type' in opts else 'Unknown'
            if map_type not in Mapping_classes_map or not Mapping_classes_map[map_type]:
                raise MadGraph5Error("Could not determine the class for the mapping of type '%s'."%str(map_type))
            target_class = Mapping_classes_map[map_type]
            return super(VirtualMapping, cls).__new__(target_class, **opts)
        else:
            return super(VirtualMapping, cls).__new__(cls, **opts)

    # TODO Review this
    def __init__(self, model = None, **opts):
        """General initialization of any mapping.
        n_legs_mapped is the number of legs 'generated' or 'removed' by the mapping.
        At NLO n_legs_mapped would be one for instance.
        The model is an optional specification that can be useful to know properties of the leg mapped.
        """
        
        # This attribute lists the name of the kinematic variables defining this mapping.
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
        """Map a given phase-space point by clustering the substructures
        and recoiling against the legs specified in singular_structure.

        :param PS_point: higher-multiplicity phase-space point,
        as a dictionary that associates integers to Lorentz5Vector's,
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
        self, PSpoint, singular_structure, momenta_dict, kinematic_variables
    ):
        """ Map the specified PS point onto another one using the kinematic variables 
        specified in the option 'kinematic_variables' and the specific splitting structure indicated.
        This function returns:
            mapped_PS_point, jacobian
        where "mapped_PS_point" is the mapped PS point and jacobian is the phase-space weight associated.
        """
        
        raise NotImplemented

    def approach_limit(
        self, singular_structure, ordering_parameter, starting_point = None
    ):
        """ Scale down starting variables given in starting_point (if specified) using the ordering_parameter
        to get closer to the limit specified by the splitting structure."""
        
        raise NotImplemented

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

        kinematic_variables = []
        # For every collinear subset of N particles,
        # a set of (N-1) z's and 2-component kt's is needed
        # plus the set's virtuality
        for substructure in singular_structure.substructures:
            # Add direct legs, last one not needed because of sum rules
            for leg in substructure.legs[:-1]:
                kinematic_variables += [
                    'z' + str(leg.n), 'kt' + str(leg.n)
                ]
            # Add virtuality of this subset
            kinematic_variables += [
                's' + str(get_parent_number(substructure, momenta_dict))
            ]

        return kinematic_variables

    def get_collinear_variables(self, PS_point, parent, children, variables):
        """Given unmapped and mapped momenta, compute the kinematic variables
        that describe the internal structure of particles going unresolved.
        Parent and children are indices that already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        """

        # WARNING This function assumes massless particles

        # Retrieve the parent's momentum
        p = PS_point[parent]
        # Compute the sum of momenta
        q = LorentzVector(4)
        for i in children:
            q += PS_point[i]
        # Pre-compute scalar products
        q2 = q.square()
        pq = p.dot(q)
        # Compute all kinematic variables
        for i in children[:-1]:
            pi = PS_point[i]
            ppi = p.dot(pi)
            qpi = q.dot(pi)
            zi = 2*qpi/q2 - ppi/pq
            kti = pi + (q2*ppi - pq*qpi)/(pq**2) * p - ppi/pq * q
            variables['z' + str(i)] = zi
            variables['kt' + str(i)] = kti
        return

    def set_collinear_variables(self, PS_point, parent, children,
                                total_momentum, variables):
        """Given the lower multiplicity parent momentum,
        the total momentum of children and collinear variables,
        compute and set the children momenta.
        Parent and children are indices that already refer to the position
        of momenta within the PS_point (no momentum dictionary used).
        """

        # WARNING This function assumes massless particles

        # Retrieve the parent's momentum
        p = PS_point[parent]
        # Rename the sum of momenta
        q = total_momentum
        # Pre-compute scalar products
        q2 = q.square()
        pq = p.dot(q)
        # Variables for sums
        z_sum = 0
        kt_sum = LorentzVector(4)
        # Set momenta for all children but the last
        for i in children[:-1]:
            zi = variables['z' + str(i)]
            kti = variables['kt' + str(i)]
            z_sum += zi
            kt_sum += kti
            kti2 = kti.square()
            PS_point[i] = kti + (zi*zi*q2+kti2)/(2*zi*pq) * p - kti2/(zi*q2) * q
        # Set momentum of the last child
        zj = 1 - z_sum
        ktj = -kt_sum
        ktj2 = ktj.square()
        PS_point[children[-1]] = ktj + (zj*zj*q2+ktj2)/(2*zj*pq) * p - ktj2/(zj*q2) * q
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
        assert self.is_valid_structure(singular_structure)

        # Build the total momentum of collinear sets and recoilers
        Q = sum(
            PS_point[leg.n]
            for leg in singular_structure.get_all_legs()
        )
        Q2 = Q.square()
        # Precompute sets and numbers...
        cluster = singular_structure.substructures[0]
        children_numbers = frozenset((leg.n for leg in cluster.legs))
        parent_number = momenta_dict.inv[children_numbers]
        # Build the collective momentum
        pC = sum(PS_point[j] for j in children_numbers)
        # Compute the parameter alpha for this collinear subset
        QpC = pC.dot(Q)
        sC = pC.square()
        alpha = (QpC - math.sqrt(QpC**2 - sC)) / Q2
        # Compute the common scale factor
        scale_factor = 1. / (1.-alpha)
        # Map all recoilers' momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler] *= scale_factor
        # Map the cluster's momentum
        PS_point[parent_number] = scale_factor * (pC - alpha * Q)
        # If needed, update the kinematic_variables dictionary
        if kinematic_variables:
            kinematic_variables['s' + str(parent_number)] = sC
            self.get_collinear_variables(
                PS_point, parent_number, children_numbers, kinematic_variables
            )

        # TODO Compute the jacobian for this mapping
        jacobian = 1.0

        return jacobian

    def map_to_higher_multiplicity(
        self, PS_point, singular_structure, momenta_dict, variables
    ):

        # Consistency checks
        assert isinstance(momenta_dict, subtraction.bidict)
        assert self.is_valid_structure(singular_structure)
        assert (
            set(variables.keys()) ==
            set(self.get_kinematic_variables_names(singular_structure,
                                                   momenta_dict))
        )

        # Find the parent's momentum
        parent = get_parent_number(singular_structure.substructures[0])
        children = singular_structure.substructures[0].get_all_legs()
        qC = PS_point[parent]
        # Compute the momentum of recoilers and the total
        qR = sum(PS_point[i] for i in singular_structure.legs)
        Q = qR + qC
        # Compute scalar products
        QqC = Q.dot(qC)
        qC2 = qC.square()
        Q2 = Q.square()
        qR2 = qR.square()
        sC = variables['s'+str(parent)]
        # Obtain parameter alpha
        alpha = 0
        if qR2 == 0:
            alpha = 0.5 * (sC - qC2) / (QqC - qC2)
        else:
            alpha = (math.sqrt(QqC**2 - Q2*qC2 + sC*qR2) - (QqC - qC2)) / qR2
        # Compute reverse-mapped momentum
        pC = (1-alpha) * qC - alpha*Q
        # Set children momenta
        self.set_collinear_variables(PS_point, parent, children, pC, variables)
        # TODO What to do with the parent's momentum? Leave as is?
        # Map recoil momenta
        for recoiler in singular_structure.legs:
            PS_point[recoiler] *= 1-alpha

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
        kinematic_variables = None
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
        p_C = [Lorentz5Vector(), ] * n_clusters
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
    print "\n%s\n"%my_PS_generator.nice_momenta_string(
                        momenta, recompute_mass=True, n_initial=my_PS_generator.n_initial)
    print "Phase-space weight : %.16e\n"%wgt,

    variables_reconstructed, wgt_reconstructed = \
                                         my_PS_generator.invertKinematics(E_cm, momenta)

    print "\n ========================="
    print " || Kinematic inversion ||"
    print " ========================="
    print "\nReconstructed random variables :\n",variables_reconstructed
    differences = [abs(variables_reconstructed[i]-random_variables[i]) 
                                    for i in range(len(variables_reconstructed))]
    print "Reconstructed weight = %.16e"%wgt_reconstructed
    if differences:
        print "\nMax. relative diff. in reconstructed variables = %.3e"%\
            max(differences[i]/random_variables[i] for i in range(len(differences)))
    print "Rel. diff. in PS weight = %.3e\n"%((wgt_reconstructed-wgt)/wgt)


# Mapping classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Mapping.
# Notice that this must be placed after all the Mapping daughter classes in this module have been declared.
Mapping_classes_map = {('single-real', 'NLO'): MappingCataniSeymourFFOne,
                       'Unknown': None}