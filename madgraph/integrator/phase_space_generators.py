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

try:
    import madgraph
except ImportError:
    MADEVENT = True
    import internal.misc as misc

else:
    MADEVENT= False
    import madgraph.various.misc as misc

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')
pjoin = os.path.join

class Lorentz5Vector(list):
    """ A convenient class to manipulate Lorentz 4-vectors while keeping
    track of their mass component.
    Format used: (E, p_x, p_y, p_z, mass)
    """

    def __init__(self, *args, **opts):
        super(Lorentz5Vector,self).__init__(*args, **opts)
        
        if len(self)==0:
            self.extend([0., 0., 0., 0., 0.])
        if len(self)==4:
            self.append(abs(math.sqrt(self[0]**2-self.rho2())))
    
    def rescaleEnergy(self):
        self[0]=math.sqrt(self[1]**2+self[2]**2+self[3]**2+self[4]**2)

    def rho2(self):
        """ Radius squared."""
        return self[1]**2+self[2]**2+self[3]**2

    def rho(self):
        """ Radius """
        return math.sqrt(self[1]**2+self[2]**2+self[3]**2)

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

    def __add__(self, y):
        newvector = Lorentz5Vector()
        newvector[0] = self[0] + y[0]
        newvector[1] = self[1] + y[1]
        newvector[2] = self[2] + y[2]
        newvector[3] = self[3] + y[3]
        return newvector 

    def __sub__(self, y):
        newvector = Lorentz5Vector()
        newvector[0] = self[0] - y[0]
        newvector[1] = self[1] - y[1]
        newvector[2] = self[2] - y[2]
        newvector[3] = self[3] - y[3]
        return newvector 

class PhaseSpaceGeneratorError(Exception):
    """Exception raised if an exception is triggered in integrators.""" 
    
class VirtualPhaseSpaceGenerator(object):

    def __init__(self, initial_masses, final_masses):
        
        self.initial_masses = initial_masses
        self.masses         = final_masses
        self.n_initial      = len(initial_masses)
        self.n_final        = len(final_masses)

    def generateKinematics(random_variables):
        raise NotImplementedError
 
    def nDimPhaseSpace(self):
        """ Return the number of random numbers required to produce a given
        multiplicity final state. """

        if self.n_final == 1:
            return 1
        return 3*self.n_final -4

    @staticmethod
    def nice_momenta_string(momenta,recompute_mass=False, n_initial=2):
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

    @classmethod
    def get_flatWeights(cls, E_cm, n):
        """ Return the phase-space volume for a n massless final states.
          Vol(E_cm, n) = (pi/2)^(n-2) *  (E_cm^2)^(n-2) / ((n-1)!*(n-2)!)
        """
        
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
            raise PhaseSpaceGeneratorError(
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
                raise PhaseSpaceGeneratorError(
                    "Phase-space of massive 2-particle collision not implemented yet (trivial though).")

        if self.n_initial == 1:
            output_momenta[0] = Lorentz5Vector([self.initial_masses[0] , 0., 0., 0., self.initial_masses[0]])
            return

        output_momenta[0] = Lorentz5Vector([E_cm/2.0 , 0., 0., E_cm/2.0, 0.])
        output_momenta[1] = Lorentz5Vector([E_cm/2.0 , 0., 0., -E_cm/2.0, 0.])
        return

    def generateKinematics(self, E_cm, random_variables):
        """ Generate a self.n_initial -> self.n_final phase-space point using 
        the random variables passed in argument."""

        # Make sure the right number of random variables are passed
        assert (len(random_variables)==self.nDimPhaseSpace())

        # The distribution weight of the generate PS point
        weight = 1.
        
        M    = [ 0. ]*(self.n_final-1)
        M[0] = E_cm

        weight *= self.generateIntermediatesMassive(M, E_cm, random_variables)
        M.append(self.masses[-1])

        output_momenta = []

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


if __name__ == '__main__':

    import random

    # Try to run the above for a 2->2.
    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [100. + 10.*i for i in range(8)])
    
    E_cm  = 5000.0
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
    print "\nMax. relative diff. in reconstructed variables = %.3e"%\
            max(differences[i]/random_variables[i] for i in range(len(differences)))
    print "Rel. diff. in PS weight = %.3e\n"%((wgt_reconstructed-wgt)/wgt)
