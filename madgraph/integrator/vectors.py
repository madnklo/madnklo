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
import copy
import numpy as np

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')


class InvalidOperation(Exception):
    pass

#=========================================================================================
# Vector
#=========================================================================================

class Vector(np.ndarray):

    def __new__(cls, *args, **opts):

        if args and isinstance(args[0], Vector):
            foo = args[0].get_copy()
        else:
            foo = np.asanyarray(*args, **opts).view(cls)
        return foo

    def eps(self):

        try: return np.finfo(self.dtype).eps
        except: return 0

    def huge(self):

        if np.issubdtype(self.dtype, np.inexact):
            return np.finfo(self.dtype).max
        elif np.issubdtype(self.dtype, np.integer):
            return np.iinfo(self.dtype).max
        else:
            raise ValueError

    def __eq__(self, other):

        eps = max(self.eps(), other.eps())
        # numpy.allclose uses abs(self-other) which would compute the norm
        # for LorentzVector, thus view self and other as numpy.ndarray's
        return np.allclose(
            self.view(type=np.ndarray), other.view(type=np.ndarray),
            math.sqrt(eps), 0.
        )

    def __ne__(self, other):

        return not self.__eq__(other)

    def __hash__(self):

        return tuple(x for x in self).__hash__()

    def get_copy(self):

        # The vector instantiated by get_copy() should be modified
        # without changing the previous instance, irrespectively of the
        # (presumably few) layers that compose entries of the vector
        # return copy.deepcopy(self)
        return copy.copy(self)

    def dot(self, v):

        assert len(self) == len(v)
        return sum(el * v[i] for i, el in enumerate(self))

    def square(self):

        return self.dot(self)

    def __abs__(self):

        return math.sqrt(self.square())

    def normalize(self):

        self.__idiv__(abs(self))
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

#=========================================================================================
# LorentzVector
#=========================================================================================

class LorentzVector(Vector):

    def __new__(cls, *args, **opts):

        if len(args) == 0:
            return super(LorentzVector, cls).__new__(cls, [0., 0., 0., 0.], **opts)
        return super(LorentzVector, cls).__new__(cls, *args, **opts)

    def space(self):
        """Return the spatial part of this LorentzVector."""

        return self[1:].view(type=Vector)

    def dot(self, v):
        """Compute the Lorentz scalar product."""
        
        pos = self[0]*v[0]
        neg = self.space().dot(v.space())
        # if pos+neg != 0 and abs(2*(pos-neg)/(pos+neg)) < 100.*self.eps(): return 0
        return pos - neg

    def square_almost_zero(self):
        """Check if the square of this LorentzVector is zero within numerical accuracy."""

        return (self.square() / self.view(Vector).square()) ** 2 < self.eps()

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
        if negative: self[0] *= -1
        return self

    def rotoboost(self, p, q):
        """Apply the Lorentz transformation that sends p in q to this vector."""

        # NOTE: when applying the same Lorentz transformation to many vectors,
        #       this function goes many times through the same checks.

        # Compute squares
        p2 = p.square()
        q2 = q.square()
        # Ratios used to decide if the two vectors are massless
        p2_sq_ratio = abs(p2)/p.view(Vector).square()
        q2_sq_ratio = abs(q2)/q.view(Vector).square()
        # Numerical tolerances
        p_eps = math.sqrt(p.eps())
        q_eps = math.sqrt(q.eps())
        # Check if both Lorentz squares are small compared to the euclidean squares,
        # in which case the alternative formula should be used
        if p2_sq_ratio < p_eps and q2_sq_ratio < q_eps:
            # Use alternative formula
            if p == self:
                for i in range(len(self)):
                    self[i] = q[i]
            else:
                logger.critical("Error in vectors.rotoboost: missing formula")
                logger.critical("Boosting %s (%.9e)" % (str(self), self.square()))
                logger.critical("p = %s (%.9e)" % (str(p), p2))
                logger.critical("q = %s (%.9e)" % (str(q), q2))
                logger.critical("Eq. (4.14) of arXiv:0706.0017v2, p. 26 not implemented")
                raise NotImplemented
            return self
        else:
            # Check that the two invariants are close,
            # else the transformation is invalid
            if abs(p2-q2)/(abs(p2)+abs(q2)) > (p_eps+q_eps):
                logger.critical("Error in vectors.rotoboost: nonzero, unequal squares")
                logger.critical("p = %s (%.9e)" % (str(p), p2))
                logger.critical("q = %s (%.9e)" % (str(q), q2))
                logger.critical("square ratios:")
                logger.critical("p: %.9e (vs %.9e)" % (p2_sq_ratio, p_eps))
                logger.critical("q: %.9e (vs %.9e)" % (q2_sq_ratio, q_eps))
                raise InvalidOperation
            # Compute scalar products
            pq = p + q
            pq2 = pq.square()
            p_s = self.dot(p)
            pq_s = self.dot(pq)
            # Assemble vector
            self.__iadd__(2 * ((p_s/q2) * q - (pq_s/pq2) * pq))
            return self

    def pt(self, axis=3):
        """Compute transverse momentum."""

        return math.sqrt(
            sum(self[i]**2 for i in range(1, len(self)) if i != axis) )

    def pseudoRap(self):
        """Compute pseudorapidity."""

        pt = self.pt()
        if pt < self.eps() and abs(self[3]) < self.eps():
            return self.huge()*(self[3]/abs(self[3]))
        th = math.atan2(pt, self[3])
        return -math.log(math.tan(th/2.))

    def rap(self):
        """Compute rapidity in the lab frame. (needs checking)"""

        if self.pt() < self.eps() and abs(self[3]) < self.eps():
            return self.huge()*(self[3]/abs(self[3]))

        return .5*math.log((self[0]+self[3])/(self[0]-self[3]))

    def getdelphi(self, p2):
        """Compute the phi-angle separation with p2."""

        pt1 = self.pt()
        pt2 = p2.pt()
        if pt1 == 0. or pt2 == 0.:
            return self.huge()
        tmp = self[1]*p2[1] + self[2]*p2[2]
        tmp /= (pt1*pt2)
        if abs(tmp) > (1.0+math.sqrt(self.eps())):
            logger.critical("Cosine larger than 1. in phase-space cuts.")
            raise ValueError
        if abs(tmp) > 1.0:
            return math.acos(tmp/abs(tmp))
        return math.acos(tmp)

    def deltaR(self, p2):
        """Compute the deltaR separation with momentum p2."""

        delta_eta = self.pseudoRap() - p2.pseudoRap()
        delta_phi = self.getdelphi(p2)
        return math.sqrt(delta_eta**2 + delta_phi**2)

    def boostVector(self):

        if self == LorentzVector():
            return Vector([0.] * 3)
        if self[0] <= 0. or self.square() < 0.:
            logger.critical("Attempting to compute a boost vector from")
            logger.critical("%s (%.9e)" % (str(self), self.square()))
            raise InvalidOperation
        return self.space()/self[0]

    def cosTheta(self):

        ptot = self.rho()
        assert (ptot > 0.)
        return self[3] / ptot

    def phi(self):

        return math.atan2(self[2], self[1])
    
    def boost(self, boost_vector, gamma=-1.):
        """Transport self into the rest frame of the boost_vector in argument.
        This means that the following command, for any vector p=(E, px, py, pz)
            p.boost(-p.boostVector())
        transforms p to (M,0,0,0).
        """

        b2 = boost_vector.square()
        if gamma < 0.:
            gamma = 1.0 / math.sqrt(1.0 - b2)

        bp = self.space().dot(boost_vector)
        gamma2 = (gamma-1.0) / b2 if b2 > 0 else 0.
        factor = gamma2*bp + gamma*self[0]
        self_space = self.space()
        self_space += factor*boost_vector
        self[0] = gamma*(self[0] + bp)

#=========================================================================================
# LorentzVectorDict
#=========================================================================================

class LorentzVectorDict(dict):
    """A simple class wrapping dictionaries that store Lorentz vectors."""

    def to_list(self):
        """Return list copy of self. Notice that the actual values of the keys
        are lost in this process."""

        return LorentzVectorList(self[k] for k in sorted(self.keys()))

    def to_dict(self):
        """Return a copy of this LorentzVectorDict """

        return LorentzVectorDict(self)

    def to_tuple(self):
        """Return a copy of this LorentzVectorDict as an immutable tuple.
        Notice that the actual values of the keys are lost in this process.
        """

        return tuple( tuple(self[k]) for k in sorted(self.keys()) )

    def __str__(self, n_initial=2):
        """Nice printout of the momenta."""

        # Use padding for minus signs
        def special_float_format(fl):
            return '%s%.16e' % ('' if fl < 0.0 else ' ', fl)

        cols_widths = [4, 25, 25, 25, 25, 25]
        template = ' '.join(
            '%%-%ds' % col_width for col_width in cols_widths
        )
        line = '-' * (sum(cols_widths) + len(cols_widths) - 1)

        out_lines = [template % ('#', ' E', ' p_x', ' p_y', ' p_z', ' M',)]
        out_lines.append(line)
        running_sum = LorentzVector()
        for i in sorted(self.keys()):
            mom = LorentzVector(self[i])
            if i <= n_initial:
                running_sum += mom
            else:
                running_sum -= mom
            out_lines.append(template % tuple(
                ['%d' % i] + [
           special_float_format(el) for el in (list(mom) + [math.sqrt(abs(mom.square()))])
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

        return type(self)((i,LorentzVector(k)) for i,k in self.items())

#=========================================================================================
# LorentzVectorList
#=========================================================================================

class LorentzVectorList(list):
    """A simple class wrapping lists that store Lorentz vectors."""

    def __str__(self, n_initial=2):
        """Nice printout of the momenta."""

        return LorentzVectorDict(
            (i + 1, v) for i, v in enumerate(self)
        ).__str__(n_initial=n_initial)

    def to_list(self):
        """Return list copy of self."""

        return LorentzVectorList(self)

    def to_tuple(self):
        """Return a copy of this LorentzVectorList as an immutable tuple."""

        return tuple( tuple(v) for v in self )

    def to_dict(self):
        """Return a copy of this LorentzVectorList as a LorentzVectorDict."""

        return LorentzVectorDict( (i+1, v) for i, v in enumerate(self) )

    def get_copy(self):
        """Return a copy that can be freely modified
        without changing the current instance.
        """

        return type(self)([LorentzVector(p) for p in self])
