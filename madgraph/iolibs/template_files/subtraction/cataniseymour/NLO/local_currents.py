##########################################################################################
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
##########################################################################################
"""Implementation of NLO colorful currents."""

import os
import math
import madgraph.integrator.mappings as mappings
import madgraph.various.misc as misc

try:
    # First try to import this in the context of the exported currents
    import SubtractionCurrents.subtraction_current_implementations_utils as utils
    import SubtractionCurrents.QCD_local_currents as currents
except ImportError:
    # If not working, then it must be within MG5_aMC context:
    import madgraph.iolibs.template_files.\
                   subtraction.subtraction_current_implementations_utils as utils
    import madgraph.iolibs.template_files.\
                   subtraction.QCD_local_currents as currents

pjoin = os.path.join

CurrentImplementationError = utils.CurrentImplementationError

#=========================================================================================
# Eikonal factor, modified by partial fractioning and without divergence
#=========================================================================================

def eikonal(pi, pj, pr):
    """Eikonal factor for soft particle with momentum pr
    emitted from the dipole with momenta pi and pj.
    """

    pipj = pi.dot(pj)
    prpi = pr.dot(pi)
    prpj = pr.dot(pj)
    return pipj / (prpi*prpj)

def mod_eikonal(pi, pj, pr):
    """Modified eikonal factor for soft particle with momentum pr
    emitted from the dipole with momenta pi and pj.
    This is obtained starting from the eikonal and:
    - ignoring 1 / sir, which is already included in the normalisation factor;
    - multiplying by the partial fraction sjr / (sir + sjr) to regulate for sjr -> 0.
    """

    pipj = pi.dot(pj)
    pijpr = pr.dot(pi+pj)
    return 2 * pipj / pijpr

def partial_fraction(pi, pj, pr, Q):

    sjr = (pj+pr).square()
    sir = (pi+pr).square()
    return sjr / (sir + sjr)
    # arj = mappings.FinalRescalingOneMapping.alpha(pr+pj, Q)
    # ari = mappings.FinalRescalingOneMapping.alpha(pr+pi, Q)
    # return arj / (ari + arj)

def pf_eikonal_pole(pi, qj, pr, Q, qir):
    """Partial-fractioned eikonal pole in the variables alpha_ir and zr,
    for soft particle with momentum pr
    emitted from a recoiler of momentum pi and spectator of (mapped) momentum qj.
    """

    pir = pi + pr
    air = mappings.FinalRescalingOneMapping.alpha(pir, Q)
    Q2 = Q.square()
    Qnorm = Q2 ** 0.5
    pirvec = pir - (Q.dot(pir)/Q2) * Q
    pirnorm = (-pirvec.square()) ** 0.5
    n = pirvec / pirnorm
    t = Q / Qnorm
    na = t + n
    nb = t - n
    zr = pr.dot(nb) / pir.dot(nb)
    bir = qir.dot(Q) / Q2
    prperp = pr - pr.dot(na) * nb/2 - pr.dot(nb) * na/2
    nperp = prperp / ((-prperp.square())**0.5)
    qjna = qj.dot(na)
    qjnb = qj.dot(nb)
    qjnperp = qj.dot(nperp)
    brk = air*(qjnb + 2*bir*Qnorm)
    brk = air*(2*qjnb + qjna)
    brk += 2*qjna*zr*bir
    brk += 2*qjnperp*((2*air*bir*zr)**0.5)
    return 2*qjna / (air*brk*Q2)

#=========================================================================================
# NLO final-collinear currents, containing the soft limits
#=========================================================================================

class QCD_final_collinear_0_qqx(currents.QCDLocalCollinearCurrent):
    """q q~ collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that the particles are a massless quark and its anti-quark in final-state
        if len(ss.legs) != 2: return None
        for leg in ss.legs:
            if not cls.is_quark(leg, model): return None
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        if not cls.are_antiparticles(ss.legs[0], ss.legs[1]): return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        return tuple(leg.n for leg in legs)

    def evaluate_kernel(self, zs, kTs, parent):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent, (kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.
        # WARNING multiplied by two because of flavor factors
        evaluation['values'][(0, 0)]['finite'] = 2. * self.TR
        evaluation['values'][(1, 0)]['finite'] = 2. * 4. * self.TR * z*(1.-z) / kT.square()
        return evaluation

class QCD_final_collinear_0_gq(currents.QCDLocalCollinearCurrent):
    """g q collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that all particles are massless final state
        for leg in ss.legs:
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        # Check that there are a quark and a gluon
        if len(ss.legs) != 2: return None
        if (cls.is_gluon(ss.legs[0], model) and cls.is_quark(ss.legs[1], model)):
            pass
        elif (cls.is_quark(ss.legs[0], model) and cls.is_gluon(ss.legs[1], model)):
            pass
        else:
            return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        if cls.is_gluon(legs[0], model): return (legs[0].n, legs[1].n)
        else: return (legs[1].n, legs[0].n)

    def evaluate_kernel(self, zs, kTs, parent, **opts):

        # Retrieve the collinear variables
        z = zs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        # We must subtract the soft-collinear (CxS *not* SxC) from this contribution:
        # P_gq           = self.CF * ((1.-z)**2 + 1.)/z
        # CxS(P_gq)      = self.CF * 2.*(1.-z) / z
        # SxC(P_gq)      = self.CF * 2. / z
        # P_gq-CxS(P_gq) = self.CF * z
        # P_gq-SxC(P_gq) = self.CF * ((1.-z)**2 - 1.)/z = self.CF * (z-2.)
        Q = opts['Q']
        pir = opts['pC']
        qir = opts['qC']
        air = mappings.FinalRescalingOneMapping.alpha(pir, Q)
        Q2 = Q.square()
        bir = qir.dot(Q) / Q2
        brk = air*(1. - 1./(2.*bir))
        evaluation['values'][(0, 0)]['finite'] = self.CF *(z - 2. + 2*brk/z)
        evaluation['values'][(0, 0)]['finite'] = self.CF * ((1.-z)**2 - 1.)/z
        return evaluation

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):

        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        parent = leg_numbers_map.inv[frozenset(children)]
        pC = sum(higher_PS_point[child] for child in children)
        qC = lower_PS_point[parent]
        ps = higher_PS_point[children[0]]
        pi = higher_PS_point[children[1]]
        pC2 = pC.square()
        alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
        Q2 = Q.square()
        beta = qC.dot(Q) / Q2
        pC2exp = 2*alpha*beta*Q2
        pC2alt = 2*alpha*Q2*(alpha*(1-2*beta)+beta)
        if self.is_cut(Q=Q, pC=pC):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        evaluation = self.evaluate_kernel(zs, kTs, parent, Q=Q, pC=pC, qC=qC, pr=ps)

        # Find all colored leg numbers except for the parent in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        color_correlation_index = 1

        # Now loop over the colored parton number pairs (parent, j)
        # and add the corresponding contributions to this current
        for j in all_colored_parton_numbers:
            # Write the eikonal for that pair
            if j == parent:
                continue
            pj = sum(higher_PS_point[child] for child in leg_numbers_map[j])
            # pj = higher_PS_point[j]
            qj = lower_PS_point[j]
            evaluation['color_correlations'].append(((parent, j),))
            # eiks = -pC2*eikonal(pi, qj, ps)*partial_fraction(pi, pj, ps, Q)
            eiks = -pC2exp*pf_eikonal_pole(pi, qj, ps, Q, qC)
            # eiks = -mod_eikonal(pi, pj, ps)
            # mod = (qj.dot(qC)) / (qj.dot(pi+ps))
            # qjmod = mod * qj
            # eiks = -mod_eikonal(pi, qjmod, ps)
            evaluation['values'][(0, color_correlation_index)] = {'finite': eiks}
            color_correlation_index += 1

        # Add the normalization factors
        norm = 8. * math.pi * alpha_s / pC2alt
        norm *= self.factor(Q=Q, pC=pC, qC=qC)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

class QCD_final_collinear_0_gg(currents.QCDLocalCollinearCurrent):
    """g g collinear tree-level current."""

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that the particles are two final-state massless gluons
        if len(ss.legs) != 2: return None
        for leg in ss.legs:
            if not cls.is_gluon(leg, model): return None
            if not cls.is_massless(leg, model): return None
            if cls.is_initial(leg): return None
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        return tuple(leg.n for leg in legs)

    def evaluate_kernel(self, zs, kTs, parent, **opts):

        # Retrieve the collinear variables
        z = zs[0]
        kT = kTs[0]
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent,( kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.

        # We must subtract the soft-collinear (CxS *not* SxC) from this contribution:
        # P_gg           = 2.*self.CA * ( (1.-z) / z + z / (1.- z) )
        # CxS(P_gg)      = 2.*self.CA * ( (1.-z) / z + z / (1.- z) )
        # SxC(P_gg)      = 2.*self.CA * ( 1 / z + 1 / (1.- z) )
        # P_gg-CxS(P_gg) = 0
        # P_gg-SxC(P_gg) = -4.*self.CA
        Q = opts['Q']
        pir = opts['pC']
        qir = opts['qC']
        air = mappings.FinalRescalingOneMapping.alpha(pir, Q)
        Q2 = Q.square()
        bir = qir.dot(Q) / Q2
        brk = (1. - 1./(2.*bir))
        brk *= air / z + air / (1-z)
        evaluation['values'][(0, 0)]['finite'] = 2.*self.CA*(-2.+brk)
        evaluation['values'][(0, 0)]['finite'] = -4.*self.CA
        evaluation['values'][(1, 0)]['finite'] = -2.*self.CA * 2.*z*(1.-z) / kT.square()
        return evaluation

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):

        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process." )
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        parent = leg_numbers_map.inv[frozenset(children)]
        pC = sum(higher_PS_point[child] for child in children)
        qC = lower_PS_point[parent]
        pC2 = pC.square()
        alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
        Q2 = Q.square()
        beta = qC.dot(Q) / Q2
        betaQ2 = qC.dot(Q)
        pC2exp = 2*alpha*beta*Q2
        pC2alt = 2*alpha*Q2*(alpha*(1-2*beta)+beta)
        if self.is_cut(Q=Q, pC=pC):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        evaluation = self.evaluate_kernel(zs, kTs, parent, Q=Q, pC=pC, qC=qC)

        # Find all colored leg numbers except for the parent in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        color_correlation_index = 1
        p0 = higher_PS_point[children[0]]
        p1 = higher_PS_point[children[1]]

        # Loop over the colored parton number pairs (parent, j)
        # and add the corresponding contributions to this current
        for j in all_colored_parton_numbers:
            # Write the eikonal for that pair
            if j == parent:
                continue
            # pj = higher_PS_point[j]
            # pj = sum(higher_PS_point[child] for child in leg_numbers_map[j])
            qj = lower_PS_point[j]
            # eik0 = -mod_eikonal(pj, p1, p0)
            # eik1 = -mod_eikonal(pj, p0, p1)
            # eik0 = -mod_eikonal(qj, qC, p0)
            # eik1 = -mod_eikonal(qj, qC, p1)
            # eik0 = -pC2*eikonal(p1, qj, p0)*partial_fraction(p1, qj, p0, Q)
            # eik1 = -pC2*eikonal(p0, qj, p1)*partial_fraction(p0, qj, p1, Q)

            eik0 = -pC2exp*pf_eikonal_pole(p0, qj, p1, Q, qC)
            eik1 = -pC2exp*pf_eikonal_pole(p1, qj, p0, Q, qC)

            evaluation['color_correlations'].append(((parent, j),))
            evaluation['values'][(0, color_correlation_index)] = {'finite': eik0 + eik1}
            color_correlation_index += 1

        # Add the normalization factors
        norm = 8. * math.pi * alpha_s / pC2alt
        norm *= self.factor(Q=Q, pC=pC, qC=qC)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

#=========================================================================================
# NLO initial-collinear currents, containing the soft limits
#=========================================================================================

class QCD_initial_collinear_0_gq(currents.QCDLocalCollinearCurrent):
    """gq collinear ISR tree-level current.
    q(initial) > g(initial_after_emission) q(final)
    """

    variables = staticmethod(currents.Q_initial_coll_variables)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')
        # Check that the particles are a massless quark and its anti-quark in final-state
        if len(ss.legs) != 2: return None
        
        n_initial_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and cls.is_initial(leg):
                n_initial_state_quarks += 1
        if n_initial_state_quarks != 1: return None
        
        n_final_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and not cls.is_initial(leg):
                n_final_state_quarks += 1
        if n_final_state_quarks != 1: return None
        
        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        kT = kTs[0]

        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent, (kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = -1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (gluon)
        initial_state_crossing_factor *= ((self.NC**2-1)/float(self.NC))
        
        z = 1./x

        # We re-use here the Altarelli-Parisi Kernel of the P_q\bar{q} final state kernel
        
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.

        norm = initial_state_crossing_factor * self.TR
        evaluation['values'][(0, 0)]['finite'] = norm
        evaluation['values'][(1, 0)]['finite'] = norm * 4. * z*(1.-z) / kT.square()

        return evaluation

class QCD_initial_collinear_0_qq(currents.QCDLocalCollinearCurrent):
    """ qq collinear ISR tree-level current.
    g(initial) > q(initial_after_emission) qx(final).
    """

    variables = staticmethod(currents.Q_initial_coll_variables)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')

        if len(ss.legs) != 2: return None

        n_initial_state_gluons = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and cls.is_initial(leg):
                n_initial_state_gluons += 1
        if n_initial_state_gluons != 1: return None
        
        n_final_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and not cls.is_initial(leg):
                n_final_state_quarks += 1
        if n_final_state_quarks != 1: return None

        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        
        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (gluon) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= (self.NC/float(self.NC**2-1))
        
        z = 1./x

        norm = initial_state_crossing_factor * self.CF
        # We re-use here the Altarelli-Parisi Kernel of the P_gq final state kernel without
        # the soft-subtractio term 2./z since the gluon is here in the initial state
        evaluation['values'][(0, 0)]['finite'] = norm * (1. + (1.-z)**2) / z

        return evaluation

class QCD_initial_collinear_0_qg(currents.QCDLocalCollinearCurrent):
    """qg collinear ISR tree-level current.
    q(initial) > q(initial_after_emission) g(final)
    """

    variables = staticmethod(currents.Q_initial_coll_variables)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')

        if len(ss.legs) != 2: return None

        n_initial_state_quarks = 0
        for leg in ss.legs:
            if cls.is_quark(leg, model) and cls.is_initial(leg):
                n_initial_state_quarks += 1
        if n_initial_state_quarks != 1: return None
        
        n_final_state_gluon = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and not cls.is_initial(leg):
                n_final_state_gluon += 1
        if n_final_state_gluon != 1: return None

        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        
        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None}}
        })
        
        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (quark) and the one of the reduced Born ME (quark)
        initial_state_crossing_factor *= 1.
        
        z = 1./x

        # We re-use here the Altarelli-Parisi Kernel of the P_qg final state kernel, including
        # its soft subtraction
        # We must subtract the soft-collinear (CxS *not* SxC) from this contribution:
        # P_qg           = self.CF * ( (1.+z**2)/(1.-z) )
        # CxS(P_qg)      = self.CF * ( 2 / (x - 1) ) = self.CF * ( 2 z / (1 - z) )
        # P_qg-CxS(P_qg) = self.CF * (1 + z**2 - 2*z) / (1 - z) = self.CF * ( 1 - z)

        norm = initial_state_crossing_factor * self.CF
        evaluation['values'][(0, 0)]['finite'] = norm * (1 - z)

        return evaluation

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):
        """Add the distributed partial fractioned soft eikonal approximation
        to this hard collinear current
        """

        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        parent = leg_numbers_map.inv[frozenset(children)]
        pC = higher_PS_point[children[0]]
        pC -= sum(higher_PS_point[child] for child in children[1:])
        if self.is_cut(Q=Q, pC=pC):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        evaluation = self.evaluate_kernel(zs, kTs, parent)

        # Find all colored leg numbers except for the parent in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        color_correlation_index = 1
        ps = higher_PS_point[children[1]]
        pi = higher_PS_point[children[0]]
        # pi = lower_PS_point[parent]

        # Loop over the colored parton number pairs (parent, j)
        # and add the corresponding contributions to this current
        for j in all_colored_parton_numbers:
            # Write the eikonal for that pair
            # (positive here since the dipole end 'children[0]' is in the initial state)
            if j == parent:
                continue
            pj = higher_PS_point[j]
            # pj = sum(higher_PS_point[child] for child in leg_numbers_map[j])
            # pj = lower_PS_point[j]
            eik1 = mod_eikonal(pi, pj, ps)
            evaluation['color_correlations'].append(((parent, j),))
            evaluation['values'][(0, color_correlation_index)] = {'finite': eik1}
            color_correlation_index += 1

        # Add the normalization factors
        pC2 = pC.square()
        norm = 8. * math.pi * alpha_s / pC2
        norm *= self.factor(Q=Q, pC=pC)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

class QCD_initial_collinear_0_gg(currents.QCDLocalCollinearCurrent):
    """gg collinear ISR tree-level current. g(initial) > g(initial_after_emission) g(final)"""

    variables = staticmethod(currents.Q_initial_coll_variables)

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD collinear tree-level currents
        init_vars = cls.common_does_implement_this_current(current, 2, 0)
        if init_vars is None: return None
        # Retrieve singular structure
        ss = current.get('singular_structure')

        if len(ss.legs) != 2: return None
        
        n_initial_state_gluons = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and cls.is_initial(leg):
                n_initial_state_gluons += 1
        if n_initial_state_gluons != 1: return None
        
        n_final_state_gluons = 0
        for leg in ss.legs:
            if cls.is_gluon(leg, model) and not cls.is_initial(leg):
                n_final_state_gluons += 1
        if n_final_state_gluons != 1: return None

        # The current is valid
        return init_vars

    @classmethod
    def get_sorted_children(cls, current, model):

        legs = current.get('singular_structure').legs
        # Always put the initial state child first
        children_numbers = [leg.n for leg in legs if leg.state == leg.INITIAL]
        # Then the final state ones
        children_numbers.extend([leg.n for leg in legs if leg.state == leg.FINAL])

        return tuple(children_numbers)

    def evaluate_kernel(self, xs, kTs, parent):

        # Retrieve the collinear variable x
        x = xs[0]
        kT = kTs[0]

        # Instantiate the structure of the result
        evaluation = utils.SubtractionCurrentEvaluation({
            'spin_correlations'  : [None, ((parent,( kT, )), ), ],
            'color_correlations' : [None],
            'values'             : {(0, 0): {'finite': None},
                                    (1, 0): {'finite': None}, }
        })

        # The factor 'x' that should be part of the initial_state_crossing_factor cancels
        # against the extra prefactor 1/x in the collinear factorisation formula
        # (see Eq. (8) of NNLO compatible NLO scheme publication arXiv:0903.1218v2)
        initial_state_crossing_factor = 1.
        # Correct for the ratio of color-averaging factor between the real ME
        # initial state flavor (gluon) and the one of the reduced Born ME (gluon)
        initial_state_crossing_factor *= 1.
        
        z = 1./x
    
        # The line below implements the g_{\mu\nu} part of the splitting kernel.
        # Notice that the extra longitudinal terms included in the spin-correlation 'None'
        # from the relation:
        #    \sum_\lambda \epsilon_\lambda^\mu \epsilon_\lambda^{\star\nu}
        #    = g^{\mu\nu} + longitudinal terms
        # are irrelevant because Ward identities evaluate them to zero anyway.

        # We re-use here the Altarelli-Parisi Kernel of the P_qg final state kernel, including
        # its soft subtraction
        # We must subtract the soft-collinear (CxS *not* SxC) from this contribution:
        # P_gg           = 2.*self.CA * ( (z/(1.-z)) + ((1.-z)/z) )
        # CxS(P_gg)      = 2.*self.CA * ( (z/(1.-z)) )
        # P_gg-CxS(P_gg) = 2.*self.CA * ((1.-z)/z)

        norm = initial_state_crossing_factor * 2. * self.CA
        evaluation['values'][(0, 0)]['finite'] =  norm * ((1.-z)/z)
        evaluation['values'][(1, 0)]['finite'] = -norm * 2.*z*(1.-z) / kT.square()
        return evaluation

    def evaluate_subtraction_current(
        self, current,
        higher_PS_point=None, lower_PS_point=None,
        leg_numbers_map=None, reduced_process=None, hel_config=None,
        Q=None, **opts ):
        """Add the distributed partial fractioned soft eikonal approximation
        to this hard collinear current
        """

        if higher_PS_point is None or lower_PS_point is None:
            raise CurrentImplementationError(
                self.name() + " needs the phase-space points before and after mapping." )
        if leg_numbers_map is None:
            raise CurrentImplementationError(
                self.name() + " requires a leg numbers map, i.e. a momentum dictionary." )
        if reduced_process is None:
            raise CurrentImplementationError(
                self.name() + " requires a reduced_process.")
        if not hel_config is None:
            raise CurrentImplementationError(
                self.name() + " does not support helicity assignment." )
        if Q is None:
            raise CurrentImplementationError(
                self.name() + " requires the total mapping momentum Q." )

        # Retrieve alpha_s and mu_r
        model_param_dict = self.model.get('parameter_dict')
        alpha_s = model_param_dict['aS']
        mu_r = model_param_dict['MU_R']

        # Include the counterterm only in a part of the phase space
        children = self.get_sorted_children(current, self.model)
        parent = leg_numbers_map.inv[frozenset(children)]
        pC = higher_PS_point[children[0]]
        pC -= sum(higher_PS_point[child] for child in children[1:])
        if self.is_cut(Q=Q, pC=pC):
            return utils.SubtractionCurrentResult.zero(
                current=current, hel_config=hel_config)

        # Evaluate collinear subtracted kernel
        zs, kTs = self.variables(higher_PS_point, lower_PS_point[parent], children, Q=Q)
        evaluation = self.evaluate_kernel(zs, kTs, parent)

        # Find all colored leg numbers except for the parent in the reduced process
        all_colored_parton_numbers = []
        for leg in reduced_process.get('legs'):
            if self.model.get_particle(leg.get('id')).get('color') == 1:
                continue
            all_colored_parton_numbers.append(leg.get('number'))

        color_correlation_index = 1
        ps = higher_PS_point[children[1]]
        pi = higher_PS_point[children[0]]
        # pi = lower_PS_point[parent]

        # Loop over the colored parton number pairs (parent, j)
        # and add the corresponding contributions to this current
        for j in all_colored_parton_numbers:
            # Write the eikonal for that pair
            # (positive here since the dipole end 'children[0]' is in the initial state)
            if j == parent:
                continue
            pj = higher_PS_point[j]
            # pj = sum(higher_PS_point[child] for child in leg_numbers_map[j])
            # pj = lower_PS_point[j]
            eik1 = mod_eikonal(pi, pj, ps)
            evaluation['color_correlations'].append(((parent, j),))
            evaluation['values'][(0, color_correlation_index)] = {'finite': eik1}
            color_correlation_index += 1

        # Add the normalization factors
        pC2 = pC.square()
        norm = 8. * math.pi * alpha_s / pC2
        norm *= self.factor(Q=Q, pC=pC)
        for k in evaluation['values']:
            evaluation['values'][k]['finite'] *= norm

        # Construct and return result
        result = utils.SubtractionCurrentResult()
        result.add_result(
            evaluation,
            hel_config=hel_config,
            squared_orders=tuple(sorted(current.get('squared_orders').items())))
        return result

#=========================================================================================
# NLO soft current
#=========================================================================================

class NoSoftCurrent(currents.QCDCurrent):
    """Trivial current returning zero for any NLO limit containing softs."""

    is_zero = True

    @classmethod
    def does_implement_this_current(cls, current, model):

        # Check the general properties common to NLO QCD soft tree-level currents
        init_vars = currents.QCDLocalSoftCurrent.\
            common_does_implement_this_current(current, 2, 0)
        if init_vars is not None:
            return init_vars
        init_vars = currents.QCDLocalSoftCollinearCurrent.\
            common_does_implement_this_current(current, 2, 0)
        return init_vars

    def evaluate_subtraction_current(
        self, current, hel_config=None, **opts ):
        """Return 0 for this current."""

        return utils.SubtractionCurrentResult.zero(current=current, hel_config=hel_config)
