def kTdiff(i, j, zs, kTs):
    """Compute the difference vector (kTi/zi - kTj/zj)."""

    return kTs[i-1]/zs[i-1] - kTs[j-1]/zs[j-1]

def sij(i, j, zs, kTs, kTdiffij=None):
    """Compute the invariant mass sij from zs and kTs."""

    if kTdiffij: kTdiffij2 = kTdiffij.square()
    else: kTdiffij2 = kTdiff(i, j, zs, kTs).square()
    return -zs[i-1]*zs[j-1]*kTdiffij2

def tijk(zi, zj, sij, sik, sjk):

    return (2*(zi*sjk-zj*sik)+(zi-zj)*sij)/(zi+zj)

def t(i, j, k, zs, kTs, s_ij=None, s_ik=None, s_jk=None):
    """Compute the invariant t_{ij,k} from zs and kTs."""

    if not s_ij: s_ij = sij(i, j, zs, kTs)
    if not s_ik: s_ik = sij(i, k, zs, kTs)
    if not s_jk: s_jk = sij(j, k, zs, kTs)
    return tijk(zs[i-1],zs[j-1],s_ij,s_ik,s_jk)


def C123_Qqqx_kernel(z1, z2, z3, s12, s13, s23):
    """Triple collinear splitting kernel"""
    s123 = s12 + s13 + s23
    t123 = tijk(z1, z2, s12, s13, s23)
    sqrbrk = -(t123 ** 2) / (s12 * s123)
    sqrbrk += (4 * z3 + (z1 - z2) ** 2) / (z1 + z2)
    sqrbrk += z1 + z2 - s12 / s123
    # return 1. / 2. * CF * TR * s123 * sqrbrk / (s12) old version
    return 1. / 2. * sqrbrk / (s12)/ s123  # now explicitly dividing by the parent mass^2L


def C123C12_Qqqx_kernel(z1, z12, k1perp, k12perp, s12, s_12_3, s123):
    """Nested kernel for qqbar collinear inside a qqbarQ triple collinear"""

    # # Retrieve momenta
    # p1 = higher_PS_point[children[0]]
    # p2 = higher_PS_point[children[1]]
    # p3 = higher_PS_point[children[2]]
    # p12 = p1 + p2
    # p123 = p12 + p3
    # p12hat = interm_PS_point[1000]
    # p3hat = interm_PS_point[children[2]]
    # p123hat = p12hat + p3hat
    # p123tilde = final__PS_point[2000]
    # Q = interm_mapping_vars['Q']
    #
    # # Compute momentum fractions and transverse momenta
    # zs_interm, kTs_interm = self.variables(
    #     higher_PS_point, p12hat, children[:2], Q=Q)
    # zs_final_, kTs_final_ = self.variables(
    #     interm_PS_point, p123tilde, (1000, children[2]), Q=Q)
    # z1 = zs_interm[0]
    # k1perp = kTs_interm[0]
    # z12 = zs_final_[0]
    # k12perp = kTs_final_[0]
    #
    # # Build scalar products
    # s12 = p12.square()
    # s12hat_3hat = 2*p3hat.dot(p12hat)
    k1perp2 = k1perp.square()
    k12perp2 = k12perp.square()
    kperpSP = 2 * k1perp.dot(k12perp)
    # Construct the iterated current C(C(1,2),3)
    perpterm = ((1 - z12) / z12) * (kperpSP ** 2) / (k1perp2 * k12perp2)
    pqg = (1 + (1 - z12) ** 2) / z12
    brk = z12 + perpterm
    C123C12_current = 4 * (pqg - 2 * z1 * (1 - z1) * brk) / (s_12_3 * s12)

    #------------ TODO FIGURE OUT THIS STUFF BELOW
    # # If jacobians are active, correct the one-step jacobian to the two-step one
    # # and correct current factors altogether
    # try:
    #     jacobian = opts['jacobian']
    #     jacobian /= (final__mapping_vars['jacobian'] * interm_mapping_vars['jacobian'])
    # except KeyError:
    #     jacobian = 1
    # # Correct current factors if they are active
    # factor_interm = self.factor(Q=Q, pC=p12, qC=p12hat)
    # factor_final_ = self.factor(Q=Q, pC=p123hat, qC=p123tilde)
    # factor_direct = self.factor(Q=Q, pC=p123, qC=p123tilde)
    # factor = factor_interm * factor_final_ / factor_direct
    # ------------ TODO FIGURE OUT THIS STUFF ABOVE

    return C123C12_current# * jacobian * factor  TODO FIGURE OUT THIS STUFF