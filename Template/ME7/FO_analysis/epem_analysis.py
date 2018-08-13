import madgraph.integrator.observables as observables
import pyjet
import numpy as np
from madgraph.integrator.vectors import LorentzVector
from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList

def ps_point_list(PS_point):

    if isinstance(PS_point, LorentzVectorList):
        return PS_point
    if isinstance(PS_point, LorentzVectorDict):
        return PS_point.to_list()
    try:
        return LorentzVectorList(LorentzVector(v) for v in PS_point)
    except:
        raise Exception("Invalid format for the phase space point: " + str(PS_point))

def is_a_parton(pdg):

    return abs(pdg) in range(1, 7) + [21]

def cpar(data_for_observables, *args, **kwargs):
    """C-parameter event shape, formula from eq. (10) of arXiv:1603.08927 [hep-ph]."""

    PS_point = ps_point_list(data_for_observables['PS_point'])
    flavors = data_for_observables['flavors']
    PS_point_FS = PS_point[len(flavors[0]):]

    num_sum = 0.
    den_sum = 0.
    for i in range(len(PS_point_FS)):
        pi = PS_point_FS[i]
        pi_space = pi.space()
        abs_pi = abs(pi_space)
        den_sum += abs_pi
        for j in range(i, len(PS_point_FS)):
            pj = PS_point_FS[j]
            pj_space = pj.space()
            abs_pj = abs(pj_space)
            costhij = pi_space.dot(pj_space) / (abs_pi*abs_pj)
            sin2thij = 1. - costhij ** 2
            # The factor 2 comes from summing for j > i
            num_sum += 2 * abs_pi * abs_pj * sin2thij
    value = 1.5 * num_sum / (den_sum ** 2)
    return ((value, 1), )

def cpar_wgt(data_for_observables, *args, **kwargs):
    """C-parameter event shape, reweighted."""

    fills = cpar(data_for_observables, *args, **kwargs)
    return tuple((v, w*v) for v, w in fills)

def eec(data_for_observables, *args, **kwargs):
    """Energy-energy correlation, formula from eq. (11) of arXiv:1603.08927 [hep-ph]."""

    PS_point = ps_point_list(data_for_observables['PS_point'])
    flavors = data_for_observables['flavors']
    PS_point_FS = PS_point[len(flavors[0]):]
    Q = sum(p for p in PS_point[:len(flavors[0])])
    Q2 = Q.square()

    ee_fills = []
    for i in range(len(PS_point_FS)):
        pi = PS_point_FS[i]
        Ei = pi[0]
        pi_space = pi.space()
        abs_pi = abs(pi_space)
        for j in range(i, len(PS_point_FS)):
            pj = PS_point_FS[j]
            Ej = pj[0]
            pj_space = pj.space()
            abs_pj = abs(pj_space)
            costhij = pi_space.dot(pj_space) / (abs_pi*abs_pj)
            # Combinatorial factor from summing for j >= i
            factor = 2
            if i == j:
                factor = 1
            ee_fills.append((-costhij, factor*Ei*Ej/Q2))

    return ee_fills

def ptj1(data_for_observables, *args, **kwargs):

    PS_point = ps_point_list(data_for_observables['PS_point'])
    flavors = data_for_observables['flavors']

    drjj_cut = .4

    # First identify all partonic jets
    jets_list = []
    for i, p in enumerate(PS_point[len(flavors[0]):]):
        if is_a_parton(flavors[1][i]):
            jets_list.append(tuple(list(p) + [i + len(flavors[0]) + 1, ]))

    # Cluster them with fastjet
    this_event = np.array(jets_list, dtype=np.dtype(
        [('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('id', 'i8')]))
    sequence = pyjet.cluster(this_event, R=drjj_cut, p=-1, ep=True)
    jets = sequence.inclusive_jets()
    all_jets = LorentzVectorList([LorentzVector(
        [a_jet.e, a_jet.px, a_jet.py, a_jet.pz]) for a_jet in jets])

    pt_max = max(j.pt() for j in all_jets)

    return ((pt_max, 1), )

ptj1_hwu = observables.HwUObservable(
    name="Leading jet p_t",
    observable_function=ptj1,
    range=[0,91.2/2],
    nbins=20
)

cpar_hwu = observables.HwUObservable(
    name="C-parameter",
    observable_function=cpar,
    range=[-1.e-16,1],
    nbins=40,
    y_axis_mode='LIN'
)

cpar_wgt_hwu = observables.HwUObservable(
    name="C-parameter, weighted",
    observable_function=cpar_wgt,
    range=[-1.e-16,1],
    nbins=40,
    y_axis_mode='LIN'
)

eec_hwu = observables.HwUObservable(
    name="Energy-energy correlation",
    observable_function=eec,
    range=[-1.-1.e-6,1.+1.e-6],
    nbins=20,
    y_axis_mode='LIN'
)

observable_list = observables.HwUObservableList([
    ptj1_hwu, cpar_hwu, cpar_wgt_hwu, eec_hwu ])
