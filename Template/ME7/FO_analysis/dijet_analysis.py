import madgraph.integrator.observables as observables
import pyjet
import numpy as np

class my_functions(object):
    @staticmethod
    def ptj1(data_for_observables,*args,**kwargs):
        from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList, LorentzVector
        def is_a_jet(pdg):
            return abs(pdg) in range(1,5)+[21]


        PS_point = data_for_observables['PS_point']
        if isinstance(PS_point, LorentzVectorDict):
            PS_point = PS_point.to_list()
        elif not isinstance(PS_point, LorentzVectorList):
            PS_point = LorentzVectorList(LorentzVector(v) for v in PS_point)
        flavors = data_for_observables['flavors']

        drjj_cut = .4
        ptj_cut = 10

        # First identify all partonic jets
        jets_list = []
        for i, p in enumerate(PS_point[len(flavors[0]):]):
            if is_a_jet(flavors[1][i]):
                jets_list.append(tuple(list(p) + [i + len(flavors[0]) + 1, ]))
        # Count partonic jets
        starting_n_jets = len(jets_list)

        # Cluster them with fastjet
        this_event = np.array(jets_list, dtype=np.dtype(
            [('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('id', 'i8')]))
        sequence = pyjet.cluster(this_event, R=drjj_cut, p=-1, ep=True)
        jets = sequence.inclusive_jets()


        if ptj_cut > 0.:
            jets = [jet for jet in jets if jet.pt >= ptj_cut]

        all_jets = LorentzVectorList([LorentzVector(
            [a_jet.e, a_jet.px, a_jet.py, a_jet.pz]) for a_jet in jets])

        return max(j.pt() for j in all_jets)

ptj1 = observables.HwUObservable(
    name="Leading jet p_t",
    observable_function=my_functions.ptj1,
    range=[0,1000],
    nbins=20
)

observable_list = observables.HwUObservableList([ptj1])