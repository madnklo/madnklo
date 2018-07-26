import madgraph.integrator.observables as observables

xsec = observables.HwUObservable(
    name                = "Inclusive cross section",
    observable_function = observables.ObservableFunctions.inclusive_xsec,
    range               = [0.,2.],
    nbins               = 3
)

h_t = observables.HwUObservable(
    name                = "H_T",
    observable_function = observables.ObservableFunctions.scalar_pt_sum,
    range               = [0,1000],
    nbins               = 20
)

observable_list = observables.HwUObservableList([
    xsec,
    h_t
])
