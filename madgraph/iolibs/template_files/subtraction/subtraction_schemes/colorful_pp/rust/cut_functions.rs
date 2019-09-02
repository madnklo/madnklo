def generalised_cuts(cut_inputs, global_variables):
    """ Function applying the correct cut for each bundles depending on the variables passed for each which can be:
        {'pA': ..., 'pC': ....} for the initial-state momenta -pA (if any, otherwise absent) and the final-state ones pC for collines
        {'pS': ...} for the soft
    """
    Q = global_variables['Q']

    for bundle_info in cut_inputs['bundles_info']:
        bundle_cut_info = bundle_info['cut_inputs']
        if 'pS' in bundle_cut_info:
            if soft_cut(pS=bundle_cut_info['pS'], Q=Q):
                return True
        elif 'pA' in bundle_cut_info:
            if initial_coll_cut(pA=bundle_cut_info['pA'], pR=bundle_cut_info['pC'], Q=Q):
                return True
        elif 'pC' in bundle_cut_info:
            if final_coll_cut(pC=bundle_cut_info['pC'], Q=Q):
                return True
        else:
            raise CurrentImplementationError("The generalised_cuts function is only applicable for mapping structure "+
                                                                                              " are collinear or soft.")
    return False


def cut_coll(**opts):

    try:
        alpha = opts['alpha']
    except KeyError:
        pC    = opts['pC']
        Q     = opts['Q']
        alpha = mappings.FinalRescalingOneMapping.alpha(pC, Q)
    # Include the counterterm only up to alpha_0
    return alpha > alpha_0

def no_cut(**opts):
    return False


def cut_initial_coll(**opts):

    pA    = opts['pA']
    pR    = opts['pR']
    Q     = opts['Q']
    y_0p  = (2.*pA.dot(pR))/Q.square()
    # Include the counterterm only up to y_0_prime
    return y_0p > y_0_prime

def cut_soft(**opts):

    try:
        y  = opts['y']
    except KeyError:
        pS = opts['pS']
        Q  = opts['Q']
        y = mappings.SoftVsFinalPureRescalingMapping.y(pS, Q)
    # Include the counterterm only up to y_0
    return y > y_0
