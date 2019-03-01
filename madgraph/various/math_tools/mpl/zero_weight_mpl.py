"""Implementation of the weight 0 MPL, which is the constant function returning 1"""


def zero_weight_mpl(entries,x):
    return 1.


def zero_weight_mpl_conditions(entries,x):
    return len(entries) == 0