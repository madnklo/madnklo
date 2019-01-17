"""Subtraction scheme without soft mappings

At NLO, the eikonals are split using partial fractionning a la FKS and shipped into the collinear counterterms that
share the same collinear divergence.

The difference with FKS is
- no reliance on sectors (can be included to improve numerics but not necessary)
- no reliance on an explicit parametrization
- independent of the choice of collinear mapping used as long as it verifies some  basic properties
"""