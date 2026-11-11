#ifndef MADNKLO_PYAMP_NATIVE_H
#define MADNKLO_PYAMP_NATIVE_H

#ifdef __cplusplus
extern "C" {
#endif

/* All functions return zero on success. Momenta are [leg][E,px,py,pz].
 * Legs and connection indices are one-based. Values include initial averages
 * and identical-particle factors. Serialize calls into one provider library. */
int mnk_pyamp_eval(int n, int ni, const int *pdgs, const double *p,
                   double alphas, int helicity, int squared_orders, double *values);
int mnk_pyamp_color(int n, int ni, const int *pdgs, int a, int b, double *value);
int mnk_pyamp_spin(int n, int ni, const int *pdgs, const double *p,
                  int leg, const double *vector, double *value);

/* Explicit-point requests need no Born call and preserve its cached state.
 * Output is [real,imaginary]. IDs may be physical or prepared "legacy/A/B".
 * Complex transverse sources use [leg][component][real,imaginary], eight
 * doubles per selected gluon. nspin=0 requests ordinary helicity sums. */
int mnk_pyamp_correlated(int n, int ni, const int *pdgs, const double *p,
                        double alphas, const char *id, double *output);
int mnk_pyamp_correlated_spin(int n, int ni, const int *pdgs, const double *p,
    double alphas, const char *id, int nspin, const int *legs,
    const double *vectors, double *output);
/* Batch named physical requests sharing one explicit point and spin assignment. */
int mnk_pyamp_correlated_many(int n, int ni, const int *pdgs, const double *p,
    double alphas, int count, const char *const *ids, int nspin,
    const int *legs, const double *vectors, double *output);
/* Batch connection pairs; output is [correlation][real,imaginary]. */
int mnk_pyamp_legacy_many(int n, int ni, const int *pdgs, const double *p,
    double alphas, int count, const int *a, const int *b, int nspin,
    const int *legs, const double *vectors, double *output);

#ifdef __cplusplus
}
#endif
#endif
