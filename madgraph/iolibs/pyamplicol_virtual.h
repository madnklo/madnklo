#ifndef MADNKLO_PYAMPLICOL_VIRTUAL_H
#define MADNKLO_PYAMPLICOL_VIRTUAL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mnk_pyamp_virtual_handle mnk_pyamp_virtual_handle;
typedef struct {
    double born, finite, single_pole, double_pole;
    /* SDK status: 0 successful/unverified, 1 verified, 2 rescued. */
    uint32_t point_status;
} mnk_pyamp_virtual_value;

typedef struct { double real, imag; } mnk_pyamp_complex;
typedef struct {
    mnk_pyamp_complex born, finite, single_pole, double_pole;
    uint32_t point_status;
} mnk_pyamp_virtual_correlated_value;

/* Exclusive, non-concurrent handles. JSON arguments are contents, not paths.
 * Only inclusive helicity sums with one Born and one virtual order are accepted.
 * Values retain the SDK's physical couplings, averages and scheme conventions;
 * metadata must be checked before using them as MadNkLO virtuals.
 */
int mnk_pyamp_virtual_open(const char *path, const char *process,
    const char *parameters_json, const char *verification_json,
    mnk_pyamp_virtual_handle **output);
void mnk_pyamp_virtual_free(mnk_pyamp_virtual_handle *handle);
int mnk_pyamp_virtual_update_parameters(mnk_pyamp_virtual_handle *handle,
    const char *changes_json);
/* p is [leg,(E,px,py,pz)], with physical incoming momenta. mu_r is a scale,
 * NOT its square. Failed calls leave *output untouched and never use a fallback.
 * NULL policy inherits the policy supplied at open (no retries by default).
 */
int mnk_pyamp_virtual_evaluate(mnk_pyamp_virtual_handle *handle,
    const double *p, size_t momentum_count, double mu_r,
    const char *verification_json, mnk_pyamp_virtual_value *output);
/* JSON size query: output=NULL, capacity=0; needed includes the final NUL. */
int mnk_pyamp_virtual_metadata_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed);
int mnk_pyamp_virtual_diagnostics_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed);
int mnk_pyamp_virtual_histories_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed);
int mnk_pyamp_virtual_spin_classes_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed);
int mnk_pyamp_virtual_history_count(const mnk_pyamp_virtual_handle *handle, size_t *count);
/* Every saved history, in catalogue order. Both interference directions and
 * complex components are retained. Sources are literal [selected_leg,component]
 * four-vectors, labelled by one-based physical legs. Zero legs selects the
 * physical helicity sum. Failed calls leave every output element untouched.
 */
int mnk_pyamp_virtual_evaluate_correlated(mnk_pyamp_virtual_handle *handle,
    const double *p, size_t momentum_count, double mu_r,
    const size_t *spin_legs, size_t spin_count, const mnk_pyamp_complex *sources,
    const char *verification_json, mnk_pyamp_virtual_correlated_value *output,
    size_t capacity);
/* Borrowed thread-local text, valid until the next adapter call on this thread.
 * All fallible functions return zero on success and one on failure.
 */
const char *mnk_pyamp_virtual_last_error(void);

/* Generated Fortran dispatch. Parameter order is compiled from the export's
 * manifest. Handles retain model state, but numerical results are not cached.
 * This entry point is non-concurrent and prints an explicit failure diagnostic.
 */
int mnk_pyamp_virtual_dispatch(int n, int ninitial, const int *pdgs,
    const double *p, int helicity, int squared_orders, double mu_r, double accuracy,
    int parameter_count, const double *parameters, double *answer,
    double *found, int *return_code);
int mnk_pyamp_virtual_correlated_dispatch(int n, int ninitial, const int *pdgs,
    const double *p, int helicity, int squared_orders, double mu_r, double accuracy,
    int parameter_count, const double *parameters, int max_vectors,
    const int *spin_counts, const mnk_pyamp_complex *spin_vectors,
    int color_count, double *answer, double *color_answer, double *found, int *return_code);
/* Diagnostics of the latest generated call for the selected process/config.
 * The same size-query convention applies. Generated calls default to binary64
 * with targeted DD/AP70/AP100 retries on failed points. Set
 * MADNKLO_PYAMP_VIRTUAL_RESCUE=off to disable automatic retries; a positive
 * requested accuracy always enables explicit cross-precision verification.
 */
int mnk_pyamp_virtual_dispatch_diagnostics_json(int n, int ninitial, const int *pdgs,
    char *output, size_t capacity, size_t *needed);

#ifdef __cplusplus
}
#endif
#endif
