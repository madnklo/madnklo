#include "pyamplicol_virtual.h"
#include <rusticol_virtual.h>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

struct mnk_pyamp_virtual_handle {
    rusticol_virtual_handle *runtime = nullptr;
    rusticol_virtual_shape_info axes = {};
    std::string diagnostics;
    ~mnk_pyamp_virtual_handle() { rusticol_virtual_free(runtime); }
};

namespace {
thread_local std::string last_error;

void check(int status) {
    if (status == RUSTICOL_VIRTUAL_OK) return;
    const char *message = rusticol_virtual_last_error_message();
    throw std::runtime_error(message && *message ? message : "virtual SDK call failed");
}

template<class Function> int guarded(Function function) {
    last_error.clear();
    try { function(); return 0; }
    catch (const std::exception &error) { last_error = error.what(); }
    catch (...) { last_error = "unknown native virtual adapter failure"; }
    return 1;
}

void require(bool valid, const char *message) {
    if (!valid) throw std::runtime_error(message);
}

template<class Function> std::string read_json(Function function) {
    size_t needed = 0;
    check(function(nullptr, 0, &needed));
    require(needed > 0, "virtual SDK returned an empty JSON buffer size");
    std::vector<char> buffer(needed);
    check(function(buffer.data(), buffer.size(), &needed));
    require(needed > 0 && needed <= buffer.size() && buffer[needed-1] == '\0',
            "virtual SDK returned invalid JSON storage");
    return std::string(buffer.data(), needed-1);
}

void copy_json(const std::string &text, char *output, size_t capacity, size_t *needed) {
    require(needed, "missing JSON size output");
    *needed = text.size()+1;
    if (!output && capacity == 0) return;
    require(output && capacity >= *needed, "JSON output buffer is too small");
    std::memcpy(output, text.c_str(), *needed);
}

struct ResultDeleter {
    void operator()(rusticol_virtual_result *result) const {
        rusticol_virtual_result_free(result);
    }
};
}

extern "C" const char *mnk_pyamp_virtual_last_error() { return last_error.c_str(); }

extern "C" int mnk_pyamp_virtual_open(const char *path, const char *process,
    const char *parameters_json, const char *verification_json,
    mnk_pyamp_virtual_handle **output) {
    return guarded([&]() {
        require(output, "missing virtual handle output");
        *output = nullptr;
        require(path && *path, "missing virtual artifact path");
        require(rusticol_virtual_abi_version() == 1u, "unsupported virtual SDK ABI");
        auto handle = std::make_unique<mnk_pyamp_virtual_handle>();
        check(rusticol_virtual_open_with_policy(path, process, parameters_json, "model",
                                               verification_json, &handle->runtime));
        check(rusticol_virtual_shape(handle->runtime, &handle->axes));
        const auto &axes = handle->axes;
        require(axes.external_count >= 3 && axes.born_orders == 1 &&
                axes.virtual_orders == 1 && axes.epsilon_count == 3,
                "virtual adapter requires one Born order, one virtual order and three Laurent coefficients");
        *output = handle.release();
    });
}

extern "C" void mnk_pyamp_virtual_free(mnk_pyamp_virtual_handle *handle) { delete handle; }

extern "C" int mnk_pyamp_virtual_update_parameters(mnk_pyamp_virtual_handle *handle,
    const char *changes_json) {
    return guarded([&]() {
        require(handle && changes_json && *changes_json, "missing handle or parameter changes");
        handle->diagnostics.clear();
        uint64_t revision = 0;
        check(rusticol_virtual_update_parameters(handle->runtime, changes_json, &revision));
    });
}

extern "C" int mnk_pyamp_virtual_evaluate(mnk_pyamp_virtual_handle *handle,
    const double *p, size_t momentum_count, double mu_r,
    const char *verification_json, mnk_pyamp_virtual_value *output) {
    return guarded([&]() {
        require(handle, "missing virtual handle");
        handle->diagnostics.clear();
        require(p && output && momentum_count == 4*handle->axes.external_count,
                "invalid virtual momenta or output");
        for (size_t i = 0; i < momentum_count; ++i)
            require(std::isfinite(p[i]), "nonfinite virtual momentum");
        const double mu_r2 = mu_r*mu_r;
        require(mu_r > 0 && std::isfinite(mu_r2) && mu_r2 > 0,
                "virtual renormalization scale must have a positive finite square");
        const uint64_t point_id = 0;
        rusticol_virtual_batch_request request = {};
        request.momenta = p;
        request.momentum_count = momentum_count;
        request.point_count = 1;
        request.point_ids = &point_id;
        request.mu_r2 = &mu_r2;
        request.scale_count = 1;
        request.verification_json = verification_json;
        rusticol_virtual_result *raw = nullptr;
        const int status = rusticol_virtual_evaluate(handle->runtime, &request, &raw);
        std::unique_ptr<rusticol_virtual_result, ResultDeleter> result(raw);
        check(status);
        require(raw, "virtual SDK returned no result");
        handle->diagnostics = read_json([&](char *buffer, size_t capacity, size_t *needed) {
            return rusticol_virtual_result_diagnostics_json(raw, buffer, capacity, needed);
        });
        rusticol_virtual_result_shape_info shape = {};
        check(rusticol_virtual_result_shape(raw, &shape));
        require(shape.requested_points == 1 && shape.resolved == 0,
                "virtual SDK returned an incompatible result layout");
        uint32_t point_status = RUSTICOL_VIRTUAL_FAILED;
        check(rusticol_virtual_result_copy_statuses(raw, &point_status, 1));
        require(point_status <= RUSTICOL_VIRTUAL_RESCUED,
                "virtual point failed or requested accuracy was not established; inspect diagnostics");
        require(shape.successful_points == 1 && shape.born_values == 1 && shape.virtual_values == 3 &&
                shape.axes.external_count == handle->axes.external_count &&
                shape.axes.born_orders == 1 && shape.axes.virtual_orders == 1 && shape.axes.epsilon_count == 3,
                "virtual SDK returned inconsistent successful-point axes");
        uint64_t successful_id = 1;
        check(rusticol_virtual_result_copy_successful_ids(raw, &successful_id, 1));
        require(successful_id == point_id, "virtual SDK returned a different point ID");
        mnk_pyamp_virtual_value value = {};
        double coefficients[3];
        check(rusticol_virtual_result_copy_born(raw, &value.born, 1));
        check(rusticol_virtual_result_copy_virtual(raw, coefficients, 3));
        require(std::isfinite(value.born) && std::isfinite(coefficients[0]) &&
                std::isfinite(coefficients[1]) && std::isfinite(coefficients[2]),
                "virtual SDK returned nonfinite coefficients");
        value.finite = coefficients[0];
        value.single_pole = coefficients[1];
        value.double_pole = coefficients[2];
        value.point_status = point_status;
        *output = value;
    });
}

extern "C" int mnk_pyamp_virtual_metadata_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed) {
    return guarded([&]() {
        require(handle, "missing virtual handle");
        const auto metadata = read_json([&](char *buffer, size_t size, size_t *count) {
            return rusticol_virtual_metadata_json(handle->runtime, buffer, size, count);
        });
        copy_json(metadata, output, capacity, needed);
    });
}

extern "C" int mnk_pyamp_virtual_diagnostics_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed) {
    return guarded([&]() {
        require(handle && !handle->diagnostics.empty(), "no virtual point diagnostics available");
        copy_json(handle->diagnostics, output, capacity, needed);
    });
}

extern "C" int mnk_pyamp_virtual_histories_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed) {
    return guarded([&]() {
        require(handle, "missing virtual handle");
        const auto text = read_json([&](char *buffer, size_t size, size_t *count) {
            return rusticol_virtual_correlation_histories_json(handle->runtime, buffer, size, count);
        });
        copy_json(text, output, capacity, needed);
    });
}

extern "C" int mnk_pyamp_virtual_spin_classes_json(const mnk_pyamp_virtual_handle *handle,
    char *output, size_t capacity, size_t *needed) {
    return guarded([&]() {
        require(handle, "missing virtual handle");
        const auto text = read_json([&](char *buffer, size_t size, size_t *count) {
            return rusticol_virtual_spin_correlation_classes_json(handle->runtime, buffer, size, count);
        });
        copy_json(text, output, capacity, needed);
    });
}

extern "C" int mnk_pyamp_virtual_history_count(const mnk_pyamp_virtual_handle *handle, size_t *count) {
    return guarded([&]() {
        require(handle && count, "missing virtual handle or history count");
        rusticol_virtual_correlated_shape_info axes = {};
        check(rusticol_virtual_correlated_shape(handle->runtime, &axes));
        require(axes.born_orders == 1 && axes.virtual_orders == 1 && axes.epsilon_count == 3,
                "incompatible correlated virtual axes");
        *count = axes.histories;
    });
}

extern "C" int mnk_pyamp_virtual_evaluate_correlated(mnk_pyamp_virtual_handle *handle,
    const double *p, size_t momentum_count, double mu_r,
    const size_t *spin_legs, size_t spin_count, const mnk_pyamp_complex *sources,
    const char *verification_json, mnk_pyamp_virtual_correlated_value *output, size_t capacity) {
    return guarded([&]() {
        require(handle, "missing virtual handle");
        handle->diagnostics.clear();
        require(p && output && momentum_count == 4*handle->axes.external_count,
                "invalid correlated virtual momenta or output");
        for (size_t i=0; i<momentum_count; ++i)
            require(std::isfinite(p[i]), "nonfinite virtual momentum");
        const double mu_r2 = mu_r*mu_r;
        require(mu_r > 0 && std::isfinite(mu_r2) && mu_r2 > 0,
                "virtual renormalization scale must have a positive finite square");
        rusticol_virtual_correlated_shape_info axes = {};
        check(rusticol_virtual_correlated_shape(handle->runtime, &axes));
        require(axes.histories && capacity >= axes.histories && axes.born_orders == 1 &&
                axes.virtual_orders == 1 && axes.epsilon_count == 3,
                "incompatible correlated virtual axes or output capacity");
        require(!spin_count || (spin_legs && sources), "missing virtual spin sources");
        std::vector<rusticol_virtual_complex64> vectors(4*spin_count);
        for (size_t i=0; i<spin_count; ++i) {
            require(spin_legs[i] > 0 && spin_legs[i] <= handle->axes.external_count,
                    "invalid virtual spin leg");
            for (size_t j=0; j<i; ++j) require(spin_legs[i] != spin_legs[j], "duplicate virtual spin leg");
            for (size_t j=0; j<4; ++j) {
                const auto &source = sources[4*i+j];
                require(std::isfinite(source.real) && std::isfinite(source.imag), "nonfinite virtual spin source");
                vectors[4*i+j] = {source.real,source.imag};
            }
        }
        const uint64_t point_id = 0;
        rusticol_virtual_batch_request request = {};
        request.momenta=p; request.momentum_count=momentum_count; request.point_count=1;
        request.point_ids=&point_id; request.mu_r2=&mu_r2; request.scale_count=1;
        request.verification_json=verification_json;
        rusticol_virtual_spin_sources spin = {spin_legs,spin_count,vectors.data(),vectors.size()};
        rusticol_virtual_correlated_result *raw = nullptr;
        const int status = spin_count ? rusticol_virtual_evaluate_correlated_with_sources(handle->runtime,&request,&spin,&raw)
                                     : rusticol_virtual_evaluate_correlated(handle->runtime,&request,&raw);
        std::unique_ptr<rusticol_virtual_correlated_result,
            decltype(&rusticol_virtual_correlated_result_free)> result(raw,rusticol_virtual_correlated_result_free);
        check(status);
        require(raw, "virtual SDK returned no correlated result");
        handle->diagnostics = read_json([&](char *buffer, size_t size, size_t *needed) {
            return rusticol_virtual_correlated_result_diagnostics_json(raw,buffer,size,needed);
        });
        rusticol_virtual_correlated_result_shape_info shape = {};
        check(rusticol_virtual_correlated_result_shape(raw,&shape));
        uint32_t point_status = RUSTICOL_VIRTUAL_FAILED;
        check(rusticol_virtual_correlated_result_copy_statuses(raw,&point_status,1));
        require(point_status <= RUSTICOL_VIRTUAL_RESCUED,
                "correlated virtual point failed or requested accuracy was not established; inspect diagnostics");
        require(shape.requested_points == 1 && shape.successful_points == 1 &&
                shape.axes.histories == axes.histories && shape.axes.born_orders == 1 &&
                shape.axes.virtual_orders == 1 && shape.axes.epsilon_count == 3 &&
                shape.born_values == axes.histories && shape.virtual_values == 3*axes.histories,
                "virtual SDK returned inconsistent correlated-point axes");
        uint64_t successful_id = 1;
        check(rusticol_virtual_correlated_result_copy_successful_ids(raw,&successful_id,1));
        require(successful_id == point_id, "virtual SDK returned a different point ID");
        std::vector<rusticol_virtual_complex64> born(axes.histories), loop(3*axes.histories);
        check(rusticol_virtual_correlated_result_copy_born(raw,born.data(),born.size()));
        check(rusticol_virtual_correlated_result_copy_virtual(raw,loop.data(),loop.size()));
        std::vector<mnk_pyamp_virtual_correlated_value> values(axes.histories);
        for (size_t i=0; i<axes.histories; ++i) {
            const rusticol_virtual_complex64 row[] = {born[i],loop[3*i],loop[3*i+1],loop[3*i+2]};
            mnk_pyamp_complex *destination[] = {&values[i].born,&values[i].finite,
                                              &values[i].single_pole,&values[i].double_pole};
            for (size_t j=0; j<4; ++j) {
                require(std::isfinite(row[j].real) && std::isfinite(row[j].imag),
                        "virtual SDK returned nonfinite correlated coefficients");
                *destination[j] = {row[j].real,row[j].imag};
            }
            values[i].point_status = point_status;
        }
        for (size_t i=0; i<axes.histories; ++i) output[i] = values[i];
    });
}
