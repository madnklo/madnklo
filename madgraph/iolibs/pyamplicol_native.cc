// Native binary64 adapter; compiled against the installed correlated Rusticol SDK.
#include <rusticol.h>
#include "pyamplicol_native.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct LegacyTerm { const char *id; double real; double imag; };
#ifdef MADNKLO_LEGACY_HEADER
#include MADNKLO_LEGACY_HEADER
#else
const std::vector<LegacyTerm> *legacy_terms(const std::string &,const std::string &) { return nullptr; }
#endif
void check(int status) {
    if (status == RUSTICOL_STATUS_OK) return;
    char message[4096] = {};
    size_t required = 0;
    rusticol_last_error_message(message, sizeof(message), &required);
    throw std::runtime_error(message);
}
struct Process {
    RusticolRuntimeHandle *handle = nullptr;
    std::vector<double> momenta;
    bool evaluated = false;
    std::string key;
    std::string rescue_path;
    bool rescue = false, rescue_loaded = false;
    ~Process() { if (handle) rusticol_runtime_free(handle); }
};
Process &process(int n, int ni, const int *pdgs) {
    if (n < 3 || ni < 1 || ni >= n) throw std::runtime_error("invalid process dimensions");
    std::string key = "p" + std::to_string(ni);
    for (int i = 0; i < n; ++i) key += "_" + std::to_string(pdgs[i]);
    static std::map<std::string, std::unique_ptr<Process>> processes;
    const char *root = std::getenv("MADNKLO_PYAMP_ROOT");
    if (!root) throw std::runtime_error("MADNKLO_PYAMP_ROOT is unset");
    auto &entry = processes[std::string(root)+"/"+key];
    if (!entry) {
        std::unique_ptr<Process> value(new Process);
        value->key = key;
        const char *precision = std::getenv("MADNKLO_PYAMP_CORRELATION_PRECISION");
        if (precision && std::string(precision) != "binary64" && std::string(precision) != "binary128")
            throw std::runtime_error("unknown correlation precision; choose binary64 or binary128");
        value->rescue = precision && std::string(precision) == "binary128";
#ifndef MADNKLO_NATIVE_RESCUE
        if (value->rescue) throw std::runtime_error("adapter was built without native rescue support");
#endif
        value->rescue_path = std::string(root)+"/"+key+".rescue.so";
        const std::string directory = std::string(root) + "/" + key;
        std::string parameters = std::string(root) + "/" + key + ".parameters.json";
        if (!std::ifstream(parameters)) parameters = std::string(root) + "/parameters.json";
        check(rusticol_runtime_load(directory.c_str(), nullptr, parameters.c_str(), &value->handle));
        entry = std::move(value);
    }
    return *entry;
}
void evaluate_requests(Process &proc,const double *p,int n,
    const RusticolCorrelatedRequest *requests,size_t count,double *output,size_t capacity) {
#ifdef MADNKLO_NATIVE_RESCUE
    if (proc.rescue) {
        if (!proc.rescue_loaded) {
            check(rusticol_runtime_load_correlated_rescue(proc.handle,proc.rescue_path.c_str()));
            proc.rescue_loaded = true;
        }
        check(rusticol_runtime_evaluate_correlated_many_rescued_f64(proc.handle,p,4*n,1,
            requests,count,nullptr,0,output,capacity));
        return;
    }
#endif
    check(rusticol_runtime_evaluate_correlated_many_f64(proc.handle,p,4*n,1,
        requests,count,nullptr,0,output,capacity));
}
void evaluate_complex(Process &proc, const double *p, int n, const std::string &id,
                      double *output, const RusticolSpinCorrelationVector *spin = nullptr,
                      size_t nspin = 0) {
    if (id.compare(0,7,"legacy/") == 0) {
        const auto *terms = legacy_terms(proc.key,id);
        if (!terms) throw std::runtime_error("unknown legacy colour ID; prepare with --legacy-order");
        std::vector<RusticolCorrelatedRequest> requests;
        requests.reserve(terms->size());
        for (const auto &term : *terms) requests.push_back({term.id,spin,nspin,0});
        std::vector<double> values(2*requests.size());
        if (!requests.empty())
            evaluate_requests(proc,p,n,requests.data(),requests.size(),values.data(),values.size());
        output[0] = output[1] = 0;
        for (size_t i = 0; i < terms->size(); ++i) {
            const auto &term = (*terms)[i];
            const double *value = values.data()+2*i;
            output[0] += term.real*value[0]-term.imag*value[1];
            output[1] += term.real*value[1]+term.imag*value[0];
        }
        if (!std::isfinite(output[0]) || !std::isfinite(output[1]))
            throw std::runtime_error("nonfinite legacy correlated result");
        return;
    }
    RusticolCorrelatedRequest request = {id.c_str(), spin, nspin, 0};
    evaluate_requests(proc,p,n,&request,1,output,2);
    if (!std::isfinite(output[0]) || !std::isfinite(output[1]))
        throw std::runtime_error("nonfinite correlated result");
}
double evaluate(Process &proc, const double *p, int n, const std::string &id,
                const RusticolSpinCorrelationVector *spin = nullptr) {
    double output[2] = {};
    evaluate_complex(proc, p, n, id, output, spin, spin ? 1 : 0);
    if (std::abs(output[1]) > 1.e-10 * std::max(std::abs(output[0]), 1.e-100))
        throw std::runtime_error("non-real or nonfinite Born result");
    return output[0];
}
void set_alpha(Process &proc, double as) {
    if (!(as > 0) || !std::isfinite(as)) throw std::runtime_error("invalid alpha_s");
    const char *parameter = std::getenv("MADNKLO_PYAMP_ALPHAS_PARAMETER");
    check(rusticol_runtime_set_model_parameter(proc.handle, parameter ? parameter : "aS", as, 0));
}
void validate_spin(int n, const int *pdgs, const double *p, int leg, const double *vector) {
    if (leg < 1 || leg > n || pdgs[leg-1] != 21)
        throw std::runtime_error("spin getter requires a declared gluon leg");
    const double *momentum = p + 4*(leg-1);
    double dot[2] = {}, pscale = 0, vscale = 0;
    for (int i = 0; i < 4; ++i) {
        if (!std::isfinite(momentum[i])) throw std::runtime_error("nonfinite spin momentum");
        pscale = std::max(pscale, std::abs(momentum[i]));
        for (int part = 0; part < 2; ++part) {
            const double value = vector[2*i+part];
            if (!std::isfinite(value)) throw std::runtime_error("nonfinite spin vector");
            dot[part] += (i == 0 ? 1 : -1)*momentum[i]*value;
            vscale = std::max(vscale, std::abs(value));
        }
    }
    if (std::max(std::abs(dot[0]),std::abs(dot[1])) > 1.e-10*pscale*vscale)
        throw std::runtime_error("spin getter requires a transverse vector (p.KT = 0)");
}
std::vector<RusticolSpinCorrelationVector> spin_vectors(int n, const int *pdgs,
    const double *p, int count, const int *legs, const double *vectors) {
    if (count < 0 || count > n || (count && (!legs || !vectors)))
        throw std::runtime_error("invalid spin vector array");
    std::vector<RusticolSpinCorrelationVector> spins;
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < i; ++j)
            if (legs[j] == legs[i]) throw std::runtime_error("duplicate spin leg");
        validate_spin(n,pdgs,p,legs[i],vectors+8*i);
        spins.push_back({size_t(legs[i]),1,vectors+8*i,8});
    }
    return spins;
}
template<class Function> void explicit_point(Process &proc, double as, Function function) {
    double previous = 0, imaginary = 0;
    if (proc.evaluated) {
        const char *parameter = std::getenv("MADNKLO_PYAMP_ALPHAS_PARAMETER");
        check(rusticol_runtime_get_model_parameter(proc.handle, parameter ? parameter : "aS",
                                                   &previous, &imaginary));
    }
    set_alpha(proc,as);
    try { function(); }
    catch (...) { if (proc.evaluated) set_alpha(proc,previous); throw; }
    if (proc.evaluated) set_alpha(proc,previous);
}
template<class Function> int guarded(Function function) {
    try { function(); return 0; }
    catch (const std::exception &error) {
        std::fprintf(stderr, "pyAmpliCol provider: %s\n", error.what());
        return 1;
    }
}
}

extern "C" int mnk_pyamp_eval(int n, int ni, const int *pdgs, const double *p,
                              double as, int hel, int nso, double *ans) {
    return guarded([&]() {
        if (hel != -1 || nso < 0 || nso > 1)
            throw std::runtime_error("only helicity sums and one squared coupling order are supported");
        if (!(as > 0) || !std::isfinite(as)) throw std::runtime_error("invalid alpha_s");
        Process &proc = process(n, ni, pdgs);
        set_alpha(proc, as);
        double value = 0;
        check(rusticol_runtime_evaluate_f64(proc.handle,p,4*n,1,&value,1));
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite Born result");
        for (int i = 0; i <= nso; ++i) ans[i] = value;
        proc.momenta.assign(p, p + 4*n);
        proc.evaluated = true;
    });
}
extern "C" int mnk_pyamp_correlated_spin(int n, int ni, const int *pdgs, const double *p,
    double as, const char *id, int nspin, const int *legs, const double *vectors, double *output) {
    return guarded([&]() {
        if (!id || !*id || !output) throw std::runtime_error("missing correlation ID or output");
        auto spins = spin_vectors(n,pdgs,p,nspin,legs,vectors);
        Process &proc = process(n, ni, pdgs);
        // Explicit-point requests must not change the legacy getters' cached state.
        explicit_point(proc,as,[&]() {
            evaluate_complex(proc, p, n, id, output, spins.data(), spins.size());
        });
    });
}
extern "C" int mnk_pyamp_legacy_many(int n, int ni, const int *pdgs, const double *p,
    double as, int count, const int *a, const int *b, int nspin,
    const int *legs, const double *vectors, double *output) {
    return guarded([&]() {
        if (count < 0 || (count && (!a || !b || !output)))
            throw std::runtime_error("invalid legacy correlation array");
        auto spins = spin_vectors(n,pdgs,p,nspin,legs,vectors);
        Process &proc = process(n,ni,pdgs);
        std::vector<const std::vector<LegacyTerm> *> aliases;
        std::map<std::string,size_t> indices;
        std::vector<RusticolCorrelatedRequest> requests;
        for (int i = 0; i < count; ++i) {
            const auto *terms = legacy_terms(proc.key,"legacy/"+std::to_string(a[i])+"/"+std::to_string(b[i]));
            if (!terms) throw std::runtime_error("unknown legacy colour ID; prepare with --legacy-order");
            aliases.push_back(terms);
            for (const auto &term : *terms) {
                if (indices.emplace(term.id,requests.size()).second)
                    requests.push_back({term.id,spins.data(),spins.size(),0});
            }
        }
        std::vector<double> values(2*requests.size());
        explicit_point(proc,as,[&]() {
            if (!requests.empty())
                evaluate_requests(proc,p,n,requests.data(),requests.size(),values.data(),values.size());
        });
        std::vector<double> result(2*count,0.);
        for (int i = 0; i < count; ++i) {
            for (const auto &term : *aliases[i]) {
                const size_t index = 2*indices.at(term.id);
                result[2*i] += term.real*values[index]-term.imag*values[index+1];
                result[2*i+1] += term.real*values[index+1]+term.imag*values[index];
            }
            if (!std::isfinite(result[2*i]) || !std::isfinite(result[2*i+1]))
                throw std::runtime_error("nonfinite legacy correlated result");
        }
        std::copy(result.begin(),result.end(),output);
    });
}
extern "C" int mnk_pyamp_correlated_many(int n, int ni, const int *pdgs, const double *p,
    double as, int count, const char *const *ids, int nspin,
    const int *legs, const double *vectors, double *output) {
    return guarded([&]() {
        if (count < 0 || (count && (!ids || !output)))
            throw std::runtime_error("invalid named correlation array");
        auto spins = spin_vectors(n,pdgs,p,nspin,legs,vectors);
        std::vector<RusticolCorrelatedRequest> requests;
        for (int i=0;i<count;++i) {
            if (!ids[i] || !*ids[i] || std::string(ids[i]).compare(0,7,"legacy/")==0)
                throw std::runtime_error("named batch requires physical colour IDs");
            requests.push_back({ids[i],spins.data(),spins.size(),0});
        }
        Process &proc = process(n,ni,pdgs);
        explicit_point(proc,as,[&]() {
            if (count) evaluate_requests(proc,p,n,requests.data(),requests.size(),output,2*count);
        });
        for (int i=0;i<2*count;++i)
            if (!std::isfinite(output[i])) throw std::runtime_error("nonfinite named correlation");
    });
}
extern "C" int mnk_pyamp_correlated(int n, int ni, const int *pdgs, const double *p,
                                     double as, const char *id, double *output) {
    return mnk_pyamp_correlated_spin(n,ni,pdgs,p,as,id,0,nullptr,nullptr,output);
}
extern "C" int mnk_pyamp_color(int n, int ni, const int *pdgs,
                               int a, int b, double *value) {
    return guarded([&]() {
        Process &proc = process(n, ni, pdgs);
        if (!proc.evaluated) throw std::runtime_error("Born hook must precede colour getter");
        if (a < 1 || b < 1 || a > n || b > n) throw std::runtime_error("invalid colour legs");
        *value = evaluate(proc, proc.momenta.data(), n,
            "T_" + std::to_string(a) + "_" + std::to_string(b));
    });
}
extern "C" int mnk_pyamp_spin(int n, int ni, const int *pdgs,
    const double *p, int leg, const double *vector, double *value) {
    return guarded([&]() {
        Process &proc = process(n, ni, pdgs);
        if (!proc.evaluated) throw std::runtime_error("Born hook must precede spin getter");
        double components[8] = {};
        for (int i = 0; i < 4; ++i) components[2*i] = vector[i];
        validate_spin(n,pdgs,p,leg,components);
        RusticolSpinCorrelationVector spin = {size_t(leg), 1, components, 8};
        *value = evaluate(proc, p, n, "born", &spin);
    });
}
