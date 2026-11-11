#include "pyamplicol_virtual.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <locale>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct VirtualParameter { const char *name; int source; bool active; double fixed; };
struct VirtualTerm { size_t history; double real,imag; };
struct VirtualAlias { int first,second,count; const VirtualTerm *terms; };
struct VirtualDefinition {
    const char *key, *process;
    int parameter_count;
    const VirtualParameter *parameters;
    int alias_count = 0;
    const VirtualAlias *aliases = nullptr;
    size_t history_count = 0, identity_history = 0;
    const char *history_json = nullptr, *spin_json = nullptr;
};
#ifdef MADNKLO_VIRTUAL_HEADER
#include MADNKLO_VIRTUAL_HEADER
#else
const VirtualDefinition *virtual_definition(const std::string &) { return nullptr; }
#endif

struct Context {
    mnk_pyamp_virtual_handle *handle = nullptr;
    std::vector<double> parameters;
    ~Context() { mnk_pyamp_virtual_free(handle); }
};

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void check(int status, Context *context = nullptr) {
    if (!status) return;
    std::string message = mnk_pyamp_virtual_last_error();
    if (context) {
        size_t size = 0;
        if (!mnk_pyamp_virtual_diagnostics_json(context->handle, nullptr, 0, &size)) {
            std::vector<char> data(size);
            if (!mnk_pyamp_virtual_diagnostics_json(context->handle, data.data(), size, &size))
                message += std::string("; diagnostics: ")+data.data();
        }
    }
    throw std::runtime_error(message);
}

std::string decimal(double value) {
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::scientific << std::setprecision(17) << value;
    return text.str();
}

std::string verification_absolute_tolerance() {
    const char *setting=std::getenv("MADNKLO_PYAMP_VIRTUAL_ABSOLUTE_TOLERANCE");
    if (!setting) return "0";
    std::istringstream input(setting);
    input.imbue(std::locale::classic());
    double value=0;
    input >> value;
    require(!input.fail() && (input >> std::ws).eof() && std::isfinite(value) && value >= 0,
            "MADNKLO_PYAMP_VIRTUAL_ABSOLUTE_TOLERANCE must be finite and nonnegative");
    return decimal(value);
}

std::string verification_policy(double accuracy, bool verify) {
    return "{\"higher\":[{\"arithmetic\":\"double-double\",\"precision\":null},"
           "{\"arithmetic\":\"arbitrary\",\"precision\":70},"
           "{\"arithmetic\":\"arbitrary\",\"precision\":100}],"
           "\"relative_tolerance\":\""+decimal(accuracy)+
           "\",\"absolute_tolerance\":\""+(verify ? verification_absolute_tolerance() : "1e-12")+
           "\",\"verify_successful\":"+(verify ? "true}" : "false}");
}

std::string default_policy() {
    const char *setting=std::getenv("MADNKLO_PYAMP_VIRTUAL_RESCUE");
    const std::string mode=setting ? setting : "on";
    require(mode == "on" || mode == "off", "MADNKLO_PYAMP_VIRTUAL_RESCUE must be on or off");
    return mode == "off" ? "" : verification_policy(1.e-8,false);
}

std::string process_key(int n,int ni,const int *pdgs) {
    require(n >= 3 && ni > 0 && ni < n && pdgs, "invalid virtual process labels");
    std::string key="p"+std::to_string(ni);
    for (int i=0; i<n; ++i) key += "_"+std::to_string(pdgs[i]);
    return key;
}

Context &context(const std::string &key, const VirtualDefinition &definition, const std::string &policy) {
    const char *explicit_root = std::getenv("MADNKLO_PYAMP_VIRTUAL_ROOT");
    const char *tree_root = std::getenv("MADNKLO_PYAMP_ROOT");
    require(explicit_root || tree_root, "set MADNKLO_PYAMP_VIRTUAL_ROOT or MADNKLO_PYAMP_ROOT");
    const std::string root = explicit_root ? explicit_root : std::string(tree_root)+"/../VirtualProvider";
    const std::string path = root+"/"+key;
    static std::map<std::pair<std::string,std::string>,std::unique_ptr<Context>> contexts;
    auto &saved = contexts[{path,policy}];
    if (!saved) {
        std::ifstream stream(path+".parameters.json");
        require(bool(stream), "missing prepared virtual parameter snapshot");
        const std::string parameters((std::istreambuf_iterator<char>(stream)), {});
        auto value = std::make_unique<Context>();
        check(mnk_pyamp_virtual_open(path.c_str(), definition.process, parameters.c_str(),
                                    policy.empty() ? nullptr : policy.c_str(), &value->handle));
        for (const auto &contract : {std::make_pair(definition.history_json,mnk_pyamp_virtual_histories_json),
                                   std::make_pair(definition.spin_json,mnk_pyamp_virtual_spin_classes_json)}) {
            if (!contract.first) continue;
            size_t size=0;
            check(contract.second(value->handle,nullptr,0,&size));
            std::vector<char> actual(size);
            check(contract.second(value->handle,actual.data(),actual.size(),&size));
            require(std::string(actual.data()) == contract.first, "virtual correlation catalogue differs from preparation");
        }
        saved = std::move(value);
    }
    return *saved;
}
}

extern "C" int mnk_pyamp_virtual_correlated_dispatch(int n, int ni, const int *pdgs,
    const double *p, int hel, int nso, double scale, double accuracy,
    int count, const double *parameters, int max_vectors, const int *spin_counts,
    const mnk_pyamp_complex *spin_vectors, int ncolor, double *answer, double *color_answer,
    double *found, int *return_code) {
    try {
        require(n >= 3 && ni > 0 && ni < n && pdgs && p && parameters && answer && found && return_code,
                "invalid generated virtual arguments");
        require(hel == -1 && nso == 1, "virtual routing requires helicity sums and one squared order");
        require(std::isfinite(accuracy) && accuracy != 0, "invalid requested virtual accuracy");
        const std::string key=process_key(n,ni,pdgs);
        const auto *definition = virtual_definition(key);
        require(definition, "virtual process was not prepared; regenerate and prepare the export");
        require(count == definition->parameter_count, "generated virtual parameter layout differs from preparation");
        require(ncolor >= 0 && (!ncolor || (color_answer && ncolor == definition->alias_count)),
                "generated virtual colour layout differs from preparation");
        std::vector<size_t> legs;
        for (int leg=0; leg<n; ++leg) {
            const int amount = spin_counts ? spin_counts[leg] : 0;
            require(amount >= 0 && amount <= max_vectors, "invalid virtual spin-vector count");
            if (amount) legs.push_back(leg+1);
        }
        require(legs.empty() || spin_vectors, "missing virtual spin vectors");
        for (int i=0; i<count; ++i) {
            require(std::isfinite(parameters[i]), "nonfinite virtual model parameter");
            const auto &binding = definition->parameters[i];
            require(binding.active || parameters[binding.source] == binding.fixed,
                    "a fixed virtual model parameter changed; regenerate the artifact");
        }
        Context &saved = context(key, *definition, default_policy());
        const std::vector<double> current(parameters,parameters+count);
        if (saved.parameters != current) {
            std::string changes = "{";
            bool first = true;
            for (int i=0; i<count; ++i) {
                const auto &binding = definition->parameters[i];
                if (!binding.active) continue;
                if (!first) changes += ",";
                first = false;
                changes += std::string("\"")+binding.name+"\":[\""+
                           decimal(parameters[binding.source])+"\",\"0\"]";
            }
            changes += "}";
            check(mnk_pyamp_virtual_update_parameters(saved.handle,changes.c_str()));
            saved.parameters = current;
        }
        std::string policy;
        if (accuracy > 0) policy=verification_policy(accuracy,true);
        std::vector<double> rows(4*(ncolor+1));
        uint32_t point_status = 1;
        const char *verification = policy.empty() ? nullptr : policy.c_str();
        if (ncolor || !legs.empty()) {
            size_t histories = 0;
            check(mnk_pyamp_virtual_history_count(saved.handle,&histories));
            require(histories == definition->history_count && definition->identity_history < histories,
                    "virtual history catalogue differs from preparation");
            std::vector<std::complex<double>> total(4*histories);
            std::vector<mnk_pyamp_complex> sources(4*legs.size());
            // Lists on selected legs are incoherent sums; each joint assignment
            // must be evaluated before summing, never by summing source vectors.
            std::function<void(size_t)> evaluate = [&](size_t selected) {
                if (selected < legs.size()) {
                    const size_t leg=legs[selected]-1;
                    for (int vector=0; vector<spin_counts[leg]; ++vector) {
                        for (size_t component=0; component<4; ++component)
                            sources[4*selected+component] = spin_vectors[leg+n*vector+n*max_vectors*component];
                        evaluate(selected+1);
                    }
                    return;
                }
                std::vector<mnk_pyamp_virtual_correlated_value> values(histories);
                check(mnk_pyamp_virtual_evaluate_correlated(saved.handle,p,4*n,scale,
                    legs.data(),legs.size(),sources.data(),verification,values.data(),values.size()), &saved);
                const uint32_t status=values[0].point_status;
                require(accuracy < 0 || status != 0, "SDK did not verify the requested virtual accuracy");
                point_status = (point_status == 0 || status == 0) ? 0 : (point_status == 2 || status == 2 ? 2 : 1);
                for (size_t history=0; history<histories; ++history) {
                    const auto &value=values[history];
                    const mnk_pyamp_complex row[]={value.born,value.finite,value.single_pole,value.double_pole};
                    for (size_t coefficient=0; coefficient<4; ++coefficient)
                        total[4*history+coefficient] += std::complex<double>(row[coefficient].real,row[coefficient].imag);
                }
            };
            evaluate(0);
            for (int row=0; row<=ncolor; ++row) for (size_t coefficient=0; coefficient<4; ++coefficient) {
                std::complex<double> value;
                if (!row) value=total[4*definition->identity_history+coefficient];
                else {
                    const auto &alias=definition->aliases[row-1];
                    for (int term=0; term<alias.count; ++term) {
                        const auto &weight=alias.terms[term];
                        require(weight.history < histories, "invalid prepared virtual colour alias");
                        value += std::complex<double>(weight.real,weight.imag)*total[4*weight.history+coefficient];
                    }
                }
                const double norm=std::abs(total[4*definition->identity_history+coefficient]);
                require(std::isfinite(value.real()) && std::isfinite(value.imag()) &&
                        std::abs(value.imag()) <= 1.e-10*std::max({std::abs(value.real()),norm,1.e-100}),
                        "legacy real virtual projection has a nonzero imaginary component");
                rows[4*row+coefficient]=value.real();
            }
        } else {
            mnk_pyamp_virtual_value value = {};
            check(mnk_pyamp_virtual_evaluate(saved.handle,p,4*n,scale,verification,&value), &saved);
            point_status=value.point_status;
            require(accuracy < 0 || point_status != 0, "SDK did not verify the requested virtual accuracy");
            rows={value.born,value.finite,value.single_pole,value.double_pole};
        }
        for (int order=0; order<=nso; ++order) {
            for (int i=0; i<4; ++i) answer[4*order+i] = rows[i];
            for (int color=0; color<ncolor; ++color) for (int i=0; i<4; ++i)
                color_answer[color+ncolor*i+4*ncolor*order]=rows[4*(color+1)+i];
            found[order] = point_status == 0 ? -1. : (accuracy > 0 ? accuracy : -1.);
        }
        *return_code = point_status == 0 ? 110 : (point_status == 1 ? 210 : 320);
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr,"pyAmpliCol virtual provider: %s\n",error.what());
        return 1;
    } catch (...) {
        std::fprintf(stderr,"pyAmpliCol virtual provider: unknown adapter failure\n");
        return 1;
    }
}

extern "C" int mnk_pyamp_virtual_dispatch(int n, int ni, const int *pdgs,
    const double *p, int hel, int nso, double scale, double accuracy,
    int count, const double *parameters, double *answer, double *found, int *return_code) {
    return mnk_pyamp_virtual_correlated_dispatch(n,ni,pdgs,p,hel,nso,scale,accuracy,count,
        parameters,0,nullptr,nullptr,0,answer,nullptr,found,return_code);
}

extern "C" int mnk_pyamp_virtual_dispatch_diagnostics_json(int n,int ni,const int *pdgs,
    char *output,size_t capacity,size_t *needed) {
    try {
        const std::string key=process_key(n,ni,pdgs);
        const auto *definition=virtual_definition(key);
        require(definition,"virtual process was not prepared");
        auto &saved=context(key,*definition,default_policy());
        return mnk_pyamp_virtual_diagnostics_json(saved.handle,output,capacity,needed);
    } catch (const std::exception &error) {
        std::fprintf(stderr,"pyAmpliCol virtual diagnostics: %s\n",error.what());
        return 1;
    }
}
