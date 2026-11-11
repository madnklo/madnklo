// Optional native provider loader. No Python or Rust SDK is needed to build it.
#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace {
void *provider_symbol(const char *name) {
    static void *library = NULL;
    if (!library) {
        const char *path = std::getenv("MADNKLO_PYAMP_LIBRARY");
        if (!path || !(library = dlopen(path, RTLD_NOW | RTLD_LOCAL))) {
            std::fprintf(stderr, "pyAmpliCol: set MADNKLO_PYAMP_LIBRARY to the prepared provider library: %s\n",
                         path ? dlerror() : "variable is unset");
            std::exit(1);
        }
    }
    void *symbol = dlsym(library, name);
    if (!symbol) {
        std::fprintf(stderr, "pyAmpliCol: missing provider symbol %s\n", name);
        std::exit(1);
    }
    return symbol;
}
void checked(int status) {
    if (status) std::exit(1);
}
}

extern "C" void mnk_pyamp_eval_(const int *n, const int *ni, const int *pdgs,
    const double *p, const double *as, const int *hel, const int *nso, double *ans) {
    typedef int (*Function)(int,int,const int*,const double*,double,int,int,double*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_eval"))(
        *n,*ni,pdgs,p,*as,*hel,*nso,ans));
}

extern "C" void mnk_pyamp_virtual_eval_(const int *n, const int *ni, const int *pdgs,
    const double *p, const int *hel, const int *nso, const double *scale,
    const double *accuracy, const int *count, const double *parameters,
    double *answer, double *found, int *return_code) {
    typedef int (*Function)(int,int,const int*,const double*,int,int,double,double,
                           int,const double*,double*,double*,int*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_virtual_dispatch"))(
        *n,*ni,pdgs,p,*hel,*nso,*scale,*accuracy,*count,parameters,answer,found,return_code));
}
extern "C" void mnk_pyamp_virtual_correlated_eval_(const int *n, const int *ni, const int *pdgs,
    const double *p, const int *hel, const int *nso, const double *scale,
    const double *accuracy, const int *count, const double *parameters, const int *max_vectors,
    const int *spin_counts, const void *spin_vectors, const int *color_count,
    double *answer, double *color_answer, double *found, int *return_code) {
    typedef int (*Function)(int,int,const int*,const double*,int,int,double,double,int,const double*,
                           int,const int*,const void*,int,double*,double*,double*,int*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_virtual_correlated_dispatch"))(
        *n,*ni,pdgs,p,*hel,*nso,*scale,*accuracy,*count,parameters,*max_vectors,
        spin_counts,spin_vectors,*color_count,answer,color_answer,found,return_code));
}
extern "C" void mnk_pyamp_color_(const int *n, const int *ni, const int *pdgs,
    const int *a, const int *b, double *value) {
    typedef int (*Function)(int,int,const int*,int,int,double*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_color"))(
        *n,*ni,pdgs,*a,*b,value));
}
extern "C" void mnk_pyamp_spin_(const int *n, const int *ni, const int *pdgs,
    const double *p, const int *leg, const double *vector, double *value) {
    typedef int (*Function)(int,int,const int*,const double*,int,const double*,double*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_spin"))(
        *n,*ni,pdgs,p,*leg,vector,value));
}

// gfortran appends the CHARACTER storage length; the explicit length trims padding.
extern "C" void mnk_pyamp_correlated_(const int *n, const int *ni, const int *pdgs,
    const double *p, const double *as, const char *id, const int *length,
    double *output, size_t storage_length) {
    if (*length <= 0 || size_t(*length) > storage_length || *length > 4096) {
        std::fprintf(stderr, "pyAmpliCol: invalid correlation ID length\n");
        std::exit(1);
    }
    const std::string name(id, *length);
    if (name.find('\0') != std::string::npos) {
        std::fprintf(stderr, "pyAmpliCol: embedded NUL in correlation ID\n");
        std::exit(1);
    }
    typedef int (*Function)(int,int,const int*,const double*,double,const char*,double*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_correlated"))(
        *n,*ni,pdgs,p,*as,name.c_str(),output));
}

extern "C" void mnk_pyamp_correlated_spin_(const int *n, const int *ni, const int *pdgs,
    const double *p, const double *as, const char *id, const int *length,
    const int *nspin, const int *legs, const double *vectors, double *output, size_t storage_length) {
    if (*length <= 0 || size_t(*length) > storage_length || *length > 4096) {
        std::fprintf(stderr,"pyAmpliCol: invalid correlation ID length\n");
        std::exit(1);
    }
    const std::string name(id,*length);
    if (name.find('\0') != std::string::npos) std::exit(1);
    typedef int (*Function)(int,int,const int*,const double*,double,const char*,int,const int*,const double*,double*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_correlated_spin"))(
        *n,*ni,pdgs,p,*as,name.c_str(),*nspin,legs,vectors,output));
}

// Legacy connection indices are process-local, not emitted-particle labels.
extern "C" void mnk_pyamp_legacy_(const int *n, const int *ni, const int *pdgs,
    const double *p, const double *as, const int *a, const int *b,
    const int *nspin, const int *legs, const double *vectors, double *output) {
    const std::string name = "legacy/"+std::to_string(*a)+"/"+std::to_string(*b);
    typedef int (*Function)(int,int,const int*,const double*,double,const char*,int,const int*,const double*,double*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_correlated_spin"))(
        *n,*ni,pdgs,p,*as,name.c_str(),*nspin,legs,vectors,output));
}

extern "C" void mnk_pyamp_legacy_many_(const int *n, const int *ni, const int *pdgs,
    const double *p, const double *as, const int *count, const int *a, const int *b,
    const int *nspin, const int *legs, const double *vectors, double *output) {
    typedef int (*Function)(int,int,const int*,const double*,double,int,const int*,const int*,int,const int*,const double*,double*);
    checked(reinterpret_cast<Function>(provider_symbol("mnk_pyamp_legacy_many"))(
        *n,*ni,pdgs,p,*as,*count,a,b,*nspin,legs,vectors,output));
}
