# Optional loader and F2PY targets; included only by provider-enabled exports.
ifndef MADNKLO_PYAMP_MK
MADNKLO_PYAMP_MK := 1
ifeq ($(.DEFAULT_GOAL),)
.DEFAULT_GOAL := all
endif
MENUM ?= 2
F2PY ?= $(DEFAULT_F2PY_COMPILER)
COMMON_FILES += pyamplicol_loader.o
PYAMP_CXXLIB = $(if $(filter Darwin,$(OS)),-lc++,-lstdc++)
STDLIB := $(STDLIB) $(filter-out $(STDLIB),$(PYAMP_CXXLIB) -ldl)
LINKLIBS := $(LINKLIBS) $(filter-out $(LINKLIBS),$(PYAMP_CXXLIB) -ldl)

pyamplicol_f2py_loader.o: ../../../Common_Files/pyamplicol_loader.cc
	$(DEFAULT_CPP_COMPILER) -std=c++11 -fPIC -c $< -o $@

ifneq ($(wildcard ../../VirtualProvider/*.json),)
PYAMP_SOURCE = $(wildcard V_user_access_subroutines.f RV_user_access_subroutines.f)
PYAMP_WRAPPER = $(wildcard V_f2py_wrapper.f RV_f2py_wrapper.f)
PYAMP_OBJECTS = $(filter-out $(notdir $(PYAMP_SOURCE:.f=.o)) born_matrix.o,$(PROCESS) $(PYAMP_VIRTUAL_RENAMED_CORE))

matrix$(MENUM)py.so: libs $(PROCESS) $(PYAMP_VIRTUAL_RENAMED_CORE) pyamplicol_f2py_loader.o
	$(MAKE) -C ../../Source/MODEL
	$(F2PY) -c $(PYAMP_SOURCE) born_matrix.f $(PYAMP_WRAPPER) \
	  $(addprefix $(OBJ)/,$(PYAMP_OBJECTS)) pyamplicol_f2py_loader.o \
	  -m matrix$(MENUM)py --f77flags="-std=legacy -ffixed-line-length-132" \
	  --f90flags="-std=legacy -ffixed-line-length-132" \
	  -I. -Iinclude -I../../Source -L../../lib -ldhelas -lmodel $(LINK_LOOP_LIBS) $(STDLIB)
else
matrix$(MENUM)py.so: $(wildcard matrix*.f) $(wildcard *.inc) pyamplicol_f2py_loader.o
	$(MAKE) -C ../../Source/DHELAS
	$(MAKE) -C ../../Source/MODEL
	$(F2PY) -c $(wildcard matrix*.f) pyamplicol_f2py_loader.o -m matrix$(MENUM)py \
	  --f77flags="-std=legacy -ffixed-line-length-132" \
	  --f90flags="-std=legacy -ffixed-line-length-132" \
	  -I. -Iinclude -I../../Source -L../../lib -ldhelas -lmodel $(STDLIB)
endif
endif
