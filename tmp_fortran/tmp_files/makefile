# FF=gfortran
# OPT=O3
# FFLAGS=-fPIC -g -fno-automatic -ffixed-line-length-132 -$(OPT)
# FDEBUG=
# FDEBUG=-fcheck=all -fbacktrace -Wall -ffpe-trap=zero
# #FDEBUG=-fcheck=bounds -fbacktrace -ffpe-trap=zero
# CC=g++
# CFLAGS=-fPIC -$(OPT)
# CDEBUG=
include ../../Source/make_opts
SHELL = /bin/bash

LIBDIR=../../lib/

LINKLIBS= -L../../lib/ -ldhelas -lmodel

LIBS=$(LIBDIR)libdhelas.$(libext) $(LIBDIR)libmodel.$(libext)

OPT=O3

OBJ=obj-$(DEFAULT_F_COMPILER)

FFLAGS+= -g -fno-automatic  -$(OPT)

FDEBUG=-fcheck=all -fbacktrace -Wall -ffpe-trap=zero

PATH_TO_COMMON_FILES=../../../../src_to_common

INCLUDE= nexternal.inc colored_partons.inc coupl.inc epem_gddx_spin_correlations.inc \
	math.inc  damping_factors.inc nsqso_born.inc leg_PDGs.inc all_sector_list.inc


#%.o: %.f $(INCLUDE) $(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@

#%.o: %.f $(INCLUDE) $(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@

%.o: %.f $(INCLUDE) $(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $^


#%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE) $(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@


SRC_COMMON_FILES= alphaS.o imap.o kinematics.o phase_space.o sectors.o vegas.o util.o


#FILES= alphaS.o imap.o kinematics.o phase_space.o sectors.o vegas.o

FILES= driver_5_3.o NLO_Rsub_5_3.o NLO_K_5_3.o NLO_IR_limits_5_3.o matrix_epem_gddx.o \
	matrix_epem_dddxdx.o 


sector_5_3: $(FILES) $(SRC_COMMON_FILES)
	$(DEFAULT_F_COMPILER)  -o $@ $^ $(LIBS)
#sector_5_3: $(FILES) $(DEFAULT_F_COMPILER) $(SRC_COMMON_FILES) -o $@ $^ $(LIBS)


.PHONY: clean
clean:
	rm -rf *.o

