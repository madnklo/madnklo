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

#INCLUDE= nexternal.inc colored_partons.inc coupl.inc epem_gddx_spin_correlations.inc \
	math.inc  damping_factors.inc nsqso_born.inc leg_PDGs.inc all_sector_list.inc

INCLUDE0=$(PWD)/include

INCLUDE =$(wildcard ./include/*.inc)

FFLAGS+= -I$(INCLUDE0)

#%.o: %.f $(INCLUDE) $(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@

#%.o: %.f $(INCLUDE) $(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@

%.o: %.f $(INCLUDE)
	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@


%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE)
	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@


COMMON_FILES= alphaS.o imap.o kinematics.o phase_space.o sectors.o vegas.o util.o


FILES= driver_5_3.o NLO_Rsub_5_3.o NLO_K_5_3.o NLO_IR_limits_5_3.o matrix_epem_gddx.o \
	matrix_epem_dddxdx.o $(COMMON_FILES)


#sector_5_3: $(FILES) $(DEFAULT_F_COMPILER) $(SRC_COMMON_FILES) -o $@ $^ $(LIBS)

sector_5_3: $(FILES)
	$(DEFAULT_F_COMPILER) -o $@ $^ $(LIBS)


.PHONY: clean
clean:
	rm -rf *.o

