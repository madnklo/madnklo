include ../../Source/make_opts
SHELL = /bin/bash

DEBUG=
#DEBUG= -ggdb -pg

LIBDIR=../../lib/

# TODO: compile libraries when generating process
LINKLIBS= -L../../lib/ -ldhelas -lmodel -lpdf -lgeneric

LIBS=$(LIBDIR)libdhelas.$(libext) $(LIBDIR)libmodel.$(libext) $(LIBDIR)libgeneric.$(libext) \
	$(LIBDIR)libpdf.$(libext)

#$(shell echo LIBS)

OPT=O3

OBJ=obj-$(DEFAULT_F_COMPILER)

FFLAGS+= -g -fno-automatic $(DEBUG)  -$(OPT)

FDEBUG=-fcheck=all -fbacktrace -Wall -ffpe-trap=zero,overflow,underflow

PATH_TO_COMMON_FILES=../../../../Template/Fortran_tmp/src_to_common


INCLUDE0=$(PWD)/include

INCLUDE =$(wildcard ./include/*.inc)

FFLAGS+= -I$(INCLUDE0)

%.o: %.f $(INCLUDE)
	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@


%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE)
	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) $< -o $@


COMMON_FILES= alphaS.o imap.o kinematics.o phase_space.o sectors.o vegas.o 


FILES= driver_3_5.o NLO_Rsub_3_5.o NLO_K_3_5.o NLO_IR_limits_3_5.o matrix_epem_gddx.o \
	matrix_epem_ddx.o $(COMMON_FILES)


sector_3_5: $(FILES)
	$(DEFAULT_F_COMPILER) -o $@ $^ $(LIBS)


.PHONY: clean
clean:
	rm -rf *.o

