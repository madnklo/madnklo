include ../../Source/make_opts
SHELL = /bin/bash



OS=$(shell uname -s)

ifeq ("$(OS)","Darwin")
STDLIB= -lc++
endif

PWD=$(shell pwd)

FJPATH=$(shell fastjet-config --prefix)
LIBSC=-L$(FJPATH)/lib \
	-lfastjettools -lfastjet -lfastjetplugins \
	-lsiscone_spherical -lsiscone $(STDLIB)
INC=-I$(FJPATH)/include


CFLAGS=-fPIC -$(OPT)
CDEBUG=

DEBUG=
#DEBUG= -ggdb -pg

LIBDIR=../../lib/

# TODO: compile libraries when generating process
LINKLIBS= -L../../lib/ -ldhelas -lmodel -lpdf -lgeneric

LIBS=$(LIBDIR)libdhelas.$(libext) $(LIBDIR)libmodel.$(libext) $(LIBDIR)libgeneric.$(libext) \
	$(LIBDIR)libpdf.$(libext)

OPT=O3

OBJ=obj-$(DEFAULT_F_COMPILER)

FFLAGS+= -g -fno-automatic $(DEBUG)  -$(OPT)

FDEBUG=-fcheck=all -fbacktrace -Wall -ffpe-trap=zero,overflow,underflow

PATH_TO_COMMON_FILES=../../../../Template/Fortran_tmp/src_to_common

PATH_TO_USR_FILES=../../../Common_Files

INCLUDE0=$(PWD)/include

INCLUDE1=$(PWD)

INCLUDE =$(wildcard ./include/*.inc ./*.inc)

SOURCE= ../../Source

FFLAGS+= -I$(INCLUDE0) -I$(INCLUDE1) -I$(SOURCE)

VPATH = $(OBJ)

%.o: %.f $(INCLUDE)
	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $< 

#%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE)
#	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

%.o: $(PATH_TO_USR_FILES)/%.f $(INCLUDE)
	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

%.o: $(PATH_TO_USR_FILES)/%.cc
	$(DEFAULT_CPP_COMPILER) -c $(CFLAGS) $(CDEBUG) $< -o $(OBJ)/$@ $(INC)

#COMMON_FILES=  gen_phase_space.o vegas.o fjcore.o fastjetfortran_core.o
COMMON_FILES=  genps.o vegas.o fjcore.o fastjetfortran_core.o invarients.o transpole.o myamp.o idenparts.o dsample.o ranmar.o

USR_FILES=alphaS.o kinematics.o LO_B.o analysis_new.o HwU.o cuts.o

PROC_FILES = matrix.o configs.o props.o decayBW.o leshouche.o

FILES = driver_n.o $(PROC_FILES) $(COMMON_FILES) $(USR_FILES)

all: libs Born

Born: $(FILES)
	$(DEFAULT_F_COMPILER) $(patsubst %,$(OBJ)/%,$(FILES)) $(LIBS) $(LIBSC) -o $@ 

.PHONY: clean libs
clean:
	rm -rf $(OBJ)/*.o ../../lib/libdhelas.a
	cd ../../../Source; make clean; cd -

libs:
	mkdir -p $(OBJ)
	cd ../../Source/DHELAS; make; cd -
	cd ../../../Source; make ; cd -

