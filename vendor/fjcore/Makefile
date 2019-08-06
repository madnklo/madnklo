CXX=g++

CXXFLAGS=-O2 -fPIC

all: libfjcore

libfjcore: fjcore.o fastjetfortran_madnklo_core.o
	ar cru libfjcore.a  fjcore.o fastjetfortran_madnklo_core.o

clean:
	rm -vf libfjcore.a
 
distclean: clean
	rm -f fjcore.o fjcorefortran.o libfjcore.a

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@

