CXX=g++

CXXFLAGS=-O2 -fPIC

all: libfjcore

libfjcore: fjcorefortran.o fjcore.o
	ar cru libfjcore.a fjcorefortran.o  fjcore.o

clean:
	rm -vf libfjcore.a
 
distclean: clean
	rm -f fjcore.o fjcorefortran.o libfjcore.a

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@

