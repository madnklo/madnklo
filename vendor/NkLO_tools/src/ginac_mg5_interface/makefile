include ../../make_opts

LIBDIR = ../../lib/

ifeq ($(origin GINAC_PATH),undefined)
GINAC_PATH =
endif


ifeq ($(origin CLN_PATH),undefined)
CLN_PATH   =	
endif

CXXFLAGS   += -fPIC --std=c++11

LDFLAGS += -shared
ifdef CLN_PATH
 LDFLAGS += -L$(CLN_PATH)
endif
LDFLAGS += -lcln

ifdef GINAC_PATH
 LDFLAGS += -L$(GINAC_PATH)
endif
LDFLAGS += -lginac

# Files for GINACINTERFACE
GINACINTERFACE_SRC = ginac_mg5_interface.cpp
GINACINTERFACE= ginac_mg5_interface.o
LIBGINACINTERFACE = libginac_mg5_interface.a
LIBGINACINTERFACE_DYNAMIC = libginac_mg5_interface.so

all: $(LIBDIR)$(LIBGINACINTERFACE) $(LIBDIR)$(LIBGINACINTERFACE_DYNAMIC)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c  $<

# Static libraries
$(LIBDIR)$(LIBGINACINTERFACE): $(GINACINTERFACE)
	$(call CREATESTATCLIB, $@, $^)
	@echo "created the static library for ginac_mg5_interface"

# Dynamic libraries
$(LIBDIR)$(LIBGINACINTERFACE_DYNAMIC): $(GINACINTERFACE_SRC)
	$(call CREATEDYNAMICLIB, $@, $^)
	@echo "created the dynamic library for ginac_mg5_interface"

clean:
	@rm -f *.o *.so $(LIBDIR)$(LIBGINACINTERFACE) $(LIBDIR)$(LIBGINACINTERFACE_DYNAMIC)


