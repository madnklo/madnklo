# ----------------------------------------------------------------------------
#
# Makefile for PDF library
# Apr. 25 2003
#
# ----------------------------------------------------------------------------

LIBRARY	      = libpdf.$(libext)
LIBDIR        = ../../lib/

include ../make_opts

# Hadron beams: always LHAPDF 6 (>= 6.2), with the path to lhapdf-config set
# in input/madnklo_config.txt (bin/run_all.py writes 'lhapdf' in ../make_opts).
# Without it (lepton beams), the built-in library below is compiled but no
# PDF is ever evaluated.
ifdef lhapdf
  CXXFLAGS+=-std=c++11
  PDF         = pdfwrap_lhapdf.o pdf_lhapdf62.o pdg2pdf_lhapdf6.o opendata.o \
                pdf_alphas_lhapdf.o
else
  PDF         = Ctq4Fn.o Ctq5Par.o Ctq5Pdf.o Partonx5.o Ctq6Pdf.o cteq3.o \
                mrs98.o mrs98lo.o mrs98ht.o mrs99.o mrst2001.o mrst2002.o jeppe02.o \
                pdfwrap.o opendata.o pdf.o PhotonFlux.o pdg2pdf.o NNPDFDriver.o \
                pdf_alphas_none.o
endif

all: $(LIBDIR)$(LIBRARY)

$(LIBDIR)$(LIBRARY): $(PDF)
	$(call CREATELIB, $@, $^)

clean:
	$(RM) *.o $(LIBDIR)$(LIBRARY)
