c -*- Fortran -*-
      integer maxdim
      parameter(maxdim=20)
      integer njet,jet(maxdim)
      double precision pjet(0:3,maxdim)
      double precision ptjet(maxdim),etajet(maxdim),yjet(maxdim)
      common/cfastjeti/njet,jet
      common/cfastjetr/pjet,ptjet,etajet,yjet
