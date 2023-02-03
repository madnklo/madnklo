c -*- Fortran -*-
c
      real * 8 pi,twopi,fourpi
      real * 8 zeta2,zeta3,eulergamma
      real * 8 dim,GeVtopb
      real * 8 CF,CA,TR,NF,gamma_g,gamma_q
      real * 8 tiny0,tiny1,tiny2
      real * 8 asMZ,nloop,alphaEM
c     number of colors	
      integer Nc
      parameter (CF=4d0/3d0, CA=3d0, TR=1d0/2d0, Nc=3)
c
c     From parameters.f in local code
c     couplings (read muR from input file in driver)
      asMZ=0.118d0
      nloop=2
      alphaEM=1/137.036d0
c
      Nf=5d0
      gamma_g=(11d0*CA-4d0*TR*Nf)/6d0
      gamma_q=3d0*CF/2d0
c
c     mathematical constants
      pi   =3.141592653589793238462643383279502884197d0
      twopi=2d0*pi
      fourpi=4d0*pi
      zeta2=1.644934066848226436472415166646025189219d0
      zeta3=1.202056903159594285399738161511449990765d0
      eulergamma=0.5772156649015328606065120900824024d0
c
c     space-time dimension
      dim=4d0
c
c     pb to GeV^-2 conversion
c     sigma(pb) = GeVtopb * sigma(GeV^-2)
      GeVtopb=0.38937966d9
c
c     small constants
      tiny0=1d-6
      tiny1=1d-8
      tiny2=1d-10
