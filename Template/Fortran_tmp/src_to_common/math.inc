c -*- Fortran -*-
c
      real * 8 pi,twopi,fourpi
      real * 8 zeta2,zeta3,eulergamma
      real * 8 dim,GeVtopb
      real * 8 CF,CA,TR,beta0
      real * 8 gamma_g,gamma_q,phi_g,phi_q
      real * 8 gamma_hc_g,gamma_hc_q,phi_hc_g,phi_hc_q
      real * 8 tiny0,tiny1,tiny2
      real * 8 asMZ
c     TODO: remove all unnecessary terms
c     number of colors
      integer Nc
      integer nloop
c
c     mathematical constants
      parameter (pi   =3.141592653589793238462643383279502884197d0)
      parameter (twopi=2d0*pi)
      parameter (fourpi=4d0*pi)
      parameter (zeta2=1.644934066848226436472415166646025189219d0)
      parameter (zeta3=1.202056903159594285399738161511449990765d0)
      parameter (eulergamma=0.5772156649015328606065120900824024d0)
      parameter (CF=4d0/3d0, CA=3d0, TR=1d0/2d0, Nc=3)
c
c     TODO: From parameters.f in local code
c     couplings (read muR from input file in driver)
      parameter (asMZ=0.118d0)
      parameter (nloop=2)
c
      parameter (beta0=(11d0*CA-4d0*TR*Nf)/3d0)
      parameter (gamma_g=beta0/2d0)
      parameter (gamma_q=CF*3d0/2d0)
      parameter (gamma_hc_g=gamma_g-2d0*CA)
      parameter (gamma_hc_q=gamma_q-2d0*CF)
      parameter (phi_g=beta0*4d0/3d0+CA*(2d0/3d0-7d0/2d0*zeta2))
      parameter (phi_q=CF*13d0/3d0+CF*(2d0/3d0-7d0/2d0*zeta2))
      parameter (phi_hc_g=beta0*4d0/3d0-16d0/3d0*CA)
      parameter (phi_hc_q=-CF)
c
c     space-time dimension
      parameter (dim=4d0)
c
c     pb to GeV^-2 conversion
c     sigma(pb) = GeVtopb * sigma(GeV^-2)
      parameter (GeVtopb=0.38937966d9)
c
c     small constants
      parameter (tiny0=1d-6)
      parameter (tiny1=1d-8)
      parameter (tiny2=1d-10)
      real * 8, parameter :: zero=0d0
