      subroutine local_counter_NLO_%(isec)d_%(jsec)d(xs,xp,xsb,xpb,wgt,xj,xjB,x,KNLO,wgt_chan,ierr)
c     local NLO counterterm for sector (isec,jsec)
      implicit none
      include 'nexternal.inc'
      integer ierr,nitR
      common/iterations/nitR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xpb(0:3,nexternal-1)
      double precision wgt,xj,xjB
      INTEGER, PARAMETER :: MXDIM = 30
      DOUBLE PRECISION X(MXDIM)
      double precision KS,KHC,KNLO,wgt_chan
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      logical default_soft
      parameter(default_soft=.true.)
      %(str_def_M2)s
c
c     initialise
      ierr=0
      KS=0d0
      KHC=0d0
      KNLO=0d0
c
c     counterterms
      %(str_M2)s
c     
c     combination
      KNLO=KS+KHC
c
      return
 999  ierr=1
      return
      end





      subroutine local_counter_chan_NLO_%(isec)d_%(jsec)d(xs,xp,wgt,xj,xjB,ichan,KNLO,KCH,ierr)
c     channel decomposition (run_card channel_decomposition = T):
c     the same local counterterm as
c     local_counter_NLO_%(isec)d_%(jsec)d, split into pieces, one
c     per momentum mapping (iU,iS,iB): each soft dipole (isec,l,m)
c     and the hard-collinear piece (isec,jsec,iref), the latter
c     evaluated on the Born of its own mapping. Returns the full
c     counterterm KNLO and the part KCH of the pieces whose mapping
c     is the one of channel ichan. Only the pieces of channel ichan
c     are filled in the histograms and in the scale variations
c     (wgt_chan = 1, 0 otherwise).
      implicit none
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      include 'colored_partons.inc'
      integer ierr,nitR,ichan
      common/iterations/nitR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision wgt,xj,xjB,KNLO,KCH,val,wch
      double precision xsbhc(nexternal-1,nexternal-1)
      double precision xpbhc(0:3,nexternal-1),xjcs
      integer l,m
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer ldip,mdip
      common/c_chan_dipole/ldip,mdip
      integer maxchan
      parameter(maxchan=100)
      integer nchan,chtrip(3,maxchan)
      common/c_nlo_channels/nchan,chtrip
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      %(str_def_M2)s
c
c     initialise
      ierr=0
      KNLO=0d0
      KCH=0d0
c
%(chan_soft_block)s
c
%(chan_hc_block)s
c
      return
 999  ierr=1
      ldip=0
      mdip=0
      return
      end
