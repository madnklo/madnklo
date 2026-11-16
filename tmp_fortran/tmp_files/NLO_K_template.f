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





      subroutine local_counter_chan_NLO_%(isec)d_%(jsec)d(xs,xp,wgt,xj,xjB,ichan,KNLO,KCH,KNC,KCHNC,ierr)
c     channel decomposition (run_card channel_decomposition = T):
c     the same local counterterm as
c     local_counter_NLO_%(isec)d_%(jsec)d, split into pieces, one
c     per momentum mapping (iU,iS,iB): each soft dipole (isec,l,m)
c     (the dipole between the two initial-state partons in its two
c     partial fractions, on (isec,2,1) and (isec,1,2)) and the
c     hard-collinear piece,
c     the latter evaluated on the Born of its own mapping, the one of
c     the sector. Returns the full
c     counterterm KNLO and the part KCH of the pieces whose mapping
c     is the one of channel ichan. Only the pieces of channel ichan
c     are filled in the histograms and in the scale variations
c     (wgt_chan = 1, 0 otherwise).
c     The pieces are evaluated twice: first without the cuts on their
c     mapped Born (chan_nocut, no fills nor scale variations), giving
c     KNC and KCHNC for the channel fractions (NLO_Rsub), then with
c     the cuts, giving KNLO and KCH.
      implicit none
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      include 'colored_partons.inc'
      integer ierr,nitR,ichan
      common/iterations/nitR
      double precision xs(nexternal,nexternal)
      double precision xp(0:3,nexternal)
      double precision wgt,xj,xjB,KNLO,KCH,KNC,KCHNC,val,wch
      double precision wsel,wfill,KT(2),KC(2)
      integer ipass
      logical chan_nocut
      common/c_chan_nocut/chan_nocut
      double precision xsbhc(nexternal-1,nexternal-1)
      double precision xpbhc(0:3,nexternal-1),xjcs
      integer l,m,ke,kr,imap,nmap
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer ldip,mdip
      common/c_chan_dipole/ldip,mdip
      integer kedip
      common/c_chan_dipmap/kedip
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
      KNC=0d0
      KCHNC=0d0
      KT=0d0
      KC=0d0
c
c     pass 1: without cuts on the mapped Born, no fills (wfill = 0);
c     pass 2: with cuts, fills of the pieces of channel ichan
      do ipass=1,2
      chan_nocut=(ipass.eq.1)
      wfill=0d0
      if(ipass.eq.2)wfill=1d0
c
%(chan_soft_block)s
c
%(chan_hc_block)s
c
      enddo
      chan_nocut=.false.
      KNC=KT(1)
      KCHNC=KC(1)
      KNLO=KT(2)
      KCH=KC(2)
c
      return
 999  ierr=1
      ldip=0
      mdip=0
      kedip=0
      chan_nocut=.false.
      return
      end
