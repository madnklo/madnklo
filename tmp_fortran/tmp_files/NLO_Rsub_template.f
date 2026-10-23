      double precision function int_real_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NLO integrand for vegas
      use sectors2_module
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'nexternal.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'cuts.inc'
      INCLUDE 'leg_PDGs.inc'
      INCLUDE 'ngraphs_%(UBgraphs)s.inc'
      INCLUDE 'qcd_power.inc'
      integer i
      integer ierr
      integer ievt,nthres,ntest
      integer ntested
      parameter(ntest=20)
      save ievt,nthres,ntested
      double precision int_real_no_cnt
      double precision sNLO(nexternal,nexternal)
      double precision sLO(nexternal-1,nexternal-1)
      double precision RNLO,KNLO
c     TODO: understand x(mxdim) definition by Vegas
      integer, parameter :: mxdim = 30
      double precision x(mxdim)
      double precision wgt,wgtpl,wgts(sv_maxwgt),wgt_chan
      logical dotechcut
      logical doplot
      common/cdoplot/doplot
      logical docut
      integer iU,iS,iB,iA
      common/cNLOmaplabels/iU,iS,iB,iA
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      double precision p(0:3,nexternal)
      double precision pb(0:3,nexternal-1)
      double precision xjac,xjacB
      integer nitr
      common/iterations/nitr
      integer %(NLO_proc_str)sfl_factor 
      common/%(NLO_proc_str)sflavour_factor/%(NLO_proc_str)sfl_factor
      double precision dummy_ans(0:1),ans(0:1) !TODO SET CORRECTLY RANGE OF ANS 
      double precision alphas, alpha_qcd
      integer, parameter :: hel=-1
      integer ich
      common/comich/ich
      integer nfail
      character*8 failure_stage
      common/cfailure_count/nfail
      common/cfailure_stage/failure_stage
      double precision  amp2(n_max_cg)
      common/to_amp2/amp2
      logical firsttime
      data firsttime/.true./
      save firsttime
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      logical sig2_called
      common/c_sig2_called/sig2_called
c     channel decomposition (run_card channel_decomposition, run.inc)
      integer maxchan
      parameter(maxchan=100)
      integer nchan,chtrip(3,maxchan)
      common/c_nlo_channels/nchan,chtrip
      integer ichan
      common/c_nlo_chan_cur/ichan
      integer jU,jS,jB
      double precision KCH,rho,Kstd,xjcs
      double precision pbs(0:3,nexternal-1)
      double precision sLOs(nexternal-1,nexternal-1)
      integer nchecked
      data nchecked/0/
      save nchecked
c
c     call initialisation function
      if(firsttime)then
         call initialise_sector()
         firsttime=.false.
      endif
c
c     initialise
      xjac = 0d0
      xjacB = 0d0
      int_real_%(isec)d_%(jsec)d=0d0
      int_real_no_cnt=0d0
      RNLO=0d0
      WGT_CHAN=1d0
      sig2_called=.false.
c
c     test phase-space singularities of matrix elements
      if(ntested.lt.ntest)then
         ntested=ntested+1
         call test_R_%(isec)d_%(jsec)d(ntested)
      endif
c     the limit tests call the counterterms: discard their contributions
c     to the renormalisation-scale variations
      call sv_point_reset()
c
c     phase space and invariants
c     (with the channel decomposition: mapping of the current channel)
      jU=iU
      jS=iS
      jB=iB
      if(channel_decomposition)then
         jU=chtrip(1,ichan)
         jS=chtrip(2,ichan)
         jB=chtrip(3,ichan)
      endif
      call phase_space_npo(x,jU,jS,jB,iA,p,pb,xjac,xjacB)
      if(xjac*xjacB.eq.0d0) then
         write(77,*) 'int_real: '
         write(77,*) 'Jacobians = 0 in phase space ', xjac, xjacB
         goto 999
      endif
      call invariants_from_p(p,nexternal,sNLO,ierr)
      if(ierr.ne.0) then
         nfail=nfail+1
         write(77,*)'ierr failure: R',ierr,ich,failure_stage
         write(77,*) 'int_real: '
         write(77,*) 'Wrong NLO invariants ', sNLO
         goto 999
      endif
      call invariants_from_p(pb,nexternal-1,sLO,ierr)  
      if(ierr.ne.0)then
         nfail=nfail+1
         write(77,*)'ierr failure: R',ierr,ich,failure_stage
         goto 999
      endif
c
c     tiny technical phase-space cut to avoid fluctuations
      if(dotechcut(snlo,nexternal,tiny1)) goto 999
c
c     possible cuts
      IF(DOCUT(P,NEXTERNAL,leg_pdgs,1))GOTO 555
c
c     renormalisation scale for the (n+1)-body kinematics
      call set_mur_from_momenta(p,nexternal,leg_pdgs)
      alphas=alphas_current
c
c     real
      call %(NLO_proc_str)sME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      RNLO = ANS(0) * %(NLO_proc_str)sfl_factor
      if(RNLO.lt.0d0.or.abs(RNLO).ge.huge(1d0).or.isnan(RNLO))then
         write(77,*) 'int_real: '
         write(77,*) 'Wrong RNLO', RNLO
         goto 999
      endif
c
c     real sector function
      call get_sig2(snlo,nexternal)
      sig2_called=.true.
      call get_W_NLO(isec,jsec)
c
c     full real in sector Wij
      int_real_no_cnt=RNLO*W_NLO*xjac
c     (with the channel decomposition: below, times the channel
c     fraction)
      if(.not.channel_decomposition)then
c     renormalisation-scale variations: R ~ alpha_s^(b+1), scale of the
c     (n+1)-body kinematics
      call sv_add(int_real_no_cnt*wgt_chan,0d0,born_qcd_power+1)
c
c     plot real
      wgtpl=int_real_no_cnt*wgt*wgt_chan
      call sv_wgts(wgtpl,0d0,born_qcd_power+1,wgts)
      if(doplot)call analysis_fill(p,sNLO,nexternal,leg_pdgs,wgts)
      endif
 555  continue
c
c     channel decomposition: integrand rho * R W - K_ch of channel
c     ichan, with K_ch the counterterm pieces whose momentum mapping
c     is the channel's one and rho = K_ch/K near the singular limits
c     (|R W - K| < channel_delta |K|), 1/nchan elsewhere. Summed
c     over the channels this is R W - K at every phase-space point.
      if(channel_decomposition)then
         call local_counter_chan_NLO_%(isec)d_%(jsec)d(sNLO,p,wgt,xjac,xjacB,ichan,KNLO,KCH,ierr)
         if(ierr.ne.0) then
            nfail=nfail+1
            write(77,*)'ierr failure: R',ierr,ichan,failure_stage
            write(77,*) 'int_real: '
            write(77,*) 'Something wrong in the counterterm', KNLO
            goto 999
         endif
         if(KNLO.ne.0d0.and.abs(int_real_no_cnt-KNLO).lt.channel_delta*abs(KNLO))then
            rho=KCH/KNLO
         else
            rho=1d0/dble(nchan)
         endif
         if(int_real_no_cnt.ne.0d0)then
            call sv_add(rho*int_real_no_cnt*wgt_chan,0d0,born_qcd_power+1)
            wgtpl=rho*int_real_no_cnt*wgt*wgt_chan
            call sv_wgts(wgtpl,0d0,born_qcd_power+1,wgts)
            if(doplot)call analysis_fill(p,sNLO,nexternal,leg_pdgs,wgts)
         endif
         int_real_%(isec)d_%(jsec)d=(rho*int_real_no_cnt-KCH)*wgt_chan
c     check, at the first points: the pieces add up to the
c     counterterm of the sector (evaluated with the Born of the
c     sector mapping, without histogram fills or scale variations)
         if(nchecked.lt.ntest.and.KNLO.ne.0d0)then
            nchecked=nchecked+1
            call phase_space_CS_inv(isec,jsec,iref,p,pbs,nexternal,leg_pdgs,xjcs,mapped_labels)
            if(xjcs.ne.0d0)then
               call invariants_from_p(pbs,nexternal-1,sLOs,ierr)
               if(ierr.eq.0)then
                  call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLOs,pbs,wgt,xjac,xjacB,x,Kstd,0d0,ierr)
                  if(ierr.eq.0.and.Kstd.ne.0d0)write(88,*)'channel check: K(sum of channel pieces)/K - 1 =',KNLO/Kstd-1d0,' channel',ichan
               endif
            endif
            ierr=0
         endif
         goto 999
      endif
c
c     counterterm
      call local_counter_NLO_%(isec)d_%(jsec)d(sNLO,p,sLO,pb,wgt,xjac,xjacB,x,KNLO,wgt_chan,ierr)
      if(ierr.ne.0) then
         nfail=nfail+1
         write(77,*)'ierr failure: R',ierr,ich,failure_stage
         write(77,*) 'int_real: '
         write(77,*) 'Something wrong in the counterterm', KNLO
         goto 999
      endif
c
c     subtraction (phase-space jacobian included in counterterm definition)
      int_real_%(isec)d_%(jsec)d=int_real_no_cnt-KNLO
      int_real_%(isec)d_%(jsec)d = int_real_%(isec)d_%(jsec)d*wgt_chan
c
 999  return
      end



      subroutine initialise_sector()
      implicit none
      include 'nexternal.inc'
      include 'leg_PDGs.inc'
      integer iU,iS,iB,iA
      common/cNLOmaplabels/iU,iS,iB,iA
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      integer ldip,mdip
      common/c_chan_dipole/ldip,mdip
      logical done
      data done/.false./
      save done
c
      if(done)return
      done=.true.
      ldip=0
      mdip=0
      isec = %(isec)d
      jsec = %(jsec)d
      ksec = 0
      lsec = 0
c
c     check we are not in the ISR case
      if(isec.le.2.or.jsec.le.2)then
         write(*,*)'ISR indices',isec,jsec
         stop
      endif
c
c     specify phase-space mapping
      %(mapping_str)s
c
c     configuration files
      call configs_%(strUB)s
      call props_%(strUB)s
      call decaybw_%(strUB)s
      call getleshouche_%(strUB)s
c
c     fill underlying pdgs, labels and flavours
      call get_underlying_pdgs(isec,jsec,ksec,lsec,nexternal-1,underlying_leg_pdgs)
      call get_mapped_labels(nexternal,isec,jsec,leg_pdgs,underlying_leg_pdgs,mapped_labels)
      call get_iref(nexternal,isec,jsec,mapped_labels,iref)
      ib = iref
c
      return
      end



      subroutine setup_channels_%(isec)d_%(jsec)d(nch)
c     channel decomposition: list of the momentum mappings
c     (iU,iS,iB) of the counterterm pieces of sector
c     (%(isec)d,%(jsec)d), without repetitions: soft dipoles
c     (isec,l,m), l > m coloured and different from isec, and the
c     hard-collinear piece (isec,jsec,iref)
      implicit none
      include 'nexternal.inc'
      include 'colored_partons.inc'
      integer nch,l,m,k
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer maxchan
      parameter(maxchan=100)
      integer nchan,chtrip(3,maxchan)
      common/c_nlo_channels/nchan,chtrip
      logical has_soft,has_hc,found
      parameter(has_soft=%(chan_has_soft)s,has_hc=%(chan_has_hc)s)
c
      call initialise_sector()
      nchan=0
      if(has_soft)then
         do m=1,nexternal-1
            if(.not.isNLOQCDparton(m))cycle
            if(m.eq.isec)cycle
            do l=m+1,nexternal
               if(.not.isNLOQCDparton(l))cycle
               if(l.eq.isec)cycle
               nchan=nchan+1
               chtrip(1,nchan)=isec
               chtrip(2,nchan)=l
               chtrip(3,nchan)=m
            enddo
         enddo
      endif
      if(has_hc)then
         found=.false.
         do k=1,nchan
            if(chtrip(2,k).eq.jsec.and.chtrip(3,k).eq.iref)found=.true.
         enddo
         if(.not.found)then
            nchan=nchan+1
            chtrip(1,nchan)=isec
            chtrip(2,nchan)=jsec
            chtrip(3,nchan)=iref
         endif
      endif
      if(nchan.eq.0.or.nchan.gt.maxchan)then
         write(*,*)'setup_channels: wrong number of channels',nchan
         stop
      endif
      nch=nchan
      return
      end
