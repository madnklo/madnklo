      double precision function int_real_%(isec)d_%(jsec)d(x,wgt)
c     (n+1)-body NLO integrand for vegas.
c     With hadron beams, the event generated at each phase-space point
c     is evaluated in its two orientations along the beam axis, related
c     by the reflection z -> -z (equal beam energies): parton 1 from
c     beam 1 (the processes of this directory, imirror = 1) and, if the
c     directory has mirror processes, parton 1 from beam 2 (the mirror
c     processes, imirror = 2), each with its parton luminosity at the
c     momentum fractions of the real emission (lumi_npo, pdf_lumi.f).
c     The integrand depends on the momenta only through invariants,
c     unchanged by the reflection, so that the sampling of the singular
c     regions is the same in the two orientations.
      implicit none
      include 'nexternal.inc'
      include 'run.inc'
      include 'initial_states.inc'
      integer, parameter :: mxdim = 30
      double precision x(mxdim),wgt
      double precision int_real_orient_%(isec)d_%(jsec)d
      external int_real_orient_%(isec)d_%(jsec)d
      integer ntest,ntested
      parameter(ntest=20)
      save ntested
      integer imirror_npo
      common/c_imirror_npo/imirror_npo
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      logical firsttime
      data firsttime/.true./
      save firsttime
c
c     initialisation of the sector
      call initialise_sector()
      if(firsttime)then
         firsttime=.false.
c        the reflection z -> -z (mirror orientation, or parton 1 in
c        the Born slot 2) needs equal beam energies
         if((lpp(1).ne.0.or.lpp(2).ne.0).and.(has_mirror.or.mapped_labels(1).eq.2).and.ebeam(1).ne.ebeam(2))then
            write(*,*)'ERROR: uneven beam energies, ebeam1 =',ebeam(1),', ebeam2 =',ebeam(2),','
            write(*,*)'are not supported yet in the real emission with two hadron beams:'
            write(*,*)'the mirror processes are evaluated by reflecting the events (z -> -z),'
            write(*,*)'which needs ebeam1 = ebeam2 (run_card.dat).'
            flush(6)
            stop 1
         endif
      endif
c
c     test phase-space singularities of matrix elements (without PDFs)
      if(ntested.lt.ntest)then
         ntested=ntested+1
         imirror_npo=0
         call test_R_%(isec)d_%(jsec)d(ntested)
      endif
c     the limit tests call the counterterms: discard their contributions
c     to the renormalisation-scale variations
      call sv_point_reset()
c
      if(lpp(1).eq.0.and.lpp(2).eq.0)then
c        lepton beams: no PDFs
         int_real_%(isec)d_%(jsec)d=int_real_orient_%(isec)d_%(jsec)d(x,wgt,0)
      else
         int_real_%(isec)d_%(jsec)d=int_real_orient_%(isec)d_%(jsec)d(x,wgt,1)
         if(has_mirror)int_real_%(isec)d_%(jsec)d=int_real_%(isec)d_%(jsec)d+int_real_orient_%(isec)d_%(jsec)d(x,wgt,2)
      endif
c
      return
      end


      double precision function int_real_orient_%(isec)d_%(jsec)d(x,wgt,imirror)
c     (n+1)-body NLO integrand for one orientation imirror of the event
c     (0: lepton beams; 1, 2: see int_real_%(isec)d_%(jsec)d)
      use sectors2_module
      implicit none
      double precision fl_factor_npo
      external fl_factor_npo
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
      integer ierr,imirror
      integer ievt,nthres,ntest
      parameter(ntest=20)
      save ievt,nthres
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
      integer chref(maxchan)
      common/c_nlo_chan_ref/chref
      integer jU,jS,jB,jA
      double precision KCH,KNC,KCHNC,rho,Kstd,xjcs,lumr
      double precision pbs(0:3,nexternal-1)
      double precision sLOs(nexternal-1,nexternal-1)
      integer nchecked
      data nchecked/0/
      save nchecked
      integer imirror_npo
      common/c_imirror_npo/imirror_npo
      double precision lumi_npo
      external lumi_npo
c
c     initialise
      imirror_npo=imirror
      xjac = 0d0
      xjacB = 0d0
      int_real_orient_%(isec)d_%(jsec)d=0d0
      int_real_no_cnt=0d0
      RNLO=0d0
      WGT_CHAN=1d0
      sig2_called=.false.
c
c     phase space and invariants
c     (with the channel decomposition: mapping of the current channel)
      jU=iU
      jS=iS
      jB=iB
      jA=iA
      if(channel_decomposition)then
         jU=chtrip(1,ichan)
         jS=chtrip(2,ichan)
         jB=chtrip(3,ichan)
         jA=chref(ichan)
      endif
      call phase_space_npo(x,jU,jS,jB,jA,p,pb,xjac,xjacB)
      if(xjac*xjacB.eq.0d0) then
         write(77,*) 'int_real: '
         write(77,*) 'Jacobians = 0 in phase space ', xjac, xjacB
         goto 999
      endif
c
c     orientation of the event (hadron beams): parton 1 along +z (from
c     beam 1) for imirror = 1, along -z (from beam 2) for imirror = 2.
c     gen_mom puts the Born slot 1 along +z, and parton 1 is in the Born
c     slot mapped_labels(1)
      if((imirror.eq.1.and.mapped_labels(1).eq.2).or.(imirror.eq.2.and.mapped_labels(1).eq.1))then
         call reflect_z(p,nexternal)
         call reflect_z(pb,nexternal-1)
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
      RNLO = ANS(0) * fl_factor_npo(%(NLO_proc_str)sfl_factor)
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
c     parton luminosity (hadron beams), at the real-emission scale
      if(imirror.ne.0)int_real_no_cnt=int_real_no_cnt*lumi_npo(p,nexternal,leg_pdgs)
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
c     (|R W - K| < channel_delta |K|), 1/nchan elsewhere, with K_ch
c     and K evaluated without the cuts on the mapped Born for rho and
c     for the criterion. Summed over the channels this is R W - K at
c     every phase-space point.
      if(channel_decomposition)then
         call local_counter_chan_NLO_%(isec)d_%(jsec)d(sNLO,p,wgt,xjac,xjacB,ichan,KNLO,KCH,KNC,KCHNC,ierr)
         if(ierr.ne.0) then
            nfail=nfail+1
            write(77,*)'ierr failure: R',ierr,ichan,failure_stage
            write(77,*) 'int_real: '
            write(77,*) 'Something wrong in the counterterm', KNLO
            goto 999
         endif
c        the fractions are computed from the pieces without the cuts
c        on their mapped Born (KNC, KCHNC): near every singular limit
c        R W -> KNC whatever the cuts, and at the cut boundaries R goes
c        to the channels that sample its singularities
         if(KNC.ne.0d0.and.abs(int_real_no_cnt-KNC).lt.channel_delta*abs(KNC))then
            rho=KCHNC/KNC
         else
            rho=1d0/dble(nchan)
         endif
         if(int_real_no_cnt.ne.0d0)then
c           PDF ratios of the muF variations for R (the counterterms
c           above have set those of their own luminosities)
            if(imirror.ne.0)lumr=lumi_npo(p,nexternal,leg_pdgs)
            call sv_add(rho*int_real_no_cnt*wgt_chan,0d0,born_qcd_power+1)
            wgtpl=rho*int_real_no_cnt*wgt*wgt_chan
            call sv_wgts(wgtpl,0d0,born_qcd_power+1,wgts)
            if(doplot)call analysis_fill(p,sNLO,nexternal,leg_pdgs,wgts)
         endif
         int_real_orient_%(isec)d_%(jsec)d=(rho*int_real_no_cnt-KCH)*wgt_chan
c     check, at the first points: the pieces add up to the
c     counterterm of the sector (evaluated with the Born of the
c     sector mapping, without histogram fills or scale variations)
         if(nchecked.lt.ntest.and.KNLO.ne.0d0)then
            nchecked=nchecked+1
            call phase_space_CS_inv(iU,iS,iB,p,pbs,nexternal,leg_pdgs,xjcs,mapped_labels)
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
      int_real_orient_%(isec)d_%(jsec)d=int_real_no_cnt-KNLO
      int_real_orient_%(isec)d_%(jsec)d = int_real_orient_%(isec)d_%(jsec)d*wgt_chan
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
      integer kedip
      common/c_chan_dipmap/kedip
      logical chan_nocut
      common/c_chan_nocut/chan_nocut
      logical done
      data done/.false./
      save done
c
      if(done)return
      done=.true.
      ldip=0
      mdip=0
      kedip=0
      chan_nocut=.false.
      isec = %(isec)d
      jsec = %(jsec)d
      ksec = 0
      lsec = 0
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
c
c     phase-space mapping (iU,iS,iB) of the sector, and reference
c     momentum iA for the azimuth of the radiation (sectors.py)
      %(mapping_str)s
c
      return
      end



      subroutine setup_channels_%(isec)d_%(jsec)d(nch)
c     channel decomposition: list of the momentum mappings
c     (iU,iS,iB) of the counterterm pieces of sector
c     (%(isec)d,%(jsec)d), without repetitions:
c       soft dipoles (isec,l,m), l > m coloured and different from
c       isec; the dipole between the two initial-state partons has
c       two mappings, (isec,2,1) and (isec,1,2), each with the partial
c       fraction of the eikonal collinear-singular along its rescaled
c       beam (see M2_S_g): one channel each;
c       the hard-collinear piece, with the mapping of the sector
c       (iU,iS,iB) (initialise_sector).
c     and, for each channel, the reference momentum chref for the
c     azimuth of the radiation, with the rule of sectors.py: 1 if no
c     initial-state parton is in the mapping, otherwise the first
c     final-state particle outside the mapping.
      implicit none
      include 'nexternal.inc'
      include 'colored_partons.inc'
      integer nch,l,m,k,imap,nmap
      integer iU,iS,iB,iA
      common/cNLOmaplabels/iU,iS,iB,iA
      integer maxchan
      parameter(maxchan=100)
      integer nchan,chtrip(3,maxchan)
      common/c_nlo_channels/nchan,chtrip
      integer chref(maxchan)
      common/c_nlo_chan_ref/chref
      logical has_soft,has_hc
      parameter(has_soft=%(chan_has_soft)s,has_hc=%(chan_has_hc)s)
c
      call initialise_sector()
      nchan=0
      if(has_soft)then
         do m=1,nexternal-1
            if(.not.isNLOQCDparton(m))cycle
            if(m.eq.iU)cycle
            do l=m+1,nexternal
               if(.not.isNLOQCDparton(l))cycle
               if(l.eq.iU)cycle
               nmap=1
               if(l.le.2)nmap=2
               do imap=1,nmap
                  if(imap.eq.1)then
                     call add_channel(iU,l,m)
                  else
                     call add_channel(iU,m,l)
                  endif
               enddo
            enddo
         enddo
      endif
      if(has_hc)call add_channel(iU,iS,iB)
      if(nchan.eq.0)then
         write(*,*)'setup_channels: no channels'
         stop 1
      endif
c     reference momentum of the azimuth
      do k=1,nchan
         chref(k)=1
         if(min(chtrip(2,k),chtrip(3,k)).le.2)then
            do l=nexternal,3,-1
               if(l.ne.chtrip(1,k).and.l.ne.chtrip(2,k).and.l.ne.chtrip(3,k))chref(k)=l
            enddo
         endif
      enddo
      nch=nchan
      return
      end


      subroutine add_channel(jU,jS,jB)
c     channel decomposition: add the mapping (jU,jS,jB) to the list of
c     channels, if not already there
      implicit none
      integer jU,jS,jB,k
      integer maxchan
      parameter(maxchan=100)
      integer nchan,chtrip(3,maxchan)
      common/c_nlo_channels/nchan,chtrip
c
      do k=1,nchan
         if(chtrip(1,k).eq.jU.and.chtrip(2,k).eq.jS.and.chtrip(3,k).eq.jB)return
      enddo
      nchan=nchan+1
      if(nchan.gt.maxchan)then
         write(*,*)'add_channel: too many channels',nchan
         stop 1
      endif
      chtrip(1,nchan)=jU
      chtrip(2,nchan)=jS
      chtrip(3,nchan)=jB
      return
      end
