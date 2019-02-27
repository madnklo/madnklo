      subroutine rotate_invar(pin,pout,cth,sth,cphi,sphi)
c Given the four momentum pin, returns the four momentum pout (in the same
c Lorentz frame) by performing a three-rotation of an angle theta 
c (cos(theta)=cth) around the y axis, followed by a three-rotation of an
c angle phi (cos(phi)=cphi) along the z axis. The components of pin
c and pout are given along these axes
      implicit none
      real*8 cth,sth,cphi,sphi,pin(0:3),pout(0:3)
      real*8 q1,q2,q3
c
      q1=pin(1)
      q2=pin(2)
      q3=pin(3)
      pout(1)=q1*cphi*cth-q2*sphi+q3*cphi*sth
      pout(2)=q1*sphi*cth+q2*cphi+q3*sphi*sth
      pout(3)=-q1*sth+q3*cth 
      pout(0)=pin(0)
      return
      end


      subroutine trp_rotate_invar(pin,pout,cth,sth,cphi,sphi)
c This subroutine performs a rotation in the three-space using a rotation
c matrix that is the transpose of that used in rotate_invar(). Thus, if
c called with the *same* angles, trp_rotate_invar() acting on the output
c of rotate_invar() will return the input of rotate_invar()
      implicit none
      real*8 cth,sth,cphi,sphi,pin(0:3),pout(0:3)
      real*8 q1,q2,q3
c
      q1=pin(1)
      q2=pin(2)
      q3=pin(3)
      pout(1)=q1*cphi*cth+q2*sphi*cth-q3*sth
      pout(2)=-q1*sphi+q2*cphi 
      pout(3)=q1*cphi*sth+q2*sphi*sth+q3*cth
      pout(0)=pin(0)
      return
      end


      subroutine getaziangles(p,cphi,sphi)
      implicit none
      real*8 p(0:3),cphi,sphi
      real*8 xlength,cth,sth
      double precision rho
      external rho
c
      xlength=rho(p)
      if(xlength.ne.0.d0)then
        cth=p(3)/xlength
        sth=sqrt(1-cth**2)
        if(sth.ne.0.d0)then
          cphi=p(1)/(xlength*sth)
          sphi=p(2)/(xlength*sth)
        else
          cphi=1.d0
          sphi=0.d0
        endif
      else
        cphi=1.d0
        sphi=0.d0
      endif
      return
      end

      subroutine phspncheck_born(ecm,xmass,xmom,pass)
c Checks four-momentum conservation.
c WARNING: works only in the partonic c.m. frame
      implicit none
      include 'nexternal.inc'
      real*8 ecm,xmass(nexternal-1),xmom(0:3,nexternal-1)
      real*8 tiny,xm,xlen4,xsum(0:3),xsuma(0:3),xrat(0:3),ptmp(0:3)
      parameter (tiny=5.d-3)
      integer jflag,npart,i,j,jj
      logical pass
c
      pass=.true.
      jflag=0
      npart=nexternal-1
      do i=0,3
        xsum(i)=0.d0
        xsuma(i)=0.d0
        do j=3,npart
          xsum(i)=xsum(i)+xmom(i,j)
          xsuma(i)=xsuma(i)+abs(xmom(i,j))
        enddo
        if(i.eq.0)xsum(i)=xsum(i)-ecm
        if(xsuma(i).lt.1.d0)then
          xrat(i)=abs(xsum(i))
        else
          xrat(i)=abs(xsum(i))/xsuma(i)
        endif
        if(xrat(i).gt.tiny.and.jflag.eq.0)then
          write(*,*)'Momentum is not conserved'
          write(*,*)'i=',i
          do j=1,npart
            write(*,'(4(d14.8,1x))') (xmom(jj,j),jj=0,3)
          enddo
          jflag=1
        endif
      enddo
      if(jflag.eq.1)then
        write(*,'(4(d14.8,1x))') (xsum(jj),jj=0,3)
        write(*,'(4(d14.8,1x))') (xrat(jj),jj=0,3)
        pass=.false.
        return
      endif
c
      do j=1,npart
        do i=0,3
          ptmp(i)=xmom(i,j)
        enddo
        xm=xlen4(ptmp)
        if(abs(xm-xmass(j))/ptmp(0).gt.tiny .and.
     &       abs(xm-xmass(j)).gt.tiny)then
          write(*,*)'Mass shell violation'
          write(*,*)'j=',j
          write(*,*)'mass=',xmass(j)
          write(*,*)'mass computed=',xm
          write(*,'(4(d14.8,1x))') (xmom(jj,j),jj=0,3)
          pass=.false.
          return
        endif
      enddo
      return
      end


      subroutine phspncheck_nocms(npart,ecm,xmass,xmom,pass)
c Checks four-momentum conservation. Derived from phspncheck;
c works in any frame
      implicit none
      integer npart,maxmom
      include "genps.inc"
      include "nexternal.inc"
      real*8 ecm,xmass(-max_branch:max_particles),
     # xmom(0:3,nexternal)
      real*8 tiny,vtiny,xm,xlen4,den,ecmtmp,xsum(0:3),xsuma(0:3),
     # xrat(0:3),ptmp(0:3)
      parameter (tiny=5.d-3)
      parameter (vtiny=1.d-6)
      integer jflag,i,j,jj
      logical pass
      double precision dot
      external dot
c
      pass=.true.
      jflag=0
      do i=0,3
        xsum(i)=-xmom(i,1)-xmom(i,2)
        xsuma(i)=abs(xmom(i,1))+abs(xmom(i,2))
        do j=3,npart
          xsum(i)=xsum(i)+xmom(i,j)
          xsuma(i)=xsuma(i)+abs(xmom(i,j))
        enddo
        if(xsuma(i).lt.1.d0)then
          xrat(i)=abs(xsum(i))
        else
          xrat(i)=abs(xsum(i))/xsuma(i)
        endif
        if(xrat(i).gt.tiny.and.jflag.eq.0)then
          write(*,*)'Momentum is not conserved [nocms]'
          write(*,*)'i=',i
          do j=1,npart
            write(*,'(4(d14.8,1x))') (xmom(jj,j),jj=0,3)
          enddo
          jflag=1
        endif
      enddo
      if(jflag.eq.1)then
        write(*,'(4(d14.8,1x))') (xsum(jj),jj=0,3)
        write(*,'(4(d14.8,1x))') (xrat(jj),jj=0,3)
        pass=.false.
        return
      endif
c
      do j=1,npart
        do i=0,3
          ptmp(i)=xmom(i,j)
        enddo
        xm=xlen4(ptmp)
        if(ptmp(0).ge.1.d0)then
          den=ptmp(0)
        else
          den=1.d0
        endif
        if(abs(xm-xmass(j))/den.gt.tiny .and.
     &       abs(xm-xmass(j)).gt.tiny)then
          write(*,*)'Mass shell violation [nocms]'
          write(*,*)'j=',j
          write(*,*)'mass=',xmass(j)
          write(*,*)'mass computed=',xm
          write(*,'(4(d14.8,1x))') (xmom(jj,j),jj=0,3)
          pass=.false.
          return
        endif
      enddo
c
      ecmtmp=sqrt(2d0*dot(xmom(0,1),xmom(0,2)))
      if(abs(ecm-ecmtmp).gt.vtiny)then
        write(*,*)'Inconsistent shat [nocms]'
        write(*,*)'ecm given=   ',ecm
        write(*,*)'ecm computed=',ecmtmp
        write(*,'(4(d14.8,1x))') (xmom(jj,1),jj=0,3)
        write(*,'(4(d14.8,1x))') (xmom(jj,2),jj=0,3)
        pass=.false.
        return
      endif

      return
      end


      function xlen4(v)
      implicit none
      real*8 xlen4,tmp,v(0:3)
c
      tmp=v(0)**2-v(1)**2-v(2)**2-v(3)**2
      xlen4=sign(1.d0,tmp)*sqrt(abs(tmp))
      return
      end


      double precision function dsig(pp,wgt,vegaswgt)
c Here are the subtraction terms, the Sij function, 
c the f-damping function, and the single diagram
c enhanced multi-channel factor included
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      include "fks_powers.inc"
      include 'coupl.inc'
      include 'q_es.inc'
      include 'run.inc'
      include 'reweight.inc'
      include 'reweightNLO.inc'

c     timing statistics
      include 'timing_variables.inc'
      real deltaTOLP,deltaTPDF,deltaTFJ

      double precision pp(0:3,nexternal),wgt,vegaswgt

      double precision fks_Sij,f_damp,dot,dlum
      external fks_Sij,f_damp,dot,dlum

      double precision x,xtot,s_ev,s_c,s_s,s_sc,ffact,fx_ev,fx_c,
     #                 fx_s,fx_sc,ev_wgt,cnt_wgt_c,cnt_wgt_s,cnt_wgt_sc,
     #                 bsv_wgt,plot_wgt,cnt_swgt_s,cnt_swgt_sc,cnt_sc,cnt_s,
     #                 prefact_cnt_ssc,prefact_cnt_ssc_c,prefact_coll,
     #                 prefact_coll_c,born_wgt,prefact_deg,prefact,prefact_c,
     #                 prefact_deg_sxi,prefact_deg_slxi,deg_wgt,deg_swgt,
     #                 deg_xi_c,deg_lxi_c,deg_xi_sc,deg_lxi_sc,
     #                 cnt_swgt,cnt_wgt,xlum_ev,xlum_c,xlum_s,xlum_sc,xsec,
     #                 bpower,cpower
      integer i,j,iplot

      integer izero,ione,itwo,mohdr,iplot_ev,iplot_cnt,iplot_born
      integer ithree,ifour,ifill2,ifill3,ifill4
      parameter (izero=0)
      parameter (ione=1)
      parameter (itwo=2)
      parameter (ithree=3)
      parameter (ifour=4)
      parameter (mohdr=-100)
      parameter (iplot_ev=11)
      parameter (iplot_cnt=12)
      parameter (iplot_born=20)

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      double precision shattmp
      double precision pi
      parameter (pi=3.1415926535897932385d0)

      logical nocntevents
      common/cnocntevents/nocntevents

      double precision xnormsv
      common/cxnormsv/xnormsv

c Multi channel stuff:
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_amps/  amp2,       jamp2

      integer          isum_hel
      logical                    multi_channel
      common/to_matrix/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig

      double precision wgt1
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn

      double precision ev_enh,enhance,rwgt,unwgtfun
      logical firsttime,passcuts
      data firsttime /.true./
      integer inoborn_ev,inoborn_cnt
      double precision xnoborn_ev,xnoborn_cnt

c FKS stuff:
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision xi_i_fks_ev,y_ij_fks_ev
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

      double precision xi_i_fks_cnt(-2:2)
      common /cxiifkscnt/xi_i_fks_cnt

      double precision xiimax_ev
      common /cxiimaxev/xiimax_ev
      double precision xiimax_cnt(-2:2)
      common /cxiimaxcnt/xiimax_cnt

      double precision xinorm_ev
      common /cxinormev/xinorm_ev
      double precision xinorm_cnt(-2:2)
      common /cxinormcnt/xinorm_cnt

      double precision p1_cnt(0:3,nexternal,-2:2)
      double precision wgt_cnt(-2:2)
      double precision pswgt_cnt(-2:2)
      double precision jac_cnt(-2:2)
      common/counterevnts/p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt

      double precision zero,one
      parameter (zero=0d0,one=1d0)
      double precision tiny
      parameter (tiny=1d-6)

      double precision xicut_used
      common /cxicut_used/xicut_used
      double precision delta_used
      common /cdelta_used/delta_used
      double precision xiScut_used,xiBSVcut_used
      common /cxiScut_used/xiScut_used,xiBSVcut_used
      double precision fkssymmetryfactor,fkssymmetryfactorBorn,
     &     fkssymmetryfactorDeg
      integer ngluons,nquarks(-6:6),nphotons
      common/numberofparticles/fkssymmetryfactor,fkssymmetryfactorBorn,
     &                  fkssymmetryfactorDeg,ngluons,nquarks,nphotons
      double precision diagramsymmetryfactor
      common /dsymfactor/diagramsymmetryfactor

      character*4 abrv
      common /to_abrv/ abrv

      logical nbody
      common/cnbody/nbody

      double precision vegas_weight
      common/cvegas_weight/vegas_weight

c For the MINT folding
      integer fold
      common /cfl/fold

c Random numbers to be used in the plotting routine
      logical needrndec
      parameter (needrndec=.true.)

      real*8 ran2
      external ran2

      real*8 rndec(10)
      common/crndec/rndec

c For tests
      real*8 fksmaxwgt,xisave,ysave
      common/cfksmaxwgt/fksmaxwgt,xisave,ysave

      logical ExceptPSpoint
      integer iminmax
      common/cExceptPSpoint/iminmax,ExceptPSpoint

      double precision central_wgt_saved
      save central_wgt_saved

      double precision dsig_max,dsig_min
      double precision total_wgt_sum,total_wgt_sum_max,
     &                 total_wgt_sum_min
      common/csum_of_wgts/total_wgt_sum,total_wgt_sum_max,
     &                 total_wgt_sum_min
      double precision virt_wgt,born_wgt_ao2pi
      common/c_fks_singular/virt_wgt,born_wgt_ao2pi

c APPLgrid stuff
      integer iappl
      common /for_applgrid/ iappl
      include "appl_common.inc" 

c FxFx merging
      logical rewgt_mohdr_calculated,rewgt_izero_calculated
      double precision rewgt_mohdr,rewgt_izero,rewgt_exp_mohdr
     $     ,rewgt_exp_izero,enhanceS,enhanceH
c$$$      logical setclscales
c$$$      double precision rewgt
c$$$      external setclscales,rewgt

      double precision pmass(nexternal)
      double precision rwgt_muR_dep_fac
      include "pmass.inc"

      vegas_weight=vegaswgt

c$$$c FxFx merging
c$$$      ktscheme=1
c$$$      rewgt_mohdr_calculated=.false.
c$$$      rewgt_izero_calculated=.false.
      if (ickkw.eq.3) then
         write (*,*) 'FKS_SINGULAR NO LONGER COMPATIBLE WITH FXFX MERGING'
         write (*,*) 'THIS NEEDS TO BE FIXED BEFORE ANY RELEASE'
         write (*,*) 'FKS_SINGULAR NO LONGER COMPATIBLE WITH FXFX MERGING'
         write (*,*) 'THIS NEEDS TO BE FIXED BEFORE ANY RELEASE'
         write (*,*) 'FKS_SINGULAR NO LONGER COMPATIBLE WITH FXFX MERGING'
         write (*,*) 'THIS NEEDS TO BE FIXED BEFORE ANY RELEASE'
         write (*,*) 'FKS_SINGULAR NO LONGER COMPATIBLE WITH FXFX MERGING'
         write (*,*) 'THIS NEEDS TO BE FIXED BEFORE ANY RELEASE'
         stop
      endif



      if (fold.eq.0) then
         calculatedBorn=.false.
         call get_helicity(i_fks,j_fks)
      endif

c Make sure that the result can be non-zero. If the jacobian from the
c PS-setup or vegas are zero, we can skip this PS point and 'return'.
c Note that all the wgts and jacs should be positive.
      dsig=0d0
      if ( (wgt.le.0d0 .and. jac_cnt(0).le.0d0 .and. jac_cnt(1).le.0d0
     &     .and. jac_cnt(2).le.0d0) .or. vegaswgt.le.0d0) return

      if (firsttime)then
         inoborn_ev=0
         xnoborn_ev=0.d0
         inoborn_cnt=0
         xnoborn_cnt=0.d0
         fksmaxwgt=0.d0
c Put here call to compute bpower
         call compute_bpower(p_born,bpower)
         wgtbpower=bpower
c Store the power of alphas of the Born events in the appl common block.
         if(iappl.ne.0) appl_bpower = wgtbpower
c Initialize histograms
         call initplot
c Compute cpower done for bottom Yukawa, routine needs to be adopted
c for other muR-dependendent factors
         call compute_cpower(p_born,cpower)
         if(dabs(cpower+1d0).lt.tiny) then
            wgtcpower=0d0
         else
            wgtcpower=cpower
         endif
c Check that things are done consistently
         if(wgtcpower.ne.cpowerinput.and.dabs(cpower+1d0).gt.tiny)then
           write(*,*)'Inconsistency in the computation of cpower',
     #               wgtcpower,cpowerinput
           write(*,*)'Check value in reweight0.inc'
           stop
         endif

         firsttime=.false.
      endif

      prefact=xinorm_ev/xi_i_fks_ev*
     #        1/(1-y_ij_fks_ev)
      if( (.not.nocntevents) .and. (.not.(abrv.eq.'born' .or. abrv.eq
     &     .'grid' .or. abrv(1:2).eq.'vi' .or. nbody))
     &     )then
        prefact_cnt_ssc=xinorm_ev/min(xiimax_ev,xiScut_used)*
     #                  log(xicut_used/min(xiimax_ev,xiScut_used))*
     #                  1/(1-y_ij_fks_ev)
        if(pmass(j_fks).eq.0.d0)then
          prefact_c=xinorm_cnt(ione)/xi_i_fks_cnt(ione)*
     #              1/(1-y_ij_fks_ev)
          prefact_cnt_ssc_c=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                      log(xicut_used/min(xiimax_cnt(ione),xiScut_used))*
     #                      1/(1-y_ij_fks_ev)
          prefact_coll=xinorm_cnt(ione)/xi_i_fks_cnt(ione)*
     #                 log(delta_used/deltaS)/deltaS
          prefact_coll_c=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                   log(xicut_used/min(xiimax_cnt(ione),xiScut_used))*
     #                   log(delta_used/deltaS)/deltaS
          prefact_deg=xinorm_cnt(ione)/xi_i_fks_cnt(ione)*
     #                1/deltaS
          prefact_deg_sxi=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                    log(xicut_used/min(xiimax_cnt(ione),xiScut_used))*
     #                    1/deltaS
          prefact_deg_slxi=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                     ( log(xicut_used)**2 -
     #                       log(min(xiimax_cnt(ione),xiScut_used))**2 )*
     #                     1/(2.d0*deltaS)
        endif
      endif

      if(needrndec)then
        do i=1,10
          rndec(i)=ran2()
        enddo
      endif

c If there was an exceptional phase-space point found for the 
c virtual corrections, at the end of this subroutine, goto 44
c and compute also the "possible" minimal and maximal weight
c these points could have gotton (based upon previous PS
c points)
      ExceptPSpoint=.false.
      iminmax=-1
 44   continue
      iminmax=iminmax+1

      ev_wgt=0.d0
      cnt_wgt=0.d0
      cnt_wgt_s=0.d0
      cnt_wgt_c=0.d0
      cnt_wgt_sc=0.d0
      bsv_wgt=0.d0
      virt_wgt=0d0
      born_wgt=0.d0
      born_wgt_ao2pi=0.d0
      cnt_swgt=0.d0
      cnt_swgt_s=0.d0
      cnt_swgt_sc=0.d0
      deg_wgt=0.d0
      deg_swgt=0.d0
      plot_wgt=0.d0
      iplot=-3

      if(doreweight)then
        call reweight_settozero()
        ifill2=0
        ifill3=0
        ifill4=0
      endif
      if (abrv.eq.'born' .or. abrv.eq.'grid' .or. abrv(1:2).eq.'vi' .or.
     &     nbody)goto 540

      call cpu_time(tBefore)
      deltaTPDF = tPDF
      deltaTFJ  = tFastJet
c Real contribution:
c Set the ybst_til_tolab before applying the cuts. 
      call set_cms_stuff(mohdr)
      if(doreweight)then
        wgtxbj(1,1)=xbk(1)
        wgtxbj(2,1)=xbk(2)
      endif
      if (passcuts(pp,rwgt)) then
c$$$c       Compute the scales and sudakov-reweighting for the FxFx merging
c$$$        if (ickkw.eq.3) then
c$$$           if (.not. setclscales(pp)) then
c$$$               write (*,*) 'ERROR in setclscales mohdr'
c$$$              stop
c$$$           endif
c$$$           rewgt_mohdr=rewgt(pp,rewgt_exp_mohdr)
c$$$           rewgt_mohdr_calculated=.true.
c$$$        endif
        call set_alphaS(pp)
        x = abs(2d0*dot(pp(0,i_fks),pp(0,j_fks))/shat)
        ffact = f_damp(x)
        s_ev = fks_Sij(pp,i_fks,j_fks,xi_i_fks_ev,y_ij_fks_ev)
        if(s_ev.gt.0.d0)then
           call sreal(pp,xi_i_fks_ev,y_ij_fks_ev,fx_ev)
           xlum_ev = dlum()
           xsec = fx_ev*s_ev*ffact*wgt*prefact*rwgt
           ev_wgt = xlum_ev*xsec * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
           if(doreweight)then
             wgtmuR2(1)=muR2_current/muR_over_ref**2
             wgtmuF12(1)=muF12_current/muF1_over_ref**2
             wgtmuF22(1)=muF22_current/muF2_over_ref**2
             call reweight_fillkin(pp,ione)
             wgtwreal(1)=xsec/g**(nint(2*wgtbpower+2.d0))
           endif
        endif
      endif
      call cpu_time(tAfter)
      tDSigR=tDSigR + (tAfter-tBefore) - (tPDF-deltaTPDF) -
     &     (tFastJet-deltaTFJ)
c
c All counterevent have the same final-state kinematics. Check that
c one of them passes the hard cuts, and they exist at all
 540  continue
      call cpu_time(tBefore)
      deltaTPDF = tPDF
      deltaTFJ  = tFastJet
      deltaTOLP = tOLP
c Set the ybst_til_tolab before applying the cuts. Update below
c for the collinear, soft and/or soft-collinear subtraction terms
      call set_cms_stuff(izero)
      if ( (.not.passcuts(p1_cnt(0,1,0),rwgt)) .or.
     #     nocntevents ) goto 547

c$$$c Compute the scales and sudakov-reweighting for the FxFx merging
c$$$      if (ickkw.eq.3) then
c$$$         if (.not. setclscales(p1_cnt(0,1,0))) then
c$$$            write (*,*) 'ERROR in setclscales izero'
c$$$            stop
c$$$         endif
c$$$         rewgt_izero=rewgt(p1_cnt(0,1,0),rewgt_exp_izero)
c$$$         rewgt_izero_calculated=.true.
c$$$      endif
      call set_alphaS(p1_cnt(0,1,0))
      if(doreweight)then
        wgtqes2(2)=QES2
        wgtqes2(3)=QES2
        wgtqes2(4)=QES2
        wgtmuR2(2)=muR2_current/muR_over_ref**2
        wgtmuF12(2)=muF12_current/muF1_over_ref**2
        wgtmuF22(2)=muF22_current/muF2_over_ref**2
        call reweight_fillkin(pp,itwo)
        ifill2=1
      endif

      if (abrv.eq.'born' .or. abrv.eq.'grid' .or. abrv(1:2).eq.'vi' .or.
     &     nbody)goto 545
c
c Collinear subtraction term:
      if( y_ij_fks_ev.gt.1d0-deltaS .and.
     #    pmass(j_fks).eq.0.d0 )then
         call set_cms_stuff(ione)
         if(doreweight)then
           wgtxbj(1,3)=xbk(1)
           wgtxbj(2,3)=xbk(2)
         endif
         s_c = fks_Sij(p1_cnt(0,1,1),i_fks,j_fks,xi_i_fks_cnt(ione),one)
         if(s_c.gt.0.d0)then
            if(abs(s_c-1.d0).gt.1.d-6.and.j_fks.le.nincoming)then
              write(*,*)'Wrong S function in dsig[c]',s_c
              stop
            endif
            call sreal(p1_cnt(0,1,1),xi_i_fks_cnt(ione),one,fx_c)
            xlum_c = dlum()
            xsec = fx_c*s_c*jac_cnt(1)*(prefact_c+prefact_coll)*rwgt
            cnt_wgt_c=cnt_wgt_c-xlum_c*xsec * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            call sreal_deg(p1_cnt(0,1,1),xi_i_fks_cnt(ione),one,
     #                     deg_xi_c,deg_lxi_c)
            deg_wgt=deg_wgt+( deg_xi_c+deg_lxi_c*log(xi_i_fks_cnt(ione)) )*
     #                      jac_cnt(1)*prefact_deg*rwgt/(shat/(32*pi**2))*
     #                      xlum_c * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            iplot=1
            if(doreweight)then
              call reweight_fillkin(pp,ithree)
              ifill3=1
              wgtwreal(3)=-xsec/g**(nint(2*wgtbpower+2.d0))
              wgtwdeg(3)=
     #              ( wgtdegrem_xi+wgtdegrem_lxi*log(xi_i_fks_cnt(ione)) )*
     #                      jac_cnt(1)*prefact_deg*rwgt/(shat/(32*pi**2))/
     #              g**(nint(2*wgtbpower+2.d0))
              wgtwdegmuf(3)=wgtdegrem_muF *
     #                      jac_cnt(1)*prefact_deg*rwgt/(shat/(32*pi**2))/
     #              g**(nint(2*wgtbpower+2.d0))
            endif
         endif
      endif
c Soft subtraction term:
 545  continue
      if (xi_i_fks_ev .lt. max(xiScut_used,xiBSVcut_used)) then
         call set_cms_stuff(izero)
         if(doreweight)then
           wgtxbj(1,2)=xbk(1)
           wgtxbj(2,2)=xbk(2)
           if(ifill2.ne.1)then
             write(*,*)'Error #1[wg] in dsig',ifill2
             stop
           endif
         endif
         s_s = fks_Sij(p1_cnt(0,1,0),i_fks,j_fks,zero,y_ij_fks_ev)
         if(nbody)s_s=1.d0
         if(s_s.gt.0.d0)then
            xlum_s = dlum()
            if (abrv.eq.'born' .or. abrv.eq.'grid' .or.
     &           abrv(1:2).eq.'vi' .or. nbody) goto 546
            if (xi_i_fks_ev .lt. xiScut_used) then
              call sreal(p1_cnt(0,1,0),zero,y_ij_fks_ev,fx_s)
              xsec=fx_s*s_s*jac_cnt(0)
              cnt_s=xlum_s*xsec
              cnt_wgt_s=cnt_wgt_s-cnt_s*prefact*rwgt
     f             * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              cnt_swgt_s=cnt_swgt_s-cnt_s*prefact_cnt_ssc*rwgt
     f             * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              if(doreweight)
     #          wgtwreal(2)=-xsec*(prefact+prefact_cnt_ssc)*rwgt/
     #                      g**(nint(2*wgtbpower+2.d0))
            endif
 546        continue
            if (abrv.eq.'real' .or. .not.nbody) goto 548
            if (xi_i_fks_ev .lt. xiBSVcut_used) then
              xsec=s_s*jac_cnt(0)*xinorm_ev/
     #             (min(xiimax_ev,xiBSVcut_used)*shat/(16*pi**2))*
     #             rwgt
              xnormsv=xlum_s*xsec
              call bornsoftvirtual(p1_cnt(0,1,0),bsv_wgt,virt_wgt
     $             ,born_wgt)
              born_wgt_ao2pi=born_wgt*g**2/(8d0*Pi**2)
c$$$c For FxFx merging, include the compensation term
c$$$              if (rewgt_izero_calculated.and.rewgt_izero.lt.1d0) then
c$$$                 bsv_wgt=bsv_wgt-g**2/(4d0*Pi)*rewgt_exp_izero*born_wgt
c$$$              endif
              if(doreweight)then
                if(wgtbpower.gt.0)then
                  wgtwborn(2)=born_wgt*xsec/g**(nint(2*wgtbpower))
                else
                  wgtwborn(2)=born_wgt*xsec
                endif
                wgtwns(2)=wgtnstmp*xsec/g**(nint(2*wgtbpower+2.d0))
c$$$                 if (rewgt_izero_calculated.and.rewgt_izero.lt.1d0) then
c$$$                    wgtwns(2)=wgtwns(2)-rewgt_exp_izero*wgtwborn(2)
c$$$     $                   /(4d0*pi)
c$$$                 endif
                wgtwnsmuf(2)=wgtwnstmpmuf*xsec/g**(nint(2*wgtbpower
     &               +2.d0))
                wgtwnsmur(2)=wgtwnstmpmur*xsec/g**(nint(2*wgtbpower
     &               +2.d0))
              endif
              bsv_wgt=bsv_wgt*xnormsv * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              virt_wgt=virt_wgt*xnormsv * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              born_wgt=born_wgt*xnormsv * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              born_wgt_ao2pi=born_wgt_ao2pi*xnormsv
     f             * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            endif
 548        continue
            iplot=0
         endif
      endif
c Soft-Collinear subtraction term:
      if (abrv.eq.'born' .or. abrv.eq.'grid' .or. abrv(1:2).eq.'vi' .or.
     &     nbody)goto 547
      if (xi_i_fks_cnt(ione) .lt. xiScut_used .and.
     #    y_ij_fks_ev .gt. 1d0-deltaS .and.
     #    pmass(j_fks).eq.0.d0 )then
         call set_cms_stuff(itwo)
         if(doreweight)then
           wgtxbj(1,4)=xbk(1)
           wgtxbj(2,4)=xbk(2)
         endif
         s_sc = fks_Sij(p1_cnt(0,1,2),i_fks,j_fks,zero,one)
         if(s_sc.gt.0.d0)then
            if(abs(s_sc-1.d0).gt.1.d-6.and.j_fks.le.nincoming)then
              write(*,*)'Wrong S function in dsig[sc]',s_sc
              stop
            endif
            xlum_sc = dlum()
            call sreal(p1_cnt(0,1,2),zero,one,fx_sc)
            xsec=fx_sc*s_sc*jac_cnt(2)
            cnt_sc=xlum_sc*xsec
            cnt_wgt_sc=cnt_wgt_sc+cnt_sc*(prefact_c+prefact_coll)*rwgt
     f           * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            cnt_swgt_sc=cnt_swgt_sc+
     &           cnt_sc*(prefact_cnt_ssc_c+prefact_coll_c)*rwgt
     f           * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            call sreal_deg(p1_cnt(0,1,2),zero,one,
     #                     deg_xi_sc,deg_lxi_sc)
            deg_wgt=deg_wgt-( deg_xi_sc+deg_lxi_sc*log(xi_i_fks_cnt(ione)) )*
     #                     jac_cnt(2)*prefact_deg*rwgt/(shat/(32*pi**2))*
     #                     xlum_sc * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            deg_swgt=deg_swgt-( deg_xi_sc*prefact_deg_sxi +
     #                     deg_lxi_sc*prefact_deg_slxi )*
     #                     jac_cnt(2)*rwgt/(shat/(32*pi**2))*
     #                     xlum_sc * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            if(iplot.ne.0)iplot=2
            if(doreweight)then
              call reweight_fillkin(pp,ifour)
              ifill4=1
              wgtwreal(4)=xsec*(prefact_c+prefact_coll+
     #                          prefact_cnt_ssc_c+prefact_coll_c)*rwgt/
     #                    g**(nint(2*wgtbpower+2.d0))
              wgtwdeg(4)=(
     # -( wgtdegrem_xi+wgtdegrem_lxi*log(xi_i_fks_cnt(ione)) )*prefact_deg
     # -( wgtdegrem_xi*prefact_deg_sxi+wgtdegrem_lxi*prefact_deg_slxi ) )*
     #                     jac_cnt(2)*rwgt/(shat/(32*pi**2))/
     #              g**(nint(2*wgtbpower+2.d0))
              wgtwdegmuf(4)=
     #          -wgtdegrem_muF*( prefact_deg+prefact_deg_sxi )*
     #                     jac_cnt(2)*rwgt/(shat/(32*pi**2))/
     #              g**(nint(2*wgtbpower+2.d0))
            endif
         endif
      endif
 547  continue
c
      call cpu_time(tAfter)
      tDSigI=tDSigI + (tAfter-tBefore) - (tPDF-deltaTPDF) -
     &     (tFastJet-deltaTFJ) - (tOLP-deltaTOLP)
c
c Enhance the one channel for multi-channel integration
c
      enhance=1.d0
      if ((ev_wgt.ne.0d0 .or. cnt_wgt_c.ne.0d0 .or. cnt_wgt_s.ne.0d0
     $     .or.cnt_wgt_sc.ne.0d0 .or. bsv_wgt.ne.0d0 .or.
     $     virt_wgt.ne.0d0 .or. deg_wgt.ne.0d0 .or.deg_swgt.ne.0d0 .or.
     $     cnt_swgt_s.ne.0d0 .or. cnt_swgt_sc.ne.0d0).and.
     $     multi_channel) then
         if (bsv_wgt.eq.0d0 .and. virt_wgt.eq.0d0 .and. deg_wgt.eq.0d0
     $        .and. deg_swgt.eq.0d0 .and. cnt_wgt_c.eq.0d0 )
     $        CalculatedBorn=.false.
         if (.not.calculatedBorn .and. p_born(0,1).gt.0d0)then
            call sborn(p_born,wgt1)
         elseif(p_born(0,1).lt.0d0)then
            enhance=0d0
         endif

         if (enhance.eq.0d0)then
            xnoborn_cnt=xnoborn_cnt+1.d0
            if(log10(xnoborn_cnt).gt.inoborn_cnt)then
               write (*,*) 
     #           'Function dsig: no Born momenta more than 10**',
     #           inoborn_cnt,'times'
               inoborn_cnt=inoborn_cnt+1
            endif
         else
            xtot=0d0
            if (mapconfig(0).eq.0) then
               write (*,*) 'Fatal error in dsig, no Born diagrams '
     &           ,mapconfig,'. Check bornfromreal.inc'
               write (*,*) 'Is fks_singular compiled correctly?'
               stop
            endif
            do i=1, mapconfig(0)
               xtot=xtot+amp2(mapconfig(i))
            enddo
            if (xtot.ne.0d0) then
               enhance=amp2(mapconfig(iconfig))/xtot
               enhance=enhance*diagramsymmetryfactor
            else
               enhance=0d0
            endif
         endif
      endif

      cnt_wgt = cnt_wgt_c + cnt_wgt_s + cnt_wgt_sc
      cnt_swgt = cnt_swgt_s + cnt_swgt_sc

c$$$c Apply the FxFx Sudakov damping on the H events
c$$$      if (ev_wgt.ne.0d0 .and. ickkw.eq.3..and.
c$$$     $     .not.rewgt_mohdr_calculated) then
c$$$         write (*,*) 'Error rewgt_mohdr_calculated'
c$$$         stop
c$$$      elseif(rewgt_mohdr_calculated) then
c$$$         if (rewgt_mohdr.gt.1d0) rewgt_mohdr=1d0
c$$$         enhanceH=enhance*rewgt_mohdr
c$$$      else
         enhanceH=enhance
c$$$      endif

      ev_wgt = ev_wgt * enhanceH

c$$$c Apply the FxFx Sudakov damping on the S events
c$$$      if(.not.(cnt_wgt.eq.0d0 .and. cnt_swgt.eq.0d0 .and. bsv_wgt.eq.0d0
c$$$     $     .and. born_wgt.eq.0d0 .and. deg_wgt.eq.0d0 .and.
c$$$     $     deg_swgt.eq.0d0 )
c$$$     $     .and. ickkw.eq.3 .and. .not.rewgt_izero_calculated) then
c$$$         write (*,*) 'Error rewgt_izero_calculated'
c$$$         stop
c$$$      elseif(rewgt_izero_calculated) then
c$$$         if (rewgt_izero.gt.1d0) rewgt_izero=1d0
c$$$         enhanceS=enhance*rewgt_izero
c$$$      else
         enhanceS=enhance
c$$$      endif


      cnt_wgt = cnt_wgt * enhanceS
      cnt_swgt = cnt_swgt * enhanceS
      bsv_wgt = bsv_wgt * enhanceS
      virt_wgt = virt_wgt * enhanceS
      born_wgt = born_wgt * enhanceS
      born_wgt_ao2pi = born_wgt_ao2pi * enhanceS
      deg_wgt = deg_wgt * enhanceS
      deg_swgt = deg_swgt * enhanceS

      if(iminmax.eq.0) then
         dsig = (ev_wgt+cnt_wgt)*fkssymmetryfactor +
     &        cnt_swgt*fkssymmetryfactor +
     &        bsv_wgt*fkssymmetryfactorBorn +
     &        virt_wgt*fkssymmetryfactorBorn +
     &        deg_wgt*fkssymmetryfactorDeg +
     &        deg_swgt*fkssymmetryfactorDeg

         if (dsig.ne.dsig) then
            write (*,*) 'ERROR #51 in dsig:',dsig,'skipping event'
            dsig=0d0
            return
         endif

         total_wgt_sum=total_wgt_sum+dsig*vegaswgt
         central_wgt_saved=dsig

         call unweight_function(p_born,unwgtfun)
         dsig=dsig*unwgtfun
         virt_wgt=virt_wgt*unwgtfun*fkssymmetryfactorBorn*vegaswgt
         born_wgt_ao2pi=born_wgt_ao2pi*unwgtfun*fkssymmetryfactorBorn
     $        *vegaswgt
         if(doreweight)then
           if(ifill2.eq.0.and.(ifill3.ne.0.or.ifill4.ne.0))then
             write(*,*)'Error #2[wg] in dsig',ifill2,ifill3,ifill4
             stop
           endif
           wgtref = dsig
           xsec = enhanceS*unwgtfun
           do i=1,4
              if (i.eq.1) then
                 wgtwreal(i)=wgtwreal(i) * xsec*fkssymmetryfactor
     $                *enhanceH/enhanceS
              else
                 wgtwreal(i)=wgtwreal(i) * xsec*fkssymmetryfactor
              endif
             wgtwdeg(i)=wgtwdeg(i) * xsec*fkssymmetryfactorDeg
             wgtwdegmuf(i)=wgtwdegmuf(i) * xsec*fkssymmetryfactorDeg
             if(i.eq.2)then
               wgtwborn(i)=wgtwborn(i) * xsec*fkssymmetryfactorBorn
               wgtwns(i)=wgtwns(i) * xsec*fkssymmetryfactorBorn
               wgtwnsmuf(i)=wgtwnsmuf(i) * xsec*fkssymmetryfactorBorn
               wgtwnsmur(i)=wgtwnsmur(i) * xsec*fkssymmetryfactorBorn
             endif
           enddo
           call reweight_fill_extra()
           if(check_reweight.and.doreweight)
     #       call check_rwgt_wgt("NLO")
           if(doreweight)call fill_rwgt_NLOplot()
c Example of reweighted cross section (scale changed)
c           dsig_new=compute_rwgt_wgt_NLO(new_muR_fact,new_muF1_fact,
c     #                                   new_muF2_fact,new_QES_fact,
c     #                                   iwgtinfo)
         endif      

c For tests
        if(abs(dsig).gt.fksmaxwgt)then
            fksmaxwgt=abs(dsig)
            xisave=xi_i_fks_ev
            ysave=y_ij_fks_ev
         endif

c Plot observables for event
         plot_wgt=ev_wgt*fkssymmetryfactor*vegaswgt
         if(abs(plot_wgt).gt.1.d-20.and.pp(0,1).ne.-99d0)
     &        call outfun(pp,ybst_til_tolab,plot_wgt,iplot_ev)
c Plot observables for counterevents and Born
         plot_wgt=( cnt_wgt*fkssymmetryfactor +
     &              cnt_swgt*fkssymmetryfactor +
     &              (bsv_wgt-born_wgt)*fkssymmetryfactorBorn +
     &              deg_wgt*fkssymmetryfactorDeg +
     &              deg_swgt*fkssymmetryfactorDeg )*vegaswgt + virt_wgt
         if(abs(plot_wgt).gt.1.d-20) then
            if(iplot.eq.-3)then
               write(*,*)'Error #1 in dsig'
               stop
            elseif (p1_cnt(0,1,iplot).ne.-99d0)then
               call outfun(p1_cnt(0,1,iplot),ybst_til_tolab,plot_wgt,
     &              iplot_cnt)
            endif
         endif
c Plot observables for Born; pass cnt momenta assuming they are
c identical to Born ones
         plot_wgt=born_wgt*fkssymmetryfactorBorn*vegaswgt
         if(abs(plot_wgt).gt.1.d-20) then
            if(iplot.eq.-3)then
               write(*,*)'Error #1 in dsig'
               stop
            elseif (p1_cnt(0,1,iplot).ne.-99d0)then
               call outfun(p1_cnt(0,1,iplot),ybst_til_tolab,plot_wgt,
     &              iplot_born)
            endif
         endif
      elseif (iminmax.eq.1 .and. ExceptPSpoint) then
c for except PS points, this is the maximal approx for the virtual         
         call unweight_function(p_born,unwgtfun)
         dsig_max = ((ev_wgt+cnt_wgt)*fkssymmetryfactor +
     &        cnt_swgt*fkssymmetryfactor +
     &        bsv_wgt*fkssymmetryfactorBorn +
     &        deg_wgt*fkssymmetryfactorDeg +
     &        deg_swgt*fkssymmetryfactorDeg)*unwgtfun+virt_wgt
         total_wgt_sum_max=total_wgt_sum_max+
     &        ((dsig_max - central_wgt_saved)*vegaswgt)**2

      elseif (iminmax.eq.2 .and. ExceptPSpoint) then
c for except PS points, this is the minimal approx for the virtual         
         call unweight_function(p_born,unwgtfun)
         dsig_min = ((ev_wgt+cnt_wgt)*fkssymmetryfactor +
     &        cnt_swgt*fkssymmetryfactor +
     &        bsv_wgt*fkssymmetryfactorBorn +
     &        deg_wgt*fkssymmetryfactorDeg +
     &        deg_swgt*fkssymmetryfactorDeg)*unwgtfun+virt_wgt
         total_wgt_sum_min=total_wgt_sum_min+
     &        ((central_wgt_saved - dsig_min)*vegaswgt)**2
      else
         write (*,*) 'Error #12 in dsig',iminmax
         stop
      endif

c If exceptional PS point found, go back to beginning recompute
c the weight for this PS point using an approximation
c based on previous PS points (done in BinothLHA.f)
      if (ExceptPSpoint .and. iminmax.le.1) goto 44

      return
      end


      subroutine dsigF(pp,wgt,vegaswgt,dsigS,dsigH)
c Here are the subtraction terms, the Sij function, 
c the f-damping function, and the single diagram
c enhanced multi-channel factor included
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      include "fks_powers.inc"
      include "madfks_mcatnlo.inc"
      include 'coupl.inc'
      include 'run.inc'
      include 'q_es.inc'
      include 'nFKSconfigs.inc'
      include 'reweight_all.inc'

c     timing statistics
      include 'timing_variables.inc'
      real deltaTOLP,deltaTPDF,deltaTFJ

      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS

      double precision pp(0:3,nexternal),wgt,vegaswgt

      double precision fks_Sij,fks_Hij,f_damp,dot,dlum
      external fks_Sij,fks_Hij,f_damp,dot,dlum

      double precision x,xtot,s_ev,s_c,s_s,s_sc,ffact,fx_c, fx_s,fx_sc
     &     ,Sxmc_wgt,Hxmc_wgt,cnt_wgt_c,cnt_wgt_s,cnt_wgt_sc, bsv_wgt
     &     ,plot_wgt ,cnt_swgt_s,cnt_swgt_sc,cnt_sc,cnt_s,
     &     prefact_cnt_ssc ,prefact_cnt_ssc_c,prefact_coll,
     &     prefact_coll_c,born_wgt ,prefact_deg,prefact,prefact_c,
     &     prefact_deg_sxi ,prefact_deg_slxi,deg_wgt,deg_swgt, deg_xi_c
     &     ,deg_lxi_c ,deg_xi_sc,deg_lxi_sc, cnt_swgt,cnt_wgt,gfactsf
     &     ,gfactcl,xmcMC ,xmcME,SxmcMC,SxmcME,HxmcMC,HxmcME, xlum_c
     &     ,xlum_s,xlum_sc ,xlum_mc,xlum_mc_save, dummy,Sev_wgt,Hev_wgt
     &     ,fx_ev,probne ,sevmc ,xlum_ev,get_ptrel, xlum_mc_fact,xnormsv
     &     ,xsec,bpower,cpower ,dsigS ,dsigH,totH_wgt,virt_wgt
      integer i,j

      integer izero,ione,itwo,mohdr,iplot_ev,iplot_cnt,iplot_born
      integer ithree,ifour
      integer ifill1H,ifill2H,ifill3H,ifill4H,ifill1S,ifill2S,ifill3S
     &     ,ifill4S,ifill2S_born
      save ifill2S_born
      parameter (izero=0)
      parameter (ione=1)
      parameter (itwo=2)
      parameter (ithree=3)
      parameter (ifour=4)
      parameter (mohdr=-100)
      parameter (iplot_ev=11)
      parameter (iplot_cnt=12)
      parameter (iplot_born=20)

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      double precision shattmp
      double precision pi
      parameter (pi=3.1415926535897932385d0)

      logical nocntevents
      common/cnocntevents/nocntevents

c Multi channel stuff:
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_amps/  amp2,       jamp2

      integer          isum_hel
      logical                    multi_channel
      common/to_matrix/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig

      double precision wgt1
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn

      double precision ev_enh,enhance,rwgt,unwgtfun,enhanceS,enhanceH
      logical firsttime,passcuts
      data firsttime /.true./
      integer inoborn_ev,inoborn_cnt
      double precision xnoborn_ev,xnoborn_cnt

c FKS stuff:
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision xi_i_fks_ev,y_ij_fks_ev
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

      double precision xi_i_fks_cnt(-2:2)
      common /cxiifkscnt/xi_i_fks_cnt

      double precision xiimax_ev
      common /cxiimaxev/xiimax_ev
      double precision xiimax_cnt(-2:2)
      common /cxiimaxcnt/xiimax_cnt

      double precision xinorm_ev
      common /cxinormev/xinorm_ev
      double precision xinorm_cnt(-2:2)
      common /cxinormcnt/xinorm_cnt

      double precision p1_cnt(0:3,nexternal,-2:2)
      double precision wgt_cnt(-2:2)
      double precision pswgt_cnt(-2:2)
      double precision jac_cnt(-2:2)
      common/counterevnts/p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt

      double precision zero,one
      parameter (zero=0d0,one=1d0)
      double precision tiny
      parameter (tiny=1d-6)

      double precision xicut_used
      common /cxicut_used/xicut_used
      double precision delta_used
      common /cdelta_used/delta_used
      double precision xiScut_used,xiBSVcut_used
      common /cxiScut_used/xiScut_used,xiBSVcut_used
      double precision fkssymmetryfactor,fkssymmetryfactorBorn,
     &     fkssymmetryfactorDeg
      integer ngluons,nquarks(-6:6),nphotons
      common/numberofparticles/fkssymmetryfactor,fkssymmetryfactorBorn,
     &                  fkssymmetryfactorDeg,ngluons,nquarks,nphotons
      double precision diagramsymmetryfactor
      common /dsymfactor/diagramsymmetryfactor

      logical nbody
      common/cnbody/nbody

      character*4 abrv
      common /to_abrv/ abrv

      double precision vegas_weight
      common/cvegas_weight/vegas_weight

      double precision zhw_used
c MC stuff
      double precision zhw(nexternal),xmcxsec(nexternal)
      integer nofpartners
      logical lzone(nexternal),flagmc
      logical MCcntcalled

      double precision emsca,scalemin,scalemax,emsca_bare
      logical emscasharp
      common/cemsca/emsca,emsca_bare,emscasharp,scalemin,scalemax

c Stuff to be written (depending on AddInfoLHE) onto the LHE file
      integer iSorH_lhe,ifks_lhe(fks_configs) ,jfks_lhe(fks_configs)
     &     ,fksfather_lhe(fks_configs) ,ipartner_lhe(fks_configs)
      double precision scale1_lhe(fks_configs),scale2_lhe(fks_configs)
      common/cto_LHE1/iSorH_lhe,ifks_lhe,jfks_lhe,
     #                fksfather_lhe,ipartner_lhe
      common/cto_LHE2/scale1_lhe,scale2_lhe

c CKWW scale, from cuts.f
      double precision scale_CKKW
      common/cscale_CKKW/scale_CKKW

c For the MINT folding
      integer fold
      common /cfl/fold

      logical unwgt
      double precision evtsgn
      common /c_unwgt/evtsgn,unwgt

c For plots
      logical plotEv,plotKin
      common/cEvKinplot/plotEv,plotKin

c For tests
      real*8 fksmaxwgt,xisave,ysave
      common/cfksmaxwgt/fksmaxwgt,xisave,ysave
      logical ExceptPSpoint
      integer iminmax
      common/cExceptPSpoint/iminmax,ExceptPSpoint


      double precision central_wgt_saved
      save central_wgt_saved

      double precision dsigS_max,dsigS_min
      double precision total_wgt_sum,total_wgt_sum_max,
     &                 total_wgt_sum_min
      common/csum_of_wgts/total_wgt_sum,total_wgt_sum_max,
     &                 total_wgt_sum_min

      double precision xm12
      integer ileg
      common/cscaleminmax/xm12,ileg

      double precision ximin
      parameter(ximin=0.05d0)

c
c This is the table that will be used to unweight. (It contains for
c arguments, 1st argument: nFKSproces; 2nd argument: S or H events; 3rd
c argument: IPROC (from parton luminosities))
c
      double precision unwgt_table(0:fks_configs,3,maxproc)
      common/c_unwgt_table/unwgt_table
      INTEGER              IPROC
      DOUBLE PRECISION PD(0:MAXPROC)
      COMMON /SUBPROC/ PD, IPROC
      DOUBLE PRECISION       CONV
      PARAMETER (CONV=389379660d0)  !CONV TO PICOBARNS
      integer i_process
      common/c_addwrite/i_process
      integer iproc_save(fks_configs),eto(maxproc,fks_configs)
     $     ,etoi(maxproc,fks_configs),maxproc_found
      common/cproc_combination/iproc_save,eto,etoi,maxproc_found
c$$$c FxFx merging
c$$$      logical rewgt_mohdr_calculated,rewgt_izero_calculated
c$$$      double precision rewgt_mohdr,rewgt_izero,rewgt_exp_mohdr
c$$$     $     ,rewgt_exp_izero
c$$$      logical setclscales
c$$$      double precision rewgt
c$$$      external setclscales,rewgt

      double precision rwgt_muR_dep_fac
      double precision ddum(6)
      logical ldum

      double precision pmass(nexternal)
      include "pmass.inc"

      vegas_weight=vegaswgt

c If there was an exceptional phase-space point found for the 
c virtual corrections, at the end of this subroutine, goto 44
c and compute also the "possible" minimal and maximal weight
c these points could have gotton (based upon previous PS
c points)
      ExceptPSpoint=.false.
      iminmax=-1
 44   continue
      iminmax=iminmax+1

      emsca=0.d0
      scalemax=0.d0

      Sev_wgt=0.d0
      Sxmc_wgt=0.d0
      cnt_wgt=0.d0
      cnt_wgt_s=0.d0
      cnt_wgt_c=0.d0
      cnt_wgt_sc=0.d0
      bsv_wgt=0.d0
      virt_wgt=0d0
      born_wgt=0.d0
      cnt_swgt=0.d0
      cnt_swgt_s=0.d0
      cnt_swgt_sc=0.d0
      deg_wgt=0.d0
      deg_swgt=0.d0
      plot_wgt=0.d0
c
      Hev_wgt=0.d0
      Hxmc_wgt=0.d0
c
      dsigS=0d0
      dsigH=0d0
      MCcntcalled=.false.
c
c$$$c FxFx merging
c$$$      ktscheme=1
c$$$      rewgt_mohdr_calculated=.false.
c$$$      rewgt_izero_calculated=.false.

      if(doreweight)then
        if(.not.AddInfoLHE)then
          write(*,*)'Error in dsigF'
          write(*,*)'  AddInfoLHE must be true when unweighting'
          stop
        endif
        call reweight_settozero()
        call reweight_settozero_all(nFKSprocess*2,nbody)
        call reweight_settozero_all(nFKSprocess*2-1,nbody)
        ifill1H=0
        ifill2H=0
        ifill3H=0
        ifill4H=0
        ifill1S=0
        ifill2S=0
        ifill3S=0
        ifill4S=0
        if (nbody) ifill2S_Born=0
      endif
      if(AddInfoLHE)then
        ifks_lhe(nFKSprocess)=i_fks
        jfks_lhe(nFKSprocess)=j_fks
        fksfather_lhe(nFKSprocess)=0
        ipartner_lhe(nFKSprocess)=0
        scale1_lhe(nFKSprocess)=0.d0
        scale2_lhe(nFKSprocess)=0.d0
      endif
c Set the upper value of the shower scale for the H and S events,
c respectively
      if (ickkw.ne.3 .and. ickkw.ne.4) then
         call set_cms_stuff(mohdr)
         call set_shower_scale_noshape(pp,nFKSprocess*2)
         call set_cms_stuff(izero)
         call set_shower_scale_noshape(pp,nFKSprocess*2-1)
      endif
c
c Make sure that the result can be non-zero. If the jacobian from the
c PS-setup or vegas are zero, we can skip this PS point and 'return'.
c Note that all the wgts and jacs should be positive.
      if ( (wgt.le.0d0 .and. jac_cnt(0).le.0d0 .and. jac_cnt(1).le.0d0
     &     .and. jac_cnt(2).le.0d0) .or. vegaswgt.le.0d0) return
c
      if (fold.eq.0) then
         calculatedBorn=.false.
         call get_helicity(i_fks,j_fks)
      endif
c
      if (firsttime)then
         inoborn_ev=0
         xnoborn_ev=0.d0
         inoborn_cnt=0
         xnoborn_cnt=0.d0
         fksmaxwgt=0.d0
         firsttime=.false.
c Put here call to compute bpower
         call compute_bpower(p_born,bpower)
         wgtbpower=bpower

c Compute cpower done for bottom Yukawa, routine needs to be adopted
c for other muR-dependendent factors
         call compute_cpower(p_born,cpower)
         if(dabs(cpower+1d0).lt.tiny) then
            wgtcpower=0d0
         else
            wgtcpower=cpower
         endif
c Check that things are done consistently
        if(wgtcpower.ne.cpowerinput.and.dabs(cpower+1d0).gt.tiny)then
           write(*,*)'Inconsistency in the computation of cpower',
     #               wgtcpower,cpowerinput
           write(*,*)'Check value in reweight0.inc'
           stop
         endif
      endif
c
      prefact=xinorm_ev/xi_i_fks_ev*
     #        1/(1-y_ij_fks_ev)
c
      if( (.not.nocntevents) .and. (.not.(abrv.eq.'born' .or. abrv.eq
     &     .'grid' .or. abrv(1:2).eq.'vi' .or. nbody))
     &     )then
        prefact_cnt_ssc=xinorm_ev/min(xiimax_ev,xiScut_used)*
     #                  log(xicut_used/min(xiimax_ev,xiScut_used))*
     #                  1/(1-y_ij_fks_ev)
        if(pmass(j_fks).eq.0.d0)then
          prefact_c=xinorm_cnt(ione)/xi_i_fks_cnt(ione)*
     #              1/(1-y_ij_fks_ev)
          prefact_cnt_ssc_c=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                      log(xicut_used/min(xiimax_cnt(ione),xiScut_used))*
     #                      1/(1-y_ij_fks_ev)
          prefact_coll=xinorm_cnt(ione)/xi_i_fks_cnt(ione)*
     #                 log(delta_used/deltaS)/deltaS
          prefact_coll_c=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                   log(xicut_used/min(xiimax_cnt(ione),xiScut_used))*
     #                   log(delta_used/deltaS)/deltaS
          prefact_deg=xinorm_cnt(ione)/xi_i_fks_cnt(ione)*
     #                1/deltaS
          prefact_deg_sxi=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                    log(xicut_used/min(xiimax_cnt(ione),xiScut_used))*
     #                    1/deltaS
          prefact_deg_slxi=xinorm_cnt(ione)/min(xiimax_cnt(ione),xiScut_used)*
     #                     ( log(xicut_used)**2 -
     #                       log(min(xiimax_cnt(ione),xiScut_used))**2 )*
     #                     1/(2.d0*deltaS)
        endif
      endif
c
      probne=1.d0
c All counterevent have the same final-state kinematics. Check that
c one of them passes the hard cuts, and they exist at all
c
c Set the ybst_til_tolab before applying the cuts. Update below
c for the collinear, soft and/or soft-collinear subtraction terms
      call set_cms_stuff(izero)
      if ( (.not.passcuts(p1_cnt(0,1,0),rwgt)) .or.
     #      nocntevents ) goto 547

c$$$c Compute the scales and sudakov-reweighting for the FxFx merging
c$$$      if (ickkw.eq.3) then
c$$$         if (.not. setclscales(p1_cnt(0,1,0))) then
c$$$            write (*,*) 'ERROR in setclscales izero'
c$$$            stop
c$$$         endif
c$$$         rewgt_izero=rewgt(p1_cnt(0,1,0),rewgt_exp_izero)
c$$$         rewgt_izero_calculated=.true.
c$$$c Set the upper value of the shower scale for the H and S events,
c$$$c respectively
c$$$         call set_cms_stuff(izero)
c$$$         call set_shower_scale_noshape(pp,nFKSprocess*2-1)
c$$$      endif

      gfactsf=1.d0
      gfactcl=1.d0
      sevmc=1.d0
      xmcMC=0.d0
      xmcME=0.d0
      SxmcMC=0.d0
      SxmcME=0.d0
      HxmcMC=0.d0
      HxmcME=0.d0
c
      deltaTOLP = tOLP
      deltaTPDF = tPDF
      deltaTFJ  = tFastJet
      call cpu_time(tBefore)
c
      if (abrv.eq.'born' .or. abrv.eq.'grid' .or.
     &     abrv(1:2).eq.'vi' .or. nbody) goto 540

c For UNLOPS, skip MC counter events
      if (ickkw.eq.4) goto 540

      call set_cms_stuff(mohdr)
c$$$c     Compute the scales and sudakov-reweighting for the FxFx merging
c$$$      if (ickkw.eq.3) then
c$$$         if (.not. setclscales(pp)) then
c$$$            write (*,*) 'ERROR in setclscales mohdr'
c$$$            stop
c$$$         endif
c$$$         rewgt_mohdr=rewgt(pp,rewgt_exp_mohdr)
c$$$         rewgt_mohdr_calculated=.true.
c$$$c Set the upper value of the shower scale for the H and S events,
c$$$c respectively
c$$$         call set_cms_stuff(mohdr)
c$$$         call set_shower_scale_noshape(pp,nFKSprocess*2)
c$$$      endif
      call set_alphaS(pp)
      if(doreweight)then
         wgtmuR2_all(1,nFKSprocess*2)=muR2_current/muR_over_ref**2
         wgtmuF12_all(1,nFKSprocess*2)=muF12_current/muF1_over_ref**2
         wgtmuF22_all(1,nFKSprocess*2)=muF22_current/muF2_over_ref**2
         call reweight_fillkin_all(pp,ione,nFKSprocess*2)
         ifill1H=1
         call reweight_fillkin_all(pp,itwo,nFKSprocess*2)
         ifill2H=1
         wgtmuR2_all(1,nFKSprocess*2-1)=muR2_current/muR_over_ref**2
         wgtmuF12_all(1,nFKSprocess*2-1)=muF12_current/muF1_over_ref**2
         wgtmuF22_all(1,nFKSprocess*2-1)=muF22_current/muF2_over_ref**2
         call reweight_fillkin_all(pp,ione,nFKSprocess*2-1)
         ifill1S=1
      endif
      if(UseSfun)then
         x = abs(2d0*dot(pp(0,i_fks),pp(0,j_fks))/shat)
         ffact = f_damp(x)
         sevmc = fks_Sij(pp,i_fks,j_fks,xi_i_fks_ev,y_ij_fks_ev)
         sevmc = sevmc*ffact
      else
         x = abs(2d0*dot(pp(0,i_fks),pp(0,j_fks))/shat)
         ffact = f_damp(x)
         sevmc = fks_Hij(pp,i_fks,j_fks)
         sevmc = sevmc*ffact
      endif

c$$$      if (ickkw.eq.3 .and. .not. (rewgt_mohdr_calculated .and.
c$$$     $     rewgt_izero_calculated) )then
c$$$         write (*,*)'Both shower scales should be set before'/
c$$$     $        /' entering the MC subtraction terms'
c$$$         stop
c$$$      endif

      call xmcsubt(pp,xi_i_fks_ev,y_ij_fks_ev,gfactsf,gfactcl,probne,
     #             dummy,nofpartners,lzone,flagmc,zhw,xmcxsec)
      MCcntcalled=.true.

      if(ileg.gt.4.or.ileg.lt.1)then
         write(*,*)'Error: unrecognized ileg in dsigF', ileg
         stop
      endif

      if(sevmc.gt.0.d0.and.flagmc)then
         if(doreweight)then
            iwgtnumpartn_all(nFKSprocess*2)=nofpartners
            iwgtnumpartn_all(nFKSprocess*2-1)=nofpartners
            xsec=sevmc*wgt*prefact*rwgt
         endif
        xlum_mc_save=-1.d8
        do i=1,nofpartners
          if(lzone(i))then
            zhw_used=zhw(i)
            call get_mc_lum(j_fks,zhw_used,xi_i_fks_ev,
     #                      xlum_mc_save,xlum_mc,xlum_mc_fact)
            xmcMC=xmcMC+xmcxsec(i)*xlum_mc
            do j=1,IPROC
               unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1
     &              ,j)+xmcxsec(i)*PD(j)*xlum_mc_fact*sevmc*wgt*prefact
     &              *rwgt*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
               unwgt_table(nFKSprocess,2,j)=unwgt_table(nFKSprocess,2
     &              ,j)-xmcxsec(i)*PD(j)*xlum_mc_fact*sevmc*wgt*prefact
     &              *rwgt*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            enddo
            if(doreweight)then
               wgtwmcxsec_all(i,nFKSprocess*2)=-xsec*xlum_mc_fact
     &              *xmcxsec(i)/g**(nint(2*wgtbpower+2.d0))
               wgtmcxbj_all(1,i,nFKSprocess*2)=xbk(1)
               wgtmcxbj_all(2,i,nFKSprocess*2)=xbk(2)
               wgtwmcxsec_all(i,nFKSprocess*2-1)=xsec*xlum_mc_fact
     &              *xmcxsec(i)/g**(nint(2*wgtbpower+2.d0))
               wgtmcxbj_all(1,i,nFKSprocess*2-1)=xbk(1)
               wgtmcxbj_all(2,i,nFKSprocess*2-1)=xbk(2)
            endif
          endif
        enddo
        SxmcMC=xmcMC*sevmc*wgt*prefact*rwgt * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
        HxmcMC=-xmcMC*sevmc*wgt*prefact*rwgt * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
      endif
c
      if( (.not.flagmc).and.gfactsf.eq.1.d0 .and.
     #   xi_i_fks_ev.lt.0.02d0  .and. particle_type(i_fks).eq.8)then
        write(*,*)'Error in dsigF: will diverge'
        stop
      endif

 540  continue

c Set scales for all counterevents, using soft kinematics as done
c in the case of parton-level NLO computations

      call set_cms_stuff(izero)
      if ( (.not.passcuts(p1_cnt(0,1,0),rwgt)) .or.
     #      nocntevents ) goto 547
c$$$c Compute the scales and sudakov-reweighting for the FxFx merging
c$$$      if (ickkw.eq.3) then
c$$$         if (.not. setclscales(p1_cnt(0,1,0))) then
c$$$            write (*,*) 'ERROR in setclscales izero'
c$$$            stop
c$$$         endif
c$$$         rewgt_izero=rewgt(p1_cnt(0,1,0),rewgt_exp_izero)
c$$$         rewgt_izero_calculated=.true.
c$$$c Set the upper value of the shower scale for the H and S events,
c$$$c respectively
c$$$         call set_cms_stuff(izero)
c$$$         call set_shower_scale_noshape(pp,nFKSprocess*2-1)
c$$$      endif
      call set_alphaS(p1_cnt(0,1,0))
      if(doreweight)then
         if (nbody) then
            wgtqes2_all(2,0)=QES2
            wgtqes2_all(3,0)=QES2
            wgtqes2_all(4,0)=QES2
            wgtmuR2_all(2,0)=muR2_current/muR_over_ref**2
            wgtmuF12_all(2,0)=muF12_current/muF1_over_ref**2
            wgtmuF22_all(2,0)=muF22_current/muF2_over_ref**2
            call reweight_fillkin_all(pp,itwo,0)
            ifill2S_born=1
         else
            wgtqes2_all(2,nFKSprocess*2-1)=QES2
            wgtqes2_all(3,nFKSprocess*2-1)=QES2
            wgtqes2_all(4,nFKSprocess*2-1)=QES2
            wgtmuR2_all(2,nFKSprocess*2-1)=muR2_current/muR_over_ref**2
            wgtmuF12_all(2,nFKSprocess*2-1)=muF12_current/muF1_over_ref**2
            wgtmuF22_all(2,nFKSprocess*2-1)=muF22_current/muF2_over_ref**2
            call reweight_fillkin_all(pp,itwo,nFKSprocess*2-1)
            ifill2S=1
         endif
      endif
      if (abrv.eq.'born' .or. abrv.eq.'grid' .or. abrv(1:2).eq.'vi' .or.
     &     nbody)goto 545
c
c Collinear subtraction term:
      if( ( y_ij_fks_ev.gt.1d0-deltaS .or. 
     #     (gfactsf.lt.1.d0.and.gfactcl.lt.1.d0 .and.
     #      probne.gt.0.d0) ) .and.
     #    pmass(j_fks).eq.0.d0 )then
         call set_cms_stuff(ione)

         if(doreweight)then
            if(gfactsf.lt.1.d0.and.probne.gt.0.d0.and.gfactcl.lt.1.d0
     &           .and.pmass(j_fks).eq.0.d0)then
               wgtxbj_all(1,3,nFKSprocess*2)=xbk(1)
               wgtxbj_all(2,3,nFKSprocess*2)=xbk(2)
               call reweight_fillkin_all(pp,ithree,nFKSprocess*2)
               ifill3H=1
            endif
            wgtxbj_all(1,3,nFKSprocess*2-1)=xbk(1)
            wgtxbj_all(2,3,nFKSprocess*2-1)=xbk(2)
            call reweight_fillkin_all(pp,ithree,nFKSprocess*2-1)
            ifill3S=1
         endif
         s_c = fks_Sij(p1_cnt(0,1,1),i_fks,j_fks,xi_i_fks_cnt(ione),one)
         if(s_c.gt.0.d0)then
            if(abs(s_c-1.d0).gt.1.d-6.and.j_fks.le.nincoming)then
               write(*,*)'Wrong S function in dsigF[c]',s_c
               stop
            endif
            call sreal(p1_cnt(0,1,1),xi_i_fks_cnt(ione),one,fx_c)
            xlum_c = dlum()
            xsec = fx_c*s_c*jac_cnt(1)*prefact_c*rwgt*(1-gfactcl)
            SxmcME=SxmcME+xlum_c*xsec
            do j=1,IPROC
               unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1
     &              ,j)+xsec*PD(j)*(1-gfactsf)*probne*CONV
     f              * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here          
            enddo
            if ((gfactsf.lt.1.d0.and.gfactcl.lt.1.d0 .and.
     &           probne.gt.0.d0) .and. pmass(j_fks).eq.0.d0) then
               HxmcME=HxmcME+xlum_c*xsec
               do j=1,IPROC
                  unwgt_table(nFKSprocess,2,j)=unwgt_table(nFKSprocess
     &                 ,2,j)-xsec*PD(j)*(1-gfactsf)*probne*CONV
     f                 * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here             
               enddo
            endif
            if(doreweight) then
               wgtwreal_all(3,nFKSprocess*2-1)=xsec*(1-gfactsf)*probne/
     &              g**(nint(2*wgtbpower+2.d0))
               if(gfactsf.lt.1.d0.and.probne.gt.0.d0.and.gfactcl.lt.1.d0
     $              .and.pmass(j_fks).eq.0.d0) wgtwreal_all(3
     $              ,nFKSprocess*2)= xsec/g**(nint(2*wgtbpower+2.d0))
            endif
            if( y_ij_fks_ev.gt.1d0-deltaS )then
               xsec = fx_c*s_c*jac_cnt(1)*(prefact_c+prefact_coll)*rwgt
               cnt_wgt_c=cnt_wgt_c-xlum_c*xsec * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
               do j=1,IPROC
                  unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess
     &                 ,1,j)-xsec*PD(j)*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
               enddo
               call sreal_deg(p1_cnt(0,1,1),xi_i_fks_cnt(ione),one,
     #                        deg_xi_c,deg_lxi_c)
               deg_wgt=deg_wgt+( deg_xi_c+deg_lxi_c*log(xi_i_fks_cnt(ione)) )*
     #                         jac_cnt(1)*prefact_deg*rwgt/(shat/(32*pi**2))*
     #                         xlum_c * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
               do j=1,IPROC
                  unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess
     &                 ,1,j)+PD(j)*( deg_xi_c+deg_lxi_c
     &                 *log(xi_i_fks_cnt(ione)) )* jac_cnt(1)
     &                 *prefact_deg*rwgt/(shat/(32*pi**2))*CONV
     f                 * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here             
               enddo
               if(doreweight)then
                  wgtwreal_all(3,nFKSprocess*2-1)=wgtwreal_all(3
     &                 ,nFKSprocess*2-1)-xsec/g**(nint(2*wgtbpower
     &                 +2.d0))
                  wgtwdeg_all(3,nFKSprocess*2-1)=( wgtdegrem_xi+
     &                 wgtdegrem_lxi*log(xi_i_fks_cnt(ione)) )*
     &                 jac_cnt(1)*prefact_deg*rwgt/(shat/(32*pi**2))/
     &                 g**(nint(2*wgtbpower+2.d0))
                  wgtwdegmuf_all(3,nFKSprocess*2-1)=wgtdegrem_muF *
     &                 jac_cnt(1)*prefact_deg*rwgt/(shat/(32*pi**2))/
     &                 g**(nint(2*wgtbpower+2.d0))
               endif
            endif
         endif
      endif
c Soft subtraction term:
 545  continue
      if ( xi_i_fks_ev .lt. max(xiScut_used,xiBSVcut_used) .or.
     &     (gfactsf.lt.1.d0.and.probne.gt.0.d0) ) then
         call set_cms_stuff(izero)
         if(doreweight)then
            if(gfactsf.lt.1.d0.and.probne.gt.0.d0)then
               wgtxbj_all(1,2,nFKSprocess*2)=xbk(1)
               wgtxbj_all(2,2,nFKSprocess*2)=xbk(2)
               wgtmuR2_all(2,nFKSprocess*2)=muR2_current/muR_over_ref**2
               wgtmuF12_all(2,nFKSprocess*2)=muF12_current/muF1_over_ref**2
               wgtmuF22_all(2,nFKSprocess*2)=muF22_current/muF2_over_ref**2
               if(ifill2H.ne.1)then
                  write(*,*)'Error #1a[wg] in dsigF',ifill2H
                  stop
               endif
            endif
            if (nbody) then
               wgtxbj_all(1,2,0)=xbk(1)
               wgtxbj_all(2,2,0)=xbk(2)
               if(ifill2S_born.ne.1)then
                  write(*,*)'Error #1b[wg] in dsigF',ifill2S
                  stop
               endif
            else
               wgtxbj_all(1,2,nFKSprocess*2-1)=xbk(1)
               wgtxbj_all(2,2,nFKSprocess*2-1)=xbk(2)
               if(ifill2S.ne.1)then
                  write(*,*)'Error #1c[wg] in dsigF',ifill2S
                  stop
               endif
            endif
         endif

         s_s = fks_Sij(p1_cnt(0,1,0),i_fks,j_fks,zero,y_ij_fks_ev)
         if(nbody)s_s=1.d0
         if(s_s.gt.0.d0)then
            xlum_s = dlum()
            if (abrv.eq.'born' .or. abrv.eq.'grid' .or.
     &           abrv(1:2).eq.'vi' .or. nbody) goto 546
            call sreal(p1_cnt(0,1,0),zero,y_ij_fks_ev,fx_s)
            xsec=fx_s*s_s*jac_cnt(0)*prefact*rwgt
            SxmcME=SxmcME+xlum_s*xsec
            do j=1,IPROC
               unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1
     &              ,j)+xsec*PD(j)*(1-gfactsf)*probne*CONV
     f              * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here          
            enddo
            if(doreweight)then
               wgtwreal_all(2,nFKSprocess*2-1)=xsec*(1-gfactsf)*probne/
     &              g**(nint(2*wgtbpower+2.d0))
             endif
            if (gfactsf.lt.1.d0.and.probne.gt.0.d0) then
               HxmcME=HxmcME+xlum_s*xsec
               do j=1,IPROC
                  unwgt_table(nFKSprocess,2,j)=unwgt_table(nFKSprocess
     &                 ,2,j)-xsec*PD(j)*(1-gfactsf)*probne*CONV
     f                 * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
               enddo
               if(doreweight)wgtwreal_all(2,nFKSprocess*2)=
     &              xsec/g**(nint(2*wgtbpower+2.d0))
            endif
            if (xi_i_fks_ev .lt. xiScut_used) then
              xsec=fx_s*s_s*jac_cnt(0)
              cnt_s=xlum_s*xsec
              cnt_wgt_s=cnt_wgt_s-cnt_s*prefact*rwgt
     f             * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              cnt_swgt_s=cnt_swgt_s-cnt_s*prefact_cnt_ssc*rwgt
     f             * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              do j=1,IPROC
                 unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1
     $                ,j)-PD(j)*xsec*(prefact+prefact_cnt_ssc)*rwgt
     $                *CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              enddo
              if(doreweight)wgtwreal_all(2,nFKSprocess*2-1)
     &             =wgtwreal_all(2,nFKSprocess*2-1)-xsec*(prefact
     &             +prefact_cnt_ssc)*rwgt/g**(nint(2*wgtbpower+2.d0))
            endif
 546        continue
            if (abrv.eq.'real' .or. .not.nbody) goto 548
            if (xi_i_fks_ev .lt. xiBSVcut_used) then
              xsec=s_s*jac_cnt(0)*xinorm_ev/
     #             (min(xiimax_ev,xiBSVcut_used)*shat/(16*pi**2))*
     #             rwgt
              xnormsv=xlum_s*xsec
              call bornsoftvirtual(p1_cnt(0,1,0),bsv_wgt,virt_wgt
     $             ,born_wgt)
c$$$c For FxFx merging, include the compensation term
c$$$              if (rewgt_izero_calculated.and.rewgt_izero.lt.1d0) then
c$$$                 bsv_wgt=bsv_wgt-g**2/(4d0*Pi)*rewgt_exp_izero*born_wgt
c$$$              endif
              if(doreweight)then
                 if(wgtbpower.gt.0)then
                    wgtwborn_all=born_wgt*xsec/g**(nint(2*wgtbpower))
                 else
                    wgtwborn_all=born_wgt*xsec
                 endif
                 wgtwns_all=wgtnstmp*xsec/g**(nint(2*wgtbpower+2.d0))
c$$$                 if (rewgt_izero_calculated.and.rewgt_izero.lt.1d0) then
c$$$                    wgtwns_all=wgtwns_all-rewgt_exp_izero*wgtwborn_all
c$$$     $                   /(4d0*pi)
c$$$                 endif
                 wgtwnsmuf_all=wgtwnstmpmuf*xsec/g**(nint(2*wgtbpower
     &                +2.d0))
                 wgtwnsmur_all=wgtwnstmpmur*xsec/g**(nint(2*wgtbpower
     &                +2.d0))
              endif
              do j=1,IPROC
                 unwgt_table(0,1,j)=unwgt_table(0,1,j)+PD(j)*bsv_wgt
     &                *xsec*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
                 unwgt_table(0,3,j)=unwgt_table(0,3,j)+PD(j)*virt_wgt
     &                *xsec*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
                 unwgt_table(1,3,j)=unwgt_table(1,3,j)+PD(j)*born_wgt
     &                *xsec*CONV*g**2/(8d0*PI**2)
     f                * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              enddo
              bsv_wgt=bsv_wgt*xnormsv * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              virt_wgt=virt_wgt*xnormsv * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              born_wgt=born_wgt*xnormsv * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            endif
 548        continue
         endif
      endif
c Soft-Collinear subtraction term:
      if (abrv.eq.'born' .or. abrv.eq.'grid' .or. abrv(1:2).eq.'vi' .or.
     &     nbody)goto 547
      if ( ( (xi_i_fks_cnt(ione) .lt. xiScut_used .and.
     #        y_ij_fks_ev .gt. 1d0-deltaS) .or.
     #        (gfactsf.lt.1.d0.and.gfactcl.lt.1.d0 .and.
     #         probne.gt.0.d0) ) .and.
     #        pmass(j_fks).eq.0.d0 )then
         call set_cms_stuff(itwo)
         if(doreweight)then
            if ((gfactsf.lt.1.d0.and.gfactcl.lt.1.d0 .and.
     &           probne.gt.0.d0) .and. pmass(j_fks).eq.0.d0) then
               wgtxbj_all(1,4,nFKSprocess*2)=xbk(1)
               wgtxbj_all(2,4,nFKSprocess*2)=xbk(2)
               call reweight_fillkin_all(pp,ifour,nFKSprocess*2)
               ifill4H=1
            endif
            wgtxbj_all(1,4,nFKSprocess*2-1)=xbk(1)
            wgtxbj_all(2,4,nFKSprocess*2-1)=xbk(2)
            call reweight_fillkin_all(pp,ifour,nFKSprocess*2-1)
            ifill4S=1
         endif
         s_sc = fks_Sij(p1_cnt(0,1,2),i_fks,j_fks,zero,one)
         if(s_sc.gt.0.d0)then
            if(abs(s_sc-1.d0).gt.1.d-6.and.j_fks.le.nincoming)then
               write(*,*)'Wrong S function in dsigF[sc]',s_sc
               stop
            endif
            call sreal(p1_cnt(0,1,2),zero,one,fx_sc)
            xlum_sc = dlum()
            xsec = fx_sc*s_sc*jac_cnt(2)*prefact_c*rwgt*(1-gfactcl)
            SxmcME=SxmcME-xlum_sc*xsec
            do j=1,IPROC
               unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1
     &              ,j)-PD(j)*xsec*(1-gfactsf)*probne*CONV
     f              * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
            enddo
            if(doreweight)then
               wgtwreal_all(4,nFKSprocess*2-1)=-xsec*(1-gfactsf)*probne/
     &              g**(nint(2*wgtbpower+2.d0))
            endif
            if ((gfactsf.lt.1.d0.and.gfactcl.lt.1.d0 .and.
     &           probne.gt.0.d0) .and. pmass(j_fks).eq.0.d0)then
               HxmcME=HxmcME-xlum_sc*xsec
               do j=1,IPROC
                  unwgt_table(nFKSprocess,2,j)=unwgt_table(nFKSprocess
     &                 ,2,j)+PD(j)*xsec*(1-gfactsf)*probne*CONV
     f                 * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here             
               enddo
               if(doreweight)wgtwreal_all(4,nFKSprocess*2)=
     &              -xsec/g**(nint(2*wgtbpower+2.d0))
            endif

            if(xi_i_fks_cnt(ione) .lt. xiScut_used .and.
     #          y_ij_fks_ev .gt. 1d0-deltaS)then
              xsec=fx_sc*s_sc*jac_cnt(2)
              cnt_sc=xlum_sc*xsec
              cnt_wgt_sc=cnt_wgt_sc+cnt_sc*(prefact_c+prefact_coll)*rwgt
     f             * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              cnt_swgt_sc=cnt_swgt_sc+
     &             cnt_sc*(prefact_cnt_ssc_c+prefact_coll_c)*rwgt
     f             * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              call sreal_deg(p1_cnt(0,1,2),zero,one,
     #                       deg_xi_sc,deg_lxi_sc)
              deg_wgt=deg_wgt-
     #                    ( deg_xi_sc+deg_lxi_sc*log(xi_i_fks_cnt(ione)) )*
     #                    jac_cnt(2)*prefact_deg*rwgt/(shat/(32*pi**2))*
     #                    xlum_sc * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              deg_swgt=deg_swgt-( deg_xi_sc*prefact_deg_sxi +
     #                       deg_lxi_sc*prefact_deg_slxi )*
     #                       jac_cnt(2)*rwgt/(shat/(32*pi**2))*
     #                       xlum_sc * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              do j=1,IPROC
                 unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1
     &                ,j)+PD(j)*(xsec*(prefact_c+prefact_coll
     &                +prefact_cnt_ssc_c+prefact_coll_c)*rwgt-
     &                (deg_xi_sc+deg_lxi_sc*log(xi_i_fks_cnt(ione)) )*
     &                jac_cnt(2)*prefact_deg*rwgt/(shat/(32*pi**2))
     &                -(deg_xi_sc*prefact_deg_sxi + deg_lxi_sc
     &                *prefact_deg_slxi )* jac_cnt(2)*rwgt/(shat/(32*pi
     &                **2)))*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
              enddo
              if(doreweight)then
                 wgtwreal_all(4,nFKSprocess*2-1)=wgtwreal_all(4
     &                ,nFKSprocess*2-1)+xsec*(prefact_c+prefact_coll
     &                +prefact_cnt_ssc_c+prefact_coll_c)*rwgt/g**(nint(2
     &                *wgtbpower+2.d0))
                 wgtwdeg_all(4,nFKSprocess*2-1)=(-( wgtdegrem_xi
     &                +wgtdegrem_lxi*log(xi_i_fks_cnt(ione)) )
     &                *prefact_deg -( wgtdegrem_xi*prefact_deg_sxi
     &                +wgtdegrem_lxi*prefact_deg_slxi ) )* jac_cnt(2)
     &                *rwgt/(shat/(32*pi**2))/ g**(nint(2*wgtbpower
     &                +2.d0))
                 wgtwdegmuf_all(4,nFKSprocess*2-1)= -wgtdegrem_muF*(
     &                prefact_deg+prefact_deg_sxi )* jac_cnt(2)*rwgt
     &                /(shat/(32*pi**2))/ g**(nint(2*wgtbpower+2.d0))
              endif
           endif
        endif
      endif
      SxmcME=SxmcME*(1-gfactsf)*probne * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
      HxmcME=-HxmcME*(1-gfactsf)*probne * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
      if(doreweight)then
         xsec=(1-gfactsf)*probne
         wgtwreal_all(2,nFKSprocess*2)=
     &        -wgtwreal_all(2,nFKSprocess*2)*xsec
         wgtwreal_all(3,nFKSprocess*2)=
     &        -wgtwreal_all(3,nFKSprocess*2)*xsec
         wgtwreal_all(4,nFKSprocess*2)=
     &        -wgtwreal_all(4,nFKSprocess*2)*xsec
      endif
      Sxmc_wgt=Sxmc_wgt+SxmcMC+SxmcME
      Hxmc_wgt=Hxmc_wgt+HxmcMC+HxmcME
 547  continue
      call cpu_time(tAfter)
      tDSigI = tDSigI + (tAfter-tBefore) - (tOLP-deltaTOLP) - 
     &     (tPDF -deltaTPDF) - (tFastJet - deltaTFJ)

c Real contribution
c
c Set the ybst_til_tolab before applying the cuts. 
      if (abrv.eq.'born' .or. abrv.eq.'grid' .or. abrv(1:2).eq.'vi' .or.
     &     nbody)goto 550
      call cpu_time(tBefore)
      deltaTPDF = tPDF
      deltaTFJ  = tFastJet
c
      call set_cms_stuff(mohdr)
      if(doreweight)then
         wgtxbj_all(1,1,nFKSprocess*2)=xbk(1)
         wgtxbj_all(2,1,nFKSprocess*2)=xbk(2)
         wgtxbj_all(1,1,nFKSprocess*2-1)=xbk(1)
         wgtxbj_all(2,1,nFKSprocess*2-1)=xbk(2)
      endif
      if (passcuts(pp,rwgt)) then
c$$$c     Compute the scales and sudakov-reweighting for the FxFx merging
c$$$        if (ickkw.eq.3) then
c$$$           if (.not. setclscales(pp)) then
c$$$              write (*,*) 'ERROR in setclscales mohdr'
c$$$              stop
c$$$           endif
c$$$           rewgt_mohdr=rewgt(pp,rewgt_exp_mohdr)
c$$$           rewgt_mohdr_calculated=.true.
c$$$c Set the upper value of the shower scale for the H and S events,
c$$$c respectively
c$$$           call set_cms_stuff(mohdr)
c$$$           call set_shower_scale_noshape(pp,nFKSprocess*2)
c$$$        endif
        call set_alphaS(pp)
        x = abs(2d0*dot(pp(0,i_fks),pp(0,j_fks))/shat)
        ffact = f_damp(x)
        s_ev = fks_Sij(pp,i_fks,j_fks,xi_i_fks_ev,y_ij_fks_ev)
        if(s_ev.gt.0.d0)then
          call sreal(pp,xi_i_fks_ev,y_ij_fks_ev,fx_ev)
          xlum_ev = dlum()
          xsec = fx_ev*s_ev*ffact*wgt*prefact*rwgt
          Sev_wgt = xlum_ev*xsec*(1-probne) * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
          Hev_wgt = xlum_ev*xsec*probne * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
          do j=1,IPROC
             unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1,j)
     &            +PD(j)*xsec*(1-probne)*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
             if (ickkw.ne.4) then
                unwgt_table(nFKSprocess,2,j)=unwgt_table(nFKSprocess,2,j)
     &               +PD(j)*xsec*probne*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
             else         ! for UNLOPS, add the H-events to the S-events
                unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1,j)
     &               +PD(j)*xsec*probne*CONV * rwgt_muR_dep_fac(scale)
com-- muR-dependent fac is reweighted here
             endif
          enddo
          if(doreweight)then
             if(ifill1H.eq.0)then
                wgtmuR2_all(1,nFKSprocess*2)=
     &               muR2_current/muR_over_ref**2
                wgtmuF12_all(1,nFKSprocess*2)=
     &               muF12_current/muF1_over_ref**2
                wgtmuF22_all(1,nFKSprocess*2)=
     &               muF22_current/muF2_over_ref**2
                call reweight_fillkin_all(pp,ione,nFKSprocess*2)
                ifill1H=1
             endif
             wgtwreal_all(1,nFKSprocess*2)=
     &            xsec/g**(nint(2*wgtbpower+2.d0))*probne
             if(ifill1S.eq.0)then
                wgtmuR2_all(1,nFKSprocess*2-1)=
     &               muR2_current/muR_over_ref**2
                wgtmuF12_all(1,nFKSprocess*2-1)=
     &               muF12_current/muF1_over_ref**2
                wgtmuF22_all(1,nFKSprocess*2-1)=
     &               muF22_current/muF2_over_ref**2
                call reweight_fillkin_all(pp,ione,nFKSprocess*2-1)
                ifill1S=1
             endif
             wgtwreal_all(1,nFKSprocess*2-1)=
     &            xsec/g**(nint(2*wgtbpower+2.d0))*(1-probne)
          endif
        endif
        if(AddInfoLHE)scale2_lhe(nFKSprocess)=get_ptrel(pp,i_fks,j_fks)
      endif
      call cpu_time(tAfter)
      tDSigR = tDSigR + (tAfter-tBefore)-(tPDF-deltaTPDF)-
     &     (tFastJet-deltaTFJ)
 550  continue

      if( (.not.MCcntcalled) .and.
     &     abrv.ne.'born'.and. abrv.ne.'grid' .and. ickkw.ne.4)then
         if(pp(0,1).ne.-99d0)then
            call set_cms_stuff(mohdr)
            call assign_emsca(pp,xi_i_fks_ev,y_ij_fks_ev)
         endif
      endif

      if(AddInfoLHE.and.UseCKKW)then
         if(scale1_lhe(nFKSprocess).eq.0.d0)
     &        scale1_lhe(nFKSprocess)=scale2_lhe(nFKSprocess)
         scale2_lhe(nFKSprocess)=scale_CKKW
      endif
c
c Enhance the one channel for multi-channel integration
c
      enhance=1.d0
      if ((Sxmc_wgt.ne.0d0 .or. Hxmc_wgt.ne.0d0 .or. cnt_wgt_c.ne.0d0
     $     .or. cnt_wgt_s.ne.0d0 .or. cnt_wgt_sc.ne.0d0 .or.
     $     bsv_wgt.ne.0d0 .or. virt_wgt.ne.0d0 .or.
     $     deg_wgt.ne.0d0.or.deg_swgt.ne.0d0 .or. cnt_swgt_s.ne.0d0 .or.
     $     cnt_swgt_sc.ne.0d0 .or.Sev_wgt.ne.0d0 .or. Hev_wgt.ne.0d0)
     $     .and. multi_channel) then
         if(bsv_wgt.eq.0d0.and.virt_wgt.eq.0d0.and.deg_wgt.eq.0d0
     $        .and.deg_swgt.eq.0d0.and.cnt_wgt_c.eq.0d0 )
     $        CalculatedBorn=.false.

         if (.not.calculatedBorn .and. p_born(0,1).gt.0d0)then
            call sborn(p_born,wgt1)
         elseif(p_born(0,1).lt.0d0)then
            enhance=0d0
         endif

         if (enhance.eq.0d0)then
            xnoborn_cnt=xnoborn_cnt+1.d0
            if(log10(xnoborn_cnt).gt.inoborn_cnt)then
               write (*,*) 
     #           'Function dsigF: no Born momenta more than 10**',
     #           inoborn_cnt,'times'
               inoborn_cnt=inoborn_cnt+1
            endif
         else
            xtot=0d0
            if (mapconfig(0).eq.0) then
               write (*,*) 'Fatal error in dsigF, no Born diagrams '
     &           ,mapconfig,'. Check bornfromreal.inc'
               write (*,*) 'Is fks_singular compiled correctly?'
               stop
            endif
            do i=1, mapconfig(0)
               xtot=xtot+amp2(mapconfig(i))
            enddo
            if (xtot.ne.0d0) then
               enhance=amp2(mapconfig(iconfig))/xtot
               enhance=enhance*diagramsymmetryfactor
            else
               enhance=0d0
            endif
         endif
      endif

      cnt_wgt = cnt_wgt_c + cnt_wgt_s + cnt_wgt_sc
      cnt_swgt = cnt_swgt_s + cnt_swgt_sc

      totH_wgt = Hev_wgt+Hxmc_wgt

c$$$c Apply the FxFx Sudakov damping on the H events
c$$$      if (totH_wgt.ne.0d0 .and. ickkw.eq.3..and.
c$$$     $     .not.rewgt_mohdr_calculated) then
c$$$         write (*,*) 'Error rewgt_mohdr_calculated',totH_wgt
c$$$         stop
c$$$      elseif(rewgt_mohdr_calculated) then
c$$$         if (rewgt_mohdr.gt.1d0) rewgt_mohdr=1d0
c$$$         enhanceH=enhance*rewgt_mohdr
c$$$      else
         enhanceH=enhance
c$$$      endif

      totH_wgt = totH_wgt * enhanceH

c$$$c Apply the FxFx Sudakov damping on the S events
c$$$      if(.not.(Sev_wgt.eq.0d0 .and. Sxmc_wgt.eq.0d0 .and. cnt_wgt.eq.0d0
c$$$     $     .and. cnt_swgt.eq.0d0 .and. bsv_wgt.eq.0d0 .and.
c$$$     $     born_wgt.eq.0d0 .and. deg_wgt.eq.0d0 .and. deg_swgt.eq.0d0)
c$$$     $     .and. ickkw.eq.3 .and. .not.rewgt_izero_calculated) then
c$$$         write (*,*) 'Error rewgt_izero_calculated'
c$$$         stop
c$$$      elseif(rewgt_izero_calculated) then
c$$$         if (rewgt_izero.gt.1d0) rewgt_izero=1d0
c$$$         enhanceS=enhance*rewgt_izero
c$$$      else
         enhanceS=enhance
c$$$      endif

      Sev_wgt = Sev_wgt * enhanceS
      Sxmc_wgt = Sxmc_wgt * enhanceS
      cnt_wgt = cnt_wgt * enhanceS
      cnt_swgt = cnt_swgt * enhanceS
      bsv_wgt = bsv_wgt * enhanceS
      virt_wgt = virt_wgt * enhanceS
      born_wgt = born_wgt * enhanceS
      deg_wgt = deg_wgt * enhanceS
      deg_swgt = deg_swgt * enhanceS

c Update the shower starting scale with the shape from montecarlocounter
      if( (.not.MCcntcalled) .and.
     &     abrv.ne.'born'.and. abrv.ne.'grid' .and. ickkw.ne.4)then
         call kinematics_driver(xi_i_fks_ev,y_ij_fks_ev,shat,pp,ileg,
     &     xm12,ddum(1),ddum(2),ddum(3),ddum(4),ddum(5),ddum(6),ldum)
      endif
      if (.not.nbody) then
         call set_cms_stuff(mohdr)
         call set_shower_scale(nFKSprocess*2,.true.)
      endif
      call set_cms_stuff(izero)
      call set_shower_scale(nFKSprocess*2-1,.false.)

      if(iminmax.eq.0) then
         dsigS = (Sev_wgt+Sxmc_wgt+cnt_wgt)*fkssymmetryfactor +
     &        cnt_swgt*fkssymmetryfactor +
     &        bsv_wgt*fkssymmetryfactorBorn +
     &        virt_wgt*fkssymmetryfactorBorn +
     &        deg_wgt*fkssymmetryfactorDeg +
     &        deg_swgt*fkssymmetryfactorDeg

c Add the H-events to the S-events for UNLOPS
         if (ickkw.eq.4) then
            dsigS=dsigS+totH_wgt*fkssymmetryfactor
         endif

         call unweight_function(p_born,unwgtfun)
         dsigS=dsigS*unwgtfun

         if (dsigS.ne.dsigS) then
            write (*,*) 'ERROR, ',dsigS,
     &           ' found for dsigS, setting dsigS to 0 for this event'
            dsigS=0
            do j=1,iproc_save(nFKSprocess)
               if (.not.nbody) then
                  unwgt_table(nFKSprocess,1,j)=0d0
               else
                  unwgt_table(0,1,j)=0d0
                  unwgt_table(0,3,j)=0d0
                  unwgt_table(1,3,j)=0d0
               endif
            enddo
         endif
         
         if (fkssymmetryfactorDeg.ne.fkssymmetryfactor) then
            write (*,*) 'FKS symmetry factors should be identical'
     $           ,fkssymmetryfactorDeg,fkssymmetryfactor
            stop
         endif
         do j=1,iproc_save(nFKSprocess)
            if (.not.nbody) then
               unwgt_table(nFKSprocess,1,j)=unwgt_table(nFKSprocess,1,j)
     $              *enhanceS*fkssymmetryfactor*unwgtfun*vegaswgt
            else
               unwgt_table(0,1,j)=unwgt_table(0,1,j)
     $              *enhanceS*fkssymmetryfactorBorn*unwgtfun*vegaswgt
               unwgt_table(0,3,j)=unwgt_table(0,3,j)
     $              *enhanceS*fkssymmetryfactorBorn*unwgtfun*vegaswgt
               unwgt_table(1,3,j)=unwgt_table(1,3,j)
     $              *enhanceS*fkssymmetryfactorBorn*unwgtfun*vegaswgt
            endif
         enddo
         if(doreweight)then
            if(ifill2S.eq.0.and.(ifill3S.ne.0.or.ifill4S.ne.0))then
               write(*,*)'Error #2[wg] in dsigF S',ifill2S ,ifill3S
     &              ,ifill4S
               stop
            endif
            if ((nbody.and.ifill2S_born.eq.1).or.
     &           (.not.nbody.and.ifill2S.eq.0.and.ifill2S_born.eq.1))
     &           then
c
c Copy the values for the nbody configuration over the nFKSprocess*2-1
c
c In some rare cases with massive j_fks the counterevents need to be
c skipped, while the Born momenta are there. So, we need to fill the
c common blocks accordingly.
c
               wgtqes2_all(2,nFKSprocess*2-1)=wgtqes2_all(2,0)
               wgtqes2_all(3,nFKSprocess*2-1)=wgtqes2_all(3,0)
               wgtqes2_all(4,nFKSprocess*2-1)=wgtqes2_all(4,0)
               wgtmuR2_all(2,nFKSprocess*2-1)=wgtmuR2_all(2,0)
               wgtmuF12_all(2,nFKSprocess*2-1)=wgtmuF12_all(2,0)
               wgtmuF22_all(2,nFKSprocess*2-1)=wgtmuF22_all(2,0)
               do i=1,nexternal
                  do j=0,3
                     if (i.lt.nexternal) then
                        wgtkin_all(j,i,2,nFKSprocess*2-1)=p_born(j,i)
                     else
                        wgtkin_all(j,i,2,nFKSprocess*2-1)=0d0
                     endif
                  enddo
               enddo
               wgtxbj_all(1,2,nFKSprocess*2-1)=wgtxbj_all(1,2,0)
               wgtxbj_all(2,2,nFKSprocess*2-1)=wgtxbj_all(2,2,0)
            endif

            if (nbody) then
               do i_process=1,iproc_save(nFKSprocess) 
                  wgtref_nbody_all(i_process)=0d0
                  do j=1,iproc_save(nFKSprocess)
                     if (eto(j,nFKSprocess).eq.i_process)
     $                    wgtref_nbody_all(i_process)
     $                    =wgtref_nbody_all(i_process)+(unwgt_table(0,1
     $                    ,j)+unwgt_table(0,3,j))/vegaswgt
                  enddo
               enddo
            endif
            do i_process=1,iproc_save(nFKSprocess) 
               wgtref_all(nFKSprocess*2-1,i_process)=0d0
               do j=1,iproc_save(nFKSprocess)
                  if (eto(j,nFKSprocess).eq.i_process)
     $                 wgtref_all(nFKSprocess*2-1,i_process)
     $                 =wgtref_all(nFKSprocess*2-1,i_process)
     $                 +unwgt_table(nFKSprocess,1,j)/vegaswgt
               enddo
            enddo
            xsec = enhanceS*unwgtfun
            do i=1,4
               if (.not.nbody) then
                  wgtwreal_all(i,nFKSprocess*2-1)=wgtwreal_all(i
     &                 ,nFKSprocess*2-1) * xsec*fkssymmetryfactor
                  wgtwdeg_all(i,nFKSprocess*2-1)=wgtwdeg_all(i
     &                 ,nFKSprocess*2-1) * xsec*fkssymmetryfactorDeg
                  wgtwdegmuf_all(i,nFKSprocess*2-1)=wgtwdegmuf_all(i
     &                 ,nFKSprocess*2-1) * xsec*fkssymmetryfactorDeg
               endif
            enddo
            if (nbody) then
               wgtwborn_all=wgtwborn_all * xsec*fkssymmetryfactorBorn
               wgtwns_all=wgtwns_all * xsec*fkssymmetryfactorBorn
               wgtwnsmuf_all=wgtwnsmuf_all * xsec*fkssymmetryfactorBorn
               wgtwnsmur_all=wgtwnsmur_all * xsec*fkssymmetryfactorBorn
            endif
            if (.not.nbody) then
               do i=1,iwgtnumpartn_all(nFKSprocess*2-1)
                  wgtwmcxsec_all(i,nFKSprocess*2-1)=wgtwmcxsec_all(i
     &                 ,nFKSprocess*2-1) * xsec*fkssymmetryfactor
               enddo
            endif
            if(check_reweight.and.doreweight .and. ickkw.ne.4) then
               do i_process=1,iproc_save(nFKSprocess)
                  if (nbody) then
                     call fill_reweight0inc_nbody(i_process)
                     call check_rwgt_wgt("nbd")
                  else
                     call fill_reweight0inc(nFKSprocess*2-1,i_process)
                     call check_rwgt_wgt("Sev")
                  endif
                  call reweight_settozero()
               enddo
            endif
c Example of reweighted cross section (scale changed)
c           dsigS_new=compute_rwgt_wgt_Sev(new_muR_fact,new_muF1_fact,
c     &                                    new_muF2_fact,new_QES_fact,
c     &                                    iwgtinfo)
         endif

         total_wgt_sum=total_wgt_sum+dsigS*vegaswgt
         central_wgt_saved=dsigS
c For tests
         if(abs(dsigS).gt.fksmaxwgt)then
            fksmaxwgt=abs(dsigS)
            xisave=xi_i_fks_ev
            ysave=y_ij_fks_ev
         endif
c Plot observables for counterevents and Born
         if (.not.unwgt) then
            plot_wgt=( (Sev_wgt+Sxmc_wgt+cnt_wgt)*fkssymmetryfactor +
     &           cnt_swgt*fkssymmetryfactor +
     &           (bsv_wgt-born_wgt)*fkssymmetryfactorBorn +
     &           virt_wgt*fkssymmetryfactorBorn +
     &           deg_wgt*fkssymmetryfactorDeg +
     &           deg_swgt*fkssymmetryfactorDeg )*vegaswgt
            if( abs(plot_wgt).gt.1.d-20.and.p1_cnt(0,1,0).ne.-99d0 .and.
     &           (plotEv.or.plotKin) )
     &           call outfun(p1_cnt(0,1,0),ybst_til_tolab,plot_wgt,iplot_cnt)
            plot_wgt=(born_wgt*fkssymmetryfactorBorn)*vegaswgt
            if( abs(plot_wgt).gt.1.d-20.and.p1_cnt(0,1,0).ne.-99d0 .and.
     &           (plotEv.or.plotKin) )
     &           call outfun(p1_cnt(0,1,0),ybst_til_tolab,plot_wgt,iplot_born)
         endif

c For UNLOPS, the H-events are already added to the S-events
         if (ickkw.eq.4) return

         dsigH = totH_wgt*fkssymmetryfactor
         call unweight_function(p_born,unwgtfun)
         dsigH=dsigH*unwgtfun

         if (dsigH.ne.dsigH) then
            write (*,*) 'ERROR, ',dsigH,
     &           ' found for dsigH, setting dsigH to 0 for this event'
            dsigH=0
            if (.not.nbody) then
               do j=1,iproc_save(nFKSprocess)
                  unwgt_table(nFKSprocess,2,j)=0d0
               enddo
         endif
         endif

         if (nbody.and.dsigH.ne.0d0) then
            write (*,*) 'When doing nbody, contribution '/
     &           /'from H-events should be zero',dsigH
            stop
         endif

         if (.not.nbody) then
            do j=1,iproc_save(nFKSprocess)
               unwgt_table(nFKSprocess,2,j)=unwgt_table(nFKSprocess,2,j)
     $              *enhanceH*fkssymmetryfactor*unwgtfun*vegaswgt
            enddo
         endif
         if(doreweight)then
            if(ifill2H.eq.0.and.(ifill3H.ne.0.or.ifill4H.ne.0))then
               write(*,*)'Error #2[wg] in dsigF H',ifill2H,ifill3H
     &              ,ifill4H
               stop
            endif
            if (.not.nbody) then
               do j=1,iproc_save(nFKSprocess)
                  wgtref_all(nFKSprocess*2,j)=unwgt_table(nFKSprocess,2
     $                 ,j)/vegaswgt
               enddo
               xsec = enhanceH*unwgtfun*fkssymmetryfactor
               do i=1,4
                  wgtwreal_all(i,nFKSprocess*2)=wgtwreal_all(i
     &                 ,nFKSprocess*2) * xsec
               enddo
               do i=1,iwgtnumpartn_all(nFKSprocess*2)
                  wgtwmcxsec_all(i,nFKSprocess*2)=wgtwmcxsec_all(i
     &                 ,nFKSprocess*2) * xsec
               enddo
               if(check_reweight.and.doreweight) then
                  do i_process=1,iproc_save(nFKSprocess)
                     call fill_reweight0inc(nFKSprocess*2,i_process)
                     call check_rwgt_wgt("Hev")
                     call reweight_settozero()
                  enddo
               endif
            else
               do j=1,iproc_save(nFKSprocess)
                  if (wgtref_all(nFKSprocess*2,j).ne.0d0) then
                     write (*,*) 'wgtref not zero',j,
     &                    wgtref_all(nFKSprocess*2,j)
                     stop
                  endif
               enddo
               do i=1,4
                  if (wgtwreal_all(i,nFKSprocess*2).ne.0d0) then
                     write (*,*) 'wgtwreal not zero',i,
     &                    wgtwreal_all(i,nFKSprocess*2)
                     stop
                  endif
               enddo
               do i=1,iwgtnumpartn_all(nFKSprocess*2)
                  if (wgtwmcxsec_all(i,nFKSprocess*2).ne.0d0) then
                     write (*,*) 'wgtwmcxsec not zero',i,
     &                    wgtwmcxsec_all(i,nFKSprocess*2)
                     stop
                  endif
               enddo
            endif
c Example of reweighted cross section (scale changed)
c           dsigH_new=compute_rwgt_wgt_Hev(new_muR_fact,new_muF1_fact,
c     &                                    new_muF2_fact,new_QES_fact,
c     &                                    iwgtinfo)
         endif      

c For tests
         if(abs(dsigH).gt.fksmaxwgt)then
            fksmaxwgt=abs(dsigH)
            xisave=xi_i_fks_ev
            ysave=y_ij_fks_ev
         endif

c Plot observables for event
         if (.not.unwgt) then
            plot_wgt=totH_wgt*fkssymmetryfactor*vegaswgt 
            if( abs(plot_wgt).gt.1.d-20.and.pp(0,1).ne.-99d0. and.
     &           (plotEv.or.plotKin) )
     &           call outfun(pp,ybst_til_tolab,plot_wgt,iplot_ev)
         endif


      elseif (iminmax.eq.1 .and. ExceptPSpoint) then
c for except PS points, this is the maximal approx for the virtual         
         call unweight_function(p_born,unwgtfun)
         dsigS_max = ((Sev_wgt+Sxmc_wgt+cnt_wgt)*fkssymmetryfactor +
     &        cnt_swgt*fkssymmetryfactor +
     &        bsv_wgt*fkssymmetryfactorBorn +
     &        virt_wgt*fkssymmetryfactorBorn +
     &        deg_wgt*fkssymmetryfactorDeg +
     &        deg_swgt*fkssymmetryfactorDeg)*unwgtfun
         total_wgt_sum_max=total_wgt_sum_max+
     &        ((dsigS_max - central_wgt_saved)*vegaswgt)**2

      elseif (iminmax.eq.2 .and. ExceptPSpoint) then
c for except PS points, this is the minimal approx for the virtual         
         call unweight_function(p_born,unwgtfun)
         dsigS_min = ((Sev_wgt+Sxmc_wgt+cnt_wgt)*fkssymmetryfactor +
     &        cnt_swgt*fkssymmetryfactor +
     &        bsv_wgt*fkssymmetryfactorBorn +
     &        virt_wgt*fkssymmetryfactorBorn +
     &        deg_wgt*fkssymmetryfactorDeg +
     &        deg_swgt*fkssymmetryfactorDeg)*unwgtfun
         total_wgt_sum_min=total_wgt_sum_min+
     &        ((central_wgt_saved - dsigS_min)*vegaswgt)**2
      else
         write (*,*) 'Error #12 in dsig',iminmax
         stop
      endif

c If exceptional PS point found, go back to beginning recompute
c the weight for this PS point using an approximation
c based on previous PS points (done in BinothLHA.f)
      if (ExceptPSpoint .and. iminmax.le.1) goto 44
      return
      end



      subroutine set_shower_scale(iFKS,Hevents)
      implicit none
      include "nexternal.inc"
      include "madfks_mcatnlo.inc"
      integer iFKS
      logical Hevents
      double precision xi_i_fks_ev,y_ij_fks_ev
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt
      double precision sqrtshat_ev,shat_ev
      common/parton_cms_ev/sqrtshat_ev,shat_ev
      double precision emsca,scalemin,scalemax,emsca_bare
      logical emscasharp
      common/cemsca/emsca,emsca_bare,emscasharp,scalemin,scalemax
      character*4 abrv
      common/to_abrv/abrv
      include 'nFKSconfigs.inc'
      double precision SCALUP(fks_configs*2)
      common /cshowerscale/SCALUP
      double precision shower_S_scale(fks_configs*2)
     &     ,shower_H_scale(fks_configs*2),ref_H_scale(fks_configs*2)
     &     ,pt_hardness
      common /cshowerscale2/shower_S_scale,shower_H_scale,ref_H_scale
     &     ,pt_hardness

      double precision xm12
      integer ileg
      common/cscaleminmax/xm12,ileg

c Initialise
      SCALUP(iFKS)=0d0
c S events
      if(.not.Hevents)then
         if(abrv.ne.'born'.and.abrv.ne.'grid'.and.
     &      dampMCsubt.and.emsca.ne.0d0)then
            SCALUP(iFKS)=min(emsca,scalemax)
         else
            call assign_scaleminmax(shat_ev,xi_i_fks_ev,scalemin,scalemax,ileg,xm12)
            SCALUP(iFKS)=scalemax
         endif
         SCALUP(iFKS)=min(SCALUP(iFKS),shower_S_scale(iFKS))
c H events
      else
         if(dampMCsubt.and.emsca.ne.0d0)then
            SCALUP(iFKS)=scalemax
         else
            call assign_scaleminmax(shat_ev,xi_i_fks_ev,scalemin,scalemax,ileg,xm12)
            SCALUP(iFKS)=scalemax
         endif
         SCALUP(iFKS)=min(SCALUP(iFKS),max(shower_H_scale(iFKS),
     &                    ref_H_scale(iFKS)-min(emsca,scalemax)))
      endif
c Minimal starting scale
      SCALUP(iFKS)=max(SCALUP(iFKS),3d0)

      return
      end


      subroutine set_shower_scale_noshape(pp,iFKS)
      implicit none
      integer iFKS,j,i,iSH,nmax
      include "nexternal.inc"
      include "madfks_mcatnlo.inc"
      include 'run.inc'
      include 'nFKSconfigs.inc'
      LOGICAL  IS_A_J(NEXTERNAL),IS_A_LP(NEXTERNAL),IS_A_LM(NEXTERNAL)
      LOGICAL  IS_A_PH(NEXTERNAL)
      COMMON /TO_SPECISA/IS_A_J,IS_A_LP,IS_A_LM,IS_A_PH
      double precision sqrtshat_ev,shat_ev
      common/parton_cms_ev/sqrtshat_ev,shat_ev
      double precision sqrtshat_cnt(-2:2),shat_cnt(-2:2)
      common/parton_cms_cnt/sqrtshat_cnt,shat_cnt
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      double precision shower_S_scale(fks_configs*2)
     &     ,shower_H_scale(fks_configs*2),ref_H_scale(fks_configs*2)
     &     ,pt_hardness
      common /cshowerscale2/shower_S_scale,shower_H_scale,ref_H_scale
     &     ,pt_hardness
      double precision ptparton,pt,pp(0:3,nexternal),ppp(0:3,nexternal)
      external pt
c jet cluster algorithm
      integer NN,NJET,JET(nexternal)
      double precision pQCD(0:3,nexternal),PJET(0:3,nexternal),rfj,sycut
     $     ,palg,amcatnlo_fastjetdmergemax,di(nexternal)
      external amcatnlo_fastjetdmergemax

c Initialise
      NN=0
      ppp=0d0
      pQCD=0d0
      pt_hardness=0d0
      do j=1,nexternal
         if (j.gt.nincoming.and.is_a_j(j)) then
            NN=NN+1
            ptparton=pt(pp(0,j))
         endif
      enddo

c Unphysical situation
      if(NN.le.0)then
         write(*,*)'Error in set_shower_scale_noshape:'
         write(*,*)'not enough QCD partons in process ',NN
         stop
c Processes without jets at the Born
      elseif(NN.eq.1)then
         shower_S_scale(iFKS)=sqrtshat_cnt(0)
         shower_H_scale(iFKS)=sqrtshat_ev-ptparton
c$$$         shower_H_scale(iFKS)=sqrtshat_cnt(0)
         ref_H_scale(iFKS)=0d0
c Processes with jets at the Born (iSH = 1 (2) means S (H) events)
      else
         do iSH=1,2
            if(iSH.eq.1)then
               nmax=nexternal-1
               do j=1,nmax
                  do i=0,3
                     ppp(i,j)=p_born(i,j)
                  enddo
               enddo
            elseif(iSH.eq.2)then
               nmax=nexternal
               do j=1,nmax
                  do i=0,3
                     ppp(i,j)=pp(i,j)
                  enddo
               enddo
            else
               write(*,*)'Wrong iSH inset_shower_scale_noshape: ',iSH
               stop
            endif
            if(ppp(0,1).gt.0d0)then
c Put all (light) QCD partons in momentum array for jet clustering.
               NN=0
               do j=nincoming+1,nmax
                  if (is_a_j(j))then
                     NN=NN+1
                     do i=0,3
                        pQCD(i,NN)=ppp(i,j)
                     enddo
                  endif
               enddo
c One MUST use kt, and no lower pt cut. The radius parameter can be changed
               palg=1d0         ! jet algorithm: 1.0=kt, 0.0=C/A, -1.0 = anti-kt
               sycut=0d0        ! minimum jet pt
               rfj=1d0          ! the radius parameter
               call amcatnlo_fastjetppgenkt_timed(pQCD,NN,rfj,sycut,palg,
     &                                            pjet,njet,jet)
               do i=1,NN
                  di(i)=sqrt(amcatnlo_fastjetdmergemax(i-1))
                  if (i.gt.1.and.di(i).gt.di(i-1))then
                     write(*,*)'Error in set_shower_scale_noshape'
                     write(*,*)NN,i,di(i),di(i-1)
                     stop
                  endif
               enddo
               if(iSH.eq.1)shower_S_scale(iFKS)=di(NN)
               if(iSH.eq.2)then
                  ref_H_scale(iFKS)=di(NN-1)
                  pt_hardness=di(NN)
c$$$                  shower_H_scale(iFKS)=ref_H_scale(iFKS)-pt_hardness
                  shower_H_scale(iFKS)=ref_H_scale(iFKS)-pt_hardness/2d0
               endif
            else
               if(iSH.eq.1)shower_S_scale(iFKS)=sqrtshat_cnt(0)
               if(iSH.eq.2)then
                  ref_H_scale(iFKS)=shower_S_scale(iFKS)
                  shower_H_scale(iFKS)=ref_H_scale(iFKS)
               endif
            endif
         enddo
      endif

      return
      end



      subroutine sreal(pp,xi_i_fks,y_ij_fks,wgt)
c Wrapper for the n+1 contribution. Returns the n+1 matrix element
c squared reduced by the FKS damping factor xi**2*(1-y).
c Close to the soft or collinear limits it calls the corresponding
c Born and multiplies with the AP splitting function or eikonal factors.
      implicit none
      include "nexternal.inc"
      include "coupl.inc"

      double precision pp(0:3,nexternal),wgt
      double precision xi_i_fks,y_ij_fks

      double precision shattmp,dot
      integer i,j

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      logical softtest,colltest
      common/sctests/softtest,colltest

      double precision zero,tiny
      parameter (zero=0d0)
      logical need_color_links, need_charge_links
      common /c_need_links/need_color_links, need_charge_links

      double precision pmass(nexternal)
      include "pmass.inc"

      if (softtest.or.colltest) then
         tiny=1d-8
      else
         tiny=1d-6
      endif

      if(pp(0,1).le.0.d0)then
c Unphysical kinematics: set matrix elements equal to zero
        wgt=0.d0
        return
      endif

c Consistency check -- call to set_cms_stuff() must be done prior to
c entering this function
      shattmp=2d0*dot(pp(0,1),pp(0,2))
      if(abs(shattmp/shat-1.d0).gt.1.d-5)then
        write(*,*)'Error in sreal: inconsistent shat'
        write(*,*)shattmp,shat
        stop
      endif

      if (1d0-y_ij_fks.lt.tiny)then
         if (pmass(j_fks).eq.zero.and.j_fks.le.nincoming)then
            call sborncol_isr(pp,xi_i_fks,y_ij_fks,wgt)
         elseif (pmass(j_fks).eq.zero.and.j_fks.gt.nincoming)then
            call sborncol_fsr(pp,xi_i_fks,y_ij_fks,wgt)
         else
            wgt=0d0
         endif
      elseif (xi_i_fks.lt.tiny)then
         if (need_color_links.or.need_charge_links)then
c has soft singularities
            call sbornsoft(pp,xi_i_fks,y_ij_fks,wgt)
         else
            wgt=0d0
         endif
      else
         call smatrix_real(pp,wgt)
         wgt=wgt*xi_i_fks**2*(1d0-y_ij_fks)
      endif

c      if(wgt.lt.0.d0)then
c         write(*,*) 'Fatal error #2 in sreal',wgt,xi_i_fks,y_ij_fks
c         do i=1,nexternal
c            write(*,*) 'particle ',i,', ',(pp(j,i),j=0,3)
c         enddo
c         stop
c      endif

      return
      end



      subroutine sborncol_fsr(p,xi_i_fks,y_ij_fks,wgt)
      implicit none
      include "nexternal.inc"
      double precision p(0:3,nexternal),wgt
      double precision xi_i_fks,y_ij_fks
C  
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      double precision xi_i_fks_ev,y_ij_fks_ev
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

      double complex xij_aor
      common/cxij_aor/xij_aor

      logical rotategranny
      common/crotategranny/rotategranny

      double precision cthbe,sthbe,cphibe,sphibe
      common/cbeangles/cthbe,sthbe,cphibe,sphibe

      double precision p_born_rot(0:3,nexternal-1)

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn

      integer i,imother_fks,iord
C ap and Q contain the QCD(1) and QED(2) Altarelli-Parisi kernel
      double precision t,z,ap(2),E_j_fks,E_i_fks,Q(2),cphi_mother,
     # sphi_mother,pi(0:3),pj(0:3),wgt_born
      double complex W1(6),W2(6),W3(6),W4(6),Wij_angle,Wij_recta
      double complex azifact

c Particle types (=color/charges) of i_fks, j_fks and fks_mother
      integer i_type,j_type,m_type
      double precision ch_i,ch_j,ch_m
      common/cparticle_types/i_type,j_type,m_type,ch_i,ch_j,ch_m

      double precision zero,vtiny
      parameter (zero=0d0)
      parameter (vtiny=1d-8)
      double complex ximag
      parameter (ximag=(0.d0,1.d0))

      include 'orders.inc'
      logical split_type(nsplitorders) 
      common /c_split_type/split_type
      complex*16 ans_cnt(2, nsplitorders), wgt1(2)
      common /c_born_cnt/ ans_cnt
      double complex ans_extra_cnt(2,nsplitorders)
      integer iextra_cnt, isplitorder_born, isplitorder_cnt
      common /c_extra_cnt/iextra_cnt, isplitorder_born, isplitorder_cnt
C  
      if(p_born(0,1).le.0.d0)then
c Unphysical kinematics: set matrix elements equal to zero
         write (*,*) "No born momenta in sborncol_fsr"
         wgt=0.d0
         return
      endif

      E_j_fks = p(0,j_fks)
      E_i_fks = p(0,i_fks)
      z = 1d0 - E_i_fks/(E_i_fks+E_j_fks)
      t = z * shat/4d0
      if(rotategranny .and. nexternal-1.ne.3)then
c Exclude 2->1 (at the Born level) processes: matrix elements are
c independent of the PS point, but non-zero helicity configurations
c might flip when rotating the momenta.
        do i=1,nexternal-1
          call trp_rotate_invar(p_born(0,i),p_born_rot(0,i),
     #                          cthbe,sthbe,cphibe,sphibe)
        enddo
        CalculatedBorn=.false.
        call sborn(p_born_rot,wgt_born)
        CalculatedBorn=.false.
        if (iextra_cnt.gt.0)
     1    call extra_cnt(p_born_rot, iextra_cnt, ans_extra_cnt)
      else
        call sborn(p_born,wgt_born)
        if (iextra_cnt.gt.0)
     1    call extra_cnt(p_born, iextra_cnt, ans_extra_cnt)
      endif
      call AP_reduced(j_type,i_type,ch_j,ch_i,t,z,ap)
      call Qterms_reduced_timelike(j_type,i_type,ch_i,ch_j,t,z,Q)
      wgt=0d0
      do iord = 1, nsplitorders
        if (.not.split_type(iord).or.(iord.ne.qed_pos.and.iord.ne.qcd_pos)) cycle
C check if any extra_cnt is needed
        if (iextra_cnt.gt.0) then
            if (iord.eq.isplitorder_born) then
            ! this is the contribution from the born ME
               wgt1(1) = ans_cnt(1,iord)
               wgt1(2) = ans_cnt(2,iord)
            else if (iord.eq.isplitorder_cnt) then
            ! this is the contribution from the extra cnt
               wgt1(1) = ans_extra_cnt(1,iord)
               wgt1(2) = ans_extra_cnt(2,iord)
            else
               write(*,*) 'ERROR in sborncol_fsr', iord
               stop
            endif
        else
           wgt1(1) = ans_cnt(1,iord)
           wgt1(2) = ans_cnt(2,iord)
        endif
        if ((abs(j_type).eq.3 .and.i_type.eq.8) .or.
     #      (dabs(ch_j).ne.0d0 .and.ch_i.eq.0d0)) then
           Q(1)=0d0
           Q(2)=0d0
           wgt1(2)=0d0
        elseif (m_type.eq.8.or.ch_m.eq.0d0) then
c Insert <ij>/[ij] which is not included by sborn()
           if (1d0-y_ij_fks.lt.vtiny)then
              azifact=xij_aor
           else
              do i=0,3
                pi(i)=p_i_fks_ev(i)
                pj(i)=p(i,j_fks)
              enddo
              if(rotategranny)then
                call trp_rotate_invar(pi,pi,cthbe,sthbe,cphibe,sphibe)
                call trp_rotate_invar(pj,pj,cthbe,sthbe,cphibe,sphibe)
              endif
              CALL IXXXSO(pi ,ZERO ,+1,+1,W1)        
              CALL OXXXSO(pj ,ZERO ,-1,+1,W2)        
              CALL IXXXSO(pi ,ZERO ,-1,+1,W3)        
              CALL OXXXSO(pj ,ZERO ,+1,+1,W4)        
              Wij_angle=(0d0,0d0)
              Wij_recta=(0d0,0d0)
              do i=1,4
                 Wij_angle = Wij_angle + W1(i)*W2(i)
                 Wij_recta = Wij_recta + W3(i)*W4(i)
              enddo
              azifact=Wij_angle/Wij_recta
           endif
c Insert the extra factor due to Madgraph convention for polarization vectors
           imother_fks=min(i_fks,j_fks)
           if(rotategranny)then
             call getaziangles(p_born_rot(0,imother_fks),
     #                       cphi_mother,sphi_mother)
           else
             call getaziangles(p_born(0,imother_fks),
     #                       cphi_mother,sphi_mother)
           endif
           wgt1(2) = -(cphi_mother-ximag*sphi_mother)**2 *
     #             wgt1(2) * azifact
        else
           write(*,*) 'FATAL ERROR in sborncol_fsr',i_type,j_type,i_fks,j_fks
           stop
        endif
        if (iord.eq.qcd_pos) wgt=wgt+dble(wgt1(1)*ap(1)+wgt1(2)*Q(1))
        if (iord.eq.qed_pos) wgt=wgt+dble(wgt1(1)*ap(2)+wgt1(2)*Q(2))
      enddo
      return
      end



      subroutine sborncol_isr(p,xi_i_fks,y_ij_fks,wgt)
      implicit none
      include "nexternal.inc"
      double precision p(0:3,nexternal),wgt
      double precision xi_i_fks,y_ij_fks
C  
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      double precision xi_i_fks_ev,y_ij_fks_ev
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

      double complex xij_aor
      common/cxij_aor/xij_aor

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn

c Particle types (=color/charges) of i_fks, j_fks and fks_mother
      integer i_type,j_type,m_type
      double precision ch_i,ch_j,ch_m
      common/cparticle_types/i_type,j_type,m_type,ch_i,ch_j,ch_m

      double precision p_born_rot(0:3,nexternal-1)

      integer i, iord
C ap and Q contain the QCD(1) and QED(2) Altarelli-Parisi kernel
      double precision t,z,ap(2),Q(2),cphi_mother,sphi_mother,
     $ pi(0:3),pj(0:3),wgt_born
      double complex W1(6),W2(6),W3(6),W4(6),Wij_angle,Wij_recta
      double complex azifact

      double precision zero,vtiny
      parameter (zero=0d0)
      parameter (vtiny=1d-8)
      double complex ximag
      parameter (ximag=(0.d0,1.d0))

      include 'orders.inc'
      logical split_type(nsplitorders) 
      common /c_split_type/split_type
      complex*16 ans_cnt(2, nsplitorders), wgt1(2)
      common /c_born_cnt/ ans_cnt
      double complex ans_extra_cnt(2,nsplitorders)
      integer iextra_cnt, isplitorder_born, isplitorder_cnt
      common /c_extra_cnt/iextra_cnt, isplitorder_born, isplitorder_cnt
C  
      if(p_born(0,1).le.0.d0)then
c Unphysical kinematics: set matrix elements equal to zero
         write (*,*) "No born momenta in sborncol_isr"
         wgt=0.d0
         return
      endif

      z = 1d0 - xi_i_fks
c sreal return {\cal M} of FKS except for the partonic flux 1/(2*s).
c Thus, an extra factor z (implicit in the flux of the reduced Born
c in FKS) has to be inserted here
      t = z*shat/4d0
      if(j_fks.eq.2 .and. nexternal-1.ne.3)then
c Rotation according to innerpin.m. Use rotate_invar() if a more 
c general rotation is needed.
c Exclude 2->1 (at the Born level) processes: matrix elements are
c independent of the PS point, but non-zero helicity configurations
c might flip when rotating the momenta.
        do i=1,nexternal-1
          p_born_rot(0,i)=p_born(0,i)
          p_born_rot(1,i)=-p_born(1,i)
          p_born_rot(2,i)=p_born(2,i)
          p_born_rot(3,i)=-p_born(3,i)
        enddo
        CalculatedBorn=.false.
        call sborn(p_born_rot,wgt_born)
        CalculatedBorn=.false.
        if (iextra_cnt.gt.0)
     1    call extra_cnt(p_born_rot, iextra_cnt, ans_extra_cnt)
      else
        call sborn(p_born,wgt_born)
        if (iextra_cnt.gt.0)
     1    call extra_cnt(p_born, iextra_cnt, ans_extra_cnt)
      endif
      call AP_reduced(m_type,i_type,ch_m,ch_i,t,z,ap)
      call Qterms_reduced_spacelike(m_type,i_type,ch_m,ch_i,t,z,Q)
      wgt=0d0
      do iord = 1, nsplitorders
        if (.not.split_type(iord).or.(iord.ne.qed_pos.and.iord.ne.qcd_pos)) cycle
C check if any extra_cnt is needed
        if (iextra_cnt.gt.0) then
            if (iord.eq.isplitorder_born) then
            ! this is the contribution from the born ME
               wgt1(1) = ans_cnt(1,iord)
               wgt1(2) = ans_cnt(2,iord)
            else if (iord.eq.isplitorder_cnt) then
            ! this is the contribution from the extra cnt
               wgt1(1) = ans_extra_cnt(1,iord)
               wgt1(2) = ans_extra_cnt(2,iord)
            else
               write(*,*) 'ERROR in sborncol_isr', iord
               stop
            endif
        else
           wgt1(1) = ans_cnt(1,iord)
           wgt1(2) = ans_cnt(2,iord)
        endif
        if (abs(m_type).eq.3.or.ch_m.ne.0d0) then
           Q(1)=0d0
           Q(2)=0d0
           wgt1(2)=0d0
        else
c Insert <ij>/[ij] which is not included by sborn()
           if (1d0-y_ij_fks.lt.vtiny)then
              azifact=xij_aor
           else
              do i=0,3
                 pi(i)=p_i_fks_ev(i)
                 pj(i)=p(i,j_fks)
              enddo
              if(j_fks.eq.2)then
c Rotation according to innerpin.m. Use rotate_invar() if a more 
c general rotation is needed
                 pi(1)=-pi(1)
                 pi(3)=-pi(3)
                 pj(1)=-pj(1)
                 pj(3)=-pj(3)
              endif
              CALL IXXXSO(pi ,ZERO ,+1,+1,W1)        
              CALL OXXXSO(pj ,ZERO ,-1,+1,W2)        
              CALL IXXXSO(pi ,ZERO ,-1,+1,W3)        
              CALL OXXXSO(pj ,ZERO ,+1,+1,W4)        
              Wij_angle=(0d0,0d0)
              Wij_recta=(0d0,0d0)
              do i=1,4
                 Wij_angle = Wij_angle + W1(i)*W2(i)
                 Wij_recta = Wij_recta + W3(i)*W4(i)
              enddo
              azifact=Wij_angle/Wij_recta
           endif
c Insert the extra factor due to Madgraph convention for polarization vectors
           if(j_fks.eq.2)then
             cphi_mother=-1.d0
             sphi_mother=0.d0
           else
             cphi_mother=1.d0
             sphi_mother=0.d0
           endif
           wgt1(2) = -(cphi_mother+ximag*sphi_mother)**2 *
     #             wgt1(2) * dconjg(azifact)
        endif
        if (iord.eq.qcd_pos) wgt=wgt+dble(wgt1(1)*ap(1)+wgt1(2)*Q(1))
        if (iord.eq.qed_pos) wgt=wgt+dble(wgt1(1)*ap(2)+wgt1(2)*Q(2))
      enddo
      return
      end


      subroutine AP_reduced(col1, col2, ch1, ch2, t, z, ap)
c Returns Altarelli-Parisi splitting function summed/averaged over helicities
c times prefactors such that |M_n+1|^2 = ap * |M_n|^2. This means
c    AP_reduced = (1-z) P_{S(col1,col2)->col1+col2}(z) * g^2/t
c Therefore, the labeling conventions for particle IDs are not as in FKS:
c col1 and col2 are the two particles emerging from the branching.
c col1 and col2 can be gluons (8) or (anti-)quarks (+-3), or photons (1).
c z is the fraction of the energy of col1 and t is the invariant mass of the
c mother.
c The first entry in AP is QCD, the second QED.
      implicit none

      integer col1, col2
      double precision ch1, ch2
      double precision z,ap(2),t

      double precision CA,TR,CF
      parameter (CA=3d0,TR=1d0/2d0,CF=4d0/3d0)
      double precision tiny
      parameter (tiny=1d-5)

      include "coupl.inc"

c sanity check
      if ((col1.eq.8.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.8.and.dabs(ch2).ne.0d0) .or.
     &    (col1.eq.1.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.1.and.dabs(ch2).ne.0d0) .or.
     &    (abs(col1).eq.3.and.dabs(dabs(ch1)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch1)-2/3d0).gt.tiny) .or.
     &    (abs(col2).eq.3.and.dabs(dabs(ch2)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch2)-2/3d0).gt.tiny)) then
         write (*,*) 'Fatal Error #0 in AP_reduced',col1,col2,ch1,ch2
         stop
      endif

c g->gg splitting
      if (col1.eq.8 .and. col2.eq.8) then
         ap(1) = 2d0 * CA * ( (1d0-z)**2/z + z + z*(1d0-z)**2 )
         ap(2) = 0d0

c g/a->qqbar splitting
      elseif (abs(col1).eq.3 .and. abs(col2).eq.3) then
         ap(1) = TR           * ( z**2 + (1d0-z)**2 )*(1d0-z)
         ap(2) = 3d0 * ch1**2 * ( z**2 + (1d0-z)**2 )*(1d0-z)

c q->gq splitting
      elseif (col1.eq.8 .and. abs(col2).eq.3) then
         ap(1) = CF     * (1d0+(1d0-z)**2)*(1d0-z)/z
         ap(2) = 0d0

c q->aq splitting
      elseif (col1.eq.1 .and. abs(col2).eq.3) then
         ap(1) = 0d0
         ap(2) = ch2**2 * (1d0+(1d0-z)**2)*(1d0-z)/z

c q->qg splitting
      elseif (abs(col1).eq.3 .and. col2.eq.8) then
         ap(1) = CF     * (1d0+z**2)
         ap(2) = 0d0

c q->qa splitting
      elseif (abs(col1).eq.3 .and. col2.eq.1) then
         ap(1) = 0d0
         ap(2) = ch1**2 * (1d0+z**2)

      else
         write (*,*) 'Fatal Error #1 in AP_reduced',col1,col2,ch1,ch2
         stop
      endif

      ap(1) = ap(1)*g**2/t
      ap(2) = ap(2)*dble(gal(1))**2/t

      return
      end


      subroutine AP_reduced_prime(col1,col2,ch1,ch2,t,z,apprime)
c Returns (1-z)*P^\prime * g^2/t, with the same conventions as AP_reduced
c The first entry in AP_prime is QCD, the second QED.
      implicit none

      integer col1, col2
      double precision ch1, ch2
      double precision z,apprime(2),t

      double precision CA,TR,CF
      parameter (CA=3d0,TR=1d0/2d0,CF=4d0/3d0)
      double precision tiny
      parameter (tiny=1d-5)

      include "coupl.inc"

c sanity check
      if ((col1.eq.8.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.8.and.dabs(ch2).ne.0d0) .or.
     &    (col1.eq.1.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.1.and.dabs(ch2).ne.0d0) .or.
     &    (abs(col1).eq.3.and.dabs(dabs(ch1)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch1)-2/3d0).gt.tiny) .or.
     &    (abs(col2).eq.3.and.dabs(dabs(ch2)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch2)-2/3d0).gt.tiny)) then
         write (*,*) 'Fatal Error #0 in AP_reduced_prime',col1,col2,ch1,ch2
         stop
      endif

c g->gg splitting
      if (col1.eq.8 .and. col2.eq.8) then
         apprime(1) = 0d0
         apprime(2) = 0d0

c g/a->qqbar splitting
      elseif (abs(col1).eq.3 .and. abs(col2).eq.3) then
         apprime(1) = -2 * TR           * z * (1d0-z)**2
         apprime(2) = -2 * 3d0 * ch1**2 * z * (1d0-z)**2

c q->gq splitting
      elseif (col1.eq.8 .and. abs(col2).eq.3) then
         apprime(1) = - CF     * z * (1d0-z)
         apprime(2) = 0d0

c q->aq splitting
      elseif (col1.eq.1 .and. abs(col2).eq.3) then
         apprime(1) = 0d0
         apprime(2) = - ch2**2 * z * (1d0-z)

c q->qg splitting
      elseif (abs(col1).eq.3 .and. col2.eq.8) then
         apprime(1) = - CF     * (1d0-z)**2
         apprime(2) = 0d0

c q->qa splitting
      elseif (abs(col1).eq.3 .and. col2.eq.1) then
         apprime(1) = 0d0
         apprime(2) = - ch1**2 * (1d0-z)**2

      else
         write (*,*) 'Fatal Error #1 in AP_reduced_prime',col1,col2,ch1,ch2
         stop
      endif

      apprime(1) = apprime(1)*g**2/t
      apprime(2) = apprime(2)*dble(gal(1))**2/t

      return
      end


      subroutine Qterms_reduced_timelike(col1,col2,ch1,ch2,t,z,Qterms)
c Eq's B.31 to B.34 of FKS paper, times (1-z)*g^2/t. The labeling
c conventions for particle IDs are the same as those in AP_reduced
c The first entry in Qterms is QCD, the second QED.
      implicit none

      integer col1, col2
      double precision ch1, ch2
      double precision z,Qterms(2),t

      double precision CA,TR,CF
      parameter (CA=3d0,TR=1d0/2d0,CF=4d0/3d0)
      double precision tiny
      parameter (tiny=1d-5)

      include "coupl.inc"

c sanity check
      if ((col1.eq.8.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.8.and.dabs(ch2).ne.0d0) .or.
     &    (col1.eq.1.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.1.and.dabs(ch2).ne.0d0) .or.
     &    (abs(col1).eq.3.and.dabs(dabs(ch1)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch1)-2/3d0).gt.tiny) .or.
     &    (abs(col2).eq.3.and.dabs(dabs(ch2)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch2)-2/3d0).gt.tiny)) then
         write (*,*) 'Fatal Error #0 in Qterms_reduced_timelike',col1,col2,ch1,ch2
         stop
      endif

c g->gg splitting
      if (col1.eq.8 .and. col2.eq.8) then
         Qterms(1) = -4d0 * CA * z*(1d0-z)**2
         Qterms(2) = 0d0

c g/a->qqbar splitting
      elseif (abs(col1).eq.3 .and. abs(col2).eq.3) then
         Qterms(1) = 4d0 * TR           * z*(1d0-z)**2
         Qterms(2) = 4d0 * 3d0 * ch1**2 * z*(1d0-z)**2

c q->gq splitting
      elseif (col1.eq.8 .and. abs(col2).eq.3) then
         Qterms(1) = 0d0
         Qterms(2) = 0d0

c q->aq splitting
      elseif (col1.eq.1 .and. abs(col2).eq.3) then
         Qterms(1) = 0d0
         Qterms(2) = 0d0

c q->qg splitting
      elseif (abs(col1).eq.3 .and. col2.eq.8) then
         Qterms(1) = 0d0
         Qterms(2) = 0d0

c q->qa splitting
      elseif (abs(col1).eq.3 .and. col2.eq.1) then
         Qterms(1) = 0d0
         Qterms(2) = 0d0

      else
         write (*,*) 'Fatal Error #1 in Qterms_reduced_timelike',col1,col2,ch1,ch2
         stop
      endif

      Qterms(1) = Qterms(1)*g**2/t
      Qterms(2) = Qterms(2)*dble(gal(1))**2/t

      return
      end


      subroutine Qterms_reduced_spacelike(col1,col2,ch1,ch2,t,z,Qterms)
c Eq's B.42 to B.45 of FKS paper, times (1-z)*gS^2/t. The labeling
c conventions for particle IDs are the same as those in AP_reduced.
c Thus, col1 has momentum fraction z, and it is the one off-shell
c (see (FKS.B.41))
c The first entry in Qterms is QCD, the second QED.
      implicit none

      integer col1, col2
      double precision ch1, ch2
      double precision z,Qterms(2),t

      double precision CA,TR,CF
      parameter (CA=3d0,TR=1d0/2d0,CF=4d0/3d0)
      double precision tiny
      parameter (tiny=1d-5)

      include "coupl.inc"

c sanity check
      if ((col1.eq.8.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.8.and.dabs(ch2).ne.0d0) .or.
     &    (col1.eq.1.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.1.and.dabs(ch2).ne.0d0) .or.
     &    (abs(col1).eq.3.and.dabs(dabs(ch1)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch1)-2/3d0).gt.tiny) .or.
     &    (abs(col2).eq.3.and.dabs(dabs(ch2)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch2)-2/3d0).gt.tiny)) then
         write (*,*) 'Fatal Error #0 in Qterms_reduced_spacelike',col1,col2,ch1,ch2
         stop
      endif

c g->gg splitting
      if (col1.eq.8 .and. col2.eq.8) then
         Qterms(1) = -4d0 * CA * (1d0-z)**2/z
         Qterms(2) = 0d0

c g/a->qqbar splitting
      elseif (abs(col1).eq.3 .and. abs(col2).eq.3) then
         Qterms(1) = 0d0
         Qterms(2) = 0d0

c q->gq splitting
      elseif (col1.eq.8 .and. abs(col2).eq.3) then
         Qterms(1) = -4d0 * CF     * (1d0-z)**2/z
         Qterms(2) = 0d0

c q->aq splitting
      elseif (col1.eq.1 .and. abs(col2).eq.3) then
         Qterms(1) = 0d0
         Qterms(2) = -4d0 * ch2**2 * (1d0-z)**2/z

c q->qg splitting
      elseif (abs(col1).eq.3 .and. col2.eq.8) then
         Qterms(1) = 0d0
         Qterms(2) = 0d0

c q->qa splitting
      elseif (abs(col1).eq.3 .and. col2.eq.1) then
         Qterms(1) = 0d0
         Qterms(2) = 0d0

      else
         write (*,*) 'Fatal Error #1 in Qterms_reduced_spacelike',col1,col2,ch1,ch2
         stop
      endif

      Qterms(1) = Qterms(1)*g**2/t
      Qterms(2) = Qterms(2)*dble(gal(1))**2/t

      return
      end


      subroutine AP_reduced_SUSY(col1, col2, ch1, ch2, t, z, ap)
c Same as AP_reduced, except for the fact that it only deals with
c   go -> go g
c   sq -> sq g
c   sq -> sq a
c splittings in SUSY. We assume this function to be called with 
c col2==colour(i_fks)
      implicit none

      integer col1, col2
      double precision ch1, ch2
      double precision z,ap(2),t

      double precision CA,TR,CF
      parameter (CA=3d0,TR=1d0/2d0,CF=4d0/3d0)
      double precision tiny
      parameter (tiny=1d-5)

      include "coupl.inc"

c sanity check
      if ((col1.eq.8.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.8.and.dabs(ch2).ne.0d0) .or.
     &    (col1.eq.1.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.1.and.dabs(ch2).ne.0d0) .or.
     &    (abs(col1).eq.3.and.dabs(dabs(ch1)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch1)-2/3d0).gt.tiny) .or.
     &    (abs(col2).eq.3.and.dabs(dabs(ch2)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch2)-2/3d0).gt.tiny)) then
         write (*,*) 'Fatal Error #0 in AP_reduced_SUSY',col1,col2,ch1,ch2
         stop
      endif
      if (col2.ne.8.and.col2.ne.1)then
         write (*,*) 'Fatal error #1 in AP_reduced_SUSY',col1,col2
         stop
      endif

c go->gog splitting
      if (col1.eq.8 .and. col2.eq.8) then
         ap(1) = CA * (1d0+z**2)
         ap(2) = 0d0

c sq->sqg splitting
      elseif (abs(col1).eq.3 .and. col2.eq.8) then
         ap(1) = 2d0 * CF     * z
         ap(2) = 0d0

c sq->sqa splitting
      elseif (abs(col1).eq.3 .and. col2.eq.1) then
         ap(1) = 0d0
         ap(2) = 2d0 * ch1**2 * z

      else
         write (*,*) 'Fatal Error #2 in AP_reduced_SUSY',col1,col2,ch1,ch2
         stop
      endif

      ap(1) = ap(1)*g**2/t
      ap(2) = ap(2)*dble(gal(1))**2/t

c check this function further!!
c check this function further!!
      return
      end


      subroutine AP_reduced_massive(col1, col2, ch1, ch2, t, z, q2, m2, ap)
c Returns massive Altarelli-Parisi splitting function summed/averaged over helicities
c times prefactors such that |M_n+1|^2 = ap * |M_n|^2. This means
c    AP_reduced = (1-z) P_{S(col1,col2)->col1+col2}(z) * g^2/t
c Therefore, the labeling conventions for particle IDs are not as in FKS:
c col1 and col2 are the two particles emerging from the branching.
c col1 and col2 can be gluons (8) or (anti-)quarks (+-3), or photons (1).
c z is the fraction of the energy of col1 and t is the invariant mass of the
c mother.
c The first entry in AP is QCD, the second QED.


      implicit none

      integer col1, col2
      double precision ch1, ch2
      double precision z,ap(2),t,q2,m2

      double precision CA,TR,CF
      parameter (CA=3d0,TR=1d0/2d0,CF=4d0/3d0)
      double precision tiny
      parameter (tiny=1d-5)

      include "coupl.inc"

c sanity check
      if ((col1.eq.8.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.8.and.dabs(ch2).ne.0d0) .or.
     &    (col1.eq.1.and.dabs(ch1).ne.0d0) .or.
     &    (col2.eq.1.and.dabs(ch2).ne.0d0) .or.
     &    (abs(col1).eq.3.and.dabs(dabs(ch1)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch1)-2/3d0).gt.tiny) .or.
     &    (abs(col2).eq.3.and.dabs(dabs(ch2)-1/3d0).gt.tiny
     &                   .and.dabs(dabs(ch2)-2/3d0).gt.tiny)) then
         write (*,*) 'Fatal Error #0 in AP_reduced_massive',col1,col2,ch1,ch2
         stop
      endif

c g->gg splitting
      if (col1.eq.8 .and. col2.eq.8) then
         ap(1) = 2d0 * CA * ( (1d0-z)**2/z + z + z*(1d0-z)**2 )
         ap(2) = 0d0

c g/a->qqbar splitting
      elseif (abs(col1).eq.3 .and. abs(col2).eq.3) then
         ap(1) = TR           * ( z**2 + (1d0-z)**2 )*(1d0-z) + TR           * 2d0*m2/(z*q2)
         ap(2) = 3d0 * ch1**2 * ( z**2 + (1d0-z)**2 )*(1d0-z) + 3d0 * ch1**2 * 2d0*m2/(z*q2)

c q->gq splitting
      elseif (col1.eq.8 .and. abs(col2).eq.3) then
         ap(1) = CF     * (1d0+(1d0-z)**2)*(1d0-z)/z - CF     * 2d0*m2/(z*q2)
         ap(2) = 0d0

c q->aq splitting
      elseif (col1.eq.1 .and. abs(col2).eq.3) then
         ap(1) = 0d0
         ap(2) = ch2**2 * (1d0+(1d0-z)**2)*(1d0-z)/z - ch2**2 * 2d0*m2/(z*q2)

c q->qg splitting
      elseif (abs(col1).eq.3 .and. col2.eq.8) then
         ap(1) = CF     * (1d0+z**2) - CF      * 2d0*m2/(z*q2)
         ap(2) = 0d0

c q->qa splitting
      elseif (abs(col1).eq.3 .and. col2.eq.1) then
         ap(1) = 0d0
         ap(2) = ch1**2 * (1d0+z**2) - ch1**2  * 2d0*m2/(z*q2)

      else
         write (*,*) 'Fatal Error #1 in AP_reduced_massive',col1,col2,ch1,ch2
         stop
      endif

      ap(1) = ap(1)*g**2/t
      ap(2) = ap(2)*dble(gal(1))**2/t

      return
      end



      subroutine sbornsoft(pp,xi_i_fks,y_ij_fks,wgt)
      implicit none

      include "nexternal.inc"
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      include "coupl.inc"

      integer m,n

      double precision softcontr,pp(0:3,nexternal),wgt,eik,xi_i_fks
     &     ,y_ij_fks
      double precision wgt1
      integer i,j

      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision zero,pmass(nexternal)
      parameter(zero=0d0)
      include "pmass.inc"
c
c Call the Born to be sure that 'CalculatedBorn' is done correctly. This
c should always be done before calling the color-correlated Borns,
c because of the caching of the diagrams.
c
      call sborn(p_born(0,1),wgt1)
c
      softcontr=0d0
      do i=1,fks_j_from_i(i_fks,0)
         do j=1,i
            m=fks_j_from_i(i_fks,i)
            n=fks_j_from_i(i_fks,j)
            if ((m.ne.n .or. (m.eq.n .and. pmass(m).ne.ZERO)) .and.
     &           n.ne.i_fks.and.m.ne.i_fks) then
C wgt includes the gs/w^2 
               call sborn_sf(p_born,m,n,wgt)
               if (wgt.ne.0d0) then
                  call eikonal_reduced(pp,m,n,i_fks,j_fks,
     #                                 xi_i_fks,y_ij_fks,eik)
                  softcontr=softcontr+wgt*eik
               endif
            endif
         enddo
      enddo

      wgt=softcontr
c Add minus sign to compensate the minus in the color factor
c of the color-linked Borns (b_sf_0??.f)
c Factor two to fix the limits.
      wgt=-2d0*wgt
      return
      end


      subroutine eikonal_reduced(pp,m,n,i_fks,j_fks,xi_i_fks,y_ij_fks,eik)
c     Returns the eikonal factor
      implicit none

      include "nexternal.inc"
      double precision eik,pp(0:3,nexternal),xi_i_fks,y_ij_fks
      double precision dot,dotnm,dotni,dotmi,fact
      integer n,m,i_fks,j_fks,i
      integer softcol

      include "coupl.inc"

      external dot
      double precision xi_i_fks_ev,y_ij_fks_ev
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      real*8 phat_i_fks(0:3)

      double precision zero,pmass(nexternal),tiny
      parameter(zero=0d0)
      parameter(tiny=1d-6)
      include "pmass.inc"

c Define the reduced momentum for i_fks
      softcol=0
      if (1d0-y_ij_fks.lt.tiny)softcol=2
      if(p_i_fks_cnt(0,softcol).lt.0d0)then
        if(xi_i_fks.eq.0.d0)then
           write (*,*) 'Error #1 in eikonal_reduced',
     #                 softcol,xi_i_fks,y_ij_fks
           stop
        endif
        if(pp(0,i_fks).ne.0.d0)then
          write(*,*)'WARNING in eikonal_reduced: no cnt momenta',
     #      softcol,xi_i_fks,y_ij_fks
          do i=0,3
            phat_i_fks(i)=pp(i,i_fks)/xi_i_fks
          enddo
        else
          write (*,*) 'Error #2 in eikonal_reduced',
     #                 softcol,xi_i_fks,y_ij_fks
          stop
        endif
      else
        do i=0,3
          phat_i_fks(i)=p_i_fks_cnt(i,softcol)
        enddo
      endif
c Calculate the eikonal factor
      dotnm=dot(pp(0,n),pp(0,m))
      if ((m.ne.j_fks .and. n.ne.j_fks) .or. pmass(j_fks).ne.ZERO) then
         dotmi=dot(pp(0,m),phat_i_fks)
         dotni=dot(pp(0,n),phat_i_fks)
         fact= 1d0-y_ij_fks
      elseif (m.eq.j_fks .and. n.ne.j_fks .and.
     &        pmass(j_fks).eq.ZERO) then
         dotni=dot(pp(0,n),phat_i_fks)
         dotmi=sqrtshat/2d0 * pp(0,j_fks)
         fact= 1d0
      elseif (m.ne.j_fks .and. n.eq.j_fks .and.
     &        pmass(j_fks).eq.ZERO) then
         dotni=sqrtshat/2d0 * pp(0,j_fks)
         dotmi=dot(pp(0,m),phat_i_fks)
         fact= 1d0
      else
         write (*,*) 'Error #3 in eikonal_reduced'
         stop
      endif

      eik = dotnm/(dotni*dotmi)*fact
      return
      end


      subroutine sreal_deg(p,xi_i_fks,y_ij_fks,
     #                     collrem_xi,collrem_lxi)
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
      include "coupl.inc"
      include 'q_es.inc'
      include "run.inc"
      include 'reweight.inc'
      include "orders.inc"

      integer iord, iap
      double precision p(0:3,nexternal),collrem_xi,collrem_lxi
      double precision xi_i_fks,y_ij_fks

      double precision p_born(0:3,nexternal-1), wgt_born
      common/pborn/p_born

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      double precision delta_used
      common /cdelta_used/delta_used

      double precision rwgt,shattmp,dot,born_wgt,oo2pi,z,t,ap(2),
     # apprime(2),xkkern,xnorm
      external dot

c Particle types (=color/charges) of i_fks, j_fks and fks_mother
      integer i_type,j_type,m_type
      double precision ch_i, ch_j, ch_m
      common/cparticle_types/i_type,j_type,m_type,ch_i,ch_j,ch_m
      complex*16 ans_cnt(2, nsplitorders), wgt1(2)
      common /c_born_cnt/ ans_cnt
      logical split_type(nsplitorders) 
      common /c_split_type/split_type
      
      double precision one,pi
      parameter (one=1.d0)
      parameter (pi=3.1415926535897932385d0)

      if(j_fks.gt.nincoming)then
c Do not include this contribution for final-state branchings
         collrem_xi=0.d0
         collrem_lxi=0.d0
         if(doreweight)then
           wgtdegrem_xi=0.d0
           wgtdegrem_lxi=0.d0
           wgtdegrem_muF=0.d0
         endif
         return
      endif

      if(p_born(0,1).le.0.d0)then
c Unphysical kinematics: set matrix elements equal to zero
         write (*,*) "No born momenta in sreal_deg"
         collrem_xi=0.d0
         collrem_lxi=0.d0
         if(doreweight)then
           wgtdegrem_xi=0.d0
           wgtdegrem_lxi=0.d0
           wgtdegrem_muF=0.d0
         endif
         return
      endif

c Consistency check -- call to set_cms_stuff() must be done prior to
c entering this function
      shattmp=2d0*dot(p(0,1),p(0,2))
      if(abs(shattmp/shat-1.d0).gt.1.d-5)then
        write(*,*)'Error in sreal: inconsistent shat'
        write(*,*)shattmp,shat
        stop
      endif

      call sborn(p_born,wgt_born)

c A factor gS^2 is included in the Altarelli-Parisi kernels
      oo2pi=one/(8d0*PI**2)

      z = 1d0 - xi_i_fks
      t = one
      call AP_reduced(m_type,i_type,ch_m,ch_i,t,z,ap)
      call AP_reduced_prime(m_type,i_type,ch_m,ch_i,t,z,apprime)

      do iord = 1, nsplitorders
        if (.not.split_type(iord).or.(iord.ne.qed_pos.and.iord.ne.qcd_pos)) cycle
        wgt1(1) = ans_cnt(1,iord)
        wgt1(2) = ans_cnt(2,iord)
        

c Insert here proper functions for PDF change of scheme. With xkkern=0.d0
c one assumes MSbar
        xkkern=0.d0

        if (iord.eq.qcd_pos) iap = 1
        if (iord.eq.qed_pos) iap = 2
        collrem_xi=ap(iap)*log(shat*delta_used/(2*q2fact(j_fks))) -
     #           apprime(iap) - xkkern 
        collrem_lxi=2*ap(iap)

c The partonic flux 1/(2*s) is inserted in genps. Thus, an extra 
c factor z (implicit in the flux of the reduced Born in FKS) 
c has to be inserted here
        xnorm=1.d0/z

        collrem_xi=oo2pi * dble(wgt1(1)) * collrem_xi * xnorm
        collrem_lxi=oo2pi * dble(wgt1(1)) * collrem_lxi * xnorm

        if(doreweight)then
          write(*,*) 'FIX REWEIGHT'
          wgtdegrem_xi=ap(iap)*log(shat*delta_used/(2*QES2)) -
     #               apprime(iap) - xkkern 
          wgtdegrem_xi=oo2pi * dble(wgt1(1)) * wgtdegrem_xi * xnorm
          wgtdegrem_lxi=collrem_lxi
          wgtdegrem_muF= - oo2pi * dble(wgt1(1)) * ap(iap) * xnorm
        endif
      enddo

      return
      end


      subroutine set_cms_stuff(icountevts)
      implicit none
      include "run.inc"

      integer icountevts

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      double precision sqrtshat_ev,shat_ev
      common/parton_cms_ev/sqrtshat_ev,shat_ev

      double precision sqrtshat_cnt(-2:2),shat_cnt(-2:2)
      common/parton_cms_cnt/sqrtshat_cnt,shat_cnt

      double precision tau_ev,ycm_ev
      common/cbjrk12_ev/tau_ev,ycm_ev

      double precision tau_cnt(-2:2),ycm_cnt(-2:2)
      common/cbjrk12_cnt/tau_cnt,ycm_cnt

      double precision xbjrk_ev(2),xbjrk_cnt(2,-2:2)
      common/cbjorkenx/xbjrk_ev,xbjrk_cnt

c rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to lab frame --
c same for event and counterevents
c This is the rapidity that enters in the arguments of the sinh() and
c cosh() of the boost, in such a way that
c       y(k)_lab = y(k)_tilde - ybst_til_tolab
c where y(k)_lab and y(k)_tilde are the rapidities computed with a generic
c four-momentum k, in the lab frame and in the \tilde{k}_1+\tilde{k}_2 
c c.m. frame respectively
      ybst_til_tolab=-ycm_cnt(0)
      if(icountevts.eq.-100)then
c set Bjorken x's in run.inc for the computation of PDFs in auto_dsig
        xbk(1)=xbjrk_ev(1)
        xbk(2)=xbjrk_ev(2)
c shat=2*k1.k2 -- consistency of this assignment with momenta checked
c in phspncheck_nocms
        shat=shat_ev
        sqrtshat=sqrtshat_ev
c rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to 
c k_1+k_2 c.m. frame
        ybst_til_tocm=ycm_ev-ycm_cnt(0)
      else
c do the same as above for the counterevents
        xbk(1)=xbjrk_cnt(1,icountevts)
        xbk(2)=xbjrk_cnt(2,icountevts)
        shat=shat_cnt(icountevts)
        sqrtshat=sqrtshat_cnt(icountevts)
        ybst_til_tocm=ycm_cnt(icountevts)-ycm_cnt(0)
      endif
      return
      end


      subroutine get_mc_lum(j_fks,zhw_used,xi_i_fks,xlum_mc_save,
     #                      xlum_mc,xlum_mc_fact)
      implicit none
      include "run.inc"
      include "nexternal.inc"
      integer j_fks
      double precision dlum
      external dlum
      double precision zhw_used,xi_i_fks,xlum_mc_save,xlum_mc,
     # xlum_mc_fact
      
      double precision xbjrk_ev(2),xbjrk_cnt(2,-2:2)
      common/cbjorkenx/xbjrk_ev,xbjrk_cnt

      if(zhw_used.lt.0.d0.or.zhw_used.gt.1.d0)then
        write(*,*)'Error #1 in get_mc_lum',zhw_used
        stop
      endif
      if(j_fks.gt.nincoming)then
        xbk(1)=xbjrk_cnt(1,0)
        xbk(2)=xbjrk_cnt(2,0)
        xlum_mc_fact=1.d0
        if(xlum_mc_save.ne.-1.d8)then
          xlum_mc=xlum_mc_save
        else
          xlum_mc = dlum()
          xlum_mc_save = xlum_mc
        endif
      elseif(j_fks.eq.1)then
        xbk(1)=xbjrk_cnt(1,0)/zhw_used
        xbk(2)=xbjrk_cnt(2,0)
c Note that this is true for Pythia since, due to event projection and to
c the definition of the shower variable x = zhw_used, the Bjorken x's for
c the event (to be used in H events) are the ones for the counterevent
c multiplied by 1/x (by 1) for the emitting (non emitting) leg 
        if(xbk(1).gt.1.d0)then
          xlum_mc = 0.d0
          xlum_mc_fact = 0.d0
        else
          xlum_mc = dlum()
          xlum_mc_fact = (1-xi_i_fks)/zhw_used
          xlum_mc = xlum_mc * xlum_mc_fact
        endif
      elseif(j_fks.eq.2)then
        xbk(1)=xbjrk_cnt(1,0)
        xbk(2)=xbjrk_cnt(2,0)/zhw_used
        if(xbk(2).gt.1.d0)then
          xlum_mc = 0.d0
          xlum_mc_fact = 0.d0
        else
          xlum_mc = dlum()
          xlum_mc_fact = (1-xi_i_fks)/zhw_used
          xlum_mc = xlum_mc * xlum_mc_fact
        endif
      else
        write(*,*)'Error in get_mc_lum: unknown j_fks',j_fks
        stop
      endif
      if( xbk(1).le.0.d0.or.xbk(2).le.0.d0.or.
     #    ( (xbk(1).gt.1.d0.or.xbk(2).gt.1.d0).and.
     #      j_fks.gt.nincoming ) .or.
     #    (xbk(2).gt.1.d0.and.j_fks.eq.1) .or.
     #    (xbk(1).gt.1.d0.and.j_fks.eq.2) )then
        write(*,*)'Error in get_mc_lum: x_i',xbk(1),xbk(2)
        stop
      endif
      return
      end


      subroutine xmom_compare(i_fks,j_fks,jac,jac_cnt,p,p1_cnt,
     #                        p_i_fks_ev,p_i_fks_cnt,
     #                        xi_i_fks_ev,y_ij_fks_ev,pass)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      integer i_fks,j_fks
      double precision p(0:3,-max_branch:max_particles)
      double precision p1_cnt(0:3,nexternal,-2:2)
      double precision jac,jac_cnt(-2:2)
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      double precision xi_i_fks_ev,y_ij_fks_ev
      integer izero,ione,itwo,iunit,isum
      logical verbose,pass,pass0
      parameter (izero=0)
      parameter (ione=1)
      parameter (itwo=2)
      parameter (iunit=6)
      parameter (verbose=.false.)
      integer i_momcmp_count
      double precision xratmax
      common/ccheckcnt/i_momcmp_count,xratmax
c
      isum=0
      if(jac_cnt(0).gt.0.d0)isum=isum+1
      if(jac_cnt(1).gt.0.d0)isum=isum+2
      if(jac_cnt(2).gt.0.d0)isum=isum+4
      pass=.true.
c
      if(isum.eq.0.or.isum.eq.1.or.isum.eq.2.or.isum.eq.4)then
c Nothing to be done: 0 or 1 configurations computed
        if(verbose)write(iunit,*)'none'
      elseif(isum.eq.3.or.isum.eq.5.or.isum.eq.7)then
c Soft is taken as reference
        if(isum.eq.7)then
          if(verbose)then
            write(iunit,*)'all'
            write(iunit,*)'    '
            write(iunit,*)'C/S'
          endif
          call xmcompare(verbose,pass0,ione,izero,i_fks,j_fks,p,p1_cnt)
          pass=pass.and.pass0
          if(verbose)then
            write(iunit,*)'    '
            write(iunit,*)'SC/S'
          endif
          call xmcompare(verbose,pass0,itwo,izero,i_fks,j_fks,p,p1_cnt)
          pass=pass.and.pass0
        elseif(isum.eq.3)then
          if(verbose)then
            write(iunit,*)'C+S'
            write(iunit,*)'    '
            write(iunit,*)'C/S'
          endif
          call xmcompare(verbose,pass0,ione,izero,i_fks,j_fks,p,p1_cnt)
          pass=pass.and.pass0
        elseif(isum.eq.5)then
          if(verbose)then
            write(iunit,*)'SC+S'
            write(iunit,*)'    '
            write(iunit,*)'SC/S'
          endif
          call xmcompare(verbose,pass0,itwo,izero,i_fks,j_fks,p,p1_cnt)
          pass=pass.and.pass0
        endif
      elseif(isum.eq.6)then
c Collinear is taken as reference
        if(verbose)then
          write(iunit,*)'SC+C'
          write(iunit,*)'    '
          write(iunit,*)'SC/C'
        endif
        call xmcompare(verbose,pass0,itwo,ione,i_fks,j_fks,p,p1_cnt)
        pass=pass.and.pass0
      else
        write(6,*)'Fatal error in xmom_compare',isum
        stop
      endif
      if(.not.pass)i_momcmp_count=i_momcmp_count +1
c
      if(jac_cnt(0).gt.0.d0.and.jac.gt.0.d0)
     #  call p_ev_vs_cnt(izero,i_fks,j_fks,p,p1_cnt,
     #                   p_i_fks_ev,p_i_fks_cnt,
     #                   xi_i_fks_ev,y_ij_fks_ev)
      if(jac_cnt(1).gt.0.d0.and.jac.gt.0.d0)
     #  call p_ev_vs_cnt(ione,i_fks,j_fks,p,p1_cnt,
     #                   p_i_fks_ev,p_i_fks_cnt,
     #                   xi_i_fks_ev,y_ij_fks_ev)
c
      return
      end


      subroutine xmcompare(verbose,pass0,inum,iden,i_fks,j_fks,p,p1_cnt)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'coupl.inc'
      logical verbose,pass0
      integer inum,iden,i_fks,j_fks,iunit,ipart,i,j,k
      double precision tiny,vtiny,xnum,xden,xrat
      double precision p(0:3,-max_branch:max_particles)
      double precision p1_cnt(0:3,nexternal,-2:2)
      parameter (iunit=6)
      parameter (tiny=1.d-4)
      parameter (vtiny=1.d-10)
      double precision pmass(nexternal),zero
      parameter (zero=0d0)
      integer i_momcmp_count
      double precision xratmax
      common/ccheckcnt/i_momcmp_count,xratmax
      include "pmass.inc"
c
      pass0=.true.
      do ipart=1,nexternal
        do i=0,3
          xnum=p1_cnt(i,ipart,inum)
          xden=p1_cnt(i,ipart,iden)
          if(verbose)then
            if(i.eq.0)then
              write(iunit,*)' '
              write(iunit,*)'part=',ipart
            endif
            call xprintout(iunit,xnum,xden)
          else
            if(ipart.ne.i_fks.and.ipart.ne.j_fks)then
              if(xden.ne.0.d0)then
                xrat=abs(1-xnum/xden)
              else
                xrat=abs(xnum)
              endif
              if(abs(xnum).eq.0d0.and.abs(xden).le.vtiny)xrat=0d0
c The following line solves some problem as well, but before putting
c it as the standard, one should think a bit about it
              if(abs(xnum).le.vtiny.and.abs(xden).le.vtiny)xrat=0d0
              if(xrat.gt.tiny .and.
     &          (pmass(ipart).eq.0d0.or.xnum/pmass(ipart).gt.vtiny))then
                 write(*,*)'Kinematics of counterevents'
                 write(*,*)inum,iden
                 write(*,*)'is different. Particle:',ipart
                 write(*,*) xrat,xnum,xden
                 do j=1,nexternal
                    write(*,*) j,(p1_cnt(k,j,inum),k=0,3)
                 enddo
                 do j=1,nexternal
                    write(*,*) j,(p1_cnt(k,j,iden),k=0,3)
                 enddo
                 xratmax=max(xratmax,xrat)
                 pass0=.false.
              endif
            endif
          endif
        enddo
      enddo
      do i=0,3
        if(j_fks.gt.2)then
          xnum=p1_cnt(i,i_fks,inum)+p1_cnt(i,j_fks,inum)
          xden=p1_cnt(i,i_fks,iden)+p1_cnt(i,j_fks,iden)
        else
          xnum=-p1_cnt(i,i_fks,inum)+p1_cnt(i,j_fks,inum)
          xden=-p1_cnt(i,i_fks,iden)+p1_cnt(i,j_fks,iden)
        endif
        if(verbose)then
          if(i.eq.0)then
            write(iunit,*)' '
            write(iunit,*)'part=i+j'
          endif
          call xprintout(iunit,xnum,xden)
        else
          if(xden.ne.0.d0)then
            xrat=abs(1-xnum/xden)
          else
            xrat=abs(xnum)
          endif
          if(xrat.gt.tiny)then
            write(*,*)'Kinematics of counterevents'
            write(*,*)inum,iden
            write(*,*)'is different. Particle i+j'
            xratmax=max(xratmax,xrat)
            pass0=.false.
          endif
        endif
      enddo
      return
      end


      subroutine xmcompare_fsr(verbose,inum,iden,i_fks,j_fks,p,p1_cnt)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      logical verbose
      integer inum,iden,i_fks,j_fks,iunit,ipart,i
      double precision tiny,xnum,xden,xrat
      double precision p(0:3,-max_branch:max_particles)
      double precision p1_cnt(0:3,nexternal,-2:2)
      parameter (iunit=6)
      parameter (tiny=1.d-4)
c
      do ipart=1,nexternal
        do i=0,3
          xnum=p1_cnt(i,ipart,inum)
          xden=p1_cnt(i,ipart,iden)
          if(verbose)then
            if(i.eq.0)then
              write(iunit,*)' '
              write(iunit,*)'part=',ipart
            endif
            call xprintout(iunit,xnum,xden)
          else
            if(ipart.ne.i_fks.and.ipart.ne.j_fks)then
              if(xden.ne.0.d0)then
                xrat=abs(1-xnum/xden)
              else
                xrat=abs(xnum)
              endif
              if(xrat.gt.tiny)then
                write(*,*)'Kinematics of counterevents'
                write(*,*)inum,iden
                write(*,*)'is different. Particle:',ipart
                stop
              endif
            endif
          endif
        enddo
      enddo
      do i=0,3
        xnum=p1_cnt(i,i_fks,inum)+p1_cnt(i,j_fks,inum)
        xden=p1_cnt(i,i_fks,iden)+p1_cnt(i,j_fks,iden)
        if(verbose)then
          if(i.eq.0)then
            write(iunit,*)' '
            write(iunit,*)'part=i+j'
          endif
          call xprintout(iunit,xnum,xden)
        else
          if(xden.ne.0.d0)then
            xrat=abs(1-xnum/xden)
          else
            xrat=abs(xnum)
          endif
          if(xrat.gt.tiny)then
            write(*,*)'Kinematics of counterevents'
            write(*,*)inum,iden
            write(*,*)'is different. Particle i+j'
            stop
          endif
        endif
      enddo
      return
      end


      subroutine xprintout(iunit,xv,xlim)
      implicit real*8(a-h,o-z)
c
      if(abs(xlim).gt.1.d-30)then
        write(iunit,*)xv/xlim,xv,xlim
      else
        write(iunit,*)xv,xlim
      endif
      return
      end


      subroutine p_ev_vs_cnt(icnt,i_fks,j_fks,p,p1_cnt,
     #                       p_i_fks_ev,p_i_fks_cnt,
     #                       xi_i_fks_ev,y_ij_fks_ev)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      integer icnt,i_fks,j_fks,ipart,i
      double precision p(0:3,-max_branch:max_particles)
      double precision p1_cnt(0:3,nexternal,-2:2)
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      double precision xi_i_fks_ev,y_ij_fks_ev,tiny
      double precision rat(0:3,nexternal+3),den(0:3,nexternal+3)
      integer maxrat
c
c This routine is obsolete; the convergence checks are done elsewhere
      return

      do ipart=1,nexternal
        do i=0,3
          den(i,ipart)=p1_cnt(i,ipart,icnt)
          if(den(i,ipart).ne.0.d0)then
            rat(i,ipart)=p(i,ipart)/den(i,ipart)
          else
            rat(i,ipart)=p(i,ipart)
          endif
        enddo
      enddo
c
      do i=0,3
        den(i,nexternal+1)=p1_cnt(i,i_fks,icnt)+p1_cnt(i,j_fks,icnt)
        if(den(i,nexternal+1).ne.0.d0)then
          rat(i,nexternal+1)=(p(i,i_fks)+p(i,j_fks))/den(i,nexternal+1)
        else
          rat(i,nexternal+1)=p(i,i_fks)+p(i,j_fks)
        endif
      enddo
c
      if(icnt.eq.0)then
        tiny=4*xi_i_fks_ev
        maxrat=nexternal+3
        do i=0,3
          den(i,nexternal+2)=p_i_fks_cnt(i,0)
          if(den(i,nexternal+2).ne.0.d0)then
            rat(i,nexternal+2)=p_i_fks_ev(i)/den(i,nexternal+2)
          else
            rat(i,nexternal+2)=p_i_fks_ev(i)
          endif
        enddo
        do i=0,3
          den(i,nexternal+3)=p_i_fks_cnt(i,0)
          if(den(i,nexternal+3).ne.0.d0)then
            rat(i,nexternal+3)=p(i,i_fks)/den(i,nexternal+3)
          else
            rat(i,nexternal+3)=p(i,i_fks)
          endif
        enddo
      else
        tiny=2*sqrt(1-y_ij_fks_ev)
        maxrat=nexternal+1
      endif
c
      return
      end


c The following has been derived with minor modifications from the
c analogous routine written for VBF
      subroutine checkres(xsecvc,xseclvc,wgt,wgtl,xp,lxp,
     #                    iflag,imax,iev,nexternal,i_fks,j_fks,iret)
c Checks that the sequence xsecvc(i), i=1,imax, converges to xseclvc.
c Due to numerical inaccuracies, the test is deemed OK if there are
c at least ithrs+1 consecutive elements in the sequence xsecvc(i)
c which are closer to xseclvc than the preceding element of the sequence.
c The counting is started when an xsecvc(i0) is encountered, which is
c such that |xsecvc(i0)/xseclvc-1|<0.1 if xseclvc#0, or such that
c |xsecvc(i0)|<0.1 if xseclvc=0. In order for xsecvc(i+1 )to be defined 
c closer to xseclvc than xsecvc(i), the condition
c   |xsecvc(i)/xseclvc-1|/|xsecvc(i+1)/xseclvc-1| > rat
c if xseclvc#0, or 
c   |xsecvc(i)|/|xsecvc(i+1)| > rat
c if xseclvc=0 must be fulfilled; the value of rat is set equal to 4 and to 2
c for soft and collinear limits respectively, since the cross section is 
c expected to scale as xii**2 and sqrt(1-yi**2), and the values of xii and yi
c are chosen as powers of 10 (thus, if scaling would be exact, rat should
c be set equal to 10 and sqrt(10)).
c If the test is passed, icount=ithrs, else icount<ithrs; in the former
c case iret=0, in the latter iret=1.
c When the test is not passed, one may choose to stop the program dead here;
c in such a case, set istop=1 below. Each time the test is not passed,
c the results are written onto fort.77; set iwrite=0 to prevent the writing
      implicit none
      real*8 xsecvc(15),xseclvc,wgt(15),wgtl,lxp(0:3,21),xp(15,0:3,21)
      real*8 ckc(15),rckc(15),rat
      integer iflag,imax,iev,nexternal,i_fks,j_fks,iret,ithrs,istop,
     # iwrite,i,k,l,imin,icount
      parameter (ithrs=3)
      parameter (istop=0)
      parameter (iwrite=1)
c
      if(imax.gt.15)then
        write(6,*)'Error in checkres: imax is too large',imax
        stop
      endif
      do i=1,imax
        if(xseclvc.eq.0.d0)then
          ckc(i)=abs(xsecvc(i))
        else
          ckc(i)=abs(xsecvc(i)/xseclvc-1.d0)
        endif
      enddo
      if(iflag.eq.0)then
        rat=4.d0
      elseif(iflag.eq.1)then
        rat=2.d0
      else
        write(6,*)'Error in checkres: iflag=',iflag
        write(6,*)' Must be 0 for soft, 1 for collinear'
        stop
      endif
c
      i=1
      do while(ckc(i).gt.0.1d0 .and. xseclvc.ne.0d0)
        i=i+1
      enddo
      imin=i
      do i=imin,imax-1
        if(ckc(i+1).ne.0.d0)then
          rckc(i)=ckc(i)/ckc(i+1)
        else
          rckc(i)=1.d8
        endif
      enddo
      icount=0
      i=imin
      do while(icount.lt.ithrs.and.i.lt.imax)
        if(rckc(i).gt.rat)then
          icount=icount+1
        else
          icount=0
        endif
        i=i+1
      enddo
c
      iret=0
      if(icount.ne.ithrs)then
        iret=1
        if(istop.eq.1)then
          write(6,*)'Test failed',iflag
          write(6,*)'Event #',iev
          stop
        endif
        if(iwrite.eq.1)then
          write(77,*)'    '
          if(iflag.eq.0)then
            write(77,*)'Soft #',iev
          elseif(iflag.eq.1)then
            write(77,*)'Collinear #',iev
          endif
          write(77,*)'ME*wgt:'
          do i=1,imax
             call xprintout(77,xsecvc(i),xseclvc)
          enddo
          write(77,*)'wgt:'
          do i=1,imax
             call xprintout(77,wgt(i),wgtl)
          enddo
c
          write(78,*)'    '
          if(iflag.eq.0)then
            write(78,*)'Soft #',iev
          elseif(iflag.eq.1)then
            write(78,*)'Collinear #',iev
          endif
          do k=1,nexternal
            write(78,*)''
            write(78,*)'part:',k
            do l=0,3
              write(78,*)'comp:',l
              do i=1,imax
                call xprintout(78,xp(i,l,k),lxp(l,k))
              enddo
            enddo
          enddo
          if(iflag.eq.0)then
            write(78,*)''
            write(78,*)'part: i_fks reduced'
            do l=0,3
              write(78,*)'comp:',l
              do i=1,imax
                call xprintout(78,xp(i,l,nexternal+1),
     #                            lxp(l,nexternal+1))
              enddo
            enddo
            write(78,*)''
            write(78,*)'part: i_fks full/reduced'
            do l=0,3
              write(78,*)'comp:',l
              do i=1,imax
                call xprintout(78,xp(i,l,i_fks),
     #                            xp(i,l,nexternal+1))
              enddo
            enddo
          elseif(iflag.eq.1)then
            write(78,*)''
            write(78,*)'part: i_fks+j_fks'
            do l=0,3
              write(78,*)'comp:',l
              do i=1,imax
                call xprintout(78,xp(i,l,i_fks)+xp(i,l,j_fks),
     #                            lxp(l,i_fks)+lxp(l,j_fks))
              enddo
            enddo
          endif
        endif
      endif
      return
      end




      subroutine checksij(xsijvc,xsijlvc,xsijlim,
     #                    xsumvc,xsumlvc,xsumlim,
     #                    check,checkl,tolerance,
     #                    iflag,imax,iev,ki,kk,ll,
     #                    i_fks,j_fks,ilim,iret)
c Analogous to checkres. Relevant to S functions
      implicit none
      real*8 xsijvc(15),xsijlvc,xsumvc(15),xsumlvc,check(15),checkl
      real*8 xsijlim,xsumlim,tolerance
      real*8 xsecvc(15),xseclvc
      real*8 ckc(15),rckc(15),rat
      logical found
      integer iflag,imax,iev,ki,kk,ll,i_fks,j_fks,ilim,iret,ithrs,
     # istop,iwrite,i,imin,icount,itype
      parameter (ithrs=3)
      parameter (istop=0)
      parameter (iwrite=1)
c
      if(imax.gt.15)then
        write(6,*)'Error in checksij: imax is too large',imax
        stop
      endif
      itype=1
      iret=0
 100  continue
      if(itype.eq.1)then
        do i=1,imax
          xsecvc(i)=xsijvc(i)
        enddo
        xseclvc=xsijlvc
      elseif(itype.eq.2)then
        do i=1,imax
          xsecvc(i)=xsumvc(i)
        enddo
        xseclvc=xsumlvc
      else
        write(6,*)'Error in checksij: itype=',itype
        stop
      endif
      do i=1,imax
        if(xseclvc.eq.0.d0)then
          ckc(i)=abs(xsecvc(i))
        else
          ckc(i)=abs(xsecvc(i)/xseclvc-1.d0)
        endif
      enddo
      if(iflag.eq.0)then
        rat=8.d0
      elseif(iflag.eq.1)then
        rat=2.d0
      else
        write(6,*)'Error in checksij: iflag=',iflag
        write(6,*)' Must be 0 for soft, 1 for collinear'
        stop
      endif
c
      i=1
      dowhile(ckc(i).gt.0.1d0)
        i=i+1
      enddo
      imin=i
      do i=imin,imax-1
        if(ckc(i+1).gt.1.d-8)then
c If this condition is replaced by .eq.0, the test will fail if the series
c is made of elements all equal to the limit
          rckc(i)=ckc(i)/ckc(i+1)
        else
c Element #i+1 of series equal to the limit, so it must pass the test
          rckc(i)=rat*1.1d0
        endif
      enddo
      icount=0
      i=imin
      dowhile(icount.lt.ithrs.and.i.lt.imax)
        if(rckc(i).gt.rat)then
          icount=icount+1
        else
          icount=0
        endif
        i=i+1
      enddo
c
      if(icount.ne.ithrs)then
        iret=iret+itype
        if(istop.eq.1)then
          write(6,*)'Test failed',iflag
          write(6,*)'Event #',iev
          stop
        endif
      endif
      if(itype.eq.1.and.ki.eq.1.and.iflag.eq.0)then
        itype=2
        goto 100
      endif
c
      if(ki.eq.1.and.ilim.eq.1)then
        found=.false.
        i=0
        do while ((.not.found).and.i.lt.imax)
          i=i+1
          if(abs(check(i)-1.d0).gt.tolerance)then
            found=.true.
            itype=4
          endif
        enddo
        if(.not.found)then
          if(abs(checkl-1.d0).gt.tolerance)itype=4
        endif
        if(itype.eq.4)iret=iret+itype
      endif
c
      if( iwrite.eq.1 .and.
     #    iret.eq.1 .or.(iret.gt.1.and.ki.eq.1) )then
        if(iret.gt.7)then
          write(6,*)'Error in checksij: iret=',iret
          stop
        endif
        write(77,*)'    '
        if(iflag.eq.0)then
          write(77,*)'Soft #',iev
        elseif(iflag.eq.1)then
          write(77,*)'Collinear #',iev
        endif
        write(77,*)'iret:',iret
        write(77,*)'i_fks,j_fks:',i_fks,j_fks
        if(iret.eq.1.or.iret.eq.3.or.iret.eq.5.or.iret.eq.7)then
          write(77,*)'S_kl'
          write(77,*)'k,kk,ll',ki,kk,ll
          do i=1,imax
             call xprintout(77,xsijvc(i),xsijlvc)
          enddo
        endif
        if(iret.eq.2.or.iret.eq.3.or.iret.eq.6.or.iret.eq.7)then
          write(77,*)'sum of S'
          do i=1,imax
             call xprintout(77,xsumvc(i),xsumlvc)
          enddo
        endif
        if(iret.eq.4.or.iret.eq.5.or.iret.eq.6.or.iret.eq.7)then
          write(77,*)'check to one'
          do i=1,imax
             call xprintout(77,check(i),checkl)
          enddo
        endif
      endif
c
      if(ilim.eq.1)then
        if( abs(xsijlvc-xsijlim).gt.1.d-6 .and. 
     #    xsijlim.ne.-1.d0 )iret=iret+10
        if( abs(xsumlvc-xsumlim).gt.1.d-6 .and.
     #    xsumlim.ne.-1.d0 .and. iflag.eq.0)iret=iret+20
        if(iwrite.eq.1.and.iret.ge.10)then
          write(77,*)'    '
          if(iflag.eq.0)then
            write(77,*)'Soft #',iev
          elseif(iflag.eq.1)then
            write(77,*)'Collinear #',iev
          endif
          write(77,*)'iret:',iret
          write(77,*)'i_fks,j_fks:',i_fks,j_fks
          if((iret.ge.10.and.iret.lt.20).or.iret.ge.30)then
            write(77,*)'limit of S_kl'
            write(77,*)'k,kk,ll',ki,kk,ll
            write(77,*)xsijlvc,xsijlim
          endif
          if(iret.ge.20)then
            write(77,*)'limit of sum_j S_ij'
            write(77,*)xsumlvc,xsumlim
          endif
        endif
      endif
      return
      end


      subroutine bornsoftvirtual(p,bsv_wgt,virt_wgt,born_wgt)
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
      include "coupl.inc"
      include 'q_es.inc'
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      double precision particle_charge(nexternal)
      common /c_charges/particle_charge
      include "run.inc"
      include "fks_powers.inc"
      include 'reweight.inc'
      double precision p(0:3,nexternal),bsv_wgt,born_wgt
      double precision pp(0:3,nexternal)
      
      double precision wgt1
      double precision rwgt,Q,Ej,wgt,contr,eikIreg,m1l_W_finite_CDR
      double precision aso2pi, aeo2pi
      double precision shattmp,dot
      integer i,j,aj,m,n,k,iord,ipos_ord

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      double precision xicut_used
      common /cxicut_used/xicut_used

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      double precision pi
      parameter (pi=3.1415926535897932385d0)

      double precision c(0:1),gamma(0:1),gammap(0:1)
      common/fks_colors/c,gamma,gammap
      double precision c_used, gamma_used, gammap_used
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      double precision double,single,xmu2
      logical ComputePoles,fksprefact
      parameter (ComputePoles=.false.)
      parameter (fksprefact=.true.)

      double precision beta0,ren_group_coeff
      common/cbeta0/beta0,ren_group_coeff

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn

      double precision virt_fraction_inc
      data virt_fraction_inc /1d0/

      integer include_virt
      double precision vol3

c For tests of virtuals
      double precision vobmax,vobmin
      common/cvirt0test/vobmax,vobmin
      double precision vNsumw,vAsumw,vSsumw,vNsumf,vAsumf,vSsumf
      common/cvirt1test/vNsumw,vAsumw,vSsumw,vNsumf,vAsumf,vSsumf
      integer nvtozero
      logical doVirtTest
      common/cvirt2test/nvtozero,doVirtTest
      double precision xnormsv
      common/cxnormsv/xnormsv
      double precision vrat

      double precision virt_wgt,ran2
      external ran2

      character*4 abrv
      common /to_abrv/ abrv

      logical ExceptPSpoint
      integer iminmax
      common/cExceptPSpoint/iminmax,ExceptPSpoint

      double precision average_virtual,virtual_fraction
      common/c_avg_virt/average_virtual,virtual_fraction
      double precision virtual_over_born
      common/c_vob/virtual_over_born

c timing statistics
      include "timing_variables.inc"

c For the MINT folding
      integer fold
      common /cfl/fold
      double precision virt_wgt_save
      save virt_wgt_save

      double precision pmass(nexternal),zero,tiny
      parameter (zero=0d0)
      parameter (tiny=1d-6)
      include 'orders.inc'
      logical firsttime
      data firsttime / .true. /
      logical need_color_links_used, need_charge_links_used
      data need_color_links_used / .false. /
      data need_charge_links_used / .false. /
      logical split_type(nsplitorders) 
      common /c_split_type/split_type
      logical split_type_used(nsplitorders)
      data split_type_used /nsplitorders * .false./
      logical need_color_links, need_charge_links
      common /c_need_links/need_color_links, need_charge_links
      complex*16 ans_cnt(2, nsplitorders)
      common /c_born_cnt/ ans_cnt
      double precision oneo8pi2
      parameter(oneo8pi2 = 1d0/(8d0*pi**2))
      include 'nFKSconfigs.inc'
      INTEGER nFKSprocess, nFKSprocess_save
      COMMON/c_nFKSprocess/nFKSprocess
      include "pmass.inc"
      
      if (firsttime) then
C check if any real emission need cahrge/color links
         nFKSprocess_save = nFKSprocess
         do nFKSprocess = 1, FKS_configs
            call fks_inc_chooser()
            need_color_links_used = need_color_links_used .or. need_color_links
            need_charge_links_used = need_charge_links_used .or. need_charge_links
            do i=1, nsplitorders
              split_type_used(i) = split_type_used(i) .or. split_type(i)
            enddo
         enddo
         if (need_charge_links_used) then
             write(*,*) 'Charge-linked born are used'
         else
             write(*,*) 'Charge-linked born are not used'
         endif
         if (need_color_links_used) then
             write(*,*) 'Color-linked born are used'
         else
             write(*,*) 'Color-linked born are not used'
         endif
         firsttime = .false.
         nFKSprocess = nFKSprocess_save
         call fks_inc_chooser()
      endif
         

      aso2pi=g**2/(8*pi**2)
      aeo2pi=dble(gal(1))**2/(8*pi**2)

      if (.not.(need_color_links_used.or.need_charge_links_used
     #    .or.abrv.eq.'grid')) then
C just return 0
         bsv_wgt=0d0
         virt_wgt=0d0
         born_wgt=0d0
         if(doreweight)then
           wgtnstmp=0d0
           wgtwnstmpmuf=0d0
           wgtwnstmpmur=0d0
         endif
         goto 999
      endif

c Consistency check -- call to set_cms_stuff() must be done prior to
c entering this function
      shattmp=2d0*dot(p(0,1),p(0,2))
      if(abs(shattmp/shat-1.d0).gt.1.d-5)then
         write(*,*)'Error in bornsoftvirtual: inconsistent shat'
         write(*,*)shattmp,shat
         stop
      endif

      call sborn(p_born,wgt1)

c Born contribution:
      bsv_wgt=wgt1
      born_wgt=wgt1
      virt_wgt=0d0

      if (abrv.eq.'born' .or. abrv.eq.'grid') goto 549
      if (abrv.eq.'virt' .or. abrv.eq.'viSC' .or.
     #       abrv.eq.'viLC') goto 547

c Q contribution eq 5.5 and 5.6 of FKS
C loop over QCD/QED (iord=1,2 respectively)
      do iord= 1,2
         Q=0d0
C skip what we don't need
         if (iord.eq.1) ipos_ord = qcd_pos
         if (iord.eq.2) ipos_ord = qed_pos
         if (.not.split_type_used(ipos_ord)) cycle
         do i=1 ,nexternal
            if (i.ne.i_fks .and. pmass(i).eq.ZERO) then
c set the various color factors according to the 
c type of the leg
               if (particle_type(i).eq.8) then
                  aj=0
               elseif(abs(particle_type(i)).eq.3) then
                  aj=1
               else
                  aj=-1
               endif
               Ej=p(0,i)
               if (pdg_type(i).eq.22) then
                   write(*,*) 'NOT IMPLEMENTED'
                   stop
               endif

               if (ipos_ord.eq.qcd_pos) then
C                 set colour factors
                  if (aj.eq.-1) cycle
                  c_used = c(aj)
                  gamma_used = gamma(aj)
                  gammap_used = gammap(aj)
               else if (ipos_ord.eq.qed_pos) then
                  if (particle_charge(i).eq.0d0) cycle
C                 set charge factors
                  c_used = particle_charge(i)**2
                  gamma_used = 3d0/2d0 * particle_charge(i)**2
                  gammap_used = (13d0/2d0 - 2d0 * pi**2 / 3d0) * particle_charge(i)**2
               endif

               if (i.gt.nincoming) then 
C Q terms for final state partons
                if(abrv.eq.'novA')then
c 2+3+4
                  Q = Q
     &             -2*dlog(shat/QES2)*dlog(xicut_used)*c_used
     &             -( dlog(deltaO/2d0)*( gamma_used-
     &                     2d0*c_used*dlog(2d0*Ej/xicut_used/sqrtshat) )
     &               +2*dlog(xicut_used)**2*c_used )
     &             +gammap_used
     &             +2d0*c_used*dlog(2d0*Ej/sqrtshat)**2
     &             -2d0*gamma_used*dlog(2d0*Ej/sqrtshat)
                elseif(abrv.eq.'novB')then
c 2+3+4_mu
                  Q = Q
     &             -2*dlog(shat/QES2)*dlog(xicut_used)*c_used
     &             -( dlog(deltaO/2d0)*( gamma_used-
     &                     2d0*c_used*dlog(2d0*Ej/xicut_used/sqrtshat) )
     &               +2*dlog(xicut_used)**2*c_used )
                elseif(abrv.eq.'viSA')then
c 1                
                  Q = Q
     &              -dlog(shat/QES2)*( gamma_used-
     &                      2d0*c_used*dlog(2d0*Ej/sqrtshat) )
                elseif(abrv.eq.'viSB')then
c 1+4_L
                  Q = Q
     &              -dlog(shat/QES2)*( gamma_used-
     &                      2d0*c_used*dlog(2d0*Ej/sqrtshat) )
     &             +gammap_used
     &             +2d0*c_used*dlog(2d0*Ej/sqrtshat)**2
     &             -2d0*gamma_used*dlog(2d0*Ej/sqrtshat)
                elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
     #                abrv.ne.'viLC')then
c 1+2+3+4
                  Q = Q+gammap_used
     &              -dlog(shat*deltaO/2d0/QES2)*( gamma_used-
     &                     2d0*c_used*dlog(2d0*Ej/xicut_used/sqrtshat) )
     &              +2d0*c_used*( dlog(2d0*Ej/sqrtshat)**2
     &              -dlog(xicut_used)**2 )
     &              -2d0*gamma_used*dlog(2d0*Ej/sqrtshat)
                else
                  write(*,*)'Error in bornsoftvirtual'
                  write(*,*)'abrv in Q:',abrv
                  stop
                endif

               else
C Q terms for initial state partons
                if(abrv.eq.'novA'.or.abrv.eq.'novB')then
c 2+3+4 or 2+3+4_mu
                    Q=Q-2*dlog(shat/QES2)*dlog(xicut_used)*c_used
     &              -dlog(q2fact(i)/shat)*(
     &              gamma_used+2d0*c_used*dlog(xicut_used) )
                elseif(abrv.eq.'viSA'.or.abrv.eq.'viSB')then
c 1 or 1+4_L
                    Q=Q-dlog(shat/QES2)*gamma_used
                elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
     #             abrv.ne.'viLC')then
c 1+2+3+4
                    Q=Q-dlog(q2fact(i)/QES2)*(
     &               gamma_used+2d0*c_used*dlog(xicut_used))
                else
                    write(*,*)'Error in bornsoftvirtual'
                    write(*,*)'abrv in Q:',abrv
                    stop
                endif
               endif
            endif
         enddo
C end of the external particle loop
         if (ipos_ord.eq.qcd_pos) bsv_wgt =
     $      bsv_wgt+aso2pi*Q*dble(ans_cnt(1,qcd_pos))
         if (ipos_ord.eq.qed_pos) bsv_wgt =
     $      bsv_wgt+aeo2pi*Q*dble(ans_cnt(1,qed_pos))

      enddo


c     If doing MC over helicities, must sum over the two
c     helicity contributions for the Q-terms of collinear limit.
 547  continue
      if (abrv.eq.'virt' .or. abrv.eq.'viSC' .or.
     #    abrv.eq.'viLC') goto 548
c
c I(reg) terms, eq 5.5 of FKS
      do iord = 1, nsplitorders
        if (iord.eq.qcd_pos) then
            if (.not. need_color_links_used) cycle
            need_color_links=need_color_links_used
            need_charge_links=.false.
        else if (iord.eq.qed_pos) then
            if (.not. need_charge_links_used) cycle
            need_charge_links=need_charge_links_used
            need_color_links=.false.
        else
            cycle
        endif
        contr=0d0
        do i=1,fks_j_from_i(i_fks,0)
           do j=1,i
              m=fks_j_from_i(i_fks,i)
              n=fks_j_from_i(i_fks,j)
              if ((m.ne.n .or. (m.eq.n .and. pmass(m).ne.ZERO)).and.
     &           n.ne.i_fks.and.m.ne.i_fks) then
c To be sure that color-correlated Borns work well, we need to have
c *always* a call to sborn(p_born,wgt) just before. This is okay,
c because there is a call above in this subroutine
C wgt includes the gs/w^2 
                 call sborn_sf(p_born,m,n,wgt)
                 if (wgt.ne.0d0) then
                    call eikonal_Ireg(p,m,n,xicut_used,eikIreg)
                    contr=contr+wgt*eikIreg
                 endif
              endif
           enddo
        enddo

C WARNING: THE FACTOR -2 BELOW COMPENSATES FOR THE MISSING -2 IN THE
C COLOUR LINKED BORN -- SEE ALSO SBORNSOFT().
C If the colour-linked Borns were normalized as reported in the paper
c we should set
c   bsv_wgt=bsv_wgt+ao2pi*contr  <-- DO NOT USE THIS LINE
c
        bsv_wgt=bsv_wgt-2*oneo8pi2*contr
      enddo

      call fks_inc_chooser()

 548  continue
c Finite part of one-loop corrections
c convert to Binoth Les Houches Accord standards
      virt_wgt=0d0
      if (fold.eq.0) then
         if ((ran2().le.virtual_fraction .and.
     $           abrv(1:3).ne.'nov').or.abrv(1:4).eq.'virt') then
               call cpu_time(tBefore)
               Call BinothLHA(p_born,born_wgt,virt_wgt)
c$$$               virt_wgt=m1l_W_finite_CDR(p_born,born_wgt)
               call cpu_time(tAfter)
               tOLP=tOLP+(tAfter-tBefore)
               virtual_over_born=virt_wgt/(born_wgt*aso2pi)
               virt_wgt=(virt_wgt-average_virtual*born_wgt*aso2pi)
               if (abrv.ne.'virt') then
                  virt_wgt=virt_wgt/virtual_fraction
               endif
               virt_wgt_save=virt_wgt
         endif
      elseif(fold.eq.1) then
         virt_wgt=virt_wgt_save
c$$$            bsv_wgt=bsv_wgt+virt_wgt_save
      endif
      if (abrv(1:4).ne.'virt')
     &        bsv_wgt=bsv_wgt+average_virtual*born_wgt*aso2pi

c eq.(MadFKS.C.13)
      if(abrv.eq.'viSA'.or.abrv.eq.'viSB')then
          write(*,*) 'FIX visA, visB'
C        bsv_wgt=bsv_wgt + 2*pi*beta0*wgtbpower*log(shat/QES2)*
C     #                       ao2pi*dble(wgt1(1))
C      elseif(abrv.eq.'novA'.or.abrv.eq.'novB')then
C        bsv_wgt=bsv_wgt + 2*pi*beta0*wgtbpower*log(q2fact(1)/shat)*
C     #                       ao2pi*dble(wgt1(1))
C      elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
C     #       abrv.ne.'viLC')then
C         bsv_wgt=bsv_wgt + 2*pi*beta0*wgtbpower*log(q2fact(1)/QES2)*
C     #                       ao2pi*dble(wgt1(1))
      endif
c eq.(MadFKS.C.14)
      if(abrv(1:2).ne.'vi')then
        
        if (dabs(log(q2fact(1)/scale**2)).gt.1d-6) then
            ! fix the next line and remove the if statement (all points)
            write(*,*) 'FIX muR!=muF'
C            bsv_wgt=bsv_wgt - 2*pi*beta0*wgtbpower*
C     #                       log(q2fact(1)/scale**2)*ao2pi*dble(wgt1(1))
        endif
      endif


 549  continue

      if(doreweight)then
        write(*,*) 'FIX REWEIGHT!!!!!'
        ! the following lines have to be uncommented and fixed
C        wgtwnstmpmuf=0.d0
C        if(abrv.ne.'born' .and. abrv.ne.'grid')then
C          if(abrv(1:2).eq.'vi')then
C            wgtwnstmpmur=0.d0
C          else
C               do i=1,nincoming
C               if (particle_type(i).ne.1)then
C                  if (particle_type(i).eq.8) then
C                    aj=0
C                   elseif(abs(particle_type(i)).eq.3) then
C                    aj=1
C                   endif
C                wgtwnstmpmuf=wgtwnstmpmuf-
C     #                   ( gamma(aj)+2d0*c(aj)*dlog(xicut_used) )
C              endif
C               enddo
C               wgtwnstmpmuf=ao2pi*wgtwnstmpmuf*dble(wgt1(1))
C               wgtwnstmpmur=2*pi*(beta0*wgtbpower
C     #      +ren_group_coeff*wgtcpower)*ao2pi*dble(wgt1(1))
C             endif
c bsv_wgt here always contains the Born; must subtract it, since 
c we need the pure NLO terms only
C             wgtnstmp=bsv_wgt+virt_wgt-born_wgt-
C     #                wgtwnstmpmuf*log(q2fact(1)/QES2)-
C     #                wgtwnstmpmur*log(scale**2/QES2)
C           else
C             wgtnstmp=0d0
C             wgtwnstmpmur=0.d0
C           endif
         endif

      if (abrv(1:2).eq.'vi') then
         bsv_wgt=bsv_wgt-born_wgt

         if(doVirtTest .and. iminmax.eq.0)then
           if(born_wgt.ne.0.d0)then
             vrat=bsv_wgt/(aso2pi*born_wgt)
             if(vrat.gt.vobmax)vobmax=vrat
             if(vrat.lt.vobmin)vobmin=vrat
             vNsumw=vNsumw+xnormsv
             vAsumw=vAsumw+vrat*xnormsv
             vSsumw=vSsumw+vrat**2*xnormsv
             vNsumf=vNsumf+1.d0
             vAsumf=vAsumf+vrat
             vSsumf=vSsumf+vrat**2
           else
             if(bsv_wgt.ne.0.d0)nvtozero=nvtozero+1
           endif
         endif

         born_wgt=0d0
      endif

      if (ComputePoles) then
         call sborn(p_born,wgt1)

         print*,"           "
         write(*,123)((p(i,j),i=0,3),j=1,nexternal)
         xmu2=q2fact(1)
         call getpoles(p,xmu2,double,single,fksprefact)
         print*,"BORN",born_wgt!/conv
         print*,"DOUBLE",double
         print*,"SINGLE",single
c         print*,"LOOP",virt_wgt!/born_wgt/ao2pi*2d0
c         print*,"LOOP2",(virtcor+born_wgt*4d0/3d0-double*pi**2/6d0)
c         stop
 123     format(4(1x,d22.16))
      endif


 999  continue
      return
      end


      subroutine compute_bpower(p_born,bpower)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"

      double precision p_born(0:3,nexternal-1)
      double precision bpower,born_wgt
      double precision wgt1

      integer           isum_hel
      logical                   multi_channel
      common/to_matrix/isum_hel, multi_channel
      integer isum_hel_orig
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn

      double precision tiny
      parameter (tiny=1d-6)

c Make sure that we sum over helicities (such that we do get a
c non-zero Born)
      isum_hel_orig = isum_hel
      isum_hel=0
      call get_helicity(i_fks,j_fks)

      calculatedBorn=.false.
      call sborn(p_born,wgt1)
c Born contribution:
      born_wgt=wgt1
      
c Multiply the strong coupling by 10
      if (g.ne.0d0) then
         g=10d0*g
      else
         write(*,*)'Error in bornsoftvirtual'
         write(*,*)'Strong coupling is zero'
         stop
      endif

c Update alphaS-dependent couplings
      call update_as_param()

c recompute the Born with the new couplings
      calculatedBorn=.false.
      call sborn(p_born,wgt1)

c Compute bpower
      bpower=Log10(wgt1/born_wgt)/2d0
      if(abs(bpower-dble(nint(bpower))) .gt. tiny) then
         write(*,*)'Error in computation of bpower:'
         write(*,*)' not an integer',bpower
         stop
      elseif (bpower.lt.-tiny) then
         write(*,*)'Error in computation of bpower:'
         write(*,*)' negative value',bpower
         stop
      else
c set it to the integer exactly
         bpower=dble(nint(bpower))
         write(*,*)'bpower is', bpower
      endif

c Change couplings back and recompute the Born to make sure that 
c nothing funny happens later on
      g=g/10d0
      call update_as_param()
      isum_hel=isum_hel_orig
      calculatedBorn=.false.
      call sborn(p_born,wgt1)

      return
      end

c       This function computes the power of a muR-dependent factor which
c       is stored in cpower. You need to modify it when you try to 
c       reweight your cross section with a muR-dependent factor
c       (runfac=1 in reweight0.inc)
c Note: The implementation below only works for the Bottom Yukawa in
c       the SM where "GC_33" contains the Yukawa, for other models
c       or general muR-dependent factors you need to change GC_33
c       to the corresponding coupling.
      subroutine compute_cpower(p_born,cpower)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      include 'reweight.inc'

      double precision p_born(0:3,nexternal-1)
      double precision cpower,born_wgt
      double complex wgt1(2)

      integer isum_hel
      logical multi_channel
      common/to_matrix/isum_hel, multi_channel
      integer isum_hel_orig
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn

      double precision tiny
      parameter (tiny=1d-6)
c comment these lines to calculate cpower
      cpower = -1d0
      return
c comment these lines to calculate cpower

c   The following is relevant for a muR-dependent bottom-mass in Yukawa.
c$$$
c$$$c Make sure that we sum over helicities (such that we do get a
c$$$c non-zero Born)
c$$$      isum_hel_orig = isum_hel
c$$$      isum_hel=0
c$$$      call get_helicity(i_fks,j_fks)
c$$$
c$$$      calculatedBorn=.false.
c$$$      call sborn(p_born,wgt1)
c$$$c Born contribution:
c$$$      born_wgt=dble(wgt1(1))
c$$$      
c$$$c Multiply the Yukawa by 10 (If you use this,
c$$$c double check that GC_33 is the yukawa! (also below))
c$$$      if (GC_33.ne.0d0) then
c$$$         GC_33 = GC_33 * 10d0
c$$$      else
c$$$         write(*,*)'Warning In Bornsoftvirtual'
c$$$         Write(*,*)'Yukawa Is Zero - Cpower Set To Zero'
c$$$         Cpower = 0d0
c$$$         Return
c$$$      Endif
c$$$
c$$$c recompute the Born with the new Yukawa
c$$$      calculatedBorn=.false.
c$$$      call sborn(p_born,wgt1)
c$$$
c$$$c Compute cpower
c$$$      cpower=Log10(dble(wgt1(1))/born_wgt)
c$$$      if(abs(cpower-dble(nint(cpower))) .gt. tiny) then
c$$$         write(*,*)'Error in computation of cpower:'
c$$$         write(*,*)' not an integer',cpower
c$$$         stop
c$$$      elseif (cpower.lt.-tiny) then
c$$$         write(*,*)'Error in computation of cpower:'
c$$$         write(*,*)' negative value',cpower
c$$$         stop
c$$$      else
c$$$c set it to the integer exactly
c$$$         cpower=dble(nint(cpower))
c$$$         write(*,*)'cpower is', cpower
c$$$c Check consistency with value used in reweighting
c$$$c$$$         if( (doreweight.or.doNLOreweight) .and.
c$$$c$$$     &        abs(cpower-wgtcpower).gt.tiny )then
c$$$c$$$            write(*,*)'Error in compute_cpower'
c$$$c$$$            write(*,*)'cpower(s) are:',cpower,wgtcpower
c$$$c$$$            stop
c$$$c$$$         endif
c$$$      endif
c$$$
c$$$c Change couplings back and recompute the Born to make sure that 
c$$$c nothing funny happens later on
c$$$      GC_33 = GC_33 / 10d0
c$$$      isum_hel=isum_hel_orig
c$$$      calculatedBorn=.false.
c$$$      call sborn(p_born,wgt1)
c$$$
c$$$      return
      end


      subroutine eikonal_Ireg(p,m,n,xicut_used,eikIreg)
      implicit none
      double precision zero,pi,pi2
      parameter (zero=0.d0)
      parameter (pi=3.1415926535897932385d0)
      parameter (pi2=pi**2)
      include "nexternal.inc"
      include 'coupl.inc'
      include 'q_es.inc'
      double precision p(0:3,nexternal),xicut_used,eikIreg
      integer m,n

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      character*4 abrv
      common /to_abrv/ abrv

      double precision Ei,Ej,kikj,rij,tmp,xmj,betaj,betai,xmi2,xmj2,
     # vij,xi0,alij,tHVvl,tHVv,arg1,arg2,arg3,arg4,xi1a,xj1a,
     # dot,ddilog
      external dot,ddilog

      double precision pmass(nexternal)
      include "pmass.inc"

      tmp=0.d0
      if(pmass(m).eq.0.d0.and.pmass(n).eq.0.d0)then
        if(m.eq.n)then
          write(*,*)'Error #2 in eikonal_Ireg',m,n
          stop
        endif
        Ei=p(0,n)
        Ej=p(0,m)
        kikj=dot(p(0,n),p(0,m))
        rij=kikj/(2*Ei*Ej)
        if(abs(rij-1.d0).gt.1.d-6)then
          if(abrv.eq.'novA')then
c 2+3+4
            tmp=2*dlog(shat/QES2)*dlog(xicut_used)+
     #          2*dlog(xicut_used)**2+
     #          2*dlog(xicut_used)*dlog(rij)-
     #          ddilog(rij)+1d0/2d0*dlog(rij)**2-
     #          dlog(1-rij)*dlog(rij)
          elseif(abrv.eq.'novB')then
c 2+3+4_mu
            tmp=2*dlog(shat/QES2)*dlog(xicut_used)+
     #          2*dlog(xicut_used)**2+
     #          2*dlog(xicut_used)*dlog(rij)
          elseif(abrv.eq.'viSA')then
c 1                
            tmp=1d0/2d0*dlog(shat/QES2)**2+
     #          dlog(shat/QES2)*dlog(rij)
          elseif(abrv.eq.'viSB')then
c 1+4_L
            tmp=1d0/2d0*dlog(shat/QES2)**2+
     #          dlog(shat/QES2)*dlog(rij)-
     #          ddilog(rij)+1d0/2d0*dlog(rij)**2-
     #          dlog(1-rij)*dlog(rij)
          elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
     #           abrv.ne.'viLC')then
c 1+2+3+4
            tmp=1d0/2d0*dlog(xicut_used**2*shat/QES2)**2+
     #          dlog(xicut_used**2*shat/QES2)*dlog(rij)-
     #          ddilog(rij)+1d0/2d0*dlog(rij)**2-
     #          dlog(1-rij)*dlog(rij)
          else
             write(*,*)'Error #11 in eikonal_Ireg',abrv
             stop
          endif
        else
          if(abrv.eq.'novA')then
c 2+3+4
            tmp=2*dlog(shat/QES2)*dlog(xicut_used)+
     #          2*dlog(xicut_used)**2-pi2/6.d0
          elseif(abrv.eq.'novB')then
c 2+3+4_mu
            tmp=2*dlog(shat/QES2)*dlog(xicut_used)+
     #          2*dlog(xicut_used)**2
          elseif(abrv.eq.'viSA')then
c 1                
            tmp=1d0/2d0*dlog(shat/QES2)**2
          elseif(abrv.eq.'viSB')then
c 1+4_L
            tmp=1d0/2d0*dlog(shat/QES2)**2-pi2/6.d0
          elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
     #           abrv.ne.'viLC')then
c 1+2+3+4
            tmp=1d0/2d0*dlog(xicut_used**2*shat/QES2)**2-pi2/6.d0
          else
             write(*,*)'Error #12 in eikonal_Ireg',abrv
             stop
          endif
        endif
      elseif( (pmass(m).ne.0.d0.and.pmass(n).eq.0.d0) .or.
     #        (pmass(m).eq.0.d0.and.pmass(n).ne.0.d0) )then
        if(m.eq.n)then
          write(*,*)'Error #3 in eikonal_Ireg',m,n
          stop
        endif
        if(pmass(m).ne.0.d0.and.pmass(n).eq.0.d0)then
          Ei=p(0,n)
          Ej=p(0,m)
          xmj=pmass(m)
          betaj=sqrt(1-xmj**2/Ej**2)
        else
          Ei=p(0,m)
          Ej=p(0,n)
          xmj=pmass(n)
          betaj=sqrt(1-xmj**2/Ej**2)
        endif
        kikj=dot(p(0,n),p(0,m))
        rij=kikj/(2*Ei*Ej)

        if(abrv.eq.'novA')then
c 2+3+4
          tmp=dlog(xicut_used)*dlog(shat/QES2)+
     #        dlog(xicut_used)**2+
     #        2*dlog(xicut_used)*dlog(kikj/(xmj*Ei))-
     #        ddilog(1-(1+betaj)/(2*rij))+ddilog(1-2*rij/(1-betaj))+
     #        1/2.d0*log(2*rij/(1-betaj))**2-pi2/12.d0-
     #        1/4.d0*dlog((1+betaj)/(1-betaj))**2
        elseif(abrv.eq.'novB')then
c 2+3+4_mu
          tmp=dlog(xicut_used)*dlog(shat/QES2)+
     #        dlog(xicut_used)**2+
     #        2*dlog(xicut_used)*dlog(kikj/(xmj*Ei))
        elseif(abrv.eq.'viSA')then
c 1                
          tmp=1/4.d0*dlog(shat/QES2)**2+
     #        dlog(shat/QES2)*dlog(kikj/(xmj*Ei))
        elseif(abrv.eq.'viSB')then
c 1+4_L
          tmp=1/4.d0*dlog(shat/QES2)**2+
     #        dlog(shat/QES2)*dlog(kikj/(xmj*Ei))-
     #        ddilog(1-(1+betaj)/(2*rij))+ddilog(1-2*rij/(1-betaj))+
     #        1/2.d0*log(2*rij/(1-betaj))**2-pi2/12.d0-
     #        1/4.d0*dlog((1+betaj)/(1-betaj))**2
        elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
     #         abrv.ne.'viLC')then
c 1+2+3+4
          tmp=dlog(xicut_used)*( dlog(xicut_used*shat/QES2)+
     #                           2*dlog(kikj/(xmj*Ei)) )-
     #        ddilog(1-(1+betaj)/(2*rij))+ddilog(1-2*rij/(1-betaj))+
     #        1/2.d0*log(2*rij/(1-betaj))**2+
     #        dlog(shat/QES2)*dlog(kikj/(xmj*Ei))-pi2/12.d0+
     #        1/4.d0*dlog(shat/QES2)**2-
     #        1/4.d0*dlog((1+betaj)/(1-betaj))**2
        else
           write(*,*)'Error #13 in eikonal_Ireg',abrv
           stop
        endif
      elseif(pmass(m).ne.0.d0.and.pmass(n).ne.0.d0)then
        if(n.eq.m)then
          Ei=p(0,n)
          betai=sqrt(1-pmass(n)**2/Ei**2)
          if(abrv.eq.'novA')then
c 2+3+4
            tmp=2*dlog(xicut_used)-
     #          1/betai*dlog((1+betai)/(1-betai))
          elseif(abrv.eq.'novB')then
c 2+3+4_mu
            tmp=2*dlog(xicut_used)
          elseif(abrv.eq.'viSA')then
c 1                
            tmp=dlog(shat/QES2)
          elseif(abrv.eq.'viSB')then
c 1+4_L
            tmp=dlog(shat/QES2)-
     #          1/betai*dlog((1+betai)/(1-betai))
          elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
     #           abrv.ne.'viLC')then
c 1+2+3+4
            tmp=dlog(xicut_used**2*shat/QES2)-
     #          1/betai*dlog((1+betai)/(1-betai))
          else
             write(*,*)'Error #14 in eikonal_Ireg',abrv
             stop
          endif
        else
          Ei=p(0,n)
          Ej=p(0,m)
          betai=sqrt(1-pmass(n)**2/Ei**2)
          betaj=sqrt(1-pmass(m)**2/Ej**2)
          xmi2=pmass(n)**2
          xmj2=pmass(m)**2
          kikj=dot(p(0,n),p(0,m))
          vij=sqrt(1-xmi2*xmj2/kikj**2)
          alij=kikj*(1+vij)/xmi2
          tHVvl=(alij**2*xmi2-xmj2)/2.d0
          tHVv=tHVvl/(alij*Ei-Ej)
          arg1=alij*Ei
          arg2=arg1*betai
          arg3=Ej
          arg4=arg3*betaj
          xi0=1/vij*log((1+vij)/(1-vij))
          xi1a=kikj**2*(1+vij)/xmi2*( xj1a(arg1,arg2,tHVv,tHVvl)-
     #                                xj1a(arg3,arg4,tHVv,tHVvl) )

          if(abrv.eq.'novA')then
c 2+3+4
            tmp=xi0*dlog(xicut_used)+1/2.d0*xi1a
          elseif(abrv.eq.'novB')then
c 2+3+4_mu
            tmp=xi0*dlog(xicut_used)
          elseif(abrv.eq.'viSA')then
c 1                
            tmp=1/2.d0*xi0*dlog(shat/QES2)
          elseif(abrv.eq.'viSB')then
c 1+4_L
            tmp=1/2.d0*xi0*dlog(shat/QES2)+1/2.d0*xi1a
          elseif(abrv.ne.'virt' .and. abrv.ne.'viSC' .and.
     #           abrv.ne.'viLC')then
c 1+2+3+4
            tmp=1/2.d0*xi0*dlog(xicut_used**2*shat/QES2)+1/2.d0*xi1a
          else
             write(*,*)'Error #15 in eikonal_Ireg',abrv
             stop
          endif
        endif
      else
        write(*,*)'Error #4 in eikonal_Ireg',m,n,pmass(m),pmass(n)
        stop
      endif
      eikIreg=tmp
      return
      end


      function xj1a(x,y,tHVv,tHVvl)
      implicit none
      real*8 xj1a,x,y,tHVv,tHVvl,ddilog
      external ddilog
c
      xj1a=1/(2*tHVvl)*( dlog((x-y)/(x+y))**2+4*ddilog(1-(x+y)/tHVv)+
     #                   4*ddilog(1-(x-y)/tHVv) )
      return
      end


      FUNCTION DDILOG(X)
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(0:19)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)
      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/
      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DDILOG=H
      RETURN
      END


      subroutine getpoles(p,xmu2,double,single,fksprefact)
c Returns the residues of double and single poles according to 
c eq.(B.1) and eq.(B.2) if fksprefact=.true.. When fksprefact=.false.,
c the prefactor (mu2/Q2)^ep in eq.(B.1) is expanded, and giving an
c extra contribution to the single pole
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      double precision particle_charge(nexternal), particle_charge_born(nexternal-1)
      common /c_charges/particle_charge
      common /c_charges_born/particle_charge_born
      include 'coupl.inc'
      include 'q_es.inc'
      double precision p(0:3,nexternal),xmu2,double,single
      logical fksprefact
      double precision c(0:1),gamma(0:1),gammap(0:1)
      common/fks_colors/c,gamma,gammap
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      double precision wgt1
      double precision born,wgt,kikj,dot,vij,aso2pi,aeo2pi
      double precision contr1, contr2
      integer aj,i,j,m,n,ilink
      double precision pmass(nexternal),zero,pi
      parameter (pi=3.1415926535897932385d0)
      parameter (zero=0d0)
      include 'orders.inc'
      complex*16 ans_cnt(2, nsplitorders)
      common /c_born_cnt/ ans_cnt
      logical need_color_links, need_charge_links
      common /c_need_links/need_color_links, need_charge_links
      double precision oneo8pi2
      parameter(oneo8pi2 = 1d0/(8d0*pi**2))
      include "pmass.inc"
c
      double=0.d0
      single=0.d0
      aso2pi=g**2/(8*pi**2)
      aeo2pi=dble(gal(1))**2/(8*pi**2)
      call sborn(p_born,wgt1)
c QCD Born terms
      contr1 = 0d0
      contr2 = 0d0
      born=dble(ans_cnt(1,qcd_pos))
      do i=1,nexternal
        if(i.ne.i_fks .and. particle_type(i).ne.1)then
          if (particle_type(i).eq.8) then
             aj=0
          elseif(abs(particle_type(i)).eq.3) then
             aj=1
          endif
          if(pmass(i).eq.ZERO)then
            contr2=contr2-c(aj)
            contr1=contr1-gamma(aj)
          else
            contr1=contr1-c(aj)
          endif
        endif
      enddo
      double=double+contr2*born*aso2pi
      single=single+contr1*born*aso2pi

c QED Born terms
      contr1 = 0d0
      contr2 = 0d0
      born=dble(ans_cnt(1,qed_pos))
      do i=1,nexternal
        if(i.ne.i_fks.and.particle_charge(i).ne.0d0)then
          if(pmass(i).eq.ZERO)then
            contr2=contr2-particle_charge(i)**2
            contr1=contr1-3d0/2d0*particle_charge(i)**2
          else
            contr1=contr1-particle_charge(i)**2
          endif
        endif
      enddo
      double=double+contr2*born*aeo2pi
      single=single+contr1*born*aeo2pi

c Colour and charge-linked Born terms
      do ilink = 1, 2
        if (ilink.eq.1) then
          need_color_links = .true.
          need_charge_links = .false.
        else
          need_color_links = .false.
          need_charge_links = .true.
        endif
        contr1=0d0
        do i=1,fks_j_from_i(i_fks,0)
          do j=1,i
            m=fks_j_from_i(i_fks,i)
            n=fks_j_from_i(i_fks,j)
            if( m.ne.n .and. n.ne.i_fks .and. m.ne.i_fks )then
C wgt includes the gs/w^2 factor
              call sborn_sf(p_born,m,n,wgt)
c The factor -2 compensate for that missing in sborn_sf
              wgt=-2*wgt
              if(wgt.ne.0.d0)then
                if(pmass(m).eq.zero.and.pmass(n).eq.zero)then
                  kikj=dot(p(0,n),p(0,m))
                  contr1=contr1+log(2*kikj/QES2)*wgt
                elseif(pmass(m).ne.zero.and.pmass(n).eq.zero)then
                  contr1=contr1-0.5d0*log(pmass(m)**2/QES2)*wgt
                  kikj=dot(p(0,n),p(0,m))
                  contr1=contr1+log(2*kikj/QES2)*wgt
                elseif(pmass(m).eq.zero.and.pmass(n).ne.zero)then
                  contr1=contr1-0.5d0*log(pmass(n)**2/QES2)*wgt
                  kikj=dot(p(0,n),p(0,m))
                  contr1=contr1+log(2*kikj/QES2)*wgt
                elseif(pmass(m).ne.zero.and.pmass(n).ne.zero)then
                  kikj=dot(p(0,n),p(0,m))
                  vij=sqrt(1-(pmass(n)*pmass(m)/kikj)**2)
                  contr1=contr1+0.5d0*1/vij*log((1+vij)/(1-vij))*wgt
                else
                  write(*,*)'Error in getpoles',i,j,n,m,pmass(n),pmass(m)
                  stop
                endif
              endif
            endif
          enddo
        enddo
        single=single+contr1*oneo8pi2
      enddo

C restore need_color/charge_links
      call fks_inc_chooser()

      if(.not.fksprefact)single=single+double*log(xmu2/QES2)
c
      return
      end


      function m1l_finite_CDR(p,born)
c Returns the finite part of virtual contribution, according to the
c definitions given in (B.1) and (B.2). This function must include
c the factor as/(2*pi)
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      include 'coupl.inc'
      include 'q_es.inc'
      double precision p(0:3,nexternal-1),m1l_finite_CDR,born
      double precision CF,pi,aso2pi,shat,dot,xlgq2os
      parameter (CF=4d0/3d0)
      parameter (pi=3.1415926535897932385d0)
c
      aso2pi=g**2/(8*pi**2)
c This is relevant to e+e- --> qqbar
      shat=2d0*dot(p(0,1),p(0,2))
      xlgq2os=log(QES2/shat)
      m1l_finite_CDR=-aso2pi*CF*(xlgq2os**2+3*xlgq2os-pi**2+8.d0)*born
      return
      end


      function m1l_W_finite_CDR(p,born)
c Returns the finite part of virtual contribution, according to the
c definitions given in (B.1) and (B.2). This function must include
c the factor as/(2*pi)
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      include 'coupl.inc'
      include 'q_es.inc'
      double precision p(0:3,nexternal-1),m1l_W_finite_CDR,born
      double precision CF,pi,aso2pi,shat,dot,xlgq2os
      parameter (CF=4d0/3d0)
      parameter (pi=3.1415926535897932385d0)
c
      aso2pi=g**2/(8*pi**2)
      shat=2d0*dot(p(0,1),p(0,2))
      xlgq2os=log(QES2/shat)

c This is relevant to qqbar -> W 
      m1l_W_finite_CDR=aso2pi*CF*(-xlgq2os**2-3d0*xlgq2os+pi**2-8d0)
      m1l_W_finite_CDR=m1l_W_finite_CDR*born

c This is relevant to gg -> H
c$$$      m1l_W_finite_CDR=aso2pi*(-3d0*xlgq2os**2+11d0+3d0*pi**2)
c$$$      m1l_W_finite_CDR=m1l_W_finite_CDR*born

c This is relevant to bbbar -> H
c$$$      m1l_W_finite_CDR=aso2pi
c$$$     f     * (-4d0/3d0*xlgq2os**2
c$$$     f        -8d0/3d0+(16d0/3d0+8d0/3d0)*pi**2/6d0)
c$$$      m1l_W_finite_CDR=m1l_W_finite_CDR*born
      return
      end


      subroutine setfksfactor(iconfig)
      implicit none

      double precision CA,CF, PI
      parameter (CA=3d0,CF=4d0/3d0)
      parameter (pi=3.1415926535897932385d0)

      double precision c(0:1),gamma(0:1),gammap(0:1)
      common/fks_colors/c,gamma,gammap

      double precision beta0,ren_group_coeff
      common/cbeta0/beta0,ren_group_coeff

      logical softtest,colltest
      common/sctests/softtest,colltest

      integer config_fks,i,j,iconfig,fac1,fac2

      double precision fkssymmetryfactor,fkssymmetryfactorBorn,
     &     fkssymmetryfactorDeg
      integer ngluons,nquarks(-6:6),nphotons
      common/numberofparticles/fkssymmetryfactor,fkssymmetryfactorBorn,
     &                  fkssymmetryfactorDeg,ngluons,nquarks,nphotons

      include 'coupl.inc'
      include 'genps.inc'
      include 'nexternal.inc'
      include 'fks_powers.inc'
      include 'nFKSconfigs.inc'
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      include 'reweight0.inc'
      include 'run.inc'
      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS

      include "born_conf.inc"

      logical firsttime,firsttime_nFKSprocess(fks_configs)
      data firsttime,firsttime_nFKSprocess/.true.,fks_configs*.true./

      double precision xicut_used
      common /cxicut_used/xicut_used
      double precision delta_used
      common /cdelta_used/delta_used
      double precision xiScut_used,xiBSVcut_used
      common /cxiScut_used/xiScut_used,xiBSVcut_used
      logical rotategranny
      common/crotategranny/rotategranny
      double precision diagramsymmetryfactor
      common /dsymfactor/diagramsymmetryfactor

      integer           isum_hel
      logical                   multi_channel
      common/to_matrix/isum_hel, multi_channel

      character*1 integrate
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      integer fac_i,fac_j,i_fks_pdg,j_fks_pdg

      integer fac_i_FKS(fks_configs),fac_j_FKS(fks_configs)
     $     ,i_type_FKS(fks_configs),j_type_FKS(fks_configs)
     $     ,m_type_FKS(fks_configs),ngluons_FKS(fks_configs)
     $     ,nphotons_FKS(fks_configs)
      double precision
     $ ch_i_FKS(fks_configs),ch_j_FKS(fks_configs),ch_m_FKS(fks_configs)
      save fac_i_FKS,fac_j_FKS,i_type_FKS,j_type_FKS,m_type_FKS
     $     ,ngluons_FKS,ch_i_FKS,ch_j_FKS,ch_m_FKS,nphotons_FKS

      character*13 filename

      character*4 abrv
      common /to_abrv/ abrv

      logical nbody
      common/cnbody/nbody

      integer fold
      common /cfl/fold

c Particle types (=color) of i_fks, j_fks and fks_mother
      integer i_type,j_type,m_type
      double precision ch_i,ch_j,ch_m
      common/cparticle_types/i_type,j_type,m_type,ch_i,ch_j,ch_m
      double precision particle_charge(nexternal), particle_charge_born(nexternal-1)
      common /c_charges/particle_charge
      common /c_charges_born/particle_charge_born
      double precision zero
      parameter (zero=0d0)

c The value of rotategranny may be superseded later if phase space
c parametrization allows it
      rotategranny=.false.

      softtest=.false.
      colltest=.false.
      fold=0

      if (j_fks.gt.nincoming)then
         delta_used=deltaO
      else
         delta_used=deltaI
      endif
      
      xicut_used=xicut
      xiScut_used=xiScut
      if( nbody .or. (abrv.eq.'born' .or. abrv.eq.'grid' .or.
     &     abrv(1:2).eq.'vi') )then
        xiBSVcut_used=1.d0
      else
        xiBSVcut_used=xiBSVcut
      endif

      c(0)=CA
      c(1)=CF
      gamma(0)=( 11d0*CA-2d0*Nf )/6d0
      gamma(1)=CF*3d0/2d0
      gammap(0)=( 67d0/9d0 - 2d0*PI**2/3d0 )*CA - 23d0/18d0*Nf
      gammap(1)=( 13/2d0 - 2d0*PI**2/3d0 )*CF
            
c Beta_0 defined according to (MadFKS.C.5)
      beta0=gamma(0)/(2*pi)
c ren_group_coeff defined accordingly
      ren_group_coeff=ren_group_coeff_in/(2*pi)

      if (firsttime_nFKSprocess(nFKSprocess)) then
         firsttime_nFKSprocess(nFKSprocess)=.false.
c---------------------------------------------------------------------
c              Symmetry Factors
c---------------------------------------------------------------------
c fkssymmetryfactor:
c Calculate the FKS symmetry factors to be able to reduce the number
c of directories to (maximum) 4 (neglecting quark flavors):
c     1. i_fks=gluon, j_fks=gluon 
c     2. i_fks=gluon, j_fks=quark
c     3. i_fks=gluon, j_fks=anti-quark
c     4. i_fks=quark, j_fks=anti-quark (or vice versa).
c This sets the fkssymmetryfactor (in which the quark flavors are taken
c into account) for the subtracted reals.
c
c fkssymmetryfactorBorn:
c Note that in the Born's included here, the final state identical
c particle factor is set equal to the identical particle factor
c for the real contribution to be able to get the correct limits for the
c subtraction terms and the approximated real contributions.
c However when we want to calculate the Born contributions only, we
c have to correct for this difference. Since we only include the Born
c related to a soft limit (this uniquely defines the Born for a given real)
c the difference is always n!/(n-1)!=n, where n is the number of final state
c gluons in the real contribution.
c
c Furthermore, because we are not integrating all the directories, we also
c have to include a fkssymmetryfactor for the Born contributions. However,
c this factor is not the same as the factor defined above, because in this
c case i_fks is fixed to the extra gluon (which goes soft and defines the
c Born contribution) and should therefore not be taken into account when
c calculating the symmetry factor. Together with the factor n above this
c sets the fkssymmetryfactorBorn equal to the fkssymmetryfactor for the
c subtracted reals.
c
c We set fkssymmetryfactorBorn to zero when i_fks not a gluon
c
         i_fks_pdg=pdg_type(i_fks)
         j_fks_pdg=pdg_type(j_fks)
      
         fac_i_FKS(nFKSprocess)=0
         fac_j_FKS(nFKSprocess)=0
         do i=nincoming+1,nexternal
            if (i_fks_pdg.eq.pdg_type(i)) fac_i_FKS(nFKSprocess) =
     $           fac_i_FKS(nFKSprocess) + 1
            if (j_fks_pdg.eq.pdg_type(i)) fac_j_FKS(nFKSprocess) =
     $           fac_j_FKS(nFKSprocess) + 1
         enddo
c Overwrite if initial state singularity
         if(j_fks.le.nincoming) fac_j_FKS(nFKSprocess)=1

c i_fks and j_fks of the same type? -> subtract 1 to avoid double counting
         if (j_fks.gt.nincoming .and. i_fks_pdg.eq.j_fks_pdg)
     $        fac_j_FKS(nFKSprocess)=fac_j_FKS(nFKSprocess)-1

c THESE TESTS WORK ONLY FOR FINAL STATE SINGULARITIES
         if (j_fks.gt.nincoming) then
            if ( i_fks_pdg.eq.j_fks_pdg .and. i_fks_pdg.ne.21) then
               write (*,*) 'ERROR, if PDG type of i_fks and j_fks '//
     &              'are equal, they MUST be gluons',
     &              i_fks,j_fks,i_fks_pdg,j_fks_pdg
               stop
            elseif(abs(particle_type(i_fks)).eq.3) then
               if ( particle_type(i_fks).ne.-particle_type(j_fks) .or.
     &              pdg_type(i_fks).ne.-pdg_type(j_fks)) then
                  write (*,*) 'ERROR, if i_fks is a color triplet,'//
     &                 ' j_fks must be its anti-particle,'//
     &                 ' or an initial state gluon.',
     &                 i_fks,j_fks,particle_type(i_fks),
     &                 particle_type(j_fks),pdg_type(i_fks),pdg_type(j_fks)
                  stop
               endif
            elseif(i_fks_pdg.ne.21.and.i_fks_pdg.ne.22) then ! if not already above, it MUST be a gluon or photon
               write (*,*) 'ERROR, i_fks is not a g/gamma and falls not'//
     $              ' in other categories',i_fks,j_fks,i_fks_pdg
     $              ,j_fks_pdg
            endif
         endif

         ngluons_FKS(nFKSprocess)=0
         nphotons_FKS(nFKSprocess)=0
         do i=nincoming+1,nexternal
            if (pdg_type(i).eq.21) ngluons_FKS(nFKSprocess)
     $           =ngluons_FKS(nFKSprocess)+1
            if (pdg_type(i).eq.22) nphotons_FKS(nFKSprocess)
     $           =nphotons_FKS(nFKSprocess)+1
         enddo



c Set color types of i_fks, j_fks and fks_mother.
         i_type=particle_type(i_fks)
         j_type=particle_type(j_fks)
         ch_i=particle_charge(i_fks)
         ch_j=particle_charge(j_fks)
         call get_mother_col_charge(i_type,ch_i,j_type,ch_j,m_type,ch_m) 
         i_type_FKS(nFKSprocess)=i_type
         j_type_FKS(nFKSprocess)=j_type
         m_type_FKS(nFKSprocess)=m_type
         ch_i_FKS(nFKSprocess)=ch_i
         ch_j_FKS(nFKSprocess)=ch_j
         ch_m_FKS(nFKSprocess)=ch_m
      endif

      i_type=i_type_FKS(nFKSprocess)
      j_type=j_type_FKS(nFKSprocess)
      m_type=m_type_FKS(nFKSprocess)
      ch_i=ch_i_FKS(nFKSprocess)
      ch_j=ch_j_FKS(nFKSprocess)
      ch_m=ch_m_FKS(nFKSprocess)

c Set matrices used by MC counterterms
c      call set_mc_matrices

      fac_i=fac_i_FKS(nFKSprocess)
      fac_j=fac_j_FKS(nFKSprocess)
      ngluons=ngluons_FKS(nFKSprocess)
      nphotons=nphotons_FKS(nFKSprocess)
c Setup the FKS symmetry factors. 
      if (nbody.and.pdg_type(i_fks).eq.21) then
         fkssymmetryfactor=dble(ngluons)
         fkssymmetryfactorDeg=dble(ngluons)
         fkssymmetryfactorBorn=dble(ngluons)
      elseif (nbody.and.pdg_type(i_fks).eq.22) then
         fkssymmetryfactor=dble(nphotons)
         fkssymmetryfactorDeg=dble(nphotons)
         fkssymmetryfactorBorn=dble(nphotons)
      else
         fkssymmetryfactor=dble(fac_i*fac_j)
         fkssymmetryfactorDeg=dble(fac_i*fac_j)
         if (pdg_type(i_fks).eq.21) then
            fkssymmetryfactorBorn=dble(fac_i*fac_j)
         else
            fkssymmetryfactorBorn=0d0
         endif
         if (abrv.eq.'grid') then
            fkssymmetryfactorBorn=1d0
            fkssymmetryfactor=0d0
            fkssymmetryfactorDeg=0d0
         endif
      endif

      if (firsttime) then
c Check to see if this channel needs to be included in the multi-channeling
         diagramsymmetryfactor=0d0
         if (multi_channel) then
            open (unit=19,file="symfact.dat",status="old",err=14)
            do i=1,mapconfig(0)
               read (19,*,err=23) fac1,fac2
               if (i.eq.iconfig) then
                  if (mapconfig(iconfig).ne.fac1) then
                     write (*,*) 'inconsistency in symfact.dat',i
     $                    ,iconfig,mapconfig(iconfig),fac1
                     stop
                  endif
                  diagramsymmetryfactor=dble(fac2)
               endif
            enddo
            close(19)
         else                   ! no multi_channel
            diagramsymmetryfactor=1d0
         endif
 12      continue
         firsttime=.false.
      endif

      return

 99   continue
      write (*,*) '"integrate.fks" or "nbodyonly.fks" not found.'
      write (*,*) 'make and run "genint_fks" first.'
      stop
 23   continue
      write (*,*) '"symfact.dat" is not of the correct format'
      stop
 14   continue
      diagramsymmetryfactor=1d0
      goto 12
      end


      subroutine get_mother_col_charge(i_type, ch_i, j_type, ch_j, m_type, ch_m)
C viven the type (color representation) and charges of i and j, return
C the type and charges of the mother particle
      implicit none
      integer i_type, j_type, m_type
      double precision ch_i, ch_j, ch_m
      include 'nexternal.inc'
      integer i_fks,j_fks

      common/fks_indices/i_fks,j_fks

      if (abs(i_type).eq.abs(j_type) .and. 
     &    abs(ch_i).eq.abs(ch_j) .and. 
     &    abs(i_type).gt.1) then
        ! neutral color octet splitting
         m_type=8
         ch_m = 0d0
         if ( (j_fks.le.nincoming .and.
     &         abs(i_type).eq.3 .and. j_type.ne.i_type) .or.
     &        (j_fks.gt.nincoming .and.
     &         abs(i_type).eq.3 .and. j_type.ne.-i_type)) then
            write(*,*)'Flavour mismatch #1col in get_mother_col_charge',
     &                 i_fks,j_fks,i_type,j_type
            stop
         endif
      elseif (abs(i_type).eq.abs(j_type) .and. 
     &        dabs(ch_i).eq.dabs(ch_j).and.abs(i_type).eq.1) then
        ! neutral color singlet splitting
         m_type=1
         ch_m = 0d0
         if ( (j_fks.le.nincoming .and.
     &         dabs(ch_i).gt.0d0 .and. ch_j.ne.ch_i) .or.
     &        (j_fks.gt.nincoming .and.
     &         dabs(ch_i).gt.0d0 .and. ch_j.ne.-ch_i)) then
            write(*,*)'Flavour mismatch #1chg in get_mother_col_charge',
     &                 i_fks,j_fks,ch_i,ch_j
            stop
         endif
      elseif ((abs(i_type).eq.3 .and. j_type.eq.8) .or.
     &        (dabs(ch_i).gt.0d0 .and. ch_j.eq.0d0) ) then
         if(j_fks.le.nincoming)then
            m_type=-i_type
            if (m_type.eq.-1) m_type=1
            ch_m = -ch_i
         else
            write(*,*) 'Error in get_mother_col_charge: (i,j)=(q,g)'
            stop
         endif
      elseif ((i_type.eq.8 .and. abs(j_type).eq.3) .or.
     &         (ch_i.eq.0d0 .and. dabs(ch_j).gt.0d0) ) then
         m_type=j_type
         ch_m= ch_j
      else
         write(*,*)'Flavour mismatch #2 in get_mother_col_charge',
     &      i_type,j_type,m_type
         stop
      endif
      return
      end


      subroutine get_helicity(i_fks,j_fks)
      implicit none
      include "nexternal.inc"
      include "born_nhel.inc"
      include "madfks_mcatnlo.inc"
      integer NHEL(nexternal,max_bhel*2),IHEL
chel  include "helicities.inc"
      include 'nFKSconfigs.inc'
      double precision hel_fac
      integer get_hel,skip(fks_configs)
      common/cBorn/hel_fac,get_hel,skip
      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn
      integer hel_wgt,hel_wgt_born,hel_wgt_real
      integer nhelreal(nexternal,4),goodhelreal(4)
      integer nhelrealall(nexternal,max_bhel*2)
      common /c_nhelreal/ nhelreal,nhelrealall,goodhelreal,hel_wgt_real
      integer nhelborn(nexternal-1,2),goodhelborn(2)
      integer nhelbornall(nexternal-1,max_bhel)
      common /c_nhelborn/ nhelborn,nhelbornall,goodhelborn,hel_wgt_born

      integer           isum_hel
      logical                   multi_channel
      common/to_matrix/isum_hel, multi_channel

      integer i,nexthel,j,i_fks,j_fks,ngood,k
      data nexthel /0/
      data ngood /0/
      logical done,firsttime,all_set,chckr
      data firsttime/.true./
      integer goodhelr(0:4,max_bhel/2),goodhelb(0:2,max_bhel/2)
      save goodhelr,goodhelb,all_set,chckr
      double precision rnd,ran2
      external ran2

      character*4 abrv
      common /to_abrv/ abrv
      logical Hevents
      common/SHevents/Hevents
      logical usexinteg,mint
      common/cusexinteg/usexinteg,mint

c Do not change these two lines, because ./bin/compile_madfks.sh might
c need to change them automatically
      logical HelSum
      parameter (HelSum=.true.)

c************
c goodhelr=2, real emission matrix element not yet calculated
c             for this helicity
c goodhelr=1, real emission matrix element calculated and non-zero
c goodhelr=0, real emission matrix element calculated and zero,
c             so can be skipped next time.
c************
      if (HelSum) return

      if (isum_hel.ne.0) then ! MC over helicities
c First, set the goodhelr and goodhelb to their starting values
      if (firsttime) then
         if ((mint .and. (.not.Hevents) .and. (abrv(1:2).eq.'vi' .or.
     &        abrv.eq.'born' .or. abrv.eq.'grid' .or.
     &        (.not.UseSudakov))) .or. (.not.mint .and. (abrv.eq.'born'
     &        .or. abrv.eq.'grid' .or. abrv(1:2).eq.'vi'))) then
c           if computing only the Born diagrams, should not
c           consider real emission helicities            
            chckr=.false.
         else
            chckr=.true.
         endif
         do i=1,fks_configs
            skip(i)=1
         enddo
c read from file if possible
         open(unit=65,file='goodhel.dat',status='old',err=532)
         all_set=.true.
         do j=0,4
            read (65,*,err=532) (goodhelr(j,i),i=1,max_bhel/2)
         enddo
         do j=0,2
            read (65,*,err=532) (goodhelb(j,i),i=1,max_bhel/2)
         enddo
         read(65,*,err=532) hel_wgt
         hel_wgt_born=hel_wgt
         hel_wgt_real=hel_wgt
         do i=1,max_bhel/2
            if ((chckr .and.
     &           (goodhelb(0,i).eq.2 .or. goodhelr(0,i).eq.2)) .or.
     &           (.not.chckr.and.goodhelb(0,i).eq.2)) all_set=.false.
         enddo
         close(65)
         goto 533
c if file does not exist or has wrong format, set all to 2
 532     close(65)
         write (*,*) 'Good helicities not found in file'
         all_set=.false.
         do j=0,4
            do i=1,max_bhel/2
               goodhelr(j,i)=2
            enddo
         enddo
         do j=0,2
            do i=1,max_bhel/2
               goodhelb(j,i)=2
            enddo
         enddo
         hel_wgt=max_bhel/2
         hel_wgt_born=hel_wgt
         hel_wgt_real=hel_wgt
 533     continue
         firsttime=.false.
         goto 534 ! no previous event, so skip to the next helicity
      endif

c From previous event, check if there is an update
      if (.not.all_set) then
c real emission
         if(goodhelr(0,ngood).eq.2) then
            if ( goodhelreal(1).eq.0 .and.
     &           goodhelreal(2).eq.0 .and.
     &           goodhelreal(3).eq.0 .and.
     &           goodhelreal(4).eq.0 ) then
               do j=0,4
                  goodhelr(j,ngood)=0
               enddo
            elseif( goodhelreal(1).le.1 .and.
     &              goodhelreal(2).le.1 .and.
     &              goodhelreal(3).le.1 .and.
     &              goodhelreal(4).le.1 ) then
               goodhelr(0,ngood)=1
               do j=1,4
                  goodhelr(j,ngood)=goodhelreal(j)
               enddo
            elseif (.not.(goodhelreal(1).eq.2 .and.
     &                    goodhelreal(2).eq.2 .and.
     &                    goodhelreal(2).eq.2 .and.
     &                    goodhelreal(2).eq.2) ) then
               write (*,*) 'Error #2 in get_helicities',
     &              ngood,(goodhelr(j,ngood),j=0,4)
               stop
            endif
         endif
c Born and counter events
         if(goodhelb(0,ngood).eq.2) then
            if ( goodhelborn(1).eq.0 .and.
     &           goodhelborn(2).eq.0 ) then
               do j=0,2
                  goodhelb(j,ngood)=0
               enddo
            elseif( goodhelborn(1).le.1 .and.
     &              goodhelborn(2).le.1 ) then
               goodhelb(0,ngood)=1
               do j=1,2
                  goodhelb(j,ngood)=goodhelborn(j)
               enddo
            elseif (.not.(goodhelborn(1).eq.2 .and.
     &                    goodhelborn(2).eq.2) ) then
               write (*,*) 'Error #3 in get_helicities',
     &              nexthel,(goodhelb(j,ngood),j=0,2)
               stop
            endif
         endif

c Calculate new hel_wgt
         hel_wgt=0
         do i=1,max_bhel/2
            if((chckr .and.
     &           (goodhelb(0,i).ge.1.or.goodhelr(0,i).ge.1)) .or.
     &           (.not.chckr .and. goodhelb(0,i).ge.1)) then
               hel_wgt=hel_wgt+1
            endif
         enddo
         hel_wgt_born=hel_wgt
         hel_wgt_real=hel_wgt

c check if all have been set, if so -> write to file
         all_set=.true.
         do i=1,max_bhel/2
            if ((chckr .and.
     &           (goodhelb(0,i).eq.2 .or. goodhelr(0,i).eq.2)) .or.
     &           (.not.chckr.and.goodhelb(0,i).eq.2)) all_set=.false.
         enddo
         if (all_set) then
            write (*,*) 'All good helicities have been found.',hel_wgt
            open(unit=65,file='goodhel.dat',status='unknown')
            do j=0,4
               write (65,*) (goodhelr(j,i),i=1,max_bhel/2)
            enddo
            do j=0,2
               write (65,*) (goodhelb(j,i),i=1,max_bhel/2)
            enddo
            write(65,*) hel_wgt
            close(65)
         endif
      else
         do i=1,4
            if (goodhelr(i,ngood).ne.goodhelreal(i)) then
               write (*,*)'Error #4 in get_helicities',i,ngood
               stop
            endif
         enddo
         do i=1,2
            if (goodhelb(i,ngood).ne.goodhelborn(i)) then
               write (*,*)'Error #5 in get_helicities',i,ngood
               stop
            endif
         enddo
      endif

c Get the next helicity
 534  continue
      done=.false.
      do while (.not.done)
         if (nexthel.eq.max_bhel*2) nexthel=0
         nexthel=nexthel+1
         if(nhel(i_fks,nexthel).eq.1.and.nhel(j_fks,nexthel).eq.1) then
            if (ngood.eq.max_bhel/2) ngood=0
            ngood=ngood+1
            if((chckr .and.
     &           (goodhelr(0,ngood).ge.1.or.goodhelb(0,ngood).ge.1)).or.
     &           (.not.chckr .and. goodhelb(0,ngood).ge.1)) then
c Using random number to see if we have to go to the next.
c Probably this is an overkill, but have to make sure that there is
c no bias considering the *semi*-random numbers from VEGAS.
               rnd=ran2()
               if (rnd.le.1d0/dble(hel_wgt)) then
                  done=.true.
               endif
            endif
         endif
      enddo

      do i=1,nexternal
         if (i.eq.i_fks) then
            nhelreal(i,1)=1
            nhelreal(i,2)=1
            nhelreal(i,3)=-1
            nhelreal(i,4)=-1
         elseif (i.eq.j_fks) then
            nhelreal(i,1)=1
            nhelreal(i,2)=-1
            nhelreal(i,3)=1
            nhelreal(i,4)=-1
         else
            nhelreal(i,1)=nhel(i,nexthel)
            nhelreal(i,2)=nhel(i,nexthel)
            nhelreal(i,3)=nhel(i,nexthel)
            nhelreal(i,4)=nhel(i,nexthel)
         endif
      enddo
      do j=1,4
         goodhelreal(j)=goodhelr(j,ngood)
      enddo

      do i=1,nexternal-1
         if (i.eq.min(i_fks,j_fks)) then
            nhelborn(i,1)=1
            nhelborn(i,2)=-1
         elseif(i.lt.max(i_fks,j_fks)) then
            nhelborn(i,1)=nhel(i,nexthel)
            nhelborn(i,2)=nhel(i,nexthel)
         else
            nhelborn(i,1)=nhel(i+1,nexthel)
            nhelborn(i,2)=nhel(i+1,nexthel)
         endif
      enddo
      do j=1,2
         goodhelborn(j)=goodhelb(j,ngood)
      enddo

      else !isum_hel is zero, sum explicitly over helicities

      do i=1,nexternal
         do j=1,max_bhel*2
            nhelrealall(i,j)=nhel(i,j)
         enddo
      enddo
      do i=1,nexternal-1
         k=0
         do j=1,max_bhel*2
            if (nhel(i_fks,j).eq.-1) then
               k=k+1
               if (i.lt.i_fks) then
                  nhelbornall(i,k)=nhel(i,j)                  
               elseif(i.ge.i_fks) then
                  nhelbornall(i,k)=nhel(i+1,j)
               endif
            endif
         enddo
      enddo

      endif
      return
      end

      function get_ptrel(pp,i_fks,j_fks)
      implicit none
      include 'nexternal.inc'
      double precision get_ptrel,pp(0:3,nexternal)
      integer i_fks,j_fks
      double precision tmp,psum(3)
      integer i
c
      if(j_fks.le.2)then
        tmp=sqrt(pp(1,i_fks)**2+pp(2,i_fks)**2)
      else
        do i=1,3
          psum(i)=pp(i,i_fks)+pp(i,j_fks)
        enddo
        tmp=( pp(2,i_fks)*psum(1)-pp(1,i_fks)*psum(2) )**2+
     #      ( pp(3,i_fks)*psum(1)-pp(1,i_fks)*psum(3) )**2+
     #      ( pp(3,i_fks)*psum(2)-pp(2,i_fks)*psum(3) )**2
        if(tmp.ne.0.d0)tmp=sqrt( tmp/
     #       (psum(1)**2+psum(2)**2+psum(3)**2) )
      endif
      get_ptrel=tmp
      return
      end



      FUNCTION FK88RANDOM(SEED)
*     -----------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
*
      IMPLICIT INTEGER(A-Z)
      REAL*8 MINV,FK88RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      FK88RANDOM = SEED*MINV
      END
