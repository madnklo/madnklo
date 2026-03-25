

      SUBROUTINE SUB_M2_HC_RV_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr,ret)
c     collinear limit C_(ia,ib) * Wcollinear for RV
c     this is meant to represent the full collinear
c     for sector (ia,ib)
      use sectors2_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      double precision ret(-2:0)
      double precision M2_C_gq(-2:0),M2_C_gq_0,M2_SC_gq(-2:0)
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,wgt,wgts(1),wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,VLO(-2:0),EIK0,EIK1(-2:0)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,y,xinit,logab,damp
      double precision wa,wb,wr,mb2,mr2
      double precision ANS(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision ddilog
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_HC_RV_gq)s_den
      common/%(proc_prefix_HC_RV_gq)s_iden/%(proc_prefix_HC_RV_gq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      ret=0d0
      M2_C_gq_0=0d0
      M2_C_gq=0d0
      M2_SC_gq=0d0
      ierr=0
      damp=0d0
c
c     checks
      if(.not.(abs(leg_pdgs(ib)).le.6.and.leg_pdgs(ia).eq.21))then
         write(*,*)'Wrong pdgs in M2_HC_RV_gq',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)))then
         write(*,*)'Wrong indices in M2_HC_RV_gq',ia,ib,isec,jsec
         stop
      endif
c
c     possible cuts
      if(docut(xpb,nexternal-1,underlying_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      mb2=pmass(ib)**2
      mr2=pmass(ir)**2
c     TODO: CHECK x
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
      logab=log(sab/scale**2)
c
      EIK0     =  SBR/(SAB*SAR) - MB2/SAB**2 - MR2/SAR**2
      EIK1(-2) =  CA*EIK0
      EIK1(-1) = -CA*EIK0*log(sab*sar/sbr/scale**2)
      EIK1( 0) =  CA*EIK0/2d0*(log(sab*sar/sbr/scale**2)**2-5d0*zeta2)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HC_RV_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
c      call %(proc_prefix_HC_RV_gq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = 0d0 !ANS(0)
c      call %(proc_prefix_HC_RV_gq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      VLO = 0d0 !ANS(0)
c
c     In the following equation the x variable is related to the quark energy
      M2_C_gq_0=BLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))
c
      M2_C_gq(-2:0) = M2_C_gq(-2:0) + VLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))
      M2_C_gq(-2) = M2_C_gq(-2) + alphas/2d0/pi*M2_C_gq_0*(-CA)
      M2_C_gq(-1) = M2_C_gq(-1) + alphas/2d0/pi*M2_C_gq_0*(CA*logab+CF*log(x*(1d0-x))-beta0/2d0)
      M2_C_gq( 0) = M2_C_gq( 0) + alphas/2d0/pi*(M2_C_gq_0*(CA*(7d0*zeta2-logab**2)/2d0+CF*(-logab*log(x*(1d0-x))+ddilog(-(1d0-x)/x)+ddilog(-x/(1d0-x))))+BLO*CF*(CA-CF))
c
      if(ia.eq.isec) then
         M2_SC_gq(-2:0) = M2_SC_gq(-2:0)+2d0*CF*(EIK0*VLO(-2:0)-alphas/2d0/pi*EIK1(-2:0)*BLO)
         M2_SC_gq(-1)   = M2_SC_gq(-1)-2d0*CF*alphas/2d0/pi*beta0/2d0*EIK0*BLO
      else
         continue
      endif
c     compute collinear limit of sector function
      call get_wc_nlo(isec,jsec,iref)
      M2_C_gq =  M2_C_gq*wc_nlo
c     account for different damping factors according to recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      ret = M2_C_gq-M2_SC_gq
c     include prefactors
      ret = ret *dble(%(proc_prefix_HC_RV_gq)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor*damp*pref/sab*xj*extra
c
c     plot
      wgtpl=-ret(0)*wgt/nit*wgt_chan
      wgts=wgtpl
c     if(doplot)call histo_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgts)
c
c     sanity check
      if(abs(ret(0)).ge.huge(1d0).or.isnan(ret(0)))then
         write(77,*)'Exception caught in M2_HC_RV_gq',ret(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      SUBROUTINE DELTA_HC_RV_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr,ret)
c     collinear limit C_(ia,ib) * Wcollinear for RV
c     this is meant to represent the full collinear
c     for sector (ia,ib)
      use sectors2_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      double precision ret(-2:0)
      double precision DELTA_C_gq_0(-2:0),DELTA_SC_gq(-2:0)
      double precision M2tmp(-2:0),M2tmp_SC(-2:0)
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,wgt,wgts,wgtpl,wgt_chan,xj,xjcs,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO,VLO(-2:0),KKVLO(-2:0),EIK0,EIK1(-2:0)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision xpb_abr(0:3,nexternal-1),xpb_bra(0:3,nexternal-1),xpb_arb(0:3,nexternal-1)
      double precision xsb_abr(nexternal-1,nexternal-1),xsb_bra(nexternal-1,nexternal-1),xsb_arb(nexternal-1,nexternal-1)
      double precision sab,sar,sbr,sabr,sal,sbl,slm,x,y,xinit,logab,damp
      double precision sb_arb_bl,sb_arb_br,sb_bl,sb_bra_al,sb_bra_ar,sb_lm
      double precision wa,wb,wr,mb2,mr2
      double precision ANS(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision %(proc_prefix_HC_RV_gq)s_GET_CCBLO
      double precision %(proc_prefix_HC_RV_gq)s_GET_KKBLO
      double precision %(proc_prefix_HC_RV_gq)s_GET_KKVLO
      double precision gamma_l,phi_l
      double precision ccblo_lm,ccblo_parent_l
      double precision blo_arb,ccblo_l_arb,ccblo_parent_l_bra,ccblo_parent_l_arb
      double precision ddilog
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      integer %(proc_prefix_real)s_den
     common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_HC_RV_gq)s_den
      common/%(proc_prefix_HC_RV_gq)s_iden/%(proc_prefix_HC_RV_gq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      integer mapped_labels_bra(nexternal),mapped_labels_arb(nexternal)
      common/c_mapped_labels/mapped_labels
      integer l,m,lb,mb
      double precision PtimesB
      logical isnloqcdparton(nexternal)
      double precision pmass(nexternal)
      include 'pmass.inc'
c     initialise
      ret=0d0
      DELTA_C_gq_0=0d0
      DELTA_SC_gq =0d0
      M2tmp=0d0
      M2tmp_SC = 0d0
      PtimesB = 0d0
      ierr=0
      damp=0d0
      BLO_arb = 0d0
c
c     checks
      if(.not.(abs(leg_pdgs(ib)).le.6.and.leg_pdgs(ia).eq.21))then
         write(*,*)'Wrong pdgs in DELTA_HC_RV_gq',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)))then
         write(*,*)'Wrong indices in DELTA_HC_RV_gq',ia,ib,isec,jsec
         stop
      endif
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      sabr = sar+sbr
      mb2=pmass(ib)**2
      mr2=pmass(ir)**2
c     Energy fraction of the quark
      x=sbr/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c     Evaluate momenta with different mappings
c     at fixed real kinematic
c     We need xbp^{(abr)}, xbp^{(bra)},xbp^{(abrprime)},
      call get_unp_mapped_labels(nexternal,ib,ir,mapped_labels_bra)
      call phase_space_CS_inv(ib,ir,ia,xp,xpb_bra,nexternal,leg_PDGs,xjCS,mapped_labels_bra)
      if(xjCS.eq.0d0)goto 999
      call invariants_from_p(xpb_bra,nexternal-1,xsb_bra,ierr)
      if(ierr.eq.1)goto 999
      call get_unp_mapped_labels(nexternal,ia,ir,mapped_labels_arb)
      call phase_space_CS_inv(ia,ir,ib,xp,xpb_arb,nexternal,leg_PDGs,xjCS,mapped_labels_arb)
      if(xjCS.eq.0d0)goto 999
      call invariants_from_p(xpb_arb,nexternal-1,xsb_arb,ierr)
      if(ierr.eq.1)goto 999
c     double sum for Eq.(5.20) (c,d) ---> (l,m)
c     Barred momenta with different mapping labels are needed
c     No double-pole
      M2tmp(-2) = 0d0
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in DELTA_HC_RV_gq',sab,sar+sbr,x
         goto 999
      endif
c     compute collinear limit of sector function
      call get_wc_nlo(isec,jsec,iref)
c     (abr) mapped Born
      ANS = 0d0
      call %(proc_prefix_HC_RV_gq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c     Eikonal E^{(i)}_{jr}
      eik0 = sbr/sab/sar - mb2/sab**2-mr2/sar**2
c     P_{gq}^{\mu\nu}B_{\mu\nu}
      PtimesB = CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*BLO
c
      do l=1,nexternal
         if(.not.isNLOQCDparton(l))cycle
         if(l.eq.ia) cycle
         if(l.eq.ib) cycle
         if(docut(xpb,nexternal-1,underlying_leg_pdgs,0)) goto 778
c     NNLO invariants
         sal = xs(ia,l)
         sbl = xs(ib,l)
c     NLO invariants (abr)
         lb = mapped_labels(l)
         sb_bl = xsb(mapped_labels(ib),lb)

C     Born matrix elements mapped (abr)
c     Colour-correlated Born_{[ab]l}
         CCBLO_parent_l=%(proc_prefix_HC_RV_gq)s_GET_CCBLO(mapped_labels(ib),lb)

         if(abs(leg_pdgs(l)).le.6) then
            gamma_l = gamma_q
            phi_l = phi_q
         elseif(leg_pdgs(l).eq.21) then
            gamma_l = gamma_g
            phi_l = phi_g
         else
            write(*,*) 'Error in evaluating gamma_c and phi_c'
            write(*,*)     'in delta_HC_RV_gq'
            write(*,*) 'c, leg_pdgs(c) = ', l, leg_pdgs(l)
            write(*,*) 'Exit...'
            stop
         endif
c
         M2tmp(-1) = M2tmp(-1) + gamma_l*PtimesB
         M2tmp(0) = M2tmp(0) + phi_l*PtimesB
         M2tmp(-1) = M2tmp(-1) + CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_parent_l*2d0*(dlog(sb_bl/sabr))
         M2tmp(0) = M2tmp(0) - CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_parent_l*2d0*1d0/2d0*dlog(sb_bl/sabr)**2
c     Soft-collinear term
         if(ia.eq.isec) then
            M2tmp_SC(-1) = M2tmp_SC(-1) + gamma_l*BLO
            M2tmp_SC(-1) = M2tmp_SC(-1) + 2d0*dlog(sb_bl/sbr)*CCBLO_parent_l
            M2tmp_SC(0) = M2tmp_SC(0) + phi_l*BLO
            M2tmp_SC(0) = M2tmp_SC(0) + 2d0*1d0/2d0*(dlog(sbr/sabr)**2-dlog(sbl/sabr)**2)*CCBLO_parent_l
         endif


         do m=1,nexternal
            if(.not.isNLOQCDparton(m))cycle
            if(m.eq.ia) cycle
            if(m.eq.ib) cycle
            if(m.eq.l) cycle
c     (ijr) ---> (ia ib ir)
c     Barred invariants (abr)
            mb = mapped_labels(m)
            sb_lm = xsb(lb,mb)
c    NNLO invariants
            slm = xs(l,m)
            sal = xs(ia,l)
            sar = xs(ia,ir)
            sbl = xs(ib,l)
            sbr = xs(ib,ir)
C     Define Born matrix elements for each pair (l,m)
C     Born matrix elements mapped (abr)
c     Colour-correlated Born_{lm}
            CCBLO_lm=%(proc_prefix_HC_RV_gq)s_GET_CCBLO(lb,mb)
c
            M2tmp(-1) = M2tmp(-1) - CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_lm*dlog(slm/sb_lm)
            M2tmp(0) = M2tmp(0) + CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_lm*1d0/2d0*dlog(slm/sb_lm)**2
c
c           if(m.gt.2.and.l.gt.2)then
c              y=sil/(sil+sim+slm)
c              z=sim/(sim+slm)
c              damp=((1d0-y)*(1d0-z))**alpha
c           elseif(m.gt.2.and.l.le.2)then
c              z=sim/(sim+slm)
c              x=1d0 - sil/(sim+slm)
c              damp=((1d0-z)*x)**alpha
c           elseif(m.le.2.and.l.le.2)then
c              x=1d0 - (sil+sim)/slm
c              damp=x**alpha
c     endif

c     Soft-Collinear
            if(ia.eq.isec) then
               M2tmp_SC(-1) = M2tmp_SC(-1)-CCBLO_lm*dlog(slm/sb_lm)
               M2tmp_SC(0) = M2tmp_SC(0)+CCBLO_lm*1d0/2d0*dlog(slm/sb_lm)**2
            endif
            ret(-2:0)=ret(-2:0)+alphas/2d0/pi/sab*M2tmp(-2:0)*wc_NLO-alphas/2d0/pi*2d0*CF*eik0*M2tmp_SC(-2:0)     
            ret = ret *dble(%(proc_prefix_HC_RV_gq)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor*damp*pref*xj*extra
c
c     plot
            wgtpl=-ret(0)*wgt/nit*wgt_chan
            wgts=wgtpl
c     if(doplot)call analysis_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgts)
c
 778        if(docut(xpb_bra,nexternal-1,underlying_leg_pdgs,0)) goto 779
c
c     Barred invariants (bra)
            sb_bra_al = xsb_bra(mapped_labels_bra(ia),mapped_labels_bra(l))
            sb_bra_ar = xsb_bra(mapped_labels_bra(ia),mapped_labels_bra(ir))
c     Mapped (bra) Born matrix element
            ANS = 0d0
            call %(proc_prefix_HC_RV_gq)s_ME_ACCESSOR_HOOK(xpb_bra,hel,alphas,ANS)
c     The mother particle for the splitting [ab] ---> a b is ib
            CCBLO_parent_l_bra=%(proc_prefix_HC_RV_gq)s_GET_CCBLO(mapped_labels(ib),lb) 
            M2tmp(-1) = M2tmp(-1) + CA/CF*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_parent_l_bra*dlog(sar/sal)
c
            M2tmp(0) = M2tmp(0) + CA/CF*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_parent_l_bra*1d0/2d0*(dlog(sb_bra_al/sb_bra_ar)**2-dlog(sar*sb_bra_al/sb_bra_ar/sal)**2)
c
            if(ia.eq.isec) then
               M2tmp_SC(-1) = M2tmp_SC(-1) + CA/CF*CCBLO_parent_l_bra*dlog(sar/sal)
               M2tmp_SC(0) = M2tmp_SC(0) + CA/CF*CCBLO_parent_l_bra*1d0/2d0*(dlog(sb_bra_al/sb_bra_ar)**2-dlog(sar*sb_bra_al/sb_bra_ar/sal)**2)
            endif
c     a <---> b
 779        if(docut(xpb_arb,nexternal-1,underlying_leg_pdgs,0)) cycle
c     Barred invariants (arb)
            sb_arb_bl = xsb_arb(mapped_labels_arb(ib),mapped_labels_arb(l))
            sb_arb_br = xsb_arb(mapped_labels_arb(ib),mapped_labels_arb(ir))
c     Mapped (arb) Born matrix element
            ANS = 0d0
            call %(proc_prefix_HC_RV_gq)s_ME_ACCESSOR_HOOK(xpb_arb,hel,alphas,ANS)
            BLO_arb = ANS(0)
c     The mother particle for the splitting [ab] ---> a b is ib
            CCBLO_parent_l_arb=%(proc_prefix_HC_RV_gq)s_GET_CCBLO(mapped_labels(ib),lb) 
c
            M2tmp(-1) = M2tmp(-1) + CA/CF*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_parent_l_arb*dlog(sbr/sbl)
c
            M2tmp(0) = M2tmp(0) + CA/CF*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*CCBLO_parent_l_arb*1d0/2d0*(dlog(sb_arb_bl/sb_arb_br)**2-dlog(sbr*sb_arb_bl/sb_arb_br/sbl)**2)
c

            if(ia.eq.isec) then      ! CHECK
               M2tmp_SC(-1) = M2tmp_SC(-1)+(2d0*CF-CA)/CF*CCBLO_parent_l_arb*dlog(sbr/sbl)
               M2tmp_SC(0) = M2tmp_SC(0)+(2d0*CF-CA)/CF*CCBLO_parent_l_arb*(dlog(sb_arb_bl/sb_arb_br)**2-dlog(sbr*sb_arb_bl/sb_arb_br/sal)**2)
            endif

            ret(-2:0)=ret(-2:0)+alphas/2d0/pi/sab*M2tmp(-2:0)*wc_NLO-alphas/2d0/pi*2d0*CF*eik0*M2tmp_SC(-2:0)     
            ret = ret *dble(%(proc_prefix_HC_RV_gq)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor*damp*pref*xj*extra
c     plot
            wgtpl=-ret(0)*wgt/nit*wgt_chan
            wgts=wgtpl
c     if(doplot)call histo_fill(xpb_arb,xsb_arb,nexternal-1,underlying_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpb_arb,xsb_arb,nexternal-1,underlying_leg_pdgs,wgts)
c
         enddo
      enddo



c     Term with rprime missing
      if(ia.eq.isec) then
         if(docut(xpb,nexternal-1,underlying_leg_pdgs,0)) goto 998
         M2tmp_SC(-1) = M2tmp_SC(-1) + alphas/2d0/pi*gamma_q*BLO
         M2tmp_SC(0) = M2tmp_SC(0) + alphas/2d0/pi*(phi_q-gamma_q*dlog(sCM/scale**2))*BLO
 998     if(docut(xpb_arb,nexternal-1,underlying_leg_pdgs,0)) return
        M2tmp_SC(-1) = M2tmp_SC(-1) - alphas/2d0/pi*gamma_q*BLO_arb
        M2tmp_SC(0) = M2tmp_SC(0) - alphas/2d0/pi*(phi_q-gamma_q*dlog(sCM/scale**2))*BLO_arb
        ret(-2:0)=ret(-2:0)-alphas/2d0/pi*2d0*CF*eik0*M2tmp_SC(-2:0)     
        ret = ret *dble(%(proc_prefix_HC_RV_gq)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor*damp*pref*xj*extra
c     plot
         wgtpl=-ret(0)*wgt/nit*wgt_chan
         wgts=wgtpl
         if(doplot)call analysis_fill(xpb_arb,xsb_arb,nexternal-1,underlying_leg_pdgs,wgts)
      endif
c
c
c     sanity check
      if(abs(ret(0)).ge.huge(1d0).or.isnan(ret(0)))then
         write(77,*)'Exception caught in DELTA_HC_RV_gq',ret(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
