

      function M2_C_RV_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
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
      double precision M2_C_RV_gq(-2:0)
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,M2tmp(-2:0),M2tmp0,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,VLO(-2:0)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,y,xinit,logab,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision alphaz
      parameter(alphaz=1d0)
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
      integer %(proc_prefix_C_RV_gq)s_den
      common/%(proc_prefix_C_RV_gq)s_iden/%(proc_prefix_C_RV_gq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      integer mapped_flavours(nexternal-1),mapped_indices_shuff(nexternal-1)
      common/c_mapped_quantities_c/mapped_labels,mapped_flavours,mapped_indices_shuff
      double precision xpb_to_ME(0:3,nexternal-1)
c
c     initialise
      M2_C_RV_gq=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpb_to_ME=0d0
c
c     checks
      if(.not.(abs(leg_pdgs(ib)).le.6.and.leg_pdgs(ia).eq.21))then
         write(*,*)'Wrong pdgs in M2_C_RV_gq',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)))then
         write(*,*)'Wrong indices in M2_C_RV_gq',ia,ib,isec,jsec
         stop
      endif
c
c     possible cuts
      IF(DOCUT(XPB,NEXTERNAL-1,MAPPED_FLAVOURS,0))RETURN
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
      logab=log(sab/scale**2)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_C_RV_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      XPB_TO_ME(0:3,MAPPED_INDICES_SHUFF(:))=XPB(0:3,:)
      call %(proc_prefix_C_RV_gg)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
      BLO = ANS(0)
      call %(proc_prefix_C_RV_gg)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
      VLO = ANS(0)
c
c     In the following equation the x variable is related to the quark energy
      M2TMP0=BLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))
c
      M2TMP(-2:0) = M2TMP(-2:0) + VLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))
      M2TMP(-2) = M2TMP(-2) + alphas/2d0/pi*M2TMP0*(-CA)
      M2TMP(-1) = M2TMP(-1) + alphas/2d0/pi*M2TMP0*(CA*logab+CF*log(x*(1d0-x))-beta0/2d0)
      M2TMP( 0) = M2TMP( 0) + alphas/2d0/pi*(M2TMP0*(CA*(7d0*zeta2-logab**2)/2d0+CF*(-logab*log(x*(1d0-x))+ddilog(-(1d0-x)/x)+ddilog(-x/(1d0-x))))+BLO*CF*(CA-CF))
c
c     Including correct multiplicity factor
      M2tmp(-2:0) = M2tmp(-2:0)*dble(%(proc_prefix_C_RV_gq)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp(-2:0)=M2tmp(-2:0)*damp

c     compute collinear limit of sector function
      call get_wc_nlo(xs,isec,jsec,iref,alphaz,nexternal)

      M2_C_RV_gq(-2:0)=M2tmp(-2:0)*pref/sab*xj*extra*wc_nlo
c     apply flavour factor
      M2_C_RV_gq(-2:0)=M2_C_RV_gq(-2:0)*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_C_RV_gq(0)*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,mapped_flavours,wgtpl)
c
c     sanity check
      if(abs(M2_C_RV_gq(0)).ge.huge(1d0).or.isnan(M2_C_RV_gq(0)))then
         write(77,*)'Exception caught in M2_C_RV_gq',M2_C_RV_gq(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      function M2_SC_RV_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     soft-collinear limit S_(ia)C_(ia,ib) for RV
c     this is meant to represent the soft-collinear
c     for sector (ia,ib)
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision M2_SC_RV_gq(-2:0)
      double precision pref,M2tmp(-2:0),wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,VLO(-2:0),EIK0,EIK1(-2:0)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision mb2,mr2,sab,sar,sbr,x,damp,y,xinit
      double precision ANS(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
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
      integer %(proc_prefix_SC_RV_gq)s_den
      common/%(proc_prefix_SC_RV_gq)s_iden/%(proc_prefix_SC_RV_gq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      integer mapped_flavours(nexternal-1),mapped_indices_shuff(nexternal-1)
      common/c_mapped_quantities_c/mapped_labels,mapped_flavours,mapped_indices_shuff
      double precision xpb_to_ME(0:3,nexternal-1)
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      M2_SC_RV_gq=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpb_to_ME=0d0
c
c     checks
      if(.not.(abs(leg_pdgs(ib)).le.6.and.leg_pdgs(ia).eq.21))then
         write(*,*)'Wrong pdgs in M2_SC_RV_gq',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)))then
         write(*,*)'Wrong indices in M2_SC_RV_gq',ia,ib,isec,jsec
         stop
      endif
c
c     possible cuts
      IF(DOCUT(XPB,NEXTERNAL-1,MAPPED_FLAVOURS,0))RETURN
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
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_SC_RV_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born and Virtual
      XPB_TO_ME(0:3,MAPPED_INDICES_SHUFF(:))=XPB(0:3,:)
      call %(proc_prefix_SC_RV_gq)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
      BLO = ANS(0)
      call %(proc_prefix_SC_RV_gq)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
      VLO = ANS(0)

      EIK0     =  SBR/(SAB*SAR) - MB2/SAB**2 - MR2/SAR**2
      EIK1(-2) =  CA*EIK0
      EIK1(-1) = -CA*EIK0*log(sab*sar/sbr/scale**2)
      EIK1( 0) =  CA*EIK0/2d0*(log(sab*sar/sbr/scale**2)**2-5d0*zeta2)
c
      M2TMP(-2:0) = M2TMP(-2:0)+2d0*CF*(EIK0*VLO(-2:0)-alphas/2d0/pi*EIK1(-2:0)*BLO)
      M2TMP(-1)   = M2TMP(-1)-2d0*CF*alphas/2d0/pi*beta0/2d0*EIK0*BLO
c     Including correct multiplicity factor
      M2tmp(-2:0) = M2tmp(-2:0)*dble(%(proc_prefix_SC_RV_gq)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp(-2:0)=M2tmp(-2:0)*damp

      M2_SC_RV_gq(-2:0)=M2tmp(-2:0)*pref*xj*extra
c     apply flavour factor
      M2_SC_RV_gq(-2:0)=M2_SC_RV_gq(-2:0)*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=+M2_SC_RV_gq(0)*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,mapped_flavours,wgtpl)
c
c     sanity check
      if(abs(M2_SC_RV_gq(0)).ge.huge(1d0).or.isnan(M2_SC_RV_gq(0)))then
         write(77,*)'Exception caught in M2_SC_RV_gq',M2_SC_RV_gq(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
