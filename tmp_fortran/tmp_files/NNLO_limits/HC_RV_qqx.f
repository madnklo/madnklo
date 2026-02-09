

      SUBROUTINE SUB_M2_HC_RV_qqx(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr,ret)
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
      double precision M2_C_qqx(-2:0),M2_C_qqx_0
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO,VLO(-2:0),KKVLO(-2:0)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,logab,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision alphaz
      parameter(alphaz=1d0)
      double precision %(proc_prefix_HC_RV_qqx)s_GET_KKBLO
      double precision %(proc_prefix_HC_RV_qqx)s_GET_KKVLO
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
      integer %(proc_prefix_HC_RV_qqx)s_den
      common/%(proc_prefix_HC_RV_qqx)s_iden/%(proc_prefix_HC_RV_qqx)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
c
c     initialise
      ret=0d0
      M2_C_qqx_0=0d0
      M2_C_qqx=0d0
      ierr=0
      damp=0d0
c
c     checks
      if(.not.(abs(leg_pdgs(ia)).le.6.and.leg_pdgs(ia)+leg_pdgs(ib).eq.0.and.underlying_leg_pdgs(mapped_labels(ib)).eq.21))then
         write(*,*)'Wrong pdgs in M2_HC_RV_qqx',leg_pdgs(ia),leg_pdgs(ib),mapped_labels(ib)
         stop
      endif
      if(.not.(ia.eq.isec.and.ib.eq.jsec))then
         write(*,*)'Wrong indices in M2_HC_RV_qqx',ia,ib,isec,jsec
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
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
      logab=log(sab/scale**2)
c
c     coefficients of kt
c     kt = wa pa + wb pb + wr pr
      wa= 1d0 - x
      wb= - x
      wr= - (1d0-2d0*x)*sab/(sar+sbr)
      kt(:) = wa*xp(:,ia) + wb*xp(:,ib) + wr*xp(:,ir)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HC_RV_qqx',sab,sar+sbr,x
         goto 999
      endif
c
      parent_leg = mapped_labels(ib)
c     TODO: improve ktmuktnuBmunu / kt^2
c     call Born and spin-correlated Born
c      call %(proc_prefix_HC_RV_qqx)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = 0d0 !ANS(0)
      KKBLO = 0d0 !%(proc_prefix_HC_RV_qqx)s_GET_KKBLO(parent_leg,xpb,kt)
c     call Virtual and spin-correlated Virtual
c      call %(proc_prefix_HC_RV_qqx)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      VLO = 0d0 !ANS(0)
      KKVLO = 0d0 !%(proc_prefix_HC_RV_qqx)s_GET_KKVLO(parent_leg,xpb,kt)
c
      M2_C_qqx_0 = TR*(BLO-4d0/sab*KKBLO)

      M2_C_qqx(-2:0) = M2_C_qqx(-2:0) + TR*(BLO-4d0/sab*KKBLO)
      M2_C_qqx(-2) = M2_C_qqx(-2) + alphas/2d0/pi*M2_C_qqx_0*(CA-2d0*CF)
      M2_C_qqx(-1) = M2_C_qqx(-1) + alphas/2d0/pi*(M2_C_qqx_0*((2d0*CF-CA)*logab+CA*log(x*(1d0-x))-beta0/2d0)+(4d0*x*(1d0-x)/sab*KKBLO-BLO)*TR*(beta0-3d0*CF))
      M2_C_qqx( 0) = M2_C_qqx( 0) + alphas/2d0/pi*(M2_C_qqx_0*((CA-2d0*CF)*(7d0*zeta2-logab**2)/2d0+CA*(-logab*log(x*(1d0-x))+ddilog(-(1d0-x)/x)+ddilog(-x/(1d0-x))))+(4d0*x*(1d0-x)/sab*KKBLO-BLO)*TR*(-logab*(beta0-3d0*CF)+7d0*CA/3d0+5d0*beta0/3d0-8d0*CF))
c
c     compute collinear limit of sector function
      call get_wc_nlo(xs,isec,jsec,iref,alphaz,nexternal)
      M2_C_qqx = M2_C_qqx*wc_nlo
c     account for different damping factors according to recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      ret = M2_C_qqx
c     include prefactors
      ret = ret *dble(%(proc_prefix_HC_RV_qqx)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor*damp*pref/sab*xj*extra
c
c     plot
      wgtpl=-ret(0)*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgtpl)
c
c     sanity check
      if(abs(ret(0)).ge.huge(1d0).or.isnan(ret(0)))then
         write(77,*)'Exception caught in M2_HC_RV_qqx',ret(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
