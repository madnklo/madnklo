

      function M2_C_RV_gg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
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
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision M2_C_RV_gg(-2:0)
      double precision pref,M2tmp(-2:0),M2TMP0,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO,VLO(-2:0),KKVLO(-2:0)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision alphaz
      parameter(alphaz=1d0)
      double precision %(proc_prefix_C_RV_gg)s_GET_KKBLO
      double precision %(proc_prefix_C_RV_gg)s_GET_KKVLO
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
      integer %(proc_prefix_C_RV_gg)s_den
      common/%(proc_prefix_C_RV_gg)s_iden/%(proc_prefix_C_RV_gg)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      double precision xpbsave(0:3,nexternal-1)

c     
c     initialise
      M2_C_RV_gg=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpbsave=xpb
c
c     possible cuts
c      call GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,UnderLying_LEG_PDGS)
      call get_collinear_mapped_labels(ia,ib,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
c     Reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal,underlying_leg_pdgs,mapped_flavours,mapped_labels,xpbsave)
      call invariants_from_p(xpbsave,nexternal-1,xsb,ierr)
      if(ierr.eq.1) goto 999

      IF(DOCUT(XPBSAVE,NEXTERNAL-1,UNDERLYING_LEG_PDGS,0))RETURN      
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
         write(77,*)'Inaccuracy 1 in M2_C_RV_gg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_C_RV_gg)s_ME_ACCESSOR_HOOK(xpbsave,hel,alphas,ANS)
      BLO = ANS(0)
c
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'M2_C_RV_gg: '
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
c     TODO: improve ktmuktnuBmunu / kt^2
      KKBLO = %(proc_prefix_C_RV_gg)s_GET_KKBLO(parent_leg,xpbsave,kt)
      KKVLO = %(proc_prefix_C_RV_gg)s_GET_KKVLO(parent_leg,xpbsave,kt)
      M2TMP0 = CA*2d0*(2d0/sab*KKBLO+x/(1d0-x)*(1d0+1d0-x**alpha)*BLO+(1d0-x)/x*(1d0+1d0-(1d0-x)**alpha)*BLO)

      M2TMP(-2:0) = M2TMP(-2:0) + CA*2d0*(2d0/sab*KKVLO+x/(1d0-x)*(1d0+1d0-x**alpha)*VLO+(1d0-x)/x*(1d0+1d0-(1d0-x)**alpha)*VLO)
      M2TMP(-2) = M2TMP(-2) - alphas/2d0/pi*CA*M2TMP0
      M2TMP(-1) = M2TMP(-1) + alphas/2d0/pi*(CA*log(sab/scale**2)+log(x*(1d0-x))-beta0/2d0)*M2TMP0
      M2TMP( 0) = M2TMP(0) + (CA/2d0*(7d0*zeta2-log(sab/scale**2)**2)-log(sab/scale**2)*log(x*(1d0-x))+ddilog(-(1d0-x)/x)+ddilog(-x/(1d0-x)))*M2TMP0+1d0/sab*KKBLO*CA*(3d0*CA-beta0)
c
c     Including correct multiplicity factor
      M2tmp(-2:0) = M2tmp(-2:0)*dble(%(proc_prefix_C_RV_gg)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir) 
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp(-2:0)=M2tmp(-2:0)*damp

c     compute collinear limit of sector function
      call get_wc_nlo(xs,ia,ib,ir,alphaz,nexternal)
      
      M2_C_RV_gg(-2:0)=M2tmp(-2:0)*pref/sab*xj*extra*wc_nlo
c     apply flavour factor
      M2_C_RV_gg(-2:0)=M2_C_RV_gg(-2:0)*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_C_RV_gg(0)*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbsave,xsb,nexternal-1,UNDERLYING_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_C_RV_gg(0)).ge.huge(1d0).or.isnan(M2_C_RV_gg(0)))then
         write(77,*)'Exception caught in M2_C_RV_gg',M2_C_RV_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      function M2_SC_RV_gg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
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
      double precision M2_SC_RV_gg(-2:0)
      double precision pref,M2tmp(-2:0),wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,VLO(-2:0),EIK0,EIK1(-2:0)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision mb2,mr2,sab,sar,sbr,x,damp,y,xinit
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
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
      integer %(proc_prefix_SC_RV_gg)s_den
      common/%(proc_prefix_SC_RV_gg)s_iden/%(proc_prefix_SC_RV_gg)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      double precision xpbsave(0:3,nexternal-1)
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INCLUDE 'pmass.inc'

c     
c     initialise
      M2_SC_RV_gg=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpbsave=xpb
c     
c     return if not gluon
      if(leg_pdgs(ia).ne.21) then
        write(77,*)'not a gluon in M2_SC_RV_gg',leg_pdgs(ia)
        goto 999
      endif
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c      
c     possible cuts
c      call GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,UnderLying_LEG_PDGS)
      call get_collinear_mapped_labels(ia,ib,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
c     Reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal,underlying_leg_pdgs,mapped_flavours,mapped_labels,xpbsave)
      call invariants_from_p(xpbsave,nexternal-1,xsb,ierr)
      if(ierr.eq.1) goto 999

      IF(DOCUT(XPBSAVE,NEXTERNAL-1,UNDERLYING_LEG_PDGS,0))RETURN
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
         write(77,*)'Inaccuracy 1 in M2_SC_RV_gg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born and Virtual
      call %(proc_prefix_SC_RV_gg)s_ME_ACCESSOR_HOOK(xpbsave,hel,alphas,ANS)
      BLO = ANS(0)
      call %(proc_prefix_SC_RV_gg)s_ME_ACCESSOR_HOOK(xpbsave,hel,alphas,ANS)
      VLO = ANS(0)

      EIK0=SBR/(SAB*SAR) - MB2/SAB**2 - MR2/SAR**2
      EIK1(-2) =  CA*EIK0
      EIK1(-1) = -CA*EIK0*log(sab*sar/sbr/scale**2)
      EIK1( 0) =  CA*EIK0/2d0*(log(sab*sar/sbr/scale**2)**2-5d0*zeta2)
c
      M2TMP(-2:0) = M2TMP(-2:0)+2d0*CA*EIK0*VLO(-2:0)-alphas/2d0/pi*EIK1(-2:0)*BLO
      M2TMP(-1)   = M2TMP(-1)-alphas/2d0/pi*beta0*EIK0*BLO
c     Including correct multiplicity factor
      M2tmp(-2:0) = M2tmp(-2:0)*dble(%(proc_prefix_SC_RV_gg)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir) 
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp(-2:0)=M2tmp(-2:0)*damp
      
      M2_SC_RV_gg(-2:0)=M2tmp(-2:0)*pref*xj*extra
c     apply flavour factor
      M2_SC_RV_gg(-2:0)=M2_SC_RV_gg(-2:0)*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=+M2_SC_RV_gg(0)*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbsave,xsb,nexternal-1,UNDERLYING_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_SC_RV_gg(0)).ge.huge(1d0).or.isnan(M2_SC_RV_gg(0)))then
         write(77,*)'Exception caught in M2_SC_RV_gg',M2_SC_RV_gg(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end  
