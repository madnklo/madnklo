      double precision function M2_S(i,xs,xp,wgt,ZSoft,xj,xjB,nit,extra,wgt_chan,ierr)
c     single-soft limit S_(i) * Zsoft
c     it returns 0 if i is not a gluon
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'      
      integer i,l,m,lb,mb,ierr,nit,idum
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,Zsoft,xj,xjB,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,ml2,mm2,y,z,x,damp
      integer mapped_labels(nexternal), mapped_flavours(NEXTERNAL)
      logical isLOQCDparton(nexternal-1)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
c     external
      integer get_color_dipole_index
      external get_color_dipole_index
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
      integer, parameter :: HEL = - 1
      double precision  %(proc_prefix_S)s_GET_CCBLO
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_S)s_den
      common/%(proc_prefix_S)s_iden/%(proc_prefix_S)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INCLUDE 'pmass.inc'
      
c
c     initialise
      M2_S=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      idum=0
c
c     return if not gluon
      if(leg_pdgs(I).ne.21)return
c
c     safety check on PDGs
      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
        WRITE(*,*) 'M2_S:'
        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
        STOP
      ENDIF
c
c     get PDGs
      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      CALL GET_SOFT_MAPPED_LABELS(I,idum,idum,NEXTERNAL,LEG_PDGS,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref=-8d0*pi*alphas
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isNLOQCDparton(m))cycle
         if(m.eq.i)cycle
         do l=1,nexternal
            if(.not.isNLOQCDparton(l))cycle
            if(l.eq.i)cycle
            if(l.eq.m)cycle
c
            lb=mapped_labels(l)
            mb=mapped_labels(m)
c
c         check labels and pdgs
          IF(.NOT.(ISLOQCDPARTON(LB).AND.ISLOQCDPARTON(MB)))THEN
            WRITE(*,*)'Wrong indices 1 in M2_S',LB,MB
            STOP
          ENDIF
          IF(leg_pdgs(l).ne.Born_leg_pdgs(lb).or.leg_pdgs(m).ne.Born_leg_pdgs(mb))THEN
            WRITE(*,*)'Wrong indices 2 in M2_S',L,M,LB,MB
            STOP
          ENDIF
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the single-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,'S',xjCS)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))CYCLE
c
c     invariant quantities
            sil=xs(i,l)
            sim=xs(i,m)
            slm=xs(l,m)
            ml2=pmass(l)**2
            mm2=pmass(m)**2
c
c     safety check
            if(sil*sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S',sil,sim
               goto 999
            endif
c
c     call colour-connected Born
            call %(proc_prefix_S)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S)s_GET_CCBLO(lb,mb)
c
c     eikonal
            M2TMP=SLM/(SIL*SIM) - ML2/SIL**2 - MM2/SIM**2
            M2TMP = CCBLO*M2TMP
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_S)s_den)/dble(%(proc_prefix_real)s_den)
c
c     damping factors
            if(m.gt.2.and.l.gt.2)then
               y=sil/(sil+sim+slm)
               z=sim/(sim+slm)
               damp=((1d0-y)*(1d0-z))**alpha
            elseif(m.gt.2.and.l.le.2)then
               z=sim/(sim+slm)
               x=1d0 - sil/(sim+slm)
               damp=((1d0-z)*x)**alpha
            elseif(m.le.2.and.l.le.2)then
               x=1d0 - (sil+sim)/slm
               damp=x**alpha
            endif
            M2tmp=M2tmp*damp*xj
            M2_S=M2_S+pref*M2tmp*Zsoft*extra
c
c     plot
            wgtpl=-pref*M2tmp*Zsoft*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_real)s_fl_factor
            if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
         enddo 
      enddo
c
c     apply flavour factor
      M2_S = M2_s * %(proc_prefix_real)s_fl_factor
c
c     sanity check
      if(abs(M2_S).ge.huge(1d0).or.isnan(M2_S))then
         write(77,*)'Exception caught in M2_S',M2_S
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end



      
      DOUBLE PRECISION FUNCTION M2_S_ALT(I,IB,IR,XS,XP,XSB,XPB,WGT,ZSOFT,XJ,XJB,NIT,EXTRA,wgt_chan,IERR)
C     single-soft limit S_(i) * Zsoft, mapped as the collinear one
C     it returns 0 if i is not a gluon
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'colored_partons.inc'
      INCLUDE 'leg_PDGs.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INTEGER I,L,M,IB,IR,LB,MB,NIT,IERR,PARENT_LEG,idum
      DOUBLE PRECISION PREF,M2TMP,WGT,WGTPL,wgt_chan,ZSOFT,XJ,XJB,EXTRA
      DOUBLE PRECISION XS(NEXTERNAL,NEXTERNAL),XSB(NEXTERNAL-1,NEXTERNAL-1)
      DOUBLE PRECISION BLO,CCBLO
      DOUBLE PRECISION XP(0:3,NEXTERNAL),XPB(0:3,NEXTERNAL-1)
      DOUBLE PRECISION SIL,SIM,SLM,X,Y,Z,DAMP
      DOUBLE PRECISION ANS(0:NSQSO_BORN)
      INTEGER MAPPED_LABELS(NEXTERNAL),MAPPED_FLAVOURS(NEXTERNAL)
      INTEGER, PARAMETER :: HEL = - 1
      DOUBLE PRECISION ALPHAS,ALPHA_QCD
      LOGICAL ISLOQCDPARTON(NEXTERNAL-1)
C     set logical doplot
      LOGICAL DOPLOT
      COMMON/CDOPLOT/DOPLOT
      DOUBLE PRECISION SCM
      COMMON/CSCM/SCM
      LOGICAL DOCUT
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      INTEGER GET_COLOR_DIPOLE_INDEX
      EXTERNAL GET_COLOR_DIPOLE_INDEX
      double precision  %(proc_prefix_S)s_GET_CCBLO
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_S)s_den
      common/%(proc_prefix_S)s_iden/%(proc_prefix_S)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
C     
C     initialise
      M2_S_ALT=0D0
      M2TMP=0D0
      IERR=0
      DAMP=0D0
      idum=0
c
c     return if not gluon
      if(leg_pdgs(I).ne.21)return
c
c     safety check on PDGs
      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
        STOP
      ENDIF
C     
C     get PDGs and possible cuts
      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      CALL GET_SOFT_MAPPED_LABELS(I,IB,IR,NEXTERNAL,LEG_PDGS,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
C     
C     possible cuts
      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN
C
C     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      PREF=-8D0*PI*ALPHAS
C
C     call colour-connected Born outside of the eikonal sum
      CALL %(proc_prefix_S)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
C     
C     eikonal double sum
      DO M=1,NEXTERNAL-1
        IF(.NOT.ISNLOQCDPARTON(M))CYCLE
        IF(M.EQ.I)CYCLE
        DO L=M+1,NEXTERNAL
          IF(.NOT.ISNLOQCDPARTON(L))CYCLE
          IF(L.EQ.I)CYCLE
c
          LB=MAPPED_LABELS(L)
          MB=MAPPED_LABELS(M)
c
c         check labels and pdgs
          IF(.NOT.(ISLOQCDPARTON(LB).AND.ISLOQCDPARTON(MB)))THEN
            WRITE(*,*)'Wrong indices 1 in M2_S_ALT',LB,MB
            STOP
          ENDIF
          IF(leg_pdgs(l).ne.Born_leg_pdgs(lb).or.leg_pdgs(m).ne.Born_leg_pdgs(mb))THEN
            WRITE(*,*)'Wrong indices 2 in M2_S_ALT',L,M,LB,MB
            STOP
          ENDIF
C         
C         invariant quantities
          SIL=XS(I,L)
          SIM=XS(I,M)
          SLM=XS(L,M)
C
C         safety check
          IF(SIL*SIM.LE.0D0)THEN
            WRITE(77,*)'Inaccuracy 1 in M2_S_ALT',SIL,SIM
            GOTO 999
          ENDIF
C
C         eikonal
          CCBLO = %(proc_prefix_S)s_GET_CCBLO(lb,mb)
          M2TMP=CCBLO*2D0*SLM/(SIL*SIM)
c
C         Including correct multiplicity factor
          M2tmp = M2tmp*dble(%(proc_prefix_S)s_den)/dble(%(proc_prefix_real)s_den)
c
c         Damping factors
          IF(M.GT.2.AND.L.GT.2)THEN
            Y=SIL/(SIL+SIM+SLM)
            Z=SIM/(SIM+SLM)
            DAMP=((1D0-Y)*(1D0-Z))**ALPHA
          ELSEIF(M.GT.2.AND.L.LE.2)THEN
            Z=SIM/(SIM+SLM)
            X=1D0 - SIL/(SIM+SLM)
            DAMP=((1D0-Z)*X)**ALPHA
          ELSEIF(M.LE.2.AND.L.LE.2)THEN
            X=1D0 - (SIL+SIM)/SLM
            DAMP=X**ALPHA
          ENDIF
          M2TMP=M2TMP*DAMP
          M2_S_ALT=M2_S_ALT+PREF*M2TMP*ZSOFT*XJ*EXTRA
        ENDDO
      ENDDO
C     apply flavour factor
      M2_S_ALT=M2_S_ALT*%(proc_prefix_real)s_fl_factor
C         
C     plot
      WGTPL=-M2_S_ALT*WGT/NIT*wgt_chan
      IF(DOPLOT)CALL HISTO_FILL(XPB,XSB,NEXTERNAL-1,WGTPL)
C     
C     sanity check
      IF(ABS(M2_S_ALT).GE.HUGE(1D0).OR.ISNAN(M2_S_ALT))THEN
        WRITE(77,*)'Exception caught in M2_S_ALT',M2_S_ALT
        GOTO 999
      ENDIF
C     
      RETURN
 999  IERR=1
      RETURN
      END


      double precision function M2_HC_gg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     hard-collinear limit C_(ia,ib) - S_(ia)C_(ia,ib) - S_(ib)C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision %(proc_prefix_HC_gg)s_GET_KKBLO
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
      integer %(proc_prefix_HC_gg)s_den
      common/%(proc_prefix_HC_gg)s_iden/%(proc_prefix_HC_gg)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)

c
c     initialise
      M2_HC_gg=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c     possible cuts
      call GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN
      
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
         write(77,*)'Inaccuracy 1 in M2_HC_gg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_HC_gg)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(ia,ib,ir,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'M2_HC_gg: '
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      KKBLO = %(proc_prefix_HC_gg)s_GET_KKBLO(parent_leg,xpb,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=CA*2d0*(2d0/sab*KKBLO+x/(1d0-x)*(1d0-x**alpha)*BLO+(1d0-x)/x*(1d0-(1d0-x)**alpha)*BLO)
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_HC_gg)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir) 
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_HC_gg=M2tmp*pref/sab*xj*extra
c     apply flavour factor
      M2_HC_gg=M2_HC_gg*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_HC_gg*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_HC_gg).ge.huge(1d0).or.isnan(M2_HC_gg))then
         write(77,*)'Exception caught in M2_HC_gg',M2_HC_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
  
                  
      double precision function M2_HC_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     hard-collinear limit C_(ia,ib) - S_(ia)C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'      
      integer ia,ib,ir,ierr,nit
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision ans(0:nsqso_born)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      double precision alphas,alpha_qcd
      integer,parameter :: HEL = - 1
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_HC_gq)s_den
      common/%(proc_prefix_HC_gq)s_iden/%(proc_prefix_HC_gq)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)

c
c     initialise
      M2_HC_gq=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c     possible cuts
      call GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      if(leg_pdgs(ia).eq.21) then
         x=sbr/(sar+sbr)
      else
         x=sar/(sar+sbr)
      endif
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HC_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_HC_gq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c     In the following equation the x variable is related to the quark energy
      M2tmp=BLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0-x**alpha))
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_HC_gq)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_HC_gq=M2tmp*pref/sab*xj*extra
c     apply flavour factor
      M2_HC_gq=M2_HC_gq*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_HC_gq*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_HC_gq).ge.huge(1d0).or.isnan(M2_HC_gq))then
         write(77,*)'Exception caught in M2_HC_gq',M2_HC_gq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_HC_qqx(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     hard-collinear limit C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      double precision alphas,alpha_qcd
      double precision %(proc_prefix_HC_qqx)s_get_kkblo
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_HC_qqx)s_den
      common/%(proc_prefix_HC_qqx)s_iden/%(proc_prefix_HC_qqx)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
c
c     initialise
      M2_HC_qqx=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c     possible cuts
      call GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN


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
         write(77,*)'Inaccuracy 1 in M2_HC_qqx',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_HC_qqx)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(ia,ib,ir,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      KKBLO = %(proc_prefix_HC_qqx)s_GET_KKBLO(parent_leg,xpb,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=TR*(BLO-4d0/sab*KKBLO)
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_HC_qqx)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_HC_qqx=M2tmp*pref/sab*xj*extra
c     apply flavour factor
      M2_HC_qqx=M2_HC_qqx*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_HC_qqx*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_HC_qqx).ge.huge(1d0).or.isnan(M2_HC_qqx))then
         write(77,*)'Exception caught in M2_HC_qqx',M2_HC_qqx
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      SUBROUTINE DUMMY_ME_ACCESSOR_HOOK(P,HEL,USER_ALPHAS,ANS)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'nsqso_born.inc'
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS(0:NSQSO_BORN)
      INTEGER HEL
      DOUBLE PRECISION USER_ALPHAS

      write(*,*) 'This subroutine should never be called!!!'
      write(*,*) 'This counterterm does not contribute to this sector.'
      write(*,*) 'Exit...'
      STOP
      END


      DOUBLE PRECISION FUNCTION DUMMY_GET_CCBLO(LB,MB)
      IMPLICIT NONE
      INTEGER LB,MB

      WRITE(*,*) 'YOU ARE IN DUMMY_GET_CCBLO'
      WRITE(*,*) 'THIS FUNCTION MUST NEVER BE CALLED!'
      WRITE(*,*) 'EXIT...'

      STOP
      END


      DOUBLE PRECISION FUNCTION DUMMY_GET_KKBLO(PARENT_LEG,XPB,KT)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INTEGER PARENT_LEG
      DOUBLE PRECISION XPB(0:3,NEXTERNAL-1),KT(0:3)

      WRITE(*,*) 'YOU ARE IN DUMMY_GET_KKBLO'
      WRITE(*,*) 'THIS FUNCTION MUST NEVER BE CALLED!'
      WRITE(*,*) 'EXIT...'

      STOP
      END







c c STUFF NOT WORKING
c c      double precision function M2_S2(i,xs,xp,wgt,ZSoft,xj,xjB,XX,nit,extra,ierr)
c cc     single-soft limit S_(i) * Zsoft
c cc     it returns 0 if i is not a gluon
c c      implicit none
c c      include 'nexternal.inc'
c c      INCLUDE 'coupl.inc'
c c      include 'math.inc'
c c      include 'damping_factors.inc'
c c      include 'colored_partons.inc'
c c      include 'leg_PDGs.inc'
c c      include 'nsqso_born.inc'
c c      INCLUDE 'input.inc'
c c      INCLUDE 'run.inc'      
c c      integer i,l,m,lb,mb,ierr,nit,idum,ia
c c      double precision pref,M2tmp,wgt,wgtpl,Zsoft,xj,xjB,xjCS,xx(3),ddum
c c      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
c c      double precision BLO,ccBLO,extra
c c      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
c c      double precision sil,sim,slm,y,z,x,damp
c c      integer mapped_labels(nexternal), mapped_flavours(NEXTERNAL)
c c      logical isLOQCDparton(nexternal-1)
c cc     set logical doplot
c c      logical doplot
c c      common/cdoplot/doplot
c c      double precision sCM
c c      common/cscm/sCM
c c      logical docut
c c      integer %(proc_prefix_real)s_fl_factor
c c      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
c cc     external
c c      integer get_color_dipole_index
c c      external get_color_dipole_index
c c      double precision alphas,ans(0:NSQSO_BORN)
c c      double precision alpha_qcd
c c      integer, parameter :: HEL = - 1
c c      double precision  %(proc_prefix_S)s_GET_CCBLO
c c      integer %(proc_prefix_real)s_den
c c      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
c c      integer %(proc_prefix_S)s_den
c c      common/%(proc_prefix_S)s_iden/%(proc_prefix_S)s_den
c c      INTEGER ISEC,JSEC
c c      COMMON/CNLOSECINDICES/ISEC,JSEC
c c      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
c cc
c cc     initialise
c c      M2_S2=0d0
c c      M2tmp=0d0
c c      ierr=0
c c      damp=0d0
c c      idum=0
c c      Zsoft=0d0
c cc
c cc     return if not gluon
c c      if(leg_pdgs(I).ne.21)return
c cc
c cc     safety check on PDGs
c c      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
c c        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
c c        STOP
c c      ENDIF
c cc
c cc     get PDGs
c c      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
c c      CALL GET_SOFT_MAPPED_LABELS(I,idum,idum,NEXTERNAL,LEG_PDGS,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
c cC
c cC     possible cuts
c c      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN
c cc
c cc     overall kernel prefix
c c      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c c      pref=-8d0*pi*alphas
c cC
c cC     call colour-connected Born outside of the eikonal sum
c c      CALL %(proc_prefix_S)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
c cc
c cc     eikonal double sum
c c      do m=1,nexternal-1
c c         if(.not.isNLOQCDparton(m))cycle
c c         if(m.eq.i)cycle
c c         do l=m+1,nexternal
c c            if(.not.isNLOQCDparton(l))cycle
c c            if(l.eq.i)cycle
c cc
c c            lb=mapped_labels(l)
c c            mb=mapped_labels(m)
c cc
c cc         check labels and pdgs
c c          IF(.NOT.(ISLOQCDPARTON(LB).AND.ISLOQCDPARTON(MB)))THEN
c c            WRITE(*,*)'Wrong indices 1 in M2_S',LB,MB
c c            STOP
c c          ENDIF
c c          IF(leg_pdgs(l).ne.Born_leg_pdgs(lb).or.leg_pdgs(m).ne.Born_leg_pdgs(mb))THEN
c c            WRITE(*,*)'Wrong indices 2 in M2_S',L,M,LB,MB
c c            STOP
c c          ENDIF
c cc
c cc     phase-space mapping according to l and m, at fixed radiation
c cc     phase-space point: the singular kernel is in the same point
c cc     as the single-real, ensuring numerical stability, while the
c cc     underlying Born configuration is remapped
c c          iA = 1                ! default azimuth for NLO
c c              call phase_space_CS(xx,i,l,m,ia,xp,xpb,nexternal,leg_PDGs,'S',xjCS)
c c            if(xjCS.eq.0d0)goto 999
c c            call invariants_from_p(xp,nexternal,xsb,ierr)
c c            if(ierr.eq.1)goto 999
c c            call get_Z_NLO(xs,ddum,ddum,ISEC,JSEC,Zsoft,'S',ierr)
c c            IF(IERR.EQ.1)GOTO 999
c cc
c cc     invariant quantities
c c            sil=xs(i,l)
c c            sim=xs(i,m)
c c            slm=xs(l,m)
c cc
c cc     safety check
c c            if(sil*sim.le.0d0)then
c c               write(77,*)'Inaccuracy 1 in M2_S',sil,sim
c c               goto 999
c c            endif
c cc
c cc     eikonal
c c            ccBLO = %(proc_prefix_S)s_GET_CCBLO(lb,mb)
c c            M2tmp=ccBLO*2d0*slm/(sil*sim) * Zsoft * xjB * xjCS
c cc
c cc     Including correct multiplicity factor
c c            M2tmp = M2tmp*dble(%(proc_prefix_S)s_den)/dble(%(proc_prefix_real)s_den)
c cc
c cc     damping factors
c c          if(alpha.ne.0d0)then
c c             write(*,*)'Implement damping factors in M2_S2!'
c c             stop
c c          endif
c c            if(m.gt.2.and.l.gt.2)then
c c               y=sil/(sil+sim+slm)
c c               z=sim/(sim+slm)
c c               damp=((1d0-y)*(1d0-z))**alpha
c c            elseif(m.gt.2.and.l.le.2)then
c c               z=sim/(sim+slm)
c c               x=1d0 - sil/(sim+slm)
c c               damp=((1d0-z)*x)**alpha
c c            elseif(m.le.2.and.l.le.2)then
c c               x=1d0 - (sil+sim)/slm
c c               damp=x**alpha
c c            endif
c c            M2tmp=M2tmp*damp
c c            M2_S2=M2_S2+pref*M2tmp*Zsoft*extra
c cc
c c         enddo 
c c      enddo
c cc
c cc     apply flavour factor
c c      M2_S2 = M2_s2 * %(proc_prefix_real)s_fl_factor
c cc
c cc     plot
c c      wgtpl=-M2_S2*wgt/nit
c c      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c cc
c cc     sanity check
c c      if(abs(M2_S2).ge.huge(1d0).or.isnan(M2_S2))then
c c         write(77,*)'Exception caught in M2_S2',M2_S2
c c         goto 999
c c      endif
c cc
c c      return
c c 999  ierr=1
c c      return
c c      end
c c
c c
c c
c c      DOUBLE PRECISION FUNCTION M2_S_DIFF(I,IB,IR,XS,XP,XSB,XPB,WGT,ZSOFT,XJ,XJB,XX,NIT,EXTRA,IERR)
c cC     difference from the soft counterterm mapped according to (ilm)
c cC     and the soft couterterm mapped according to (ijr), multiplied
c cC     by Zsoft
c cC     it returns 0 if i is not a gluon
c c      IMPLICIT NONE
c c      INCLUDE 'nexternal.inc'
c c      INCLUDE 'coupl.inc'
c c      INCLUDE 'math.inc'
c c      INCLUDE 'damping_factors.inc'
c c      INCLUDE 'colored_partons.inc'
c c      INCLUDE 'leg_PDGs.inc'
c c      INCLUDE 'nsqso_born.inc'
c c      INCLUDE 'input.inc'
c c      INCLUDE 'run.inc'
c c      INTEGER I,L,M,iB,iR,iA,LB,MB,IERR,NIT,idum
c c      DOUBLE PRECISION PREF,M2TMP,WGT,WGTPL,ZSOFT,ZSOFT_ilm,ddum
c c      DOUBLE PRECISION XJ,XJB,XJCS_ILM
c c      DOUBLE PRECISION XS(NEXTERNAL,NEXTERNAL),XSB(NEXTERNAL-1,NEXTERNAL-1),XS_ilm(NEXTERNAL,NEXTERNAL),XX(3)
c c      DOUBLE PRECISION BLO,CCBLO,EXTRA
c c      DOUBLE PRECISION XP(0:3,NEXTERNAL),XPB(0:3,NEXTERNAL-1)
c c      DOUBLE PRECISION XP_ilm(0:3,NEXTERNAL)
c c      DOUBLE PRECISION SIL,SIM,SLM,Y,Z,X,DAMP
c c      DOUBLE PRECISION SIL_ilm,SIM_ilm,SLM_ilm
c c      INTEGER MAPPED_LABELS(NEXTERNAL), MAPPED_FLAVOURS(NEXTERNAL)
c c      LOGICAL ISLOQCDPARTON(NEXTERNAL-1)
c cC     set logical doplot
c c      LOGICAL DOPLOT
c c      COMMON/CDOPLOT/DOPLOT
c c      DOUBLE PRECISION SCM
c c      COMMON/CSCM/SCM
c c      LOGICAL DOCUT
c c      integer %(proc_prefix_real)s_fl_factor
c c      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
c cC     external
c c      INTEGER GET_COLOR_DIPOLE_INDEX
c c      EXTERNAL GET_COLOR_DIPOLE_INDEX
c c      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
c c      DOUBLE PRECISION ALPHA_QCD
c c      INTEGER, PARAMETER :: HEL = - 1
c c      double precision  %(proc_prefix_S)s_GET_CCBLO
c c      integer %(proc_prefix_real)s_den
c c      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
c c      integer %(proc_prefix_S)s_den
c c      common/%(proc_prefix_S)s_iden/%(proc_prefix_S)s_den
c c      INTEGER ISEC,JSEC
c c      COMMON/CNLOSECINDICES/ISEC,JSEC
c c      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
c c      integer iseed
c cC
c cC     initialise
c c      M2_S_DIFF=0D0
c c      M2TMP=0D0
c c      IERR=0
c c      DAMP=0D0
c c      idum=0
c c      ddum=1d0
c cc
c cc     return if not gluon
c c      if(leg_pdgs(I).ne.21)return
c cc
c cc     safety check on PDGs
c c      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
c c        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
c c        STOP
c c      ENDIF
c cc
c cc     get PDGs
c c      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
c c      CALL GET_SOFT_MAPPED_LABELS(I,idum,idum,NEXTERNAL,LEG_PDGS,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
c cC     
c cC     possible cuts
c c      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN
c cC
c cC     overall kernel prefix
c c      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
c c      PREF=-8D0*PI*ALPHAS
c cC
c cC     call colour-connected Born outside of the eikonal sum
c c      call %(proc_prefix_S)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
c cC
c cC     eikonal double sum
c c      DO M=1,NEXTERNAL-1
c c        IF(.NOT.ISNLOQCDPARTON(M))CYCLE
c c        IF(M.EQ.I)CYCLE
c c        DO L=M+1,NEXTERNAL
c c          IF(.NOT.ISNLOQCDPARTON(L))CYCLE
c c          IF(L.EQ.I)CYCLE
c cc
c c          LB=MAPPED_LABELS(L)
c c          MB=MAPPED_LABELS(M)
c cc
c cc         check labels and pdgs 
c c          IF(.NOT.(ISLOQCDPARTON(LB).AND.ISLOQCDPARTON(MB)))THEN
c c            WRITE(*,*)'Wrong indices 1 in M2_S_DIFF',LB,MB
c c            STOP
c c          ENDIF
c c          IF(leg_pdgs(l).ne.Born_leg_pdgs(lb).or.leg_pdgs(m).ne.Born_leg_pdgs(mb))THEN
c c            WRITE(*,*)'Wrong indices 2 in M2_S_DIFF',L,M,LB,MB
c c            STOP
c c          ENDIF
c cc
c cc         invariant quantities with mapping (ijr), at fixed Born kinematics
c c          SIL=XS(I,L)
c c          SIM=XS(I,M)
c c          SLM=XS(L,M)
c cc         invariant quantities with mapping (ilm), at fixed Born kinematics
c c          iA = 1  ! default azimuth for NLO
c cc TODO: understand below I,M,L assignment
c c          CALL PHASE_SPACE_CS(XX,I,M,L,IA,XP_ILM,XPB,NEXTERNAL,LEG_PDGS,'S',XJCS_ILM)
c c          IF(XJCS_ILM.EQ.0D0)GOTO 999
c c          CALL INVARIANTS_FROM_P(XP_ilm,NEXTERNAL,XS_ilm,IERR)
c c          IF(IERR.EQ.1)GOTO 999
c c          call get_Z_NLO(xs_ilm,ddum,ddum,ISEC,JSEC,Zsoft_ilm,'S',ierr)
c c          IF(IERR.EQ.1)GOTO 999
c c          SIL_ilm=XS_ilm(I,L)
c c          SIM_ilm=XS_ilm(I,M)
c c          SLM_ilm=XS_ilm(L,M)
c cC         
c cC         safety check
c c          IF(SIL*SIM.LE.0D0.or.SIL_ilm*SIM_ilm.LE.0D0)THEN
c c            WRITE(77,*)'Inaccuracy 1 in M2_S_DIFF'
c c            WRITE(77,*)SIL,SIM,SIL_ilm,SIM_ilm
c c            GOTO 999
c c          ENDIF
c cC         
c cC         eikonal difference
c c          ccBLO = %(proc_prefix_S)s_GET_CCBLO(lb,mb)
c c          M2TMP = SLM_ilm/(SIL_ilm*SIM_ilm) * ZSOFT_ilm * XJB * XJCS_ILM
c cc          M2TMP = SLM_ilm/(SIL_ilm*SIM_ilm) * ZSOFT * XJB * XJCS_ILM
c c          M2TMP = M2TMP - SLM/(SIL*SIM) * ZSOFT * XJ
c c          M2TMP = M2TMP * CCBLO*2D0
c cc
c cC         Including correct multiplicity factor
c c          M2tmp = M2tmp*dble(%(proc_prefix_S)s_den)/dble(%(proc_prefix_real)s_den)
c cc
c cc         Damping factors
c c          if(alpha.ne.0d0)then
c c             write(*,*)'Implement damping factors in M2_S_DIFF!'
c c             stop
c c          endif
c c          IF(M.GT.2.AND.L.GT.2)THEN
c c            Y=SIL/(SIL+SIM+SLM)
c c            Z=SIM/(SIM+SLM)
c c            DAMP=((1D0-Y)*(1D0-Z))**ALPHA
c c          ELSEIF(M.GT.2.AND.L.LE.2)THEN
c c            Z=SIM/(SIM+SLM)
c c            X=1D0 - SIL/(SIM+SLM)
c c            DAMP=((1D0-Z)*X)**ALPHA
c c          ELSEIF(M.LE.2.AND.L.LE.2)THEN
c c            X=1D0 - (SIL+SIM)/SLM
c c            DAMP=X**ALPHA
c c          ENDIF
c c          M2TMP=M2TMP*DAMP
c c          M2_S_DIFF=M2_S_DIFF+PREF*M2TMP*EXTRA
c c        ENDDO
c c      ENDDO
c cC     apply flavour factor
c c      M2_S_DIFF = M2_S_DIFF * %(proc_prefix_real)s_fl_factor
c cC
c cC     plot
c c      WGTPL=-M2_S_DIFF*WGT/NIT
c c      IF(DOPLOT)CALL HISTO_FILL(XPB,XSB,NEXTERNAL-1,WGTPL)
c cC     
c cC     sanity check
c c      IF(ABS(M2_S_DIFF).GE.HUGE(1D0).OR.ISNAN(M2_S_DIFF))THEN
c c        WRITE(77,*)'Exception caught in M2_S_DIFF',M2_S_DIFF
c c        GOTO 999
c c      ENDIF
c cC     
c c      RETURN
c c 999  IERR=1
c c      RETURN
c c      END

