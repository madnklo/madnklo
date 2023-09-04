      double precision function M2_S(i,xs,xp,wgt,Wsoft,xj,nit,extra,ierr)
c     single-soft limit S_(i) * Wsoft
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
      double precision pref,M2tmp,wgt,wgtpl,Wsoft,xj,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,y,z,x,damp
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
c
c     initialise
      M2_S=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      idum=0
c
c     get PDGs
      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      CALL GET_SOFT_MAPPED_LABELS(I,idum,idum,NEXTERNAL,LEG_PDGS
     $,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
         WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),
     &   NEXTERNAL
         STOP
      ENDIF
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref=-8d0*pi*alphas
c
c     eikonal double sum
      do m=1,nexternal-1
         if(.not.isNLOQCDparton(m))cycle
         if(m.eq.i)cycle
         do l=m+1,nexternal
            if(.not.isNLOQCDparton(l))cycle
            if(l.eq.i)cycle
c
            lb=mapped_labels(l)
            mb=mapped_labels(m)
c     check LO color labels 
            if(.not.(isLOQCDparton(lb).and.isLOQCDparton(mb)))then
               write(*,*)'Wrong LO indices in soft kernel',lb,mb
               stop
            endif
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
            M2tmp=ccBLO*2d0*slm/(sil*sim)
c
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
            M2tmp=M2tmp*damp
            M2_S=M2_S+pref*M2tmp*Wsoft*extra
c
c     plot
            wgtpl=-pref*M2tmp*Wsoft*extra*xj*wgt/nit
            wgtpl = wgtpl*%(proc_prefix_real)s_fl_factor
            if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
         enddo 
      enddo

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



      DOUBLE PRECISION FUNCTION M2_S_ALT(I,IB,IR,XS,XP,XSB,XPB,WGT
     $ ,WSOFT,XJ,NIT,EXTRA,IERR)
C     single-soft limit S_(i) * Wsoft, mapped as the collinear one
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
      DOUBLE PRECISION PREF,M2TMP,WGT,WGTPL,WSOFT,XJ,EXTRA
      DOUBLE PRECISION XS(NEXTERNAL,NEXTERNAL),XSB(NEXTERNAL-1
     $ ,NEXTERNAL-1)
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
c     possible cuts
      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN
C     
C     get PDGs and possible cuts
      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      CALL GET_SOFT_MAPPED_LABELS(I,IB,IR,NEXTERNAL,LEG_PDGS
     $ ,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
         WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),
     &   NEXTERNAL
         STOP
      ENDIF
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
C         check LO color labels 
          IF(.NOT.(ISLOQCDPARTON(LB).AND.ISLOQCDPARTON(MB)))THEN
            WRITE(*,*)'Wrong LO indices in soft kernel',LB,MB
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
c         TODO: add check that flavour(L)=flavour(LB) and so for M, MB
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
          M2_S_ALT=M2_S_ALT+PREF*M2TMP*WSOFT*EXTRA
        ENDDO
      ENDDO
C     apply flavour factor
      M2_S_ALT=M2_S_ALT*%(proc_prefix_real)s_fl_factor
C         
C     plot
      WGTPL=-M2_S_ALT*XJ*WGT/NIT
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



      DOUBLE PRECISION FUNCTION M2_S_DIFF(I,IB,IR,XS,XP,XSB,XPB,WGT
     $ ,WSOFT,XJ,XJB,XX,NIT,EXTRA,IERR)
C     difference from the soft counterterm mapped according to (ilm)
C     and the soft couterterm mapped according to (ijr), multiplied
C     by Wsoft
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
      INTEGER I,L,M,iB,iR,iA,LB,MB,IERR,NIT,idum
      DOUBLE PRECISION PREF,M2TMP,WGT,WGTPL,WSOFT
      DOUBLE PRECISION XJ,XJB,XJCS_ILM
      DOUBLE PRECISION XS(NEXTERNAL,NEXTERNAL),XSB(NEXTERNAL-1
     $ ,NEXTERNAL-1),XS_ilm(NEXTERNAL,NEXTERNAL),XX(3)
      DOUBLE PRECISION BLO,CCBLO,EXTRA
      DOUBLE PRECISION XP(0:3,NEXTERNAL),XPB(0:3,NEXTERNAL-1)
      DOUBLE PRECISION XP_ilm(0:3,NEXTERNAL)
      DOUBLE PRECISION SIL,SIM,SLM,Y,Z,X,DAMP
      DOUBLE PRECISION SIL_ilm,SIM_ilm,SLM_ilm
      INTEGER MAPPED_LABELS(NEXTERNAL), MAPPED_FLAVOURS(NEXTERNAL)
      LOGICAL ISLOQCDPARTON(NEXTERNAL-1)
C     set logical doplot
      LOGICAL DOPLOT
      COMMON/CDOPLOT/DOPLOT
      DOUBLE PRECISION SCM
      COMMON/CSCM/SCM
      LOGICAL DOCUT
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
C     external
      INTEGER GET_COLOR_DIPOLE_INDEX
      EXTERNAL GET_COLOR_DIPOLE_INDEX
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      double precision  %(proc_prefix_S)s_GET_CCBLO
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_S)s_den
      common/%(proc_prefix_S)s_iden/%(proc_prefix_S)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      integer iseed
C
C     initialise
      M2_S_DIFF=0D0
      M2TMP=0D0
      IERR=0
      DAMP=0D0
      idum=0
C
C     possible cuts
      IF(DOCUT(XPB,NEXTERNAL-1,BORN_LEG_PDGS,0))RETURN
c
c     get PDGs
      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      CALL GET_SOFT_MAPPED_LABELS(I,idum,idum,NEXTERNAL,LEG_PDGS
     $,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
         WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),
     &   NEXTERNAL
         STOP
      ENDIF
C
C     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      PREF=-8D0*PI*ALPHAS
C
C     call colour-connected Born outside of the eikonal sum
      call %(proc_prefix_S)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
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
C         check LO color labels 
          IF(.NOT.(ISLOQCDPARTON(LB).AND.ISLOQCDPARTON(MB)))THEN
            WRITE(*,*)'Wrong LO indices in DIFF soft kernel',LB,MB
            STOP
          ENDIF
c
c         invariant quantities with mapping (ijr), at fixed Born kinematics
          SIL=XS(I,L)
          SIM=XS(I,M)
          SLM=XS(L,M)
c         invariant quantities with mapping (ilm), at fixed Born kinematics
          iA = 1  ! default azimuth for NLO
          CALL PHASE_SPACE_CS(XX,I,L,M,IA,XP_ILM,XPB,NEXTERNAL,
     $    LEG_PDGS,'S',XJCS_ILM)
          IF(XJCS_ILM.EQ.0D0)GOTO 999
          CALL INVARIANTS_FROM_P(XP_ilm,NEXTERNAL,XS_ilm,IERR)
          IF(IERR.EQ.1)GOTO 999
          SIL_ilm=XS_ilm(I,L)
          SIM_ilm=XS_ilm(I,M)
          SLM_ilm=XS_ilm(L,M)
C         
C         safety check
          IF(SIL*SIM.LE.0D0.or.SIL_ilm*SIM_ilm.LE.0D0)THEN
            WRITE(77,*)'Inaccuracy 1 in M2_S_DIFF'
            WRITE(77,*)SIL,SIM,SIL_ilm,SIM_ilm
            GOTO 999
          ENDIF
C         
C         eikonal difference
          ccBLO = %(proc_prefix_S)s_GET_CCBLO(lb,mb)
          M2TMP = SLM_ilm/(SIL_ilm*SIM_ilm) * XJB * XJCS_ILM
          M2TMP = M2TMP - SLM/(SIL*SIM) * XJ
          M2TMP = M2TMP * CCBLO*2D0
c
C         Including correct multiplicity factor
          M2tmp = M2tmp*dble(%(proc_prefix_S)s_den)/dble(%(proc_prefix_real)s_den)
c
c         Damping factors
          if(alpha.ne.0d0)then
             write(*,*)'Implement damping factors in M2_S_DIFF!'
             stop
          endif
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
          M2_S_DIFF=M2_S_DIFF+PREF*M2TMP*WSOFT*EXTRA
        ENDDO
      ENDDO
C     apply flavour factor
      M2_S_DIFF = M2_S_DIFF * %(proc_prefix_real)s_fl_factor
C
C     plot
      WGTPL=-M2_S_DIFF*WGT/NIT
      IF(DOPLOT)CALL HISTO_FILL(XPB,XSB,NEXTERNAL-1,WGTPL)
C     
C     sanity check
      IF(ABS(M2_S_DIFF).GE.HUGE(1D0).OR.ISNAN(M2_S_DIFF))THEN
        WRITE(77,*)'Exception caught in M2_S_DIFF',M2_S_DIFF
        GOTO 999
      ENDIF
C     
      RETURN
 999  IERR=1
      RETURN
      END


      double precision function M2_H_C_FgFg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
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
      double precision pref,M2tmp,wgt,wgtpl,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision %(proc_prefix_H_C_FgFg)s_GET_KKBLO
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
      integer %(proc_prefix_H_C_FgFg)s_den
      common/%(proc_prefix_H_C_FgFg)s_iden/%(proc_prefix_H_C_FgFg)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)

c
c     initialise
      M2_H_C_FgFg=0d0
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
         write(77,*)'Inaccuracy 1 in M2_H_C_FgFg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_H_C_FgFg)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(ia,ib,ir,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      KKBLO = %(proc_prefix_H_C_FgFg)s_GET_KKBLO(parent_leg,xpb,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=CA*2d0*(2d0/sab*KKBLO+x/(1d0-x)*(1d0-x**alpha)*BLO+(1d0-x)/x*(1d0-(1d0-x)**alpha)*BLO)
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_H_C_FgFg)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir) 
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_H_C_FgFg=M2tmp*pref/sab*extra
c     apply flavour factor
      M2_H_C_FgFg=M2_H_C_FgFg*%(proc_prefix_real)s_fl_factor
c
c     plot
c      wgtpl=-M2_H_C_FgFg*xj*wgt/nit/2d0/sCM
      wgtpl=-M2_H_C_FgFg*xj*wgt/nit
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C_FgFg).ge.huge(1d0).or.isnan(M2_H_C_FgFg))then
         write(77,*)'Exception caught in M2_H_C_FgFg',M2_H_C_FgFg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
  
                  
      double precision function M2_H_C_FgFq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
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
      double precision pref,M2tmp,wgt,wgtpl,xj,extra
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
      integer %(proc_prefix_H_C_FgFq)s_den
      common/%(proc_prefix_H_C_FgFq)s_iden/%(proc_prefix_H_C_FgFq)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)

c
c     initialise
      M2_H_C_FgFq=0d0
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
         write(77,*)'Inaccuracy 1 in M2_H_C_FgFq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_H_C_FgFq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c     In the following equation the x variable is related to the quark energy
      M2tmp=BLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0-x**alpha))
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_H_C_FgFq)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_H_C_FgFq=M2tmp*pref/sab*extra
c     apply flavour factor
      M2_H_C_FgFq=M2_H_C_FgFq*%(proc_prefix_real)s_fl_factor
c
c     plot
c      wgtpl=-M2_H_C_FgFq*xj*wgt/nit/2d0/sCM
      wgtpl=-M2_H_C_FgFq*xj*wgt/nit
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C_FgFq).ge.huge(1d0).or.isnan(M2_H_C_FgFq))then
         write(77,*)'Exception caught in M2_H_C_FgFq',M2_H_C_FgFq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_H_C_FqFqx(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
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
      double precision pref,M2tmp,wgt,wgtpl,xj,extra
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
      double precision %(proc_prefix_H_C_FqFqx)s_get_kkblo
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_H_C_FqFqx)s_den
      common/%(proc_prefix_H_C_FqFqx)s_iden/%(proc_prefix_H_C_FqFqx)s_den
      INTEGER ISEC,JSEC
      COMMON/CNLOSECINDICES/ISEC,JSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
c
c     initialise
      M2_H_C_FqFqx=0d0
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
         write(77,*)'Inaccuracy 1 in M2_H_C_FqFqx',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_H_C_FqFqx)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(ia,ib,ir,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      KKBLO = %(proc_prefix_H_C_FqFqx)s_GET_KKBLO(parent_leg,xpb,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=TR*(BLO-4d0/sab*KKBLO)
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_H_C_FqFqx)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_H_C_FqFqx=M2tmp*pref/sab*extra
c     apply flavour factor
      M2_H_C_FqFqx=M2_H_C_FqFqx*%(proc_prefix_real)s_fl_factor
c
c     plot
c      wgtpl=-M2_H_C_FqFqx*xj*wgt/nit/2d0/sCM
      wgtpl=-M2_H_C_FqFqx*xj*wgt/nit
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C_FqFqx).ge.huge(1d0).or.isnan(M2_H_C_FqFqx))then
         write(77,*)'Exception caught in M2_H_C_FqFqx',M2_H_C_FqFqx
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







