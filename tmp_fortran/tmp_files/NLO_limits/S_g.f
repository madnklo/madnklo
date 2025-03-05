
      
      double precision function M2_S_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     single-soft limit S_(i) * Zsoft
c     it returns 0 if i is not a gluon
      use sectors2_module
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjB,xjCS
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
      double precision alphaZ
      parameter(alphaZ=1d0)
      integer, parameter :: HEL = - 1
      double precision  %(proc_prefix_S_g)s_GET_CCBLO
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_S_g)s_den
      common/%(proc_prefix_S_g)s_iden/%(proc_prefix_S_g)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INCLUDE 'pmass.inc'
      
c
c     initialise
      M2_S_g=0d0
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
        WRITE(*,*) 'M2_S_g:'
        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
        STOP
      ENDIF
c
c     get PDGs
c      CALL GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,UNDERLYING_LEG_PDGS)

      CALL GET_SOFT_MAPPED_LABELS(I,NEXTERNAL,LEG_PDGS,MAPPED_LABELS,MAPPED_FLAVOURS,ISLOQCDPARTON)
c
c     call Z soft
      if(i.eq.isec) then
         CALL GET_ZS_NLO(ISEC,JSEC)
      elseif(i.eq.jsec) then
         CALL GET_ZS_NLO(JSEC,ISEC)
      else
         write(*,*)'In M2_S_g i should be = isec or = jsec',i,isec,jsec
         stop
      endif
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
            WRITE(*,*)'Wrong indices 1 in M2_S_g',LB,MB
            STOP
          ENDIF
          IF(leg_pdgs(l).ne.UnderLying_leg_pdgs(lb).or.leg_pdgs(m).ne.UnderLying_leg_pdgs(mb))THEN
            WRITE(*,*)'Wrong indices 2 in M2_S_g',L,M,LB,MB
            STOP
          ENDIF
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the single-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            IF(DOCUT(XPB,NEXTERNAL-1,UNDERLYING_LEG_PDGS,0))CYCLE
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
               write(77,*)'Inaccuracy 1 in M2_S_g',sil,sim
               goto 999
            endif
c
c     call colour-connected Born
            call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S_g)s_GET_CCBLO(lb,mb)
c
c     eikonal
            M2TMP=SLM/(SIL*SIM) - ML2/SIL**2 - MM2/SIM**2
            M2TMP = CCBLO*M2TMP
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_real)s_den)
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
            M2_S_g=M2_S_g+pref*M2tmp*ZS_NLO*extra
c
c     plot
            wgtpl=-pref*M2tmp*ZS_NLO*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_real)s_fl_factor
            if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
         enddo 
      enddo
c
c     apply flavour factor
      M2_S_g = M2_S_g * %(proc_prefix_real)s_fl_factor
c
c     sanity check
      if(abs(M2_S_g).ge.huge(1d0).or.isnan(M2_S_g))then
         write(77,*)'Exception caught in M2_S_g',M2_S_g
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

