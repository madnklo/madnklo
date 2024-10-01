
      
      double precision function M2_SS_qqx(i,j,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     double-soft limit S_(i,j) * ZSS_NNLO
c     it returns 0 if i is not a gluon
      use sectors3_module
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
      integer i,j,l,m,ierr,nit,idum
      integer jb,lb,mb
      integer jbb,lbb,mbb
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sil,sim,slm,sij,sjl,sjm,ml2,mm2,y,z,x,damp
      integer NLO_mapped_labels(nexternal), NLO_mapped_flavours(NEXTERNAL)
      integer LO_mapped_labels(nexternal), LO_mapped_flavours(NEXTERNAL)
      logical isNLOmappedQCDparton(nexternal-1)
      logical isLOmappedQCDparton(nexternal-2)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_rr)s_fl_factor
      common/%(proc_prefix_rr)s_flavour_factor/%(proc_prefix_rr)s_fl_factor
c     external
      integer get_color_dipole_index
      external get_color_dipole_index
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
      double precision alphaZ
      parameter(alphaZ=2d0)
      integer, parameter :: HEL = - 1
      double precision  %(proc_prefix_S_g)s_GET_CCBLO
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      integer %(proc_prefix_S_g)s_den
      common/%(proc_prefix_S_g)s_iden/%(proc_prefix_S_g)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER REAL_LEG_PDGS(NEXTERNAL-1)
      INTEGER BORN_LEG_PDGS(NEXTERNAL-2)
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INCLUDE 'pmass.inc'
c
c     initialise
      M2_SS_qqx=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      idum=0
c
c     return if not a qqb pair
      if((leg_pdgs(i) + leg_pdgs(j)).ne.0)return
c
c     safety check on PDGs
      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
        WRITE(*,*) 'M2_SS_qqx:'
        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
        STOP
      ENDIF
c
c     get PDGs
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,REAL_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)
      CALL GET_COLLINEAR_MAPPED_LABELS(ISEC,JSEC,NEXTERNAL,LEG_PDGS,NLO_MAPPED_LABELS,NLO_MAPPED_FLAVOURS)
      JB = NLO_MAPPED_LABELS(J)
      do l=1,nexternal
         if(l.eq.isec) cycle
          if(abs(NLO_mapped_flavours(l)).le.6.or.NLO_mapped_flavours(l).eq.21)isNLOmappedQCDparton(NLO_mapped_labels(l)) = .true.
      enddo
      CALL GET_COLLINEAR_MAPPED_LABELS(JB,NLO_MAPPED_LABELS(KSEC),NEXTERNAL-1,REAL_LEG_PDGS,LO_MAPPED_LABELS,LO_MAPPED_FLAVOURS)
      do l=1,nexternal-1
         if(l.eq.jb) cycle
          if(abs(LO_mapped_flavours(l)).le.6.or.LO_mapped_flavours(l).eq.21)isLOmappedQCDparton(LO_mapped_labels(l)) = .true.
      enddo
c
c     call Z double-soft
      call get_sigNNLO(XS,alphaz,nexternal)
      CALL GET_ZSS_NNLO(I,KSEC,J,LSEC)
      if(ierr.eq.1)goto 999
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref=32d0*pi**2*alphas**2
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(M))cycle
         if(m.eq.i.or.m.eq.j)cycle
         do l=1,nexternal
            if(.not.ISNNLOQCDPARTON(L))cycle
            if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c
            lb = NLO_mapped_labels(l)
            mb = NLO_mapped_labels(m)
            lbb = LO_mapped_labels(lb)
            mbb = LO_mapped_labels(mb)
c     
c         check labels and pdgs
            IF(.NOT.(ISNLOMAPPEDQCDPARTON(LB).AND.ISNLOMAPPEDQCDPARTON(MB)))THEN
               WRITE(*,*)'Wrong indices 1 in M2_SS_qqx',LB,MB
               STOP
            ENDIF
            IF(.NOT.(ISLOMAPPEDQCDPARTON(LBB).AND.ISLOMAPPEDQCDPARTON(MBB)))THEN
               WRITE(*,*)'Wrong indices 2 in M2_SS_qqx',LBB,MBB
               STOP
            ENDIF
          IF(leg_pdgs(l).ne.born_leg_pdgs(lbb).or.leg_pdgs(m).ne.born_leg_pdgs(mbb))THEN
            WRITE(*,*)'Wrong indices 3 in M2_SS_qqx',L,M,LBB,MBB
            STOP
          ENDIF
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the double-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
          call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS1)
          call phase_space_CS_inv(jb,lb,mb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2)
          if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
          call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
          if(ierr.eq.1)goto 999
c
c     possible cuts
            IF(DOCUT(XPBB,NEXTERNAL-2,BORN_LEG_PDGS,0))CYCLE
c
c     invariant quantities
c            (c,d) in the paper --> (m,l)
            sij = xs(i,j)
            sil = xs(i,l)
            sim = xs(i,m)
            sjl = xs(j,l)
            sjm = xs(j,m)
            slm = xs(l,m)
c
c     safety check
          IF(SIJ.LE.0D0.or.(SIL+SJL).le.0d0.or.(SIM+SJM).le.0d0)THEN
            WRITE(77,*)'Inaccuracy 1 in M2_SS_qqx',SIJ, SIL+SJL, SIM+SJM
            GOTO 999
          ENDIF
c
c     call colour-connected Born
c     TODO: fix strings for the associated underlying Born
            call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S_g)s_GET_CCBLO(lbb,mbb)
c
c     eikonal
c     See file K2_I2_G_v2.pdf in the DropBox directory
c     (c,d) -> (m,l)      
            M2tmp = 2d0*TR*(((sil*sjm-sim*sjl)**2-slm*sij*(sil+sjl)*(sim+sjm))/(sij**2*(sil+sjl)**2*(sim+sjm)**2))
            M2TMP = CCBLO*M2TMP
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_SS_qqx=M2_SS_qqx+pref*M2tmp*ZSS_NNLO*extra
c
c     plot
            wgtpl=-pref*M2tmp*ZSS_NNLO*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,wgtpl)
         enddo 
      enddo
c
c     apply flavour factor
      M2_SS_qqx = M2_SS_qqx * %(proc_prefix_rr)s_fl_factor
c
c     sanity check
      if(abs(M2_SS_qqx).ge.huge(1d0).or.isnan(M2_SS_qqx))then
         write(77,*)'Exception caught in M2_SS_qqx',M2_SS_qqx
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

