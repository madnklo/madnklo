
      
      double precision function M2_SS_qq(i,j,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sil,sim,slm,ml2,mm2,y,z,x,damp
      integer NLO_mapped_labels(nexternal), NLO_mapped_flavours(NEXTERNAL)
      integer LO_mapped_labels(nexternal), LO_mapped_flavours(NEXTERNAL)
      logical isNLOmappedQCDparton(nexternal-1)
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
      parameter(alphaZ=1d0)
      integer, parameter :: HEL = - 1
      double precision  %(proc_prefix_S_g)s_GET_CCBLO
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      integer %(proc_prefix_S_g)s_den
      common/%(proc_prefix_S_g)s_iden/%(proc_prefix_S_g)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INCLUDE 'pmass.inc'
      
c
c     initialise
      M2_SS_qq=0d0
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
        WRITE(*,*) 'M2_SS_qq:'
        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
        STOP
      ENDIF
c
c     get PDGs
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,REAL_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)

      CALL GET_SOFT_MAPPED_LABELS(I,NEXTERNAL,LEG_PDGS,NLO_MAPPED_LABELS,NLO_MAPPED_FLAVOURS,ISNLOMAPPEDQCDPARTON)
      jb = NLO_mapped_labels(j)
      CALL GET_SOFT_MAPPED_LABELS(JB,NEXTERNAL-1,REAL_LEG_PDGS,LO_MAPPED_LABELS,LO_MAPPED_FLAVOURS,ISLOMAPPEDQCDPARTON)


c
c     call Z double-soft
      call get_ZSS_NNLO(i,ksec,j,lsec)
c      else
c         write(*,*)'In M2_S_g i should be = isec or = jsec',i,isec,jsec
c         stop
c      endif
      if(ierr.eq.1)goto 999
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref=32d0*pi**2*alphas**2
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.ismappedQCDparton(m))cycle
         if(m.eq.i.or.m.eq.j)cycle
         do l=1,nexternal
            if(.not.ismappedQCDparton(l))cycle
            if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c     
            
            lb = NLO_mapped_labels(l)
            mb = NLO_mapped_labels(m)
            lbb = LO_mapped_labels(lb)
            mbb = LO_mapped_labels(mb)
c     
c         check labels and pdgs
            IF(.NOT.(ISNLOMAPPEDQCDPARTON(LB).AND.ISNLOMAPPEDQCDPARTON(MB)))THEN
               WRITE(*,*)'Wrong indices 1 in M2_SS_qq',LB,MB
               STOP
            ENDIF
         
            IF(.NOT.(ISLOMAPPEDQCDPARTON(LBB).AND.ISLOMAPPEDQCDPARTON(MBB)))THEN

               WRITE(*,*)'Wrong indices 2 in M2_SS_qq',LBB,MBB
               STOP
            ENDIF
          IF(leg_pdgs(l).ne.UnderLying_leg_pdgs(lb).or.leg_pdgs(m).ne.UnderLying_leg_pdgs(mb))THEN
            WRITE(*,*)'Wrong indices 2 in M2_SS_qq',L,M,LB,MB
            STOP
          ENDIF
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the single-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
          call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS1)
          call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS2)
          if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
          call invariants_from_p(xpb,nexternal-1,xsb,ierr)
          if(ierr.eq.1)goto 999
c
c     possible cuts
            IF(DOCUT(XPBB,NEXTERNAL-2,BORN_LEG_PDGS,0))CYCLE
c
c     invariant quantities
c            (c,d) in the paper --> (m,l)

            sij = xs(i,j)
            sim = xs(i,m)
            sjl = xs(j,l)
            sil = xs(i,l)
            sjm = xs(j,m)
c
c     safety check
            if(sil*sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S_g',sil,sim
               goto 999
            endif
c
c     call colour-connected Born
            call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S_g)s_GET_CCBLO(lbb,mbb)
c
c     eikonal

c            See file K2_I2_G_v2.pdf in the DropBox directory
c            (c,d) -> (m,l)
            
            M2tmp = 2d0*TR*(((sil*sjm-sim*sjl)**2-slm*sij*(sil+sjl)*(sim+sjm))/(sij**2*(sil+sjl)**2*(sim+sjm)**2))

            M2TMP = CCBLO*M2TMP
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_SS_g=M2_SS_g+pref*M2tmp*ZSS_NNLO*extra
c
c     plot
            wgtpl=-pref*M2tmp*ZS_NLO*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
         enddo 
      enddo
c
c     apply flavour factor
      M2_S_g = M2_S_g * %(proc_prefix_rr)s_fl_factor
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

