
      
      double precision function M2_HC_SS_QQX(ia,ib,ik,ir,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     double-soft limit S_(i,j) * ZSS_NNLO
c     it returns 0 if i is not a gluon
      use sectors2_module
      use sectors4_module
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
      integer i,j,k,r
      integer ia,ib,ik,ir,l,m,ierr,nit,idum,parent,sec_index(2)
      integer jb,lb,mb
      integer jbb,lbb,mbb
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2), kt(0:3)
      double precision sij,sir,sjr
      double precision wa,wb,wr
      double precision sblm,sbjl,sbjm,ktkl,ktkm,kt2
      double precision x,y,xinit,damp
      double precision dot
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
      double precision   %(proc_prefix_Born)s_GET_CCBLO
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
c      integer (proc_prefix_S_g)s_den
c      common/(proc_prefix_S_g)s_iden/(proc_prefix_S_g)s_den
      integer %(proc_prefix_Born)s_den
      common/%(proc_prefix_Born)s_iden/%(proc_prefix_Born)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER REAL_LEG_PDGS(NEXTERNAL-1)
      INTEGER BORN_LEG_PDGS(NEXTERNAL-2)
      DOUBLE PRECISION PMASS(NEXTERNAL)
      integer mapped_sec(2,nexternal)
      integer ic,id
      INCLUDE 'pmass.inc'
c
c     initialise
      M2_HC_SS_QQX=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      idum=0
      wa = 0d0
      wb = 0d0
      wr = 0d0
      sij = 0d0
      sir = 0d0
      sjr = 0d0
      x   = 0d0
      y   = 0d0
      xinit = 0d0
      kt = 0d0
      kt2 = 0d0
      sblm = 0d0
      sbjl = 0d0
      sbjm = 0d0
      ktkl = 0d0
      ktkm = 0d0
      kt2 = 0d0
      ic = 0
      id = 0
c
c     return if not a qqb pair
      if((leg_pdgs(ia) + leg_pdgs(ib)).ne.0)return

c     initial checks and label assignment
      if(lsec.eq.0)then
         if((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)) then
            ic = ksec
         elseif((ia.eq.isec.and.ib.eq.ksec).or.(ia.eq.ksec.and.ib.eq.isec)) then
            ic = jsec
         elseif((ia.eq.jsec.and.ib.eq.ksec).or.(ia.eq.ksec.and.ib.eq.jsec)) then
            ic = isec
         else
            write(*,*)'Wrong indices 1 in M2_HC_qqx'
            write(*,*)ia,ib,isec,jsec,ksec
            stop
         endif
      else
         if((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)) then
            ic = ksec
            id = lsec
         elseif((ia.eq.ksec.and.ib.eq.lsec).or.(ia.eq.lsec.and.ib.eq.ksec)) then
            ic = isec
            id = jsec
         else
            write(*,*)'Wrong indices 2 in M2_HC_qqx'
            write(*,*)ia,ib,isec,jsec,ksec,lsec
            stop
         endif
      endif



      
c
c     safety check on PDGs
      IF(SIZE(LEG_PDGS).NE.NEXTERNAL)THEN
        WRITE(*,*) 'M2_HC_SS_QQX:'
        WRITE(*,*) 'Wrong dimension for leg_PDGs',SIZE(LEG_PDGS),NEXTERNAL
        STOP
      ENDIF
c
c     get PDGs
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,REAL_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)
      CALL GET_COLLINEAR_MAPPED_LABELS(ISEC,JSEC,NEXTERNAL,LEG_PDGS,NLO_MAPPED_LABELS,NLO_MAPPED_FLAVOURS)
      call reshuffle_momenta(nexternal,real_leg_pdgs,nlo_mapped_flavours,nlo_mapped_labels,xpb)

      JB = NLO_MAPPED_LABELS(j)
      PARENT = JB
      do l=1,nexternal
         if(l.eq.isec) cycle
          if(abs(NLO_mapped_flavours(l)).le.6.or.NLO_mapped_flavours(l).eq.21)isNLOmappedQCDparton(NLO_mapped_labels(l)) = .true.
      enddo
      CALL GET_COLLINEAR_MAPPED_LABELS(JB,NLO_MAPPED_LABELS(KSEC),NEXTERNAL-1,REAL_LEG_PDGS,LO_MAPPED_LABELS,LO_MAPPED_FLAVOURS)
      call reshuffle_momenta(nexternal-1,born_leg_pdgs,lo_mapped_flavours,lo_mapped_labels,xpbb)
      do l=1,nexternal-1
         if(l.eq.jb) cycle
          if(abs(LO_mapped_flavours(l)).le.6.or.LO_mapped_flavours(l).eq.21)isLOmappedQCDparton(LO_mapped_labels(l)) = .true.
       enddo
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
         write(77,*)'Inaccuracy 1 in M2_HC_SS_QQX',sab,sar+sbr,x
         goto 999
      endif
c
c
      call get_sig2(xsb,1d0,nexternal-1)
      if(lsec.eq.0)then
         sec_index(1) = parent
         sec_index(2) = nlo_mapped_labels(ic)
      else
         sec_index(1) = nlo_mapped_labels(ic)
         sec_index(2) = nlo_mapped_labels(id)
      endif
c     Fill  mapped_sec_list(2,nexternal) with the pairs of all the
c     final state particles after mapping n+2 --> n+1

      k = 0
      do i=3,nexternal-1
         do j=i+1,nexternal
            k=k+1
            mapped_sec(1,k) = nlo_mapped_labels(i)
            mapped_sec(2,k) = nlo_mapped_labels(j)
         enddo
      enddo

     call get_ZS_NNLO(sec_index(1),sec_index(2),mapped_sec)
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref = -32d0*pi**2*alphas**2
      call phase_space_CS_inv(ia,ib,ir,xp,xpb,nexternal,leg_PDGs,xjCS1)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(M))cycle
         if(m.eq.i.or.m.eq.j)cycle
         do l=1,nexternal
            if(.not.ISNNLOQCDPARTON(L))cycle
            if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c
            lb  = NLO_mapped_labels(l)
            mb  = NLO_mapped_labels(m)
            lbb = LO_mapped_labels(lb)
            mbb = LO_mapped_labels(mb)
c     
c         check labels and pdgs
            IF(.NOT.(ISNLOMAPPEDQCDPARTON(LB).AND.ISNLOMAPPEDQCDPARTON(MB)))THEN
               WRITE(*,*)'Wrong indices 1 in M2_HC_SS_QQX',LB,MB
               STOP
            ENDIF
            IF(.NOT.(ISLOMAPPEDQCDPARTON(LBB).AND.ISLOMAPPEDQCDPARTON(MBB)))THEN
               WRITE(*,*)'Wrong indices 2 in M2_HC_SS_QQX',LBB,MBB
               STOP
            ENDIF
          IF(leg_pdgs(l).ne.born_leg_pdgs(lbb).or.leg_pdgs(m).ne.born_leg_pdgs(mbb))THEN
            WRITE(*,*)'Wrong indices 3 in M2_HC_SS_QQX',L,M,LBB,MBB
            STOP
          ENDIF
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the double-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped

          call phase_space_CS_inv(jb,lb,mb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2)
          if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
          call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
          if(ierr.eq.1)goto 999
c
c     possible cuts
            IF(DOCUT(XPBB,NEXTERNAL-2,BORN_LEG_PDGS,0))CYCLE
c
c     invariant quantities
c     (c,d) in the paper --> (m,l)
            sblm = xsb(lb,mb)
            sbjl = xsb(jb,lb)
            sbjm = xsb(jb,mb)
            ktkl = dot(kt(:),xpb(:,lb))
            ktkm = dot(kt(:),xpb(:,mb))
            kt2=dot(kt(:),kt(:))
            
c
c     safety check
          IF(SIJ.LE.0D0.or.SBJL.le.0d0.or.SBJM.le.0d0.or.kt2.eq.0d0)THEN
            WRITE(77,*)'Inaccuracy 1 in M2_HC_SS_QQX',SIJ, SBJL, SBJM, KT2
            GOTO 999
          ENDIF
c
c     call colour-connected Born
c     TODO: fix strings for the associated underlying Born
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_Born)s_GET_CCBLO(lbb,mbb)
c
c     eikonal
c     See eq.1618 in file K2_I2_G_v2.pdf in the DropBox directory
c     (c,d) -> (m,l)
            M2tmp = TR*(sblm/(sbjl*sbjm)+4d0*x*(1d0-x)/kt2*(ktkl/sbjl-ktkm/sbjm)**2)
            M2TMP = CCBLO*M2TMP
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_HC_SS_QQX=M2_HC_SS_QQX+pref*M2tmp*ZS_NNLO*extra
c
c     plot
            wgtpl=-pref*M2tmp*ZSS_NNLO*extra*wgt/sab/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,wgtpl)
         enddo 
      enddo
c
c     apply flavour factor
      M2_HC_SS_QQX = M2_HC_SS_QQX * %(proc_prefix_rr)s_fl_factor
c
c     sanity check
      if(abs(M2_HC_SS_QQX).ge.huge(1d0).or.isnan(M2_HC_SS_QQX))then
         write(77,*)'Exception caught in M2_HC_SS_QQX',M2_HC_SS_QQX
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

