      double precision function M2_C_SS_QQX(ia,ib,ir,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j) S(i,j) kernel times WC_SS: i, j are a q-qb pair
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
      integer ia,ib,ik,ir,l,m,ierr,nit,parent,map1,map2
      integer jb,lb,mb
      integer jbb,lbb,mbb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal)
      double precision xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO,extra,damp
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2),kt(0:3),kt2
      double precision sab,sar,sbr,wa,wb,wr,x
      double precision sblm,sbjl,sbjm,Ebjlm,ktkl,ktkm
      double precision dot
      logical flavourmatch
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
      double precision pij,qij,sij,zi,zj
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
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer asec,bsec,csec,dsec
      common/csecindices/asec,bsec,csec,dsec
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      integer mapped_sec(2,nexternal)
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      M2_C_SS_QQX=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology (only appears in ijjk)
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_C_SS_qqx',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = (leg_PDGs(ia).eq.-leg_PDGs(ib)).and.(abs(leg_PDGs(ia)).le.5)
      if(.not.(flavourmatch))then
       write(*,*) 'Flavour mismatch in M2_C_SS_qqx', leg_PDGs(ia),leg_PDGs(ib)
       stop 1
      endif
c
      call phase_space_cs_inv(ia,ib,ir,xp,xpb,nexternal,leg_pdgs,xjcs1,real_mapped_labels)
      if(xjcs1.eq.0d0)goto 999
c
c     parent leg
      jb = real_mapped_labels(ib)
      parent = jb
      do l=1,nexternal
         if(l.eq.isec) cycle
          if(abs(leg_pdgs(l)).le.6.or.leg_pdgs(l).eq.21) isNLOmappedQCDparton(real_mapped_labels(l)) = .true.
      enddo
      do l=1,nexternal-1
         if(l.eq.jb) cycle
          if(abs(real_leg_pdgs(l)).le.6.or.real_leg_pdgs(l).eq.21) isLOmappedQCDparton(Born_mapped_labels(l)) = .true.
       enddo
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sar/(sar+sbr)
c
c     coefficients of kt
c     kt = wa pa + wb pb + wr pr
      wa = 1d0-x
      wb = -x
      wr = -(1d0-2d0*x)*sab/(sar+sbr)
      kt(:) = wa*xp(:,ia) + wb*xp(:,ib) + wr*xp(:,ir)
      kt2 = dot(kt(:),kt(:))
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_C_SS_QQX',sab,sar+sbr,x
         goto 999
      endif
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref = -64d0*pi**2*alphas**2
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
c
c     Eikonal double sum starts here
c
      do mb=1,nexternal-1
         if(.not.ISNLOMAPPEDQCDPARTON(MB))cycle
         if(mb.eq.jb)cycle
         do lb=1,nexternal-1
            if(.not.ISNLOMAPPEDQCDPARTON(LB))cycle
            if(lb.eq.jb.or.lb.eq.mb)cycle
            lbb = Born_mapped_labels(lb)
            mbb = Born_mapped_labels(mb)
c
c        check labels and pdgs
            IF(.NOT.(ISLOMAPPEDQCDPARTON(LBB).AND.ISLOMAPPEDQCDPARTON(MBB)))THEN
               WRITE(*,*)'Wrong indices 1 in M2_C_SS_QQX',LBB,MBB
               STOP
            ENDIF
            IF(real_leg_pdgs(lb).ne.Born_leg_pdgs(lbb).or.real_leg_pdgs(mb).ne.Born_leg_pdgs(mbb)) THEN
               WRITE(*,*)'Wrong indices 2 in M2_C_SS_QQX',LB,MB,LBB,MBB
               STOP
            ENDIF
c
c        phase-space mapping according to lb and mb, at fixed radiation
c        phase-space point: the singular kernel is in the same point
c        as the double-real, ensuring numerical stability, while the
c        underlying Born configuration is remapped
            call phase_space_CS_inv(jb,lb,mb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            if(docut(xpbb,nexternal-2,born_leg_pdgs,2))cycle
c
c     invariant quantities
c     (c,d) in the paper --> (m,l)
            sblm = xsb(lb,mb)
            sbjl = xsb(jb,lb)
            sbjm = xsb(jb,mb)
            ktkl = dot(kt(:),xpb(:,lb))
            ktkm = dot(kt(:),xpb(:,mb))
c
c     safety check
            IF(SAB.LE.0D0.or.SBJL.le.0d0.or.SBJM.le.0d0.or.kt2.eq.0d0)THEN
               WRITE(77,*)'Inaccuracy 2 in M2_C_SS_QQX',SAB, SBJL, SBJM, KT2
               GOTO 999
            ENDIF
c
c     call colour-connected Born
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_Born)s_GET_CCBLO(lbb,mbb)
c
c     collinear double-soft kernel, eq. (C.36) of 2212.11190v2
            Ebjlm = sblm/sbjl/sbjm
            M2tmp = TR*(Ebjlm+4d0*x*(1d0-x)/kt2*(ktkl/sbjl-ktkm/sbjm)**2)
            m2tmp = m2tmp/sab*ccBLO
c     Include collinear double-soft sector functions, eq. (C.80) of 2212.11190v2
            call get_sig2(xs,alphaz,nexternal)
            call get_wc_nlo(ia,ib,ksec,ir)
            call get_sig2(xsb,1d0,nexternal-1)
            map1=real_mapped_labels(csec)
            map2=real_mapped_labels(dsec)
            call get_wsbar_nlo(map1,map2)
            M2TMP=M2TMP*wc_nlo*wsbar_nlo
c
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_C_SS_QQX=M2_C_SS_QQX+pref*M2tmp*extra
c
c     plot
            wgtpl=-pref*M2tmp*extra*wgt/nit*wgt_chan
            wgtpl=wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo
c
c     Double sum ends here
c
c     apply flavour factor
      M2_C_SS_QQX = M2_C_SS_QQX * %(proc_prefix_rr)s_fl_factor
c
c     sanity check
      if(abs(M2_C_SS_QQX).ge.huge(1d0).or.isnan(M2_C_SS_QQX))then
         write(77,*)'Exception caught in M2_C_SS_QQX',M2_C_SS_QQX
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
