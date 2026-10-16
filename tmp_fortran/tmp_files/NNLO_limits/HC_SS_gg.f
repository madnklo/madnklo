
      double precision function M2_HC_SS_gg(i,j,r,xs,xp,xsb,xpb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j) S(i,j) *  WC_SS - S(i) C(i,j) S(i,j) * WS_C_SS
c     where  i, j are a g-g pair
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'mapped_labels_common.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,r
      integer l,m,ierr,nit,map1,map2
      integer jb,lb,mb,jbb,lbb,mbb
      double precision pref,M2_C_SS_gg,M2_SC_SS_gg,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS2,extra,damp
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1),xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),xpbb(0:3,nexternal-2)
      double precision sij,sir,sjr,zi,zj,sbml,sbjm,sbjl
      double precision kt(0:3),kt2,ktkm,ktkl,kmkl
      double precision Ei_jr,Ebj_ml,Pij,Qij
      double precision dot
      logical flavourmatch
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
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      M2_HC_SS_gg=0d0
      M2_C_SS_gg=0d0
      M2_SC_SS_gg=0d0
      damp=1d0
      ierr=0
c
c     check sector topology (only appears in ijjk)
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_HC_SS_gg',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = (leg_PDGs(i).eq.leg_PDGs(j)).and.(leg_PDGs(i).eq.21)
      if(.not.(flavourmatch))then
       write(*,*) 'Flavour mismatch in M2_HC_SS_gg', leg_PDGs(i),leg_PDGs(j)
       stop 1
      endif
c
c     fill mapped labels
      call fill_real_mapped_labels(i,j,leg_pdgs,real_leg_pdgs)
c     m=c & l=d in (76) of dropbox
      do m=1,nexternal
         if(.not.isNNLOQCDparton(m)) cycle
         if(m.eq.i.or.m.eq.j) cycle
         call fill_born_mapped_labels(real_labels(j,i,m),real_labels(m,i,m),real_leg_pdgs,born_leg_pdgs)
      enddo
c
c     parent leg
      jb = real_labels(j,i,j)
c
c     invariant quantities
      sij=xs(i,j)
      sir=xs(i,r)
      sjr=xs(j,r)
      zi=sir/(sir+sjr)
      zj=1d0-zi
c
      Ei_jr = sjr/sij/sir
      Pij = 2d0*CA*(zi/zj+zj/zi+zi*zj)
      Qij = -2d0*CA*zi*zj
      kt = zj*xp(:,i) - zi*xp(:,j) - (1d0-2d0*zi)*sij*xp(:,r)/(sir+sjr)
      kt2 = dot(kt,kt)
c
c     safety check
      if(sij.le.0d0.or.sir+sjr.le.0d0.or.zi.le.0d0.or.zi.ge.1d0)then
        write(77,*)'Inaccuracy 1 in M2_HC_SS_gg',sij,sir+sjr
        goto 999
      endif
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref = -(8d0*pi*alphas)**2
c
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
c
c     compute wc_nlo
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(i,j,ksec,r)
c
c     call wsbar
      call get_sig2(xsb,1d0,nexternal-1)
      map1=real_labels(csec,i,j)
      map2=real_labels(dsec,i,j)
      call get_wsbar_nlo(map1,map2)
c
c     Eikonal double sum starts here
c
      do m=1,nexternal
        if(.not.isnnloqcdparton(m))cycle
        if(m.eq.i.or.m.eq.j)cycle
        do l=1,nexternal
          if(.not.isnnloqcdparton(l))cycle
          if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle

          lb = real_labels(l,i,j)
          mb = real_labels(m,i,j)
          lbb = born_labels(lb,jb,mb)
          mbb = born_labels(mb,jb,mb)
c
          sbml = xsb(mb,lb)
          sbjm = xsb(jb,mb)
          sbjl = xsb(jb,lb)
          Ebj_ml = sbml/sbjm/sbjl
          ktkm = dot(kt(:),xpb(:,mb))
          ktkl = dot(kt(:),xpb(:,lb))
          kmkl = dot(xpb(:,mb),xpb(:,lb))
c
c         safety check
          if(sbjm.le.0d0.or.sbjm.le.0d0.or.kt2.eq.0d0) then
            write(77,*)'Inaccuracy 2 in M2_C_SS_gg',sbjm,sbjl,kt2
            goto 999
          endif
c
          call phase_space_cs_inv(jb,mb,lb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,mb))
          if(xjcs2.eq.0d0)goto 999
          if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c         call colour-connected born
          call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
          ccBLO = %(proc_prefix_Born)s_get_ccBLO(mbb,lbb)
c
c         collinear double-soft kernel, eq. (76) on dropbox
          M2_C_SS_gg = M2_C_SS_gg + (Pij/sij*Ebj_ml + Qij/sij*(+2d0/kt2*(ktkm/sbjm-ktkl/sbjl)**2+2d0*kmkl/sbjm/sbjl))*ccBLO * wc_nlo*wsbar_nlo
c
c         soft-collinear double-soft kernel, eq. (78) on dropbox
          M2_SC_SS_gg = M2_SC_SS_gg+2d0*CA*Ei_jr*Ebj_ml*ccBLO * wsbar_nlo
c
c         plot
          wgtpl = -(M2_C_SS_gg - M2_SC_SS_gg)
          wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
          wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
          wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
          wgts=wgtpl
c         if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs
c         ,wgtpl)
          if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
        enddo
      enddo
c
c     Double sum ends here
c
      M2_HC_SS_gg = M2_C_SS_gg - M2_SC_SS_gg
c
c     apply flavour factor
      M2_HC_SS_gg = M2_HC_SS_gg*pref*xj*damp*extra*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_HC_SS_gg = M2_HC_SS_gg * %(proc_prefix_rr)s_fl_factor
c
      if(test_sector_function) M2_HC_SS_gg = wc_nlo*wsbar_nlo-wsbar_nlo
c
      call ct_log('KC_SS                  ',M2_C_SS_gg*pref*xj*damp*extra*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)* %(proc_prefix_rr)s_fl_factor)
      call ct_log('KS_C_SS                ',M2_SC_SS_gg*pref*xj*damp*extra*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)* %(proc_prefix_rr)s_fl_factor)
c
c     sanity check
      if(abs(M2_HC_SS_gg).ge.huge(1d0).or.isnan(M2_HC_SS_gg))then
         write(77,*)'Exception caught in M2_HC_SS_gg',M2_HC_SS_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
