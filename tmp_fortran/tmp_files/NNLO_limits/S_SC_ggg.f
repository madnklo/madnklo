      double precision function M2_S_SC_ggg(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) SC(i,j,k), gluon pair: Eq. (C.27) of 2212.11190v2.
c     Every splitting tensor and collinear sector weight uses the first
c     soft-mapped point, not the original double-real momenta.
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
      integer i,j,k,r,nit,ierr,c,d
      integer isec,jsec,ksec,lsec,iref,asec,bsec,csec,dsec
      common/cpartindices/isec,jsec,ksec,lsec,iref
      common/csecindices/asec,bsec,csec,dsec
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer %(proc_prefix_rr)s_fl_factor,%(proc_prefix_rr)s_den,%(proc_prefix_Born)s_den
      common/%(proc_prefix_rr)s_flavour_factor/%(proc_prefix_rr)s_fl_factor
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      common/%(proc_prefix_Born)s_iden/%(proc_prefix_Born)s_den
      double precision xs(nexternal,nexternal),xp(0:3,nexternal)
      double precision pmass(nexternal),wgt,xj,xjb,extra,wgt_chan,ws,total
      logical docut,doplot,test_sector_function
      common/cdoplot/doplot
      common/ctestsecfun/test_sector_function
      include 'pmass.inc'

      M2_S_SC_ggg=0d0
      total=0d0
      ierr=0
      if(min(i,j,k,r).le.nincoming.or.max(i,j,k,r).gt.nexternal)then
        write(*,*) 'S_SC_ggg requires four distinct final-state legs'
        stop 1
      endif
      if(i.eq.j.or.i.eq.k.or.i.eq.r.or.j.eq.k.or.j.eq.r.or.k.eq.r)then
        write(*,*) 'S_SC_ggg requires four distinct final-state legs'
        stop 1
      endif
      if(leg_pdgs(i).ne.21.or.leg_pdgs(j).ne.21.or.leg_pdgs(k).ne.21.or..not.isnnloqcdparton(r))then
        write(*,*) 'S_SC_ggg requires three gluons and a coloured recoiler'
        stop 1
      endif
      if(i.ne.asec)then
        write(*,*) 'S_SC_ggg: soft leg does not match the ordered sector'
        stop 1
      endif
      do c=1,nexternal
        if(.not.isnnloqcdparton(c))cycle
        if(c.le.nincoming.or.pmass(c).ne.0d0)then
          write(*,*) 'S_SC_ggg: only massless final-state QCD is implemented'
          stop 1
        endif
      enddo
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
      ws=ws_nlo
      call %(proc_prefix_Born)s_RESET_SPIN_CORRELATION_VECTORS()
      do c=1,nexternal
        if(.not.isnnloqcdparton(c))cycle
        if(c.eq.i.or.c.eq.j.or.c.eq.k)cycle
c       Spectator pairs, with the canonical second map (j,k,r).
        do d=1,nexternal
          if(.not.isnnloqcdparton(d))cycle
          if(d.eq.i.or.d.eq.j.or.d.eq.k.or.d.eq.c)cycle
          call add_term(c,d,1,1d0)
          if(ierr.ne.0)goto 999
        enddo
c       Each parent-spectator eikonal has both soft-map orientations.
c       rho^(C)_jk = rho^(C)_kj = 1, and the d_ABC term vanishes for gg.
        call add_term(j,c,2,.5d0)
        if(ierr.ne.0)goto 999
        call add_term(c,j,2,.5d0)
        if(ierr.ne.0)goto 999
        call add_term(k,c,3,.5d0)
        if(ierr.ne.0)goto 999
        call add_term(c,k,3,.5d0)
        if(ierr.ne.0)goto 999
      enddo
      M2_S_SC_ggg=total
      if(abs(total).ge.huge(1d0).or.isnan(total))goto 999
      call ct_log('KS_SC                  ',M2_S_SC_ggg)
      return
 999  ierr=1
      M2_S_SC_ggg=0d0
      return

      contains

      subroutine add_term(a,b,mapping,multiplicity)
      integer a,b,mapping,ab,bb,jb,kb,rb,parent,cb,db
      double precision multiplicity,xjcs1,xjcs2,zj,zk,sjk,sjr,skr,Pjk,Qjk,kt2
      double precision xpb(0:3,nexternal-1),xsb(nexternal-1,nexternal-1)
      double precision xpbb(0:3,nexternal-2),xsbb(nexternal-2,nexternal-2),kt(0:3)
      double precision alphas,ccBLO,kkccBLO,ans(0:NSQSO_BORN),eik,contraction,term,wc,wgts(1)
      double precision %(proc_prefix_Born)s_GET_CCBLO
## if (pyamplicol) {
      double precision %(proc_prefix_Born)s_GET_KKCCBLO
## } else {
      double precision %(proc_prefix_Born)s_GET_KKBLO,unused
## }
      if(xs(i,a).le.0d0.or.xs(i,b).le.0d0)then
        ierr=1
        return
      endif
      call fill_real_mapped_labels(i,a,leg_pdgs,real_leg_pdgs)
      ab=real_labels(a,i,a)
      bb=real_labels(b,i,a)
      jb=real_labels(j,i,a)
      kb=real_labels(k,i,a)
      rb=real_labels(r,i,a)
      call phase_space_CS_inv(i,a,b,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,i,a))
      if(xjcs1.eq.0d0)then
        ierr=1
        return
      endif
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.ne.0)return
      sjk=xsb(jb,kb)
      sjr=xsb(jb,rb)
      skr=xsb(kb,rb)
      if(min(sjk,sjr,skr).le.0d0)then
        ierr=1
        return
      endif
      zj=sjr/(sjr+skr)
      zk=skr/(sjr+skr)
      kt2=-zj*zk*sjk
      if(min(zj,zk).le.0d0.or.kt2.eq.0d0)then
        ierr=1
        return
      endif
      Pjk=2d0*CA*(zj/zk+zk/zj+zj*zk)
      Qjk=-2d0*CA*zj*zk
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      call get_wcbar_nlo(real_labels(csec,i,a),real_labels(dsec,i,a),rb)
      wc=wcbar_nlo
      if(test_sector_function)then
        total=ws*wc
        return
      endif
      call fill_born_mapped_labels(jb,kb,real_leg_pdgs,Born_leg_pdgs)
      if(mapping.eq.1)then
        call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,kb))
        kt=zj*xpb(:,kb)-zk*xpb(:,jb)-(zj-zk)*sjk/(sjr+skr)*xpb(:,rb)
      elseif(mapping.eq.2)then
        call phase_space_CS_inv(kb,rb,jb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,kb))
c       Same finite transverse-plane prescription as SC_ggg.
        kt=zj*xpb(:,kb)+zk*(sjk-sjr)/(sjk+sjr)*xpb(:,jb)-sjk/(sjr+skr)*xpb(:,rb)
      else
        call phase_space_CS_inv(jb,rb,kb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,kb))
        kt=-zk*xpb(:,jb)+zj*(skr-sjk)/(sjk+skr)*xpb(:,kb)+sjk/(sjr+skr)*xpb(:,rb)
      endif
      if(xjcs2.eq.0d0)then
        ierr=1
        return
      endif
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.ne.0)return
      parent=born_labels(kb,jb,kb)
      cb=born_labels(ab,jb,kb)
      db=born_labels(bb,jb,kb)
      call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
      alphas=alphas_current
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,-1,alphas,ans)
      ccBLO=%(proc_prefix_Born)s_GET_CCBLO(cb,db)
## if (pyamplicol) {
      kkccBLO=%(proc_prefix_Born)s_GET_KKCCBLO(xpbb,alphas,parent,kt,cb,db)
## } else {
      unused=%(proc_prefix_Born)s_GET_KKBLO(parent,xpbb,kt)
      kkccBLO=%(proc_prefix_Born)s_GET_CCBLO(cb,db)
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,-1,alphas,ans)
## }
      contraction=Pjk*ccBLO+Qjk*(ccBLO+2d0*kkccBLO/kt2)
      eik=xs(a,b)/xs(i,a)/xs(i,b)
      term=-(8d0*pi*alphas)**2*multiplicity*eik*contraction/sjk*ws*wc
      term=term*extra*xj*%(proc_prefix_rr)s_fl_factor*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      total=total+term
      if(doplot)then
        wgts(1)=-term*wgt/nit*wgt_chan
        call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      endif
      end subroutine add_term
      end function M2_S_SC_ggg
