      double precision function M2_HC_SC_gluon(k,i,j,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     Ordered-sector C(i,j) SC(k,i,j) - S(i) C(i,j) SC(k,i,j).
c     Eqs. (C.33),(C.34),(C.78),(C.79) of 2212.11190v2. The S(j)
c     overlap vanishes with this ordered sector, not with the matrix element.
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
      integer k,i,j,r,nit,ierr,c,d,rp,ib,jb,kb,cb,rb
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
      double precision xpb(0:3,nexternal-1),xsb(nexternal-1,nexternal-1)
      double precision xpbb(0:3,nexternal-2),xsbb(nexternal-2,nexternal-2)
      double precision pmass(nexternal),kt(0:3),kt2,sij,sir,sjr,zi,zj,Pij,Qij,SiCij
      double precision wgt,xj,xjb,extra,wgt_chan,wc,ws,xjcs1,xjcs2
      double precision alphas,ccBLO,kkccBLO,ans(0:NSQSO_BORN)
      double precision collinear,soft_collinear,wgts(1)
      double precision %(proc_prefix_Born)s_GET_CCBLO
## if (pyamplicol) {
      double precision %(proc_prefix_Born)s_GET_KKCCBLO
## } else {
      double precision %(proc_prefix_Born)s_GET_KKBLO,unused
## }
      logical gg,qqx,docut,doplot,test_sector_function
      common/cdoplot/doplot
      common/ctestsecfun/test_sector_function
      include 'pmass.inc'

      M2_HC_SC_gluon=0d0
      collinear=0d0
      soft_collinear=0d0
      ierr=0
      if(min(i,j,k,r).le.nincoming.or.max(i,j,k,r).gt.nexternal)then
        write(*,*) 'HC_SC: requires four distinct final-state legs'
        stop 1
      endif
      if(i.eq.j.or.i.eq.k.or.i.eq.r.or.j.eq.k.or.j.eq.r.or.k.eq.r)then
        write(*,*) 'HC_SC: requires four distinct final-state legs'
        stop 1
      endif
      gg=leg_pdgs(i).eq.21.and.leg_pdgs(j).eq.21
      qqx=abs(leg_pdgs(i)).ge.1.and.abs(leg_pdgs(i)).le.6.and.leg_pdgs(i).eq.-leg_pdgs(j)
      if((.not.gg.and..not.qqx).or.leg_pdgs(k).ne.21.or..not.isnnloqcdparton(r))then
        write(*,*) 'HC_SC: requires a gluon parent, soft gluon and coloured recoiler'
        stop 1
      endif
c     Match the symmetric three-leg recoiler rule in get_recoiler.
      rp=0
      do c=1,nexternal
        if(.not.isnnloqcdparton(c))cycle
        if(c.le.nincoming.or.pmass(c).ne.0d0)then
          write(*,*) 'HC_SC: only massless final-state QCD is implemented'
          stop 1
        endif
        if(c.eq.i.or.c.eq.j.or.c.eq.k)cycle
        if(rp.eq.0)then
          rp=c
        elseif(leg_pdgs(c).lt.leg_pdgs(rp))then
          rp=c
        endif
      enddo
      sij=xs(i,j)
      sir=xs(i,r)
      sjr=xs(j,r)
      if(min(sij,sir,sjr).le.0d0)goto 999
      zi=sir/(sir+sjr)
      zj=sjr/(sir+sjr)
      kt2=-zi*zj*sij
      if(min(zi,zj).le.0d0.or.kt2.eq.0d0)goto 999
      kt=zj*xp(:,i)-zi*xp(:,j)-(zj-zi)*sij/(sir+sjr)*xp(:,r)
      SiCij=0d0
      if(gg)then
        Pij=2d0*CA*(zi/zj+zj/zi+zi*zj)
        Qij=-2d0*CA*zi*zj
        SiCij=2d0*CA*sjr/sij/sir
      else
        Pij=TR*(1d0-2d0*zi*zj)
        Qij=2d0*TR*zi*zj
      endif
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(i,j,k,r)
      wc=wc_nlo
      call fill_real_mapped_labels(i,j,leg_pdgs,real_leg_pdgs)
      ib=real_labels(i,i,j)
      jb=real_labels(j,i,j)
      kb=real_labels(k,i,j)
      rb=real_labels(rp,i,j)
      call phase_space_CS_inv(i,j,r,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,i,j))
      if(xjcs1.eq.0d0)goto 999
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.ne.0)goto 999
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      call get_wsbar_nlo(real_labels(csec,i,j),real_labels(dsec,i,j))
      ws=wsbar_nlo
      if(test_sector_function)then
        M2_HC_SC_gluon=wc*ws
        if(gg)M2_HC_SC_gluon=M2_HC_SC_gluon-ws
        return
      endif
      call %(proc_prefix_Born)s_RESET_SPIN_CORRELATION_VECTORS()
      do c=1,nexternal
        if(.not.isnnloqcdparton(c))cycle
        if(c.eq.i.or.c.eq.j.or.c.eq.k)cycle
        cb=real_labels(c,i,j)
        if(c.ne.rp)then
          do d=1,nexternal
            if(.not.isnnloqcdparton(d))cycle
            if(d.eq.i.or.d.eq.j.or.d.eq.k.or.d.eq.rp.or.d.eq.c)cycle
            call add_pair(cb,real_labels(d,i,j),1d0)
            if(ierr.ne.0)goto 999
          enddo
          call add_pair(cb,rb,2d0)
          if(ierr.ne.0)goto 999
        endif
        call add_pair(cb,jb,2d0)
        if(ierr.ne.0)goto 999
      enddo
      M2_HC_SC_gluon=collinear-soft_collinear
      if(abs(M2_HC_SC_gluon).ge.huge(1d0).or.isnan(M2_HC_SC_gluon))goto 999
      call ct_log('KC_SC                  ',collinear)
      call ct_log('KS_C_SC                ',soft_collinear)
      return
 999  ierr=1
      M2_HC_SC_gluon=0d0
      return

      contains

      subroutine add_pair(a,b,multiplicity)
      integer a,b,ab,bb,parent
      double precision multiplicity,eik,contraction,pref,cterm,sterm
      if(xsb(kb,a).le.0d0.or.xsb(kb,b).le.0d0)then
        ierr=1
        return
      endif
      call fill_born_mapped_labels(kb,a,real_leg_pdgs,Born_leg_pdgs)
      call phase_space_CS_inv(kb,a,b,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,kb,a))
      if(xjcs2.eq.0d0)then
        ierr=1
        return
      endif
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.ne.0)return
      ab=born_labels(a,kb,a)
      bb=born_labels(b,kb,a)
      parent=born_labels(jb,kb,a)
      call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
      alphas=alphas_current
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,-1,alphas,ans)
      ccBLO=%(proc_prefix_Born)s_GET_CCBLO(ab,bb)
## if (pyamplicol) {
      kkccBLO=%(proc_prefix_Born)s_GET_KKCCBLO(xpbb,alphas,parent,kt,ab,bb)
## } else {
      unused=%(proc_prefix_Born)s_GET_KKBLO(parent,xpbb,kt)
      kkccBLO=%(proc_prefix_Born)s_GET_CCBLO(ab,bb)
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,-1,alphas,ans)
## }
      contraction=Pij*ccBLO+Qij*(ccBLO+2d0*kkccBLO/kt2)
      eik=xsb(a,b)/xsb(kb,a)/xsb(kb,b)
      pref=-(8d0*pi*alphas)**2*multiplicity*eik*extra*xj*%(proc_prefix_rr)s_fl_factor
      pref=pref*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      cterm=pref*contraction/sij*wc*ws
c     The additional soft limit is spin independent; do not use KKCCBLO.
      sterm=pref*SiCij*ccBLO*ws
      collinear=collinear+cterm
      soft_collinear=soft_collinear+sterm
      if(doplot)then
        wgts(1)=-(cterm-sterm)*wgt/nit*wgt_chan
        call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      endif
      end subroutine add_pair
      end function M2_HC_SC_gluon
