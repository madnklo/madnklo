      double precision function M2_SC_ggg(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) C(j,k), with i,j,k gluons: Eq. (C.13) of 2212.11190v2.
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'colored_partons.inc'
      include 'mapped_labels_common.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,r,nit,ierr,c,d,mapping,ib,kb,parent
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
      double precision xpb(0:3,nexternal-1),xpbb(0:3,nexternal-2)
      double precision xsbb(nexternal-2,nexternal-2),pmass(nexternal)
      double precision sjk,sjr,skr,zj,zk,Pjk,Qjk,kt(0:3),kt2
      double precision wgt,wgts(1),wgt_chan,xj,xjb,extra,wsc_nnlo,total
      double precision xjcs1,xjcs2,alphas,ccBLO,kkccBLO,ans(0:NSQSO_BORN)
      double precision %(proc_prefix_Born)s_GET_CCBLO
## if (pyamplicol) {
      double precision %(proc_prefix_Born)s_GET_KKCCBLO
## } else {
      double precision %(proc_prefix_Born)s_GET_KKBLO,unused
## }
      logical docut,doplot,test_sector_function
      common/cdoplot/doplot
      common/ctestsecfun/test_sector_function
      include 'pmass.inc'

      M2_SC_ggg=0d0
      total=0d0
      ierr=0
      if(min(i,j,k,r).le.nincoming.or.max(i,j,k,r).gt.nexternal)then
        write(*,*) 'SC_ggg requires four distinct final-state legs'
        stop 1
      endif
      if(i.eq.j.or.i.eq.k.or.i.eq.r.or.j.eq.k.or.j.eq.r.or.k.eq.r)then
        write(*,*) 'SC_ggg requires four distinct final-state legs'
        stop 1
      endif
      if(leg_pdgs(i).ne.21.or.leg_pdgs(j).ne.21.or.leg_pdgs(k).ne.21)then
        write(*,*) 'SC_ggg requires three gluons'
        stop 1
      endif
      if(.not.isnnloqcdparton(r))then
        write(*,*) 'SC_ggg requires a coloured recoiler'
        stop 1
      endif
c     The existing CS inverse map is final-final; do not cross it silently.
      do c=1,nexternal
        if(.not.isnnloqcdparton(c))cycle
        if(c.le.nincoming)then
          write(*,*) 'SC_ggg: coloured initial states are not implemented'
          stop 1
        endif
        if(pmass(c).ne.0d0)then
          write(*,*) 'SC_ggg: massive coloured legs are not implemented'
          stop 1
        endif
      enddo
      sjk=xs(j,k)
      sjr=xs(j,r)
      skr=xs(k,r)
      if(min(sjk,sjr,skr).le.0d0)goto 999
      zj=sjr/(sjr+skr)
      zk=skr/(sjr+skr)
      kt2=-zj*zk*sjk
      if(min(zj,zk).le.0d0.or.kt2.eq.0d0)goto 999
      Pjk=2d0*CA*(zj/zk+zk/zj+zj*zk)
      Qjk=-2d0*CA*zj*zk

      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      if(i.eq.asec)then
        call get_wsc_ikl_nnlo(asec,bsec,csec,dsec,iref)
        wsc_nnlo=wsc_ikl_nnlo
      elseif(i.eq.csec)then
        call get_wsc_kij_nnlo(asec,bsec,csec,dsec,iref)
        wsc_nnlo=wsc_kij_nnlo
      else
        write(*,*) 'SC_ggg: soft leg does not match the sector'
        stop 1
      endif
      if(test_sector_function)then
        M2_SC_ggg=wsc_nnlo
        return
      endif
      call fill_real_mapped_labels(j,k,leg_pdgs,real_leg_pdgs)
      ib=real_labels(i,j,k)
      kb=real_labels(k,j,k)
      parent=kb
      call %(proc_prefix_Born)s_RESET_SPIN_CORRELATION_VECTORS()

c     All three maps use the same merged (j,k) flavour slot. In the
c     alternate maps this slot contains the rescaled spectator j or k.
      do mapping=1,3
        if(mapping.eq.1)then
          call phase_space_CS_inv(j,k,r,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,j,k))
          kt=zj*xp(:,k)-zk*xp(:,j)-(zj-zk)*sjk/(sjr+skr)*xp(:,r)
        elseif(mapping.eq.2)then
          call phase_space_CS_inv(k,r,j,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,j,k))
c         Finite azimuthal continuation in the mapped transverse plane:
c         orthogonal to j and mapped r, norm -zj*zk*sjk. The collinear
c         limit and the vanishing azimuthal average of Q are unchanged.
          kt=zj*xp(:,k)+zk*(sjk-sjr)/(sjk+sjr)*xp(:,j)-sjk/(sjr+skr)*xp(:,r)
        else
          call phase_space_CS_inv(j,r,k,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,j,k))
          kt=-zk*xp(:,j)+zj*(skr-sjk)/(sjk+skr)*xp(:,k)+sjk/(sjr+skr)*xp(:,r)
        endif
        if(xjcs1.eq.0d0)goto 999
        do c=1,nexternal
          if(.not.isnnloqcdparton(c))cycle
          if(c.eq.i.or.c.eq.j.or.c.eq.k)cycle
          if(mapping.eq.1)then
            if(c.eq.r)cycle
c           Spectator pairs, including the two equal (c,r)/(r,c) terms.
            do d=1,nexternal
              if(.not.isnnloqcdparton(d))cycle
              if(d.eq.i.or.d.eq.j.or.d.eq.k.or.d.eq.r.or.d.eq.c)cycle
              call add_pair(c,d,1d0)
              if(ierr.ne.0)goto 999
            enddo
            call add_pair(c,r,2d0)
          elseif(mapping.eq.2)then
c           rho^(C)_jk = rho^(C)_kj = 1 for g -> g g; no d_ABC term.
            call add_pair(c,j,1d0)
          else
            call add_pair(c,k,1d0)
          endif
          if(ierr.ne.0)goto 999
        enddo
      enddo
      M2_SC_ggg=total
      if(abs(total).ge.huge(1d0).or.isnan(total))goto 999
      call ct_log('KSC                    ',M2_SC_ggg)
      return
 999  ierr=1
c     Do not expose a partial sum when a mapped point fails.
      M2_SC_ggg=0d0
      return

      contains

      subroutine add_pair(a,b,multiplicity)
      integer a,b,ab,bb,cb,db,pb_parent
      double precision multiplicity,eik,contraction,term
      if(xs(i,a).le.0d0.or.xs(i,b).le.0d0)then
        ierr=1
        return
      endif
      ab=real_labels(a,j,k)
      bb=real_labels(b,j,k)
      call fill_born_mapped_labels(ib,ab,real_leg_pdgs,Born_leg_pdgs)
      call phase_space_CS_inv(ib,ab,bb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,ib,ab))
      if(xjcs2.eq.0d0)then
        ierr=1
        return
      endif
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.ne.0)return
      cb=born_labels(ab,ib,ab)
      db=born_labels(bb,ib,ab)
      pb_parent=born_labels(parent,ib,ab)
      call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
      alphas=alphas_current
c     Keep the spin-summed dipole and the rank-one contraction distinct.
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,-1,alphas,ans)
      ccBLO=%(proc_prefix_Born)s_GET_CCBLO(cb,db)
## if (pyamplicol) {
      kkccBLO=%(proc_prefix_Born)s_GET_KKCCBLO(xpbb,alphas,pb_parent,kt,cb,db)
## } else {
c     The MadGraph spin hook also fills its colour-correlated table.
      unused=%(proc_prefix_Born)s_GET_KKBLO(pb_parent,xpbb,kt)
      kkccBLO=%(proc_prefix_Born)s_GET_CCBLO(cb,db)
c     Restore the spin-summed table after the legacy spin hook.
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,-1,alphas,ans)
## }
      contraction=Pjk*ccBLO+Qjk*(ccBLO+2d0*kkccBLO/kt2)
      eik=xs(a,b)/(xs(i,a)*xs(i,b))
      term=-(8d0*pi*alphas)**2*multiplicity*eik*contraction/sjk
      term=term*wsc_nnlo*extra*xj*%(proc_prefix_rr)s_fl_factor
      term=term*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      total=total+term
      if(doplot)then
        wgts(1)=-term*wgt/nit*wgt_chan
        call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      endif
      end subroutine add_pair
      end function M2_SC_ggg
