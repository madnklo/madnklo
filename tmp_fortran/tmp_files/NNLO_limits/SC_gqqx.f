## if (pyamplicol) {
      double precision function M2_SC_gqqx(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) C(j,k), with soft gluon i and a q-qbar pair j,k: Eq. (C.13) of 2212.11190v2.
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
      double precision sjk,sjr,skr,zj,zk,Pjk,Qjk,kt(0:3),kt2,dsign
      double precision dccBLO,kkdccBLO
      double precision wgt,wgts(1),wgt_chan,xj,xjb,extra,wsc_nnlo,total
      double precision xjcs1,xjcs2,alphas,ccBLO,kkccBLO,ans(0:NSQSO_BORN)
      double precision %(proc_prefix_Born)s_GET_CCBLO
      double precision %(proc_prefix_Born)s_GET_KKCCBLO
      double precision %(proc_prefix_Born)s_GET_DCCBLO
      double precision %(proc_prefix_Born)s_GET_KKDCCBLO
      logical docut,doplot,test_sector_function
      common/cdoplot/doplot
      common/ctestsecfun/test_sector_function
      include 'pmass.inc'

      M2_SC_gqqx=0d0
      total=0d0
      ierr=0
      if(min(i,j,k,r).le.nincoming.or.max(i,j,k,r).gt.nexternal)then
        write(*,*) 'SC_gqqx requires four distinct final-state legs'
        stop 1
      endif
      if(i.eq.j.or.i.eq.k.or.i.eq.r.or.j.eq.k.or.j.eq.r.or.k.eq.r)then
        write(*,*) 'SC_gqqx requires four distinct final-state legs'
        stop 1
      endif
      if(leg_pdgs(i).ne.21.or.abs(leg_pdgs(j)).lt.1.or.abs(leg_pdgs(j)).gt.5.or.leg_pdgs(j).ne.-leg_pdgs(k))then
        write(*,*) 'SC_gqqx requires a soft gluon and a q-qbar pair'
        stop 1
      endif
      if(.not.isnnloqcdparton(r))then
        write(*,*) 'SC_gqqx requires a coloured recoiler'
        stop 1
      endif
c     The existing CS inverse map is final-final; do not cross it silently.
      do c=1,nexternal
        if(.not.isnnloqcdparton(c))cycle
        if(c.le.nincoming)then
          write(*,*) 'SC_gqqx: coloured initial states are not implemented'
          stop 1
        endif
        if(pmass(c).ne.0d0)then
          write(*,*) 'SC_gqqx: massive coloured legs are not implemented'
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
      Pjk=TR*(1d0-2d0*zj*zk)
      Qjk=2d0*TR*zj*zk
c     Directed flavour factor: +1 for (q,qbar), -1 for (qbar,q).
      dsign=sign(1d0,dble(leg_pdgs(j)))

      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      if(i.eq.asec)then
        call get_wsc_ikl_nnlo(asec,bsec,csec,dsec,iref)
        wsc_nnlo=wsc_ikl_nnlo
      elseif(i.eq.csec)then
        call get_wsc_kij_nnlo(asec,bsec,csec,dsec,iref)
        wsc_nnlo=wsc_kij_nnlo
      else
        write(*,*) 'SC_gqqx: soft leg does not match the sector'
        stop 1
      endif
      if(test_sector_function)then
        M2_SC_gqqx=wsc_nnlo
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
c           rho^(C)_jk = rho^(C)_kj = 1; retain the directed d_ABC term.
            call add_pair(c,j,1d0)
          else
            call add_pair(c,k,1d0)
          endif
          if(ierr.ne.0)goto 999
        enddo
      enddo
      M2_SC_gqqx=total
      if(abs(total).ge.huge(1d0).or.isnan(total))goto 999
      call ct_log('KSC                    ',M2_SC_gqqx)
      return
 999  ierr=1
c     Do not expose a partial sum when a mapped point fails.
      M2_SC_gqqx=0d0
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
      ccBLO=%(proc_prefix_Born)s_GET_CCBLO(xpbb,alphas,cb,db)
      kkccBLO=%(proc_prefix_Born)s_GET_KKCCBLO(xpbb,alphas,pb_parent,kt,cb,db)
      contraction=Pjk*ccBLO+Qjk*(ccBLO+2d0*kkccBLO/kt2)
      if(mapping.ne.1)then
c       The adjoint parent is the d_ABC root, not the ordinary spectator.
        dccBLO=%(proc_prefix_Born)s_GET_DCCBLO(xpbb,alphas,pb_parent,cb)
        kkdccBLO=%(proc_prefix_Born)s_GET_KKDCCBLO(xpbb,alphas,pb_parent,kt,pb_parent,cb)
        if(mapping.eq.2)then
          contraction=contraction+dsign*(Pjk*dccBLO+Qjk*(dccBLO+2d0*kkdccBLO/kt2))
        else
          contraction=contraction-dsign*(Pjk*dccBLO+Qjk*(dccBLO+2d0*kkdccBLO/kt2))
        endif
      endif
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
      end function M2_SC_gqqx

## } else {

      double precision function M2_SC_gqqx(itmp,jtmp,ktmp,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     SC(i,j,k) kernel times WSC;
c     i is a gluon, j and k is a q qb pair)
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
      integer i,j,k,r,m,ierr,nit,itmp,jtmp,ktmp
      integer ib,jb,kb,mb,kbb,mbb
      double precision M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision ccBLO_krj_imj,ccBLO_jrk_imk
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sik,sjk,sjr,skr,sim,sjm,skm
      double precision zj,zk,Pjkr,Ei_jm,Ei_km
      double precision pref,extra,wsc_nnlo
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
c     set logical doplot
      logical flavourmatch_ijjk
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
      integer, parameter :: HEL = - 1
      double precision   %(proc_prefix_Born)s_GET_CCBLO
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      integer %(proc_prefix_Born)s_den
      common/%(proc_prefix_Born)s_iden/%(proc_prefix_Born)s_den
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer asec,bsec,csec,dsec
      common/csecindices/asec,bsec,csec,dsec
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
c
c     initialise
      M2_SC_gqqx=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_SC_gqqx',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch_ijjk = leg_PDGs(jtmp).eq.leg_PDGs(ktmp).and.abs(leg_PDGs(ktmp)).le.5.and.leg_PDGs(itmp).ne.leg_PDGs(ktmp)
      if(leg_pdgs(itmp).eq.0 .or. .not.(flavourmatch_ijjk)) then
        write(*,*) 'Flavour mismatch in M2_SC_gqqx',leg_PDGs(itmp),leg_PDGs(jtmp),leg_PDGs(ktmp)
        stop 1
      endif
c
c     rearrange the flavours to always land at gqqx TODO: to modify
      if(abs(leg_PDGs(itmp)).le.5)then
              i = ktmp
              j = jtmp
              k = itmp
      elseif(abs(leg_PDGs(jtmp)).le.5)then
              i = itmp
              j = ktmp
              k = jtmp
      else
              i = itmp
              j = jtmp
              k = ktmp
      endif
c
c     overall kernel prefix
c     renormalisation scale and prefactor: set right before each matrix
c     element below, on the remapped kinematics passed to it
c
c     invariant quantities
      sij = xs(i,j)
      sik = xs(i,k)
      sjk = xs(j,k)
      sjr = xs(j,r)
      skr = xs(k,r)
      zj = sjr/(sjr+skr)
      zk = skr/(sjr+skr)
      Pjkr = CF*(2d0*zk/zj+zj)
c
c     soft-collinear sector function, (C.55) or (C.56) of 2212.11190
      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      if(i.eq.asec)then
        call get_wsc_ikl_nnlo(asec,bsec,csec,dsec,iref)
        wsc_nnlo = wsc_ikl_nnlo
      elseif(i.eq.csec)then
        call get_wsc_kij_nnlo(asec,bsec,csec,dsec,iref)
        wsc_nnlo = wsc_kij_nnlo
      endif
c
c     mapping 1: [krj,imj]
c
c     fill mapped labels
      call fill_real_mapped_labels(j,k,leg_pdgs,real_leg_pdgs)
      ib = real_labels(i,j,k)
      jb = real_labels(j,j,k)
      kb = real_labels(k,j,k)
      call fill_born_mapped_labels(ib,jb,real_leg_pdgs,born_leg_pdgs)
      call phase_space_CS_inv(k,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,j,k))
c
      kbb = born_labels(kb,ib,jb)
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_labels(m,j,k)
         kbb = born_labels(kb,ib,jb)
         mbb = born_labels(mb,ib,jb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(ib,mb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,ib,jb))
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
c     renormalisation scale for the remapped kinematics xpbb
         call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
         alphas=alphas_current
         pref=-(8d0*pi*alphas)**2
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_krj_imj = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c     invariant quantities
         sim = xs(i,m)
         sjm = xs(j,m)
         Ei_jm = sjm/sij/sim
c
c     soft-collinear kernel, (C.13)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = Pjkr/sjk*Ei_jm*CA/CF*ccBLO_krj_imj
         M2tmp = M2tmp*pref*wsc_nnlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_SC_gqqx = M2_SC_gqqx + M2tmp
c
c     plot
         wgtpl=-M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     mapping 2: [jrk,imk]
      call fill_born_mapped_labels(ib,kb,real_leg_pdgs,born_leg_pdgs)
      call phase_space_CS_inv(j,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,j,k))
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_labels(m,j,k)
         kbb = born_labels(kb,ib,kb)
         mbb = born_labels(mb,ib,kb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(ib,mb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,ib,kb))
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
c     renormalisation scale for the remapped kinematics xpbb
         call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
         alphas=alphas_current
         pref=-(8d0*pi*alphas)**2
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
         ccBLO_jrk_imk = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c     invariant quantities
         sim = xs(i,m)
         skm = xs(k,m)
         Ei_km = skm/sik/sim
c
c     soft-collinear kernel according to (C.13)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = Pjkr/sjk*Ei_km*(2d0*CF-CA)/CF*ccBLO_jrk_imk
         M2tmp = M2tmp*pref*wsc_nnlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_SC_gqqx = M2_SC_gqqx + M2tmp
c
c     plot
         wgtpl=-M2tmp*wgt/nit*wgt_chan
         wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      if(test_sector_function) M2_SC_gqqx = wsc_nnlo
c
      call ct_log('KSC                    ',M2_SC_gqqx)
c
c     sanity check
      if(abs(M2_SC_gqqx).ge.huge(1d0).or.isnan(M2_SC_gqqx))then
         write(77,*)'Exception caught in M2_SC_gqqx',M2_SC_gqqx
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

## }
