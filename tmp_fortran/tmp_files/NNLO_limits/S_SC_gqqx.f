## if (pyamplicol) {
      double precision function M2_S_SC_gqqx(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) SC(i,j,k), quark-antiquark pair: Eq. (C.27) of 2212.11190v2.
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
      double precision pmass(nexternal),wgt,xj,xjb,extra,wgt_chan,ws,total,dsign
      logical docut,doplot,test_sector_function
      common/cdoplot/doplot
      common/ctestsecfun/test_sector_function
      include 'pmass.inc'

      M2_S_SC_gqqx=0d0
      total=0d0
      ierr=0
      if(min(i,j,k,r).le.nincoming.or.max(i,j,k,r).gt.nexternal)then
        write(*,*) 'S_SC_gqqx requires four distinct final-state legs'
        stop 1
      endif
      if(i.eq.j.or.i.eq.k.or.i.eq.r.or.j.eq.k.or.j.eq.r.or.k.eq.r)then
        write(*,*) 'S_SC_gqqx requires four distinct final-state legs'
        stop 1
      endif
      if(leg_pdgs(i).ne.21.or.abs(leg_pdgs(j)).lt.1.or.abs(leg_pdgs(j)).gt.5.or.leg_pdgs(j).ne.-leg_pdgs(k).or..not.isnnloqcdparton(r))then
        write(*,*) 'S_SC_gqqx requires a gluon, q-qbar pair and coloured recoiler'
        stop 1
      endif
      if(i.ne.asec)then
        write(*,*) 'S_SC_gqqx: soft leg does not match the ordered sector'
        stop 1
      endif
      do c=1,nexternal
        if(.not.isnnloqcdparton(c))cycle
        if(c.le.nincoming.or.pmass(c).ne.0d0)then
          write(*,*) 'S_SC_gqqx: only massless final-state QCD is implemented'
          stop 1
        endif
      enddo
c     Directed flavour factor for the two parent-spectator families.
      dsign=sign(1d0,dble(leg_pdgs(j)))
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
c       rho^(C)_jk = rho^(C)_kj = 1; retain the directed d_ABC term.
        call add_term(j,c,2,.5d0)
        if(ierr.ne.0)goto 999
        call add_term(c,j,2,.5d0)
        if(ierr.ne.0)goto 999
        call add_term(k,c,3,.5d0)
        if(ierr.ne.0)goto 999
        call add_term(c,k,3,.5d0)
        if(ierr.ne.0)goto 999
      enddo
      M2_S_SC_gqqx=total
      if(abs(total).ge.huge(1d0).or.isnan(total))goto 999
      call ct_log('KS_SC                  ',M2_S_SC_gqqx)
      return
 999  ierr=1
      M2_S_SC_gqqx=0d0
      return

      contains

      subroutine add_term(a,b,mapping,multiplicity)
      integer a,b,mapping,ab,bb,jb,kb,rb,parent,cb,db,spectator
      double precision multiplicity,xjcs1,xjcs2,zj,zk,sjk,sjr,skr,Pjk,Qjk,kt2
      double precision dccBLO,kkdccBLO
      double precision xpb(0:3,nexternal-1),xsb(nexternal-1,nexternal-1)
      double precision xpbb(0:3,nexternal-2),xsbb(nexternal-2,nexternal-2),kt(0:3)
      double precision alphas,ccBLO,kkccBLO,ans(0:NSQSO_BORN),eik,contraction,term,wc,wgts(1)
      double precision %(proc_prefix_Born)s_GET_CCBLO
      double precision %(proc_prefix_Born)s_GET_KKCCBLO
      double precision %(proc_prefix_Born)s_GET_DCCBLO
      double precision %(proc_prefix_Born)s_GET_KKDCCBLO
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
      Pjk=TR*(1d0-2d0*zj*zk)
      Qjk=2d0*TR*zj*zk
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
c       Same finite transverse-plane prescription as SC_gqqx.
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
      kkccBLO=%(proc_prefix_Born)s_GET_KKCCBLO(xpbb,alphas,parent,kt,cb,db)
      contraction=Pjk*ccBLO+Qjk*(ccBLO+2d0*kkccBLO/kt2)
      if(mapping.ne.1)then
c       In either soft-map orientation the d_ABC root is the parent;
c       the ordinary colour charge belongs to the other spectator.
        spectator=cb
        if(a.eq.j.or.a.eq.k)spectator=db
        dccBLO=%(proc_prefix_Born)s_GET_DCCBLO(xpbb,alphas,parent,spectator)
        kkdccBLO=%(proc_prefix_Born)s_GET_KKDCCBLO(xpbb,alphas,parent,kt,parent,spectator)
        if(mapping.eq.2)then
          contraction=contraction+dsign*(Pjk*dccBLO+Qjk*(dccBLO+2d0*kkdccBLO/kt2))
        else
          contraction=contraction-dsign*(Pjk*dccBLO+Qjk*(dccBLO+2d0*kkdccBLO/kt2))
        endif
      endif
      eik=xs(a,b)/xs(i,a)/xs(i,b)
      term=-(8d0*pi*alphas)**2*multiplicity*eik*contraction/sjk*ws*wc
      term=term*extra*xj*%(proc_prefix_rr)s_fl_factor*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      total=total+term
      if(doplot)then
        wgts(1)=-term*wgt/nit*wgt_chan
        call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      endif
      end subroutine add_term
      end function M2_S_SC_gqqx

## } else {

      double precision function M2_S_SC_gqqx(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) SC(i,j,k) kernel times WS_SC;
c     j, k is a q-qb pair, i is a g
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
      integer i,j,k,r,a,b,m,ierr,nit
      integer ib,ab,bb,rb,mb,bbb,mbb
      double precision M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1),xsbb(nexternal-2,nexternal-2)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),xpbb(0:3,nexternal-2)
      double precision sia,sib,sab,sim,sam,sbm,sbar,sbbr,sbab,zq
      double precision Ei_am,Ei_bm,Pab_ima,Pab_imb,rho_ab,rho_ba
      double precision pref,extra,damp
      double precision alphas,ans(0:NSQSO_BORN),ccblo
      double precision alpha_qcd
      logical flavourmatch_ijjk,flavourmatch_ijkj
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
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
c
c     initialise
      M2_S_SC_gqqx=0d0
      M2tmp=0d0
      ierr=0
      damp=1d0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_S_SC_gqqx',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch_ijjk = leg_PDGs(i).eq.leg_PDGs(j).and.abs(leg_PDGs(k)).le.5.and.leg_PDGs(i).ne.leg_PDGs(k)
      flavourmatch_ijkj = leg_PDGs(i).eq.leg_PDGs(k).and.abs(leg_PDGs(j)).le.5.and.leg_PDGs(i).ne.leg_PDGs(j)
      if(leg_pdgs(i).eq.0 .or. .not.(flavourmatch_ijjk.or.flavourmatch_ijkj)) then
        write(*,*) 'Flavour mismatch in M2_S_SC_gqqx',leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     overall kernel prefix
c     renormalisation scale and prefactor: set right before each matrix
c     element below, on the remapped kinematics passed to it
c
      a = csec
      b = dsec
      sia = xs(i,a)
      sib = xs(i,b)
c
      if(leg_PDGs(a).eq.21.and.abs(leg_PDGs(b)).le.6)then
        rho_ab = rho_gq
        rho_ba = rho_qg
      elseif(leg_PDGs(b).eq.21.and.abs(leg_PDGs(a)).le.6)then
        rho_ab = rho_qg
        rho_ba = rho_gq
      endif
c
c     compute soft soft-collinear sector function, (C.73) of 2212.11190
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
c
c     soft soft-collinear kernel, (126) on Dropbox
c     ToDo: some contributions are 0 for ee->jj

c     mapping 1: [ima,bra]
c
c     fill mapped labels
      call fill_real_mapped_labels(i,a,leg_pdgs,real_leg_pdgs)
      call fill_born_mapped_labels(real_labels(b,i,a),real_labels(a,i,a),real_leg_pdgs,born_leg_pdgs)
c
      ab = real_labels(a,i,a)
      bb = real_labels(b,i,a)
      rb = real_labels(r,i,a)
      bbb = Born_labels(bb,bb,ab)
c
c     eikonal sum
      do m=1,nexternal
        if(.not.isnnloqcdparton(m))cycle
        if(m.eq.i.or.m.eq.a.or.m.eq.b)cycle
c
        sim = xs(i,m)
        sam = xs(a,m)
        if(sia*sim.lt.0) then
          write(77,*)'Inaccuracy 2 in M2_S_SC_gqqx',sim,sia
          goto 999
        endif
        Ei_am = sam/sia/sim
c
        mb = real_labels(m,i,a)
        mbb = Born_labels(mb,bb,ab)
        call phase_space_cs_inv(i,m,a,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,i,a))
        call invariants_from_p(xpb,nexternal-1,xsb,ierr)
        if(ierr.eq.1)goto 999
        sbar = xsb(ab,rb)
        sbbr = xsb(bb,rb)
        sbab = xsb(ab,bb)
        if(sbar+sbbr.lt.0d0.and.sbab.lt.0) then
          write(77,*)'Inaccuracy 2 in M2_S_SC_gqqx',sbar,sbbr,sbab
          goto 999
        endif
        if(leg_PDGs(a).eq.21.and.abs(leg_PDGs(b)).le.6)then
          zq = sbbr/(sbar+sbbr)
        elseif(leg_PDGs(b).eq.21.and.abs(leg_PDGs(a)).le.6)then
          zq = sbar/(sbar+sbbr)
        endif
        Pab_ima = CF*(2d0*zq/(1d0-zq)+(1d0-zq))
c
        call phase_space_cs_inv(bb,rb,ab,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_labels(:,bb,ab))
        if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c       possible cuts
        if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c       call wcbar
        call get_sig2(xsb,alpha_mod_bar,nexternal-1)
        call get_wcbar_nlo(ab,bb,rb)
c
c       call colour-connected born matrix element
c     renormalisation scale for the remapped kinematics xpbb
        call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
        alphas=alphas_current
        pref=-(8d0*pi*alphas)**2
        call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
        ccBLO = %(proc_prefix_Born)s_get_ccBLO(bbb,mbb)
c
        M2tmp = M2tmp+Ei_am*Pab_ima/sbab*rho_ab*ccblo*wcbar_nlo
c
c       plot
        wgtpl = -M2tmp*ws_nlo
        wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
        wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
        wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
        wgts=wgtpl
c       if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,born_leg_pdgs,w
c       gtpl)
        if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,born_leg_pdgs,wgts)
      enddo
c
c     mapping 2: [imb,arb]
c
c     fill mapped labels
      call fill_real_mapped_labels(i,b,leg_pdgs,real_leg_pdgs)
      call fill_born_mapped_labels(real_labels(a,i,b),real_labels(b,i,b),real_leg_pdgs,born_leg_pdgs)
c
c     eikonal sum
      do m=1,nexternal
        if(.not.isnnloqcdparton(m))cycle
        if(m.eq.i.or.m.eq.a.or.m.eq.b)cycle
c
        sim = xs(i,m)
        sbm = xs(b,m)
        if(sib*sim.lt.0) then
          write(77,*)'Inaccuracy 2 in M2_S_SC_gqqx',sim,sib
          goto 999
        endif
        Ei_bm = sbm/sib/sim
c
        mb = real_labels(m,i,b)
        mbb = Born_labels(mb,ab,bb)
        call phase_space_cs_inv(i,m,b,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,i,b))
        call invariants_from_p(xpb,nexternal-1,xsb,ierr)
        if(ierr.eq.1)goto 999
        sbar = xsb(ab,rb)
        sbbr = xsb(bb,rb)
        sbab = xsb(ab,bb)
        if(sbar+sbbr.lt.0d0.or.sbab.lt.0d0) then
          write(77,*)'Inaccuracy 2 in M2_S_SC_gqqx',sbar,sbbr,sbab
          goto 999
        endif
        if(leg_PDGs(a).eq.21.and.abs(leg_PDGs(b)).le.6)then
          zq = sbbr/(sbar+sbbr)
        elseif(leg_PDGs(b).eq.21.and.abs(leg_PDGs(a)).le.6)then
          zq = sbar/(sbar+sbbr)
        endif
        Pab_imb = CF*(2d0*zq/(1d0-zq)+(1d0-zq))/sbab
c
        call phase_space_cs_inv(ab,rb,bb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_labels(:,ab,bb))
        if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c       possible cuts
        if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c       call wcbar
        call get_sig2(xsb,alpha_mod_bar,nexternal-1)
        call get_wcbar_nlo(ab,bb,rb)
c
c       call colour-connected born matrix element
c     renormalisation scale for the remapped kinematics xpbb
        call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
        alphas=alphas_current
        pref=-(8d0*pi*alphas)**2
        call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
        ccblo = %(proc_prefix_Born)s_get_ccBLO(bbb,mbb)
c
        M2tmp = M2tmp+Ei_bm*Pab_imb*rho_ba*ccblo*wcbar_nlo
c
c       plot
        wgtpl = -M2tmp*ws_nlo
        wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
        wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
        wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
        wgts=wgtpl
c       if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,w
c       gtpl)
        if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      M2_S_SC_gqqx = M2tmp*pref*ws_nlo*xj*damp*extra*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_S_SC_gqqx = M2_S_SC_gqqx * %(proc_prefix_rr)s_fl_factor
c
      if(test_sector_function) M2_S_SC_gqqx = wcbar_nlo*ws_nlo
c
      call ct_log('KS_SC                  ',M2_S_SC_gqqx)
c
c     sanity check
      if(abs(M2_S_SC_gqqx).ge.huge(1d0).or.isnan(M2_S_SC_gqqx))then
         write(77,*)'Exception caught in M2_S_SC_gqqx',M2_S_SC_gqqx
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

## }
