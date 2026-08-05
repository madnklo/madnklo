
      double precision function M2_S_SS_gg_SC_ggq(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) S(i,k) SC(i,k,l) kernel times WS_SS_SC
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,r,a,b,m,ierr,nit
      integer ib,ab,bb,rb,mb,bbb,mbb
      double precision M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1),xsbb(nexternal-2,nexternal-2)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),xpbb(0:3,nexternal-2)
      double precision ccBLO_ima_bra,ccBLO_imb_arb
      double precision sia,sim,sam,sib,sbm,sbab,sbar,sbbr
      double precision Ei_am,Ei_bm,Eba_br
      double precision pref,extra,damp
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
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
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_s_sc_1_mapped_labels(nexternal),real_s_sc_2_mapped_labels(nexternal),Born_s_sc_1_mapped_labels(nexternal-1),Born_s_sc_2_mapped_labels(nexternal-1)
      common/c_nnlo_s_sc_mapped_labels/real_s_sc_1_mapped_labels,real_s_sc_2_mapped_labels,Born_s_sc_1_mapped_labels,Born_s_sc_2_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_S_SS_gg_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
      damp=1d0
c
c     check sector topology
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_S_SS_gg_SC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_S_SS_gg_SC_ggq', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=-2d0*(8d0*pi*alphas)**2
c
      a = csec
      b = dsec
c
c     compute soft soft-collinear sector function, (c.73) of 2212.11190
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
c
c     mapping 1: [ima,bra]
c
c     get mapped labels
      ab = real_s_sc_1_mapped_labels(a)
      bb = real_s_sc_1_mapped_labels(b)
      rb = real_s_sc_1_mapped_labels(r)
c
c     eikonal sum
      do m=1,nexternal
        if(.not.isnnloqcdparton(m))cycle
        if(m.eq.i.or.m.eq.a.or.m.eq.b)cycle
        mb = real_s_sc_1_mapped_labels(m)
        bbb = Born_s_sc_1_mapped_labels(bb)
        mbb = Born_s_sc_1_mapped_labels(mb)
c
c       underlying born configuration is remapped
        call phase_space_cs_inv(i,m,a,xp,xpb,nexternal,leg_pdgs,xjcs1,real_s_sc_1_mapped_labels)
        call invariants_from_p(xpb,nexternal-1,xsb,ierr)
        if(ierr.eq.1)goto 999
        sbab = xsb(ab,bb)
        sbar = xsb(ab,rb)
        sbbr = xsb(bb,rb)
        if(sbab.lt.0d0.or.sbar.lt.0d0.or.sbbr.lt.0d0) then
          write(77,*)'inaccuracy 2 in m2_s_sc_ggq',sbab,sbar,sbbr
          goto 999
        endif
        eba_br = sbbr/sbab/sbar
c
        call phase_space_cs_inv(bb,rb,ab,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_s_sc_1_mapped_labels)
        if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c       possible cuts
        if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c       call colour-connected born matrix element
        call epem_ccx_me_accessor_hook(xpbb,hel,alphas,ans)
        ccblo_ima_bra = epem_ccx_get_ccblo(bbb,mbb)
c
c       invariant quantities
        sim = xs(i,m)
        sia = xs(i,a)
        sam = xs(a,m)
        Ei_am = sam/sia/sim
c
c       soft soft-soft soft-collinear kernel, (143) on Dropbox
c       todo: some contributions are 0 for ee->jj
        m2tmp = m2tmp+ca*Ei_am*Eba_br*ccblo_ima_bra
c
c       plot
        wgtpl = -m2tmp*ws_nlo
        wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
        wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
        wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
        wgts=wgtpl
c       if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,w
c       gtpl)
        if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     mapping 2: [imb,arb]
c
c     get mapped labels
      ab = real_s_sc_2_mapped_labels(a)
      bb = real_s_sc_2_mapped_labels(b)
      rb = real_s_sc_2_mapped_labels(r)
c
c     eikonal sum
      do m=1,nexternal
        if(.not.isnnloqcdparton(m))cycle
        if(m.eq.i.or.m.eq.a.or.m.eq.b)cycle
        mb = real_s_sc_2_mapped_labels(m)
        bbb = Born_s_sc_2_mapped_labels(bb)
        mbb = Born_s_sc_2_mapped_labels(mb)
c
c       underlying born configuration is remapped
        call phase_space_cs_inv(i,m,b,xp,xpb,nexternal,leg_pdgs,xjcs1,real_s_sc_2_mapped_labels)
        call invariants_from_p(xpb,nexternal-1,xsb,ierr)
        if(ierr.eq.1)goto 999
        sbab = xsb(ab,bb)
        sbar = xsb(ab,rb)
        sbbr = xsb(bb,rb)
        if(sbab.lt.0d0.or.sbar.lt.0d0.or.sbbr.lt.0d0) then
          write(77,*)'inaccuracy 2 in m2_s_sc_ggq',sbab,sbar,sbbr
          goto 999
        endif
        eba_br = sbbr/sbab/sbar
c
        call phase_space_cs_inv(ab,rb,bb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_s_sc_2_mapped_labels)
        if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c       possible cuts
        if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c       call colour-connected born matrix element
        call epem_ccx_me_accessor_hook(xpbb,hel,alphas,ans)
        ccblo_imb_arb = epem_ccx_get_ccblo(bbb,mbb)
c
c       invariant quantities
        sim = xs(i,m)
        sib = xs(i,b)
        sbm = xs(b,m)
        Ei_bm = sbm/sib/sim
c
c       soft soft-soft soft-collinear kernel, (143) on Dropbox
c       todo: some contributions are 0 for ee->jj
        m2tmp = m2tmp+(2d0*cf-ca)*Ei_bm*Eba_br*ccblo_imb_arb
c
c       plot
        wgtpl = -m2tmp*ws_nlo
        wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
        wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
        wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
        wgts=wgtpl
c       if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,w
c       gtpl)
        if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      m2_s_ss_gg_sc_ggq = m2tmp*pref*ws_nlo*xj*damp*extra*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      m2_s_ss_gg_sc_ggq = m2_s_ss_gg_sc_ggq * %(proc_prefix_rr)s_fl_factor
c
      if(test_sector_function) M2_S_SS_gg_SC_ggq = ws_nlo
c
      call ct_log('KS_SS_SC             ',M2_S_SS_gg_SC_ggq)
c
c     sanity check
      if(abs(M2_S_SS_gg_SC_ggq).ge.huge(1d0).or.isnan(M2_S_SS_gg_SC_ggq))then
         write(77,*)'Exception caught in M2_S_SS_gg_SC_ggq',M2_S_SS_gg_SC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
