
      double precision function M2_S_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     single-soft limit S_(i) * Wsoft
c     it returns 0 if i is not a gluon
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
      integer i,nit,ierr
      double precision xs(nexternal,nexternal),xp(0:3,nexternal),wgt,xj,xjB,extra,wgt_chan
      integer a,b,c,ab,bb,cb
      double precision xpb(0:3,nexternal-1),xsb(nexternal-1,nexternal-1)
      double precision Ei_ca,Ei_cb,Ei_ab,ccRNLO
      double precision M2tmp,pref,damp,xjCS1
      double precision wgts(1),wgtpl
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
      double precision  %(proc_prefix_S_g)s_GET_CCBLO
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      integer %(proc_prefix_S_g)s_den
      common/%(proc_prefix_S_g)s_iden/%(proc_prefix_S_g)s_den
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer asec,bsec,csec,dsec
      common/csecindices/asec,bsec,csec,dsec
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      M2_S_g=0d0
      M2tmp=0d0
      ierr=0
      damp=1d0
c
c     return if not gluon
      if(leg_pdgs(i).ne.21)return
c
c     safety check on PDGs
      if(size(leg_pdgs).ne.nexternal)then
        write(*,*) 'M2_S_g:'
        write(*,*) 'Wrong dimension for leg_pdgs',size(leg_pdgs),nexternal
        stop
      endif
c
c     call soft limit of sector function according to eq. (C.51)
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=-(8d0*pi*alphas)
c
      a = csec
      b = dsec
      ab = real_mapped_labels(a)
      bb = real_mapped_labels(b)
c
c     mapping 1: [ica]
      do c=1,nexternal
         if(.not.isnnloQCDparton(c))cycle
         if(c.eq.i.or.c.eq.a.or.c.eq.b)cycle
c
         call phase_space_CS_inv(i,c,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
         if(xjCS1.eq.0d0)goto 999
c        possible cuts
         if(docut(xpb,nexternal-1,real_leg_pdgs,1))cycle
         call invariants_from_p(xpb,nexternal-1,xsb,ierr)
         if(ierr.eq.1)goto 999
c
c        call Wbar
         call get_sig2(xsb,alpha_mod_bar,nexternal-1)
         map1=real_mapped_labels(csec)
         map2=real_mapped_labels(dsec)
         call get_wbar_nlo(map1,map2)
c
c        call colour-connected Born
         cb=real_mapped_labels(c)
         call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
         ccRNLO = %(proc_prefix_S_g)s_GET_CCBLO(cb,ab)
c
c        safety check
         if(xs(i,c)*xs(i,a).le.0d0)then
            write(77,*)'Inaccuracy 1 in M2_S_g',xs(i,c)*xs(i,a)
            goto 999
         endif
         Ei_ca = xs(c,a)/xs(i,c)/xs(i,a)
         M2tmp = 2d0*Ei_ca*ccRNLO*wbar_nlo*ws_nlo
         M2_s_g = M2_s_g + M2tmp
c
c     plot
         wgtpl = -m2tmp
         wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
         wgtpl = wgtpl*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_rr)s_den)
         wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
         wgts = wgtpl
c     if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
         if(doplot)call analysis_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgts)
      enddo
c
c     mapping 2: [icb]
      do c=1,nexternal
         if(.not.isnnloQCDparton(c))cycle
         if(c.eq.i.or.c.eq.a.or.c.eq.b)cycle
c
         call phase_space_CS_inv(i,c,b,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
         if(xjCS1.eq.0d0)goto 999
c        possible cuts
         if(docut(xpb,nexternal-1,real_leg_pdgs,1))cycle
         call invariants_from_p(xpb,nexternal-1,xsb,ierr)
         if(ierr.eq.1)goto 999
c
c        call Wbar
         call get_sig2(xsb,alpha_mod_bar,nexternal-1)
         map1=real_mapped_labels(csec)
         map2=real_mapped_labels(dsec)
         call get_wbar_nlo(map1,map2)
c
c        call colour-connected Born
         cb=real_mapped_labels(c)
         call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
         ccRNLO = %(proc_prefix_S_g)s_GET_CCBLO(cb,bb)
c
c        safety check
         if(xs(i,c)*xs(i,b).le.0d0)then
            write(77,*)'Inaccuracy 1 in M2_S_g',xs(i,c),xs(i,b)
            goto 999
         endif
         Ei_cb = xs(c,b)/xs(i,c)/xs(i,b)
         M2tmp = 2d0*Ei_cb*ccRNLO*wbar_nlo*ws_nlo
         M2_s_g = M2_s_g + M2tmp
c
c     plot
         wgtpl = -m2tmp
         wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
         wgtpl = wgtpl*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_rr)s_den)
         wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
         wgts = wgtpl
c     if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
         if(doplot)call analysis_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgts)
      enddo
c
c     mapping 3: [iba]
      call phase_space_CS_inv(i,b,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
      if(xjCS1.eq.0d0)goto 999
c     possible cuts
      if(docut(xpb,nexternal-1,real_leg_pdgs,1))return
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Wbar
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wbar_nlo(map1,map2)
c
c     call colour-connected Born
      call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      ccRNLO = %(proc_prefix_S_g)s_GET_CCBLO(ab,bb)
c     safety check
      if(xs(i,a)*xs(i,b).le.0d0)then
         write(77,*)'Inaccuracy 1 in M2_S_g',xs(i,a),xs(i,b)
         goto 999
      endif
      Ei_ab = xs(a,b)/xs(i,a)/xs(i,b)
      M2tmp = 2d0*Ei_ab*ccRNLO*wbar_nlo*ws_nlo
      M2_s_g = M2_s_g + M2tmp
c
c     plot
      wgtpl = -m2tmp
      wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
      wgtpl = wgtpl*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_rr)s_den)
      wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
      wgts = wgtpl
c     if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
      if(doplot)call analysis_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgts)
c
c     Include correct pre-factors
      M2_s_g = M2_s_g*pref*extra*damp*xj
      M2_s_g = M2_s_g*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_S_g = M2_S_g * %(proc_prefix_rr)s_fl_factor
c
      if(test_sector_function) M2_S_g = ws_nlo*wbar_nlo
c
      call ct_log('KS                   ',M2_S_g)
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
