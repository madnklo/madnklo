
      double precision function M2_SC_GGQ(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     SC(i,j,k) kernel times WSC;
c     i, j is a g-g pair, k is a q (or qb)
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
      integer i,j,k,r,m,ierr,nit
      integer ib,jb,kb,mb,kbb,mbb
      double precision M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision ccBLO_krj_imj,ccBLO_jrk_imk
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sik,sjk,sjr,skr,sim,sjm,skm
      double precision zj,zk,Pjkr,Ei_jm,Ei_km
      double precision pref,extra
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
      integer real_sc_mapped_labels(nexternal),Born_sc1_mapped_labels(nexternal-1),Born_sc2_mapped_labels(nexternal-1)
      common/c_NNLO_sc_mapped_labels/real_sc_mapped_labels,Born_sc1_mapped_labels,Born_sc2_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_SC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_SC_ggq', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=-(8d0*pi*alphas)**2
c
c     get mapped labels
      ib = real_sc_mapped_labels(i)
      jb = real_sc_mapped_labels(j)
      kb = real_sc_mapped_labels(k)
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
c     soft-collinear sector function, (C.55) of 2212.11190
      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      call get_wsc_nnlo(asec,bsec,csec,dsec,iref)
c
c     mapping 1: [krj,imj]
      call phase_space_CS_inv(k,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_sc_mapped_labels(m)
         kbb = Born_sc1_mapped_labels(kb)
         mbb = Born_sc1_mapped_labels(mb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(ib,mb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc1_mapped_labels)
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
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
         M2_SC_ggq = M2_SC_ggq + M2tmp
c
c     plot
         wgtpl=-M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     mapping 2: [jrk,imk]
      call phase_space_CS_inv(j,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_sc_mapped_labels(m)
         kbb = Born_sc2_mapped_labels(kb)
         mbb = Born_sc2_mapped_labels(mb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(ib,mb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
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
         M2_SC_ggq = M2_SC_ggq + M2tmp
c
c     plot
         wgtpl=-M2tmp*wgt/nit*wgt_chan
         wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      if(test_sector_function) M2_SC_ggq = wsc_nnlo
c
c     sanity check
      if(abs(M2_SC_ggq).ge.huge(1d0).or.isnan(M2_SC_ggq))then
         write(77,*)'Exception caught in M2_SC_ggq',M2_SC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
