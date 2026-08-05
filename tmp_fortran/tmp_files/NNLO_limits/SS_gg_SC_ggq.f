
      double precision function M2_SS_gg_SC_ggq(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i,j)SC(i,j,k) kernel times WSS_SC
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
      integer ib,jb,kb,rb,mb
      integer jbb,kbb,rbb,mbb
      double precision M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision ccBLO_krj_imj,ccBLO_jrk_imk
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sik,sjk,sir,sjr,skr,sim,sjm,skm
      double precision zj,zk,Ej_kr,Ei_jm,Ei_km
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
      logical consistency_check
      common/cconscheck/consistency_check
c
c     initialise
      M2_SS_gg_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_SS_gg_SC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_SS_gg_SC_ggq', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=-2d0*(8d0*pi*alphas)**2
c
c     get PDGs
      ib = real_sc_mapped_labels(i)
      jb = real_sc_mapped_labels(j)
      kb = real_sc_mapped_labels(k)
      rb = real_sc_mapped_labels(r)
c
c     Invariant quantities
      sij = xs(i,j)
      sik = xs(i,k)
      sjk = xs(j,k)
      sir = xs(i,r)
      sjr = xs(j,r)
      skr = xs(k,r)
      zj = sjr/(sjr+skr)
      zk = skr/(sjr+skr)
c
c     call W double-soft soft-collinear, (C.57-C.58) of 2212.11190
      call get_sig2(xs,1d0,nexternal)
      call get_wss_sc_nnlo(asec,bsec,csec,dsec)
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
c
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
         Ej_kr = skr/sjk/sjr
         Ei_jm = sjm/sij/sim
c
c     double-soft soft-collinear kernel, (C.14)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = Ej_kr*CA*Ei_jm*ccBLO_krj_imj
         M2tmp = M2tmp*pref*wss_sc_nnlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_SS_gg_SC_ggq = M2_SS_gg_SC_ggq + M2tmp
c
c     plot
         wgtpl=+M2tmp*wgt/nit*wgt_chan
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
         Ej_kr = skr/sjk/sjr
         Ei_km = skm/sik/sim
c
c     double-soft soft-collinear kernel, (C.14)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = Ej_kr*(2d0*CF-CA)*Ei_km*ccBLO_jrk_imk
         M2tmp = M2tmp*pref*wss_sc_nnlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_SS_gg_SC_ggq = M2_SS_gg_SC_ggq + M2tmp
c
c     plot
         wgtpl=+M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      if(test_sector_function) M2_SS_gg_SC_ggq = wss_sc_nnlo
c
      call ct_log('KSS_SC               ',M2_SS_gg_SC_ggq)
c
c     sanity check
      if(abs(M2_SS_gg_SC_ggq).ge.huge(1d0).or.isnan(M2_SS_gg_SC_ggq))then
         write(77,*)'Exception caught in M2_SS_gg_SC_ggq',M2_SS_gg_SC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
