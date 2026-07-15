
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
      integer i,j,k,r,m,ierr,nit
      integer ib,jb,kb,rb,mb
      integer jbb,kbb,rbb,mbb
      double precision M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision ccBLO_ijm_krj,ccBLO_imj_krj,ccBLO_ikm_jrk,ccBLO_imk_jrk
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sik,sjk,sir,sjr,skr,sim,sjm,skm,sbjr,sbkr,sbjk
      double precision Ei_jm,Ei_km,Ebj_kr_ijm,Ebj_kr_imj,Ebj_kr_ikm,Ebj_kr_imk
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
      integer real_s_sc_mapped_labels(nexternal),Born_s_sc_mapped_labels(nexternal-1)
      common/c_NNLO_s_sc_mapped_labels/real_s_sc_mapped_labels,Born_s_sc_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_S_SS_gg_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
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
      pref=-(8d0*pi*alphas)**2
c
c     get PDGs
      jb = real_s_sc_mapped_labels(j)
      kb = real_s_sc_mapped_labels(k)
      rb = real_s_sc_mapped_labels(r)
c
c     Invariant quantities
      sij = xs(i,j)
      sik = xs(i,k)
      sjk = xs(j,k)
      sir = xs(i,r)
      sjr = xs(j,r)
      skr = xs(k,r)
c
c     compute soft double-soft soft-collinear sector function, (C.74) of 2212.11190
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
c
c     mapping 1: [ijm,krj]
c     eikonal sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_s_sc_mapped_labels(m)
         kbb = Born_s_sc_mapped_labels(kb)
         mbb = Born_s_sc_mapped_labels(mb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(i,j,m,xp,xpb,nexternal,leg_PDGs,xjCS1,real_s_sc_mapped_labels)
         call invariants_from_p(xpb,nexternal-1,xsb,ierr)
         if(ierr.eq.1)goto 999
         sbjr = xsb(jb,rb)
         sbkr = xsb(kb,rb)
         sbjk = xsb(jb,kb)
         if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
           write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_SC_ggq',sbjk,sbjr,sbkr
           goto 999
         endif
c
         Ebj_kr_ijm = sbkr/sbjk/sbjr
c
         call phase_space_CS_inv(kb,rb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_s_sc_mapped_labels)
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_ijm_krj = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c     invariant quantities
         sim = xs(i,m)
         sjm = xs(j,m)
         Ei_jm = sjm/sij/sim
c
c     soft double-soft soft-collinear kernel, (C.28)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = CA*Ei_jm*Ebj_kr_ijm*ccBLO_ijm_krj
         M2tmp = M2tmp*pref*ws_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_S_SS_gg_SC_ggq = M2_S_SS_gg_SC_ggq + M2tmp
c
c     plot
         wgtpl=+M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     mapping 2: [imj,krj]
c     eikonal sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_s_sc_mapped_labels(m)
         kbb = Born_s_sc_mapped_labels(kb)
         mbb = Born_s_sc_mapped_labels(mb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(i,m,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_s_sc_mapped_labels)
         call invariants_from_p(xpb,nexternal-1,xsb,ierr)
         if(ierr.eq.1)goto 999
         sbjr = xsb(jb,rb)
         sbkr = xsb(kb,rb)
         sbjk = xsb(jb,kb)
         if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
           write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_SC_ggq',sbjk,sbjr,sbkr
           goto 999
         endif
c
         Ebj_kr_imj = sbkr/sbjk/sbjr
c
         call phase_space_CS_inv(kb,rb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_s_sc_mapped_labels)
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_imj_krj = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c     invariant quantities
         sim = xs(i,m)
         sjm = xs(j,m)
         Ei_jm = sjm/sij/sim
c
c     soft double-soft soft-collinear kernel, (C.28)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = CA*Ei_jm*Ebj_kr_imj*ccBLO_imj_krj
         M2tmp = M2tmp*pref*ws_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_S_SS_gg_SC_ggq = M2_S_SS_gg_SC_ggq + M2tmp
c
c     plot
         wgtpl=+M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     mapping 3: [ikm,jrk]
c     eikonal sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_s_sc_mapped_labels(m)
         kbb = Born_s_sc_mapped_labels(kb)
         mbb = Born_s_sc_mapped_labels(mb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(i,k,m,xp,xpb,nexternal,leg_PDGs,xjCS1,real_s_sc_mapped_labels)
         call invariants_from_p(xpb,nexternal-1,xsb,ierr)
         if(ierr.eq.1)goto 999
         sbjr = xsb(jb,rb)
         sbkr = xsb(kb,rb)
         sbjk = xsb(jb,kb)
         if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
           write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_SC_ggq',sbjk,sbjr,sbkr
           goto 999
         endif
c
         Ebj_kr_ikm = sbkr/sbjk/sbjr
c
         call phase_space_CS_inv(jb,rb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_s_sc_mapped_labels)
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_ikm_jrk = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c     invariant quantities
         sim = xs(i,m)
         skm = xs(k,m)
         Ei_km = skm/sik/sim
c
c     soft double-soft soft-collinear kernel, (C.28)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = (2d0*CF-CA)*Ei_km*Ebj_kr_ikm*ccBLO_ikm_jrk
         M2tmp = M2tmp*pref*ws_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_S_SS_gg_SC_ggq = M2_S_SS_gg_SC_ggq + M2tmp
c
c     plot
         wgtpl=+M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     mapping 4: [imk,jrk]
c     eikonal sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         mb = real_s_sc_mapped_labels(m)
         kbb = Born_s_sc_mapped_labels(kb)
         mbb = Born_s_sc_mapped_labels(mb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(i,m,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_s_sc_mapped_labels)
         call invariants_from_p(xpb,nexternal-1,xsb,ierr)
         if(ierr.eq.1)goto 999
         sbjr = xsb(jb,rb)
         sbkr = xsb(kb,rb)
         sbjk = xsb(jb,kb)
         if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
           write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_SC_ggq',sbjk,sbjr,sbkr
           goto 999
         endif
c
         Ebj_kr_imk = sbkr/sbjk/sbjr
c
         call phase_space_CS_inv(jb,rb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_s_sc_mapped_labels)
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_imk_jrk = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c     invariant quantities
         sim = xs(i,m)
         sik = xs(i,k)
         Ei_km = skm/sik/sim
c
c     soft double-soft soft-collinear kernel, (C.28)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = (2d0*CF-CA)*Ei_km*Ebj_kr_imk*ccBLO_imk_jrk
         M2tmp = M2tmp*pref*ws_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_S_SS_gg_SC_ggq = M2_S_SS_gg_SC_ggq + M2tmp
c
c     plot
         wgtpl=+M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      if(test_sector_function) M2_S_SS_gg_SC_ggq = ws_nlo
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
