
      double precision function M2_SS_GG_SC_GGQ(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
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
      integer ib,jb,mb,rb,kb
      integer jbb,mbb,rbb,kbb
      logical isNLOmappedQCDparton(nexternal-1)
      logical isLOmappedQCDparton(nexternal-2)
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO_lrkimk,ccBLO_krliml,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sjk,sjm,sik,sim,sjr,skr,skm,zk,zj,ml2,mm2,y,z,x,damp
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
c      integer (proc_prefix_S_g)s_den
c      common/(proc_prefix_S_g)s_iden/(proc_prefix_S_g)s_den
      integer %(proc_prefix_Born)s_den
      common/%(proc_prefix_Born)s_iden/%(proc_prefix_Born)s_den
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
c
c     initialise
      M2_SS_gg_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_SS_gg_SC_ggq', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     check sector topology
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_SS_gg_SC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     call W double-soft soft-collinear, eq. (C.57-C.58) of 2212.11190
c     a small detail is that sig2 is always called with alpha=1 in the limit
c     the necessary sig2's are raised to the respective alpha in the soft-collinear sector functions
      call get_sig2(xs,1d0,nexternal)
      call get_wss_sc_nnlo(asec,bsec,csec,dsec)
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=32d0*pi**2*alphas**2
c
c     get PDGs
      ib = real_mapped_labels(i)
      jb = real_mapped_labels(j)
      rb = real_mapped_labels(r)
      kb = real_mapped_labels(k)
      rbb = Born_mapped_labels(rb)
      kbb = Born_mapped_labels(kb)
c
c     Invariant quantities
      sij = xs(i,j)
      sjk = xs(j,k)
      sik = xs(i,k)
      skr = xs(k,r)
      sjr = xs(j,r)
      zk = skr/(skr+sjr)
      zj = sjr/(sjr+skr)
c
c     Mapping 1 for B[lrk,imk]
      call phase_space_CS_inv(j,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(m))cycle
         if(m.eq.i.or.m.eq.k.or.m.eq.j)cycle
c
            mb = real_mapped_labels(m)
            mbb = Born_mapped_labels(mb)
c
c           Mapping 1 for B[lrk,imk]
c
c           underlying Born configuration is remapped
            call phase_space_CS_inv(ib,mb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c           call colour-connected Born matrix element with the mapping [lrk,imk]
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO_lrkimk = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c           possible cuts
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c           invariant quantities: c --> m
            skm = xs(k,m)
            sim = xs(i,m)
c
c           Double-soft soft-collinear kernel according to the eq.(C.14)
            M2tmp = -2d0*sjr/skr/sjk*CA*skm/sik/sim*ccBLO_lrkimk
c           Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_SS_gg_SC_ggq=M2_SS_gg_SC_ggq+pref*M2tmp*wss_sc_nnlo*extra
      enddo
c
c     Mapping 2 for B[krl,iml]
      call phase_space_CS_inv(k,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(m))cycle
         if(m.eq.i.or.m.eq.k.or.m.eq.j)cycle
c
            mb = real_mapped_labels(m)
            mbb = Born_mapped_labels(mb)
c
c           Mapping 2 for B[krl,iml]
c
c           underlying Born configuration is remapped
            call phase_space_CS_inv(ib,mb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c           call colour-connected Born matrix element with the mapping [krl,iml]
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
            ccBLO_krliml = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c           possible cuts
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c           invariant quantities: (c --> m)
            sim = xs(i,m)
            sjm = xs(j,m)
c
c           safety check
            if(sij.le.0d0.or.skr.le.0d0.or.sjk.le.0d0.or.sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_SS_gg_SC_ggq',sij, skr, sjk, sim
               goto 999
            endif
c
c           Double-soft soft-collinear kernel according to the eq.(C.14)
            M2tmp = -2d0*sjr/skr/sjk*(2d0*CF-CA)*sjm/sim/sij*ccBLO_krliml
c           Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_SS_gg_SC_ggq=M2_SS_gg_SC_ggq+pref*M2tmp*wss_sc_nnlo*extra
c
c           plot
            wgtpl=-pref*M2tmp*wss_sc_nnlo*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     apply flavour factor
      M2_SS_gg_SC_ggq = M2_SS_gg_SC_ggq * %(proc_prefix_rr)s_fl_factor
      if(test_sector_function) M2_SS_gg_SC_ggq = WSS_SC_NNLO
c
c     sanity check
      if(abs(M2_SS_gg_SC_ggq).ge.huge(1d0).or.isnan(M2_SS_gg_SC_ggq))then
         write(77,*)'Exception caught in M2_SS_SC_ggq',M2_SS_gg_SC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
