
      double precision function M2_S_SS_gg(i,j,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     S(i) S(i,k) kernel times sector function
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'leg_pdgs.inc'
      include 'nsqso_born.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,r,l,m,ierr,nit
      integer ib,jb,lb,mb
      integer jbb,lbb,mbb
      logical isNLOmappedQCDparton(nexternal-1)
      logical isLOmappedQCDparton(nexternal-2)
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision sblj,sbmj,sbml,Ei_jl,Ei_jm,Eb_ijl,Eb_jlm
      double precision BLO,ccBLO_imj_jml,ccBLO_ijl_jml,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sir,sjr,sim,sjm,sil,sjl,sml,ml2,mm2,y,z,x,damp
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
      M2_S_SS_gg=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_S_SS_gg',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_S_SS_gg', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=16d0*pi**2*alphas**2
c
c     get PDGs
      ib = real_mapped_labels(i)
      jb = real_mapped_labels(j)
c
c     Invariant quantities
      sij = xs(i,j)
c
c     compute soft double-soft sector function from eq. (C.72)
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
c
c     mapping 1: [imj,jml]
c     eikonal double sum (c,d) -> (m,l)
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(m))cycle
         if(m.eq.i.or.m.eq.j)cycle
         do l=1,nexternal
            if(.not.ISNNLOQCDPARTON(l))cycle
            if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c
            lb = real_mapped_labels(l)
            mb = real_mapped_labels(m)
            lbb = Born_mapped_labels(lb)
            mbb = Born_mapped_labels(mb)
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the double-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
c
            call phase_space_CS_inv(i,m,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
            sbml = xsb(mb,lb)
            sbmj = xsb(mb,jb)
            sblj = xsb(lb,jb)
            Eb_jlm = sml/sbmj/sblj
            call phase_space_CS_inv(jb,mb,lb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c     call wsbar
            call get_sig2(xsb,1d0,nexternal-1)
            map1=real_mapped_labels(csec)
            map2=real_mapped_labels(dsec)
            call get_wsbar_nlo(map1,map2)
c
c     invariant quantities: (c,d) in the paper --> (m,l)
            sim = xs(i,m)
            sml = xs(m,l)
            sjm = xs(j,m)
            Ei_jm = sjm/sij/sim
c
c     safety check
            if(sij.le.0d0.or.sml.le.0d0.or.sjm.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S_SS_gg',sij, sml, sjm
               goto 999
            endif
c
c     call colour-connected Born matrix element
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_imj_jml = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
c
c     eikonal (c,d) -> (m,l) eq. (C.26)
            M2tmp = M2tmp - 2*CA*Ei_jm*Eb_jlm*ccBLO_imj_jml
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_S_SS_gg=M2_S_SS_gg+pref*M2tmp*extra*ws_nlo*wsbar_nlo
c
c     plot
            wgtpl=-pref*M2tmp*extra*wgt/nit*wgt_chan*wsbar_nlo*ws_nlo
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo
c
c     mapping 2: [ijl,jml]
c     eikonal double sum (c,d) -> (m,l)
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(m))cycle
         if(m.eq.i.or.m.eq.j)cycle
         do l=1,nexternal
            if(.not.ISNNLOQCDPARTON(l))cycle
            if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c
            lb = real_mapped_labels(l)
            mb = real_mapped_labels(m)
            lbb = Born_mapped_labels(lb)
            mbb = Born_mapped_labels(mb)
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the double-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
c
            call phase_space_CS_inv(i,j,l,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
            sbml = xsb(mb,lb)
            sbmj = xsb(mb,jb)
            sblj = xsb(lb,jb)
            Eb_ijl = sml/sbmj/sblj
            call phase_space_CS_inv(jb,mb,lb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c     call wsbar
            call get_sig2(xsb,1d0,nexternal-1)
            map1=real_mapped_labels(csec)
            map2=real_mapped_labels(dsec)
            call get_wsbar_nlo(map1,map2)
c
c     invariant quantities: (c,d) in the paper --> (m,l)
            sil = xs(i,l)
            sjl = xs(j,l)
            Ei_jl = sjl/sij/sil
c
c     safety check
            if(sij.le.0d0.or.sjl.le.0d0.or.sil.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S_SS_gg',sij, sil, sjl
               goto 999
            endif
c
c     call colour-connected Born matrix element
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_ijl_jml = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
c
c     eikonal (c,d) -> (m,l) eq. (C.26)
            M2tmp = M2tmp - 2*CA*Ei_jl*Eb_ijl*ccBLO_ijl_jml
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_S_SS_gg=M2_S_SS_gg+pref*M2tmp*extra*ws_nlo*wsbar_nlo
c
c     plot
            wgtpl=-pref*M2tmp*extra*wgt/nit*wgt_chan*wsbar_nlo*ws_nlo
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo
c
c     apply flavour factor
      M2_S_SS_gg = M2_S_SS_gg * %(proc_prefix_rr)s_fl_factor
      if(test_sector_function) M2_S_SS_gg = wsbar_nlo*ws_nlo
c
c     sanity check
      if(abs(M2_S_SS_gg).ge.huge(1d0).or.isnan(M2_S_SS_gg))then
         write(77,*)'Exception caught in M2_S_SS_gg',M2_S_SS_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
