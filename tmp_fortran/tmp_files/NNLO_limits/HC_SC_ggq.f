
      double precision function M2_HC_SC_GGQ(itmp,jtmp,ktmp,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,k) SC(i,j,k) - S(i) C(i,k) SC(i,j,k) kernel;
c     i, j is a g-g pair, k is a q (or qb)
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
      double precision M2tmp,M2_C_SC_ggq,M2_SC_SC_ggq,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision ccBLO_ijr_kmj
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sik,sjk,sjr,skr,sbim,sbjm,sbkm,sbik
      double precision zj,zk,Pjkr,Ebi_km,Ej_kr
      double precision pref,extra
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
      logical flavourmatch_ijkj,flavourmatch_ijjk
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
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_HC_SC_ggq=0d0
      M2_C_SC_ggq=0d0
      M2_SC_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology (only appears in ijkj)
      if(bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_HC_SC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch_ijjk = leg_PDGs(itmp).eq.leg_PDGs(jtmp).and.abs(leg_PDGs(ktmp)).le.5.and.leg_PDGs(itmp).ne.leg_PDGs(ktmp)
      flavourmatch_ijkj = leg_PDGs(itmp).eq.leg_PDGs(ktmp).and.abs(leg_PDGs(jtmp)).le.5.and.leg_PDGs(itmp).ne.leg_PDGs(jtmp)
      if(leg_pdgs(itmp).eq.0 .or. .not.(flavourmatch_ijjk.or.flavourmatch_ijkj)) then
        write(*,*) 'Flavour mismatch in M2_HC_SC_ggq',leg_PDGs(itmp),leg_PDGs(jtmp),leg_PDGs(ktmp)
        stop 1
      endif
c
c     rearrange the flavours to always land at ggq
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
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=-(8d0*pi*alphas)**2
c
c     compute wc_nlo
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(i,j,k,r)
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
      Ej_kr = skr/sjk/sjr
c
c     fill mapped labels [ijr,kmj]
      call fill_real_mapped_labels(i,j,leg_pdgs,real_leg_pdgs)
      do m=1,nexternal
         if(.not.isNNLOQCDparton(m)) cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k) cycle
         call fill_born_mapped_labels(real_labels(k,i,m),real_labels(m,i,m),real_leg_pdgs,born_leg_pdgs)
      enddo
c
c     mapping: [ijr,kmj]
      call phase_space_CS_inv(i,j,r,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,j))
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isnnloqcdparton(m))cycle
         if(m.eq.i.or.m.eq.j.or.m.eq.k)cycle
         ib = real_labels(i,i,j)
         jb = real_labels(j,i,j)
         kb = real_labels(k,i,j)
         mb = real_labels(m,i,j)
         kbb = born_labels(kb,kb,mb)
         mbb = born_labels(mb,kb,mb)
c
c     underlying Born configuration is remapped
         call phase_space_CS_inv(kb,mb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,kb,mb))
         if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
         if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
         call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
         if(ierr.eq.1)goto 999
c
c     call colour-connected Born matrix element
         call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
         ccBLO_ijr_kmj = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c     invariant quantities
         sbim = xsb(i,m)
         sbkm = xsb(k,m)
         sbik = xsb(i,k)
         Ebi_km = sbkm/sbik/sbim
c
c     collinear soft-collinear kernel, eq. (C.33)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = 2d0*Pjkr/sjk*Ebi_km*ccBLO_ijr_kmj
         M2tmp = M2tmp*pref*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_C_SC_ggq = M2_C_SC_ggq + M2tmp
c
c     collinear soft-collinear sector, eq. (C.78)
         call get_sig2(xsb,alpha_mod_bar,nexternal-1)
         map1=real_labels(csec,i,j)
         map2=real_labels(dsec,i,j)
         call get_wsbar_nlo(map1,map2)
         M2_C_SC_ggq=M2_C_SC_ggq*wc_nlo*wsbar_nlo
c
c     soft-collinear soft-collinear kernel, eq. (C.34)
c     TODO: some contributions are 0 for ee->jj
         M2tmp = 4d0*CF*Ej_kr*Ebi_km*ccBLO_ijr_kmj
         M2tmp = M2tmp*pref*wsbar_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj
         M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
         M2_SC_SC_ggq = M2_SC_SC_ggq + M2tmp
c
c     plot
         wgtpl=-M2tmp*wgt/nit*wgt_chan
         wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
         if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      M2_HC_SC_ggq = M2_C_SC_ggq - M2_SC_SC_ggq
c
      if(test_sector_function) M2_HC_SC_ggq = wsbar_nlo*wc_nlo-wsbar_nlo
c
c     sanity check
      if(abs(M2_HC_SC_ggq).ge.huge(1d0).or.isnan(M2_HC_SC_ggq))then
         write(77,*)'Exception caught in M2_HC_SC_ggq',M2_HC_SC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
