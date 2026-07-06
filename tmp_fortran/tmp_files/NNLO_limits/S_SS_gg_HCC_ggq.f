
      double precision function M2_S_SS_GGQ_HCC_GGQ(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) S(i,j) C(i,j,k) * (1 - SC(i,j,k)) * WS_SS_HCC
c     where  i, j is a g-g pair while k is a q (or qb)
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,r,ierr,nit,parent_leg
      integer ib,jb,kb,rb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO_ijr_jkr,BLO_ijr_krj,BLO_irj_jkr,BLO_irj_krj,BLO_ikr_jkr
      double precision BLO_ikr_jrk,BLO_irk_jkr,BLO_irk_jrk,BLO_ijk_jkr,BLO_ikj_jkr
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ans(0:NSQSO_BORN)
      double precision sij,sik,sjk,sir,sjr,skr,sbjk,sbjr,sbkr
      double precision Ei_jr,Ei_kr,Ei_jk
      double precision Ebj_kr_ijr,Ebj_kr_irj,Ebj_kr_ikr,Ebj_kr_irk,Ebj_kr_ijk,Ebj_kr_ikj
      integer, parameter :: hel = - 1
      logical flavourmatch
      double precision alphas,alpha_qcd
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_rr)s_fl_factor
      common/%(proc_prefix_rr)s_flavour_factor/%(proc_prefix_rr)s_fl_factor
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
      integer real_sc_mapped_labels(nexternal),Born_sc1_mapped_labels(nexternal-1),Born_sc2_mapped_labels(nexternal-1)
      common/c_NNLO_sc_mapped_labels/real_sc_mapped_labels,Born_sc1_mapped_labels,Born_sc2_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_S_SS_ggq_HCC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_S_SS_ggq_HCC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.leg_PDGs(j).and.abs(leg_PDGs(k)).le.5.and.leg_PDGs(i).ne.leg_PDGs(k)
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_S_SS_ggq_HCC_ggq',leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=(8d0*pi*alphas)**2
c
c     get mapped labels   ! TODO: fix mappings!!
      ib = real_mapped_labels(i)
      jb = real_mapped_labels(j)
      kb = real_mapped_labels(k)
      rb = real_mapped_labels(r)
c
c     invariant quantities
      sij  = xs(i,j)
      sjk  = xs(j,k)
      sik  = xs(i,k)
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
c
c     Global Eikonals
      Ei_jr = sjr/sij/sir
      Ei_kr = skr/sik/sir
      Ei_jk = sjk/sij/sik
c
c     safety check
      if(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0.or.zj.lt.0d0.or.zk.lt.0d0)then
        write(77,*)'Inaccuracy 1 in M2_S_SS_ggq_HCC_ggq',sij,sik,sjk,zj,zk
        goto 999
      endif
c
c     soft double-soft hard-doublecollinear sector function, (C.76-77) of 2212.11190
      call get_sig2(xs,1d0,nexternal)
      call get_ws_ss_hcc_nnlo(asec,bsec,csec,dsec)
c
c     mapping 1: [ijr,jkr]
      call phase_space_CS_inv(i,j,r,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
c
      Ebj_kr_ijr = sbkr/sbjk/sbjr
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc1_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ijr_kjr = ans(0)
c
c     mapping 2: [ijr,krj]
      call phase_space_CS_inv(i,j,r,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call phase_space_CS_inv(kb,rb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc1_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ijr_krj = ans(0)
c
c     mapping 3: [irj,jkr]
      call phase_space_CS_inv(i,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
c
      Ebj_kr_irj = sbkr/sbjk/sbjr
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irj_jkr = ans(0)
c
c     mapping 4: [irj,krj]
      call phase_space_CS_inv(i,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call phase_space_CS_inv(kb,rb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irj_krj = ans(0)
c
c     mapping 5: [ikr,jkr]
      call phase_space_CS_inv(i,k,r,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
c
      Ebj_kr_ikr = sbkr/sbjk/sbjr
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ikr_jkr = ans(0)
c
c     mapping 6: [ikr,jrk]
      call phase_space_CS_inv(i,k,r,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call phase_space_CS_inv(jb,rb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ikr_jrk = ans(0)
c
c     mapping 7: [irk,jkr]
      call phase_space_CS_inv(i,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
c
      Ebj_kr_irk = sbkr/sbjk/sbjr
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irk_jkr = ans(0)
c
c     mapping 8: [irk,jrk]
      call phase_space_CS_inv(i,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call phase_space_CS_inv(jb,rb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irk_jrk = ans(0)
c
c     mapping 9: [ijk,jkr]
      call phase_space_CS_inv(i,j,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
c
      Ebj_kr_ijk = sbkr/sbjk/sbjr
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ijk_jkr = ans(0)
c
c     mapping 10: [ikj,jkr]
      call phase_space_CS_inv(i,k,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_ggq_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
c
      Ebj_kr_ikj = sbkr/sbjk/sbjr
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ikj_jkr = ans(0)
c
c     soft double-soft hard-doublecollinear kernel, eq. (C.31)
      M2tmp = CA*Ei_jr*(Ebj_kr_ijr*(BLO_irj_jkr-BLO_ijr_krj)+Ebj_kr_irj*(BLO_irj_jkr-BLO_irj_krj))
      M2tmp = M2tmp + (2d0*CF-CA)*Ei_kr*(Ebj_kr_ikr*(BLO_ikr_jkr-BLO_ikr_jrk)+Ebj_kr_irk*(BLO_irk_jkr-BLO_irk_jrk))
      M2tmp = M2tmp + CA*Ei_jk*(Ebj_kr_ijk*BLO_ijk_jkr + Ebj_kr_ikj*BLO_ikj_jkr)
      M2tmp = CF*M2tmp*pref*ws_ss_hcc_nnlo*extra*%(proc_prefix_rr)s_fl_factor*xj
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_S_SS_ggq_HCC_ggq = M2tmp
c
      if(test_sector_function) M2_S_SS_ggq_HCC_ggq = ws_ss_hcc_nnlo
c
c     plot
      wgtpl=+M2_S_SS_ggq_HCC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_S_SS_ggq_HCC_ggq).ge.huge(1d0).or.isnan(M2_S_SS_ggq_HCC_ggq))then
         write(77,*)'Exception caught in M2_S_SS_ggq_HCC_ggq',M2_S_SS_ggq_HCC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
