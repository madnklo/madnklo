
      double precision function M2_SS_GG_SC_GGQ_CC_GGQ(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i,j)C(i,j,k)SC(i,j,k) kernel times WSS_CC_SC: i, j is a g-g pair
c     while k is a q (or qb)
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
      integer jbb,kbb,rbb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO_krj_irj,BLO_jrk_irk
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ans(0:NSQSO_BORN)
      double precision sij,sik,sjk,sir,sjr,skr,Ei_jr,Ei_kr,Ej_kr
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
      M2_SS_gg_SC_ggq_CC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_SS_gg_SC_ggq_CC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.leg_PDGs(j).and.abs(leg_PDGs(k)).le.5.and.leg_PDGs(i).ne.leg_PDGs(k)
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_SS_gg_SC_ggq_CC_ggq', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=(8d0*pi*alphas)**2
c
c     get mapped labels
      ib = real_sc_mapped_labels(i)
      jb = real_sc_mapped_labels(j)
      kb = real_sc_mapped_labels(k)
      rb = real_sc_mapped_labels(r)
c
c     invariant quantities
      sij  = xs(i,j)
      sjk  = xs(j,k)
      sik  = xs(i,k)
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
c
c     safety check
      if(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0)then
        write(77,*)'Inaccuracy 1 in M2_SS_gg_SC_ggq_CC_ggq',sij,sik,sjk
        goto 999
      endif
c
      Ei_jr = sjr/sij/sir
      Ei_kr = skr/sik/sir
      Ej_kr = skr/sjk/sjr
c
c     double-soft soft-collinear double-collinear sector function, (C.65-C.67) of 2212.11190
      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      call get_wss_sc_cc_nnlo(asec,bsec,csec,dsec,iref)
c
c     mapping 1: [krj,irj]
      call phase_space_CS_inv(k,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call phase_space_CS_inv(ib,rb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc1_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_krj_irj = ans(0)
c
c     mapping 2: [jrk,irk]
      call phase_space_CS_inv(j,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_sc_mapped_labels)
      call phase_space_CS_inv(ib,rb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_sc2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_jrk_irk = ans(0)
c
c     double-soft soft-collinear double-colinear kernel, (C.19)
      M2tmp = 2d0*CF*Ej_kr*(CA*Ei_jr*BLO_krj_irj+(2*CF-CA)*Ei_kr*BLO_jrk_irk)
      M2tmp = M2tmp*pref*wss_sc_cc_nnlo*extra*%(proc_prefix_rr)s_fl_factor*xj
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_SS_gg_SC_ggq_CC_ggq = M2tmp
c
      if(test_sector_function) M2_SS_gg_SC_ggq_CC_ggq = wss_sc_cc_nnlo
c
c     plot
      wgtpl=-M2_SS_gg_SC_ggq_CC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_SS_gg_SC_ggq_CC_ggq).ge.huge(1d0).or.isnan(M2_SS_gg_SC_ggq_CC_ggq))then
         write(77,*)'Exception caught in M2_SS_gg_SC_ggq_CC_ggq',M2_SS_gg_SC_ggq_CC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
