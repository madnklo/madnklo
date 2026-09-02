
      double precision function M2_S_SS_GG_HCC_GGQ(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
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
      integer i,j,k,a,b,r,ierr,nit,parent_leg
      integer ib,jb,kb,rb,ab,bb
      double precision pref,M2tmp,M2_S_SS_gg_CC_ggq,M2_S_SS_gg_CC_ggq_SC_ggq,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO_ira_jkr,BLO_irb_jkr,BLO_iba_jkr,BLO_ira_bra,BLO_irb_arb
      double precision sia,sab,sib,sir,sbr,sar,sbbr,sbar,sbab
      double precision Ei_ar,Ei_br,Ei_ab
      double precision Eba_br_ira,Eba_br_irb,Eba_br_iba
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ans(0:NSQSO_BORN)
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
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      integer real_s_sc_1_mapped_labels(nexternal),real_s_sc_2_mapped_labels(nexternal),Born_s_sc_1_mapped_labels(nexternal-1),Born_s_sc_2_mapped_labels(nexternal-1)
      common/c_nnlo_s_sc_mapped_labels/real_s_sc_1_mapped_labels,real_s_sc_2_mapped_labels,Born_s_sc_1_mapped_labels,Born_s_sc_2_mapped_labels
      integer real_hcc_ia_mapped_labels(nexternal),Born_hcc_ia_mapped_labels(nexternal-1)
      common/c_NNLO_hcc_mapped_labels/real_hcc_ia_mapped_labels,Born_hcc_ia_mapped_labels
      integer real_hcc_ib_mapped_labels(nexternal),Born_hcc_ib_mapped_labels(nexternal-1)
      common/c_NNLO_hcc_mapped_labels/real_hcc_ib_mapped_labels,Born_hcc_ib_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
c
c     initialise
      M2_S_SS_gg_HCC_ggq=0d0
      M2_S_SS_gg_CC_ggq=0d0
      M2_S_SS_gg_CC_ggq_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_S_SS_gg_HCC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.leg_PDGs(j).and.abs(leg_PDGs(k)).le.5.and.leg_PDGs(i).ne.leg_PDGs(k)
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_S_SS_gg_HCC_ggq',leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=(8d0*pi*alphas)**2
c
c
      a = csec
      b = dsec
c
c     invariant quantities
      sia  = xs(i,a)
      sab  = xs(a,b)
      sib  = xs(i,b)
      sir  = xs(i,r)
      sbr  = xs(b,r)
      sar  = xs(a,r)
c
c     Global Eikonals
      Ei_ar = sar/sia/sir
      Ei_br = sbr/sib/sir
      Ei_ab = sab/sib/sia
c
c     safety check
      if(sia.lt.0d0.or.sir.lt.0d0.or.sib.lt.0d0)then
        write(77,*)'Inaccuracy 1 in M2_S_SS_gg_HCC_ggq',sia,sir,sib
        goto 999
      endif
c
c     compute soft double-soft double-hardcollinear sector function, (C.75) of 2212.11190
      call get_sig2(xs,1d0,nexternal)
      call get_ws12_nlo(isec,jsec,ksec)
c
c     building S(i)S(i,j)C(i,j,k) according to eq.(17) on dropbox
c     mapping 1: [ira,jkr]
      call phase_space_CS_inv(i,r,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_hcc_ia_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_hcc_ia_mapped_labels(a)
      bb = real_hcc_ia_mapped_labels(b)
      jb = real_hcc_ia_mapped_labels(j)
      kb = real_hcc_ia_mapped_labels(k)
      rb = real_hcc_ia_mapped_labels(r)
      sbbr = xsb(bb,rb)
      sbar = xsb(ab,rb)
      sbab = xsb(ab,bb)
      if(sbbr.lt.0d0.or.sbar.lt.0d0.or.sbab.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_gg_HCC_ggq',sbbr,sbar,sbab
        goto 999
      endif
      Eba_br_ira = sbbr/sbab/sbar
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_hcc_ia_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ira_jkr = ans(0)
c
c     mapping 2: [irb,jkr]
      call phase_space_CS_inv(i,r,b,xp,xpb,nexternal,leg_PDGs,xjCS1,real_hcc_ib_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_hcc_ib_mapped_labels(a)
      bb = real_hcc_ib_mapped_labels(b)
      jb = real_hcc_ib_mapped_labels(j)
      kb = real_hcc_ib_mapped_labels(k)
      rb = real_hcc_ib_mapped_labels(r)
      sbbr = xsb(bb,rb)
      sbar = xsb(ab,rb)
      sbab = xsb(ab,bb)
      if(sbbr.lt.0d0.or.sbar.lt.0d0.or.sbab.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_gg_HCC_ggq',sbbr,sbar,sbab
        goto 999
      endif
      Eba_br_irb = sbbr/sbab/sbar
c
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_hcc_ib_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irb_jkr = ans(0)
c
c     mapping 3: [iba,jkr]
      call phase_space_CS_inv(i,b,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_hcc_ia_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_hcc_ia_mapped_labels(a)
      bb = real_hcc_ia_mapped_labels(b)
      jb = real_hcc_ia_mapped_labels(j)
      kb = real_hcc_ia_mapped_labels(k)
      rb = real_hcc_ia_mapped_labels(r)
      sbbr = xsb(bb,rb)
      sbar = xsb(ab,rb)
      sbab = xsb(ab,bb)
      if(sbbr.lt.0d0.or.sbar.lt.0d0.or.sbab.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_SS_gg_HCC_ggq',sbbr,sbar,sbab
        goto 999
      endif
      Eba_br_iba = sbbr/sbab/sbar
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_hcc_ia_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_iba_jkr = ans(0)
c
      M2_S_SS_gg_CC_ggq = 2d0*CF*(CA*Ei_ar*Eba_br_ira*BLO_ira_jkr+(2d0*CF-CA)*Ei_br*Eba_br_irb*BLO_irb_jkr+CA*Ei_ab*Eba_br_iba*BLO_iba_jkr)
c
c     building S(i)S(i,a)SC(i,a,b)C(i,j,k) according to eq.(19) on dropbox
c     mapping 1: [ira,bra]
      call phase_space_CS_inv(i,r,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_s_sc_1_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_s_sc_1_mapped_labels(a)
      bb = real_s_sc_1_mapped_labels(b)
      rb = real_s_sc_1_mapped_labels(r)
c
      call phase_space_CS_inv(bb,rb,ab,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_s_sc_1_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ira_bra = ans(0)
c
c     mapping 2: [irb,arb]
      call phase_space_CS_inv(i,r,b,xp,xpb,nexternal,leg_PDGs,xjCS1,real_s_sc_2_mapped_labels)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_s_sc_2_mapped_labels(a)
      bb = real_s_sc_2_mapped_labels(b)
      rb = real_s_sc_2_mapped_labels(r)
c
      call phase_space_CS_inv(ab,rb,bb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_s_sc_2_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irb_arb = ans(0)
c
      M2_S_SS_gg_CC_ggq_SC_ggq = 2d0*CF*(CA*Ei_ar*Eba_br_ira*BLO_ira_bra+(2d0*CF-CA)*Ei_br*Eba_br_irb*BLO_irb_arb)
c
c     soft double-soft double-hardcollinear kernel, TODO: Q_ij contribution is zero for ee->jj
      M2_S_SS_gg_HCC_ggq = M2_S_SS_gg_CC_ggq - M2_S_SS_gg_CC_ggq_SC_ggq
      M2_S_SS_gg_HCC_ggq = M2_S_SS_gg_HCC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj
      M2_S_SS_gg_HCC_ggq = M2_S_SS_gg_HCC_ggq*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
      if(test_sector_function) M2_S_SS_gg_HCC_ggq = ws12_nlo
c
      call ct_log('KS_SS_CC               ',M2_S_SS_gg_CC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den))
      call ct_log('KS_SS_CC_SC            ',M2_S_SS_gg_CC_ggq_SC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den))
c
c     plot
      wgtpl=+M2_S_SS_gg_HCC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_S_SS_gg_HCC_ggq).ge.huge(1d0).or.isnan(M2_S_SS_gg_HCC_ggq))then
         write(77,*)'Exception caught in M2_S_SS_gg_HCC_ggq',M2_S_SS_gg_HCC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
