
      double precision function M2_S_HCC_ggq(itmp,jtmp,ktmp,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i) C(i,j,k) * (1 - SC(i,j,k)) * WS_HCC
c     where  i, j is a g-g pair while k is a q (or qb)
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'nsqso_born.inc'
      include 'mapped_labels_common.inc'
      include 'leg_PDGs.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,r,a,b,ierr,nit,parent_leg,map1,map2
      integer ib,jb,kb,ab,bb,rb,itmp,jtmp,ktmp,parent
      double precision pref,M2tmp,M2_S_CC_ggq,M2_S_SC_CC_ggq,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO_ira_jkr,BLO_irb_jkr,BLO_iba_jkr,BLO_ira_bra,BLO_irb_arb
      double precision sia,sab,sib,sir,sbr,sar,sij,sik,sjk,sjr,skr,sbjk,sbjr,sbkr,zbj,zbk
      double precision Ei_ar,Ei_br,Ei_jk,Ei_ab,r1,r2,r3
      double precision Pb_jkr_ira,Pb_jkr_irb,Pb_jkr_iba
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ans(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      logical flavourmatch_ijjk,flavourmatch_ijkj
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
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
c
c     initialise
      M2_S_HCC_ggq=0d0
      M2_S_CC_ggq=0d0
      M2_S_SC_CC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_S_HCC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch_ijjk = leg_PDGs(itmp).eq.leg_PDGs(jtmp).and.abs(leg_PDGs(ktmp)).le.5.and.leg_PDGs(itmp).ne.leg_PDGs(ktmp)
      flavourmatch_ijkj = leg_PDGs(itmp).eq.leg_PDGs(ktmp).and.abs(leg_PDGs(jtmp)).le.5.and.leg_PDGs(itmp).ne.leg_PDGs(jtmp)
      if(leg_pdgs(itmp).eq.0 .or. .not.(flavourmatch_ijjk.or.flavourmatch_ijkj)) then
        write(*,*) 'Flavour mismatch in M2_SC_ggq',leg_PDGs(itmp),leg_PDGs(jtmp),leg_PDGs(ktmp)
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
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=(8d0*pi*alphas)**2
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
      Ei_jk = sab/sib/sia
c
c     safety check
      if(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0)then
        write(77,*)'Inaccuracy 1 in M2_S_HCC_ggq',sij,sik,sjk
        goto 999
      endif
c
c     compute soft double-hardcollinear sector function, (C.75) of 2212.11190
      call get_sig2(xs,1d0,nexternal)
      call get_ws12_nlo(isec,jsec,ksec)
c
c     call rho structures
c     r1 is rho_ab, r2 is rho_ba, r3 is rho_[ab] = rho_[ba]
      if(abs((leg_PDGs(a))).le.5.and.abs(leg_PDGs(b)).eq.21)then
        parent = a
      elseif(abs((leg_PDGs(a))).eq.21.and.abs(leg_PDGs(b)).le.5)then
        parent = b
      endif
      call get_rhos(a,b,parent,real_leg_pdgs,r1,r2,r3)
c
c     building S(i)C(i,j,k) according to eq.(15) on dropbox
c
c     mapping 1: [ira,jkr]
c
c     fill mapped labels
      call fill_real_mapped_labels(i,a,leg_pdgs,real_leg_pdgs)
      call fill_born_mapped_labels(real_labels(b,i,a),real_labels(a,i,a),real_leg_pdgs,born_leg_pdgs)
c
      call phase_space_CS_inv(i,r,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,a))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_labels(a,i,a)
      bb = real_labels(b,i,a)
      jb = real_labels(j,i,a)
      kb = real_labels(k,i,a)
      rb = real_labels(r,i,a)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
      Pb_jkr_ira = CF*(2d0*zbk/zbj+zbj)/sbjk
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,i,a)
      map2=real_labels(dsec,i,a)
      call get_wcbar_nlo(map1,map2,rb)
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,bb,ab))
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ira_jkr = ans(0)*wcbar_nlo
c
c     mapping 2: [irb,jkr]
c
c     fill mapped labels
      call fill_real_mapped_labels(i,b,leg_pdgs,real_leg_pdgs)
      call fill_born_mapped_labels(real_labels(a,i,b),real_labels(b,i,b),real_leg_pdgs,born_leg_pdgs)
c
      call phase_space_CS_inv(i,r,b,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,b))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_labels(a,i,b)
      bb = real_labels(b,i,b)
      jb = real_labels(j,i,b)
      kb = real_labels(k,i,b)
      rb = real_labels(r,i,b)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
      Pb_jkr_irb = CF*(2d0*zbk/zbj+zbj)/sbjk
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,i,b)
      map2=real_labels(dsec,i,b)
      call get_wcbar_nlo(map1,map2,rb)
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,ab,bb))
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irb_jkr = ans(0)*wcbar_nlo
c
c     mapping 3: [iba,jkr]
      call phase_space_CS_inv(i,b,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,a))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      jb = real_labels(j,i,a)
      kb = real_labels(k,i,a)
      rb = real_labels(r,i,a)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
      Pb_jkr_iba = CF*(2d0*zbk/zbj+zbj)/sbjk
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,i,a)
      map2=real_labels(dsec,i,a)
      call get_wcbar_nlo(map1,map2,rb)
c
      call phase_space_CS_inv(jb,kb,rb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,bb,ab))
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_iba_jkr = ans(0)*wcbar_nlo
c
      M2_S_CC_ggq = CF*(r1*Ei_ar*Pb_jkr_ira*BLO_ira_jkr+r2*Ei_br*Pb_jkr_irb*BLO_irb_jkr-r3*Ei_jk*Pb_jkr_iba*BLO_iba_jkr)
c
c     building S(i)SC(i,a,b)C(i,j,k) according to eq.(18) on dropbox
c
c     mapping 1: [ira,bra]
c
      call phase_space_CS_inv(i,r,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,a))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_labels(a,i,a)
      bb = real_labels(b,i,a)
      jb = real_labels(j,i,a)
      kb = real_labels(k,i,a)
      rb = real_labels(r,i,a)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
      Pb_jkr_ira = CF*(2d0*zbk/zbj+zbj)/sbjk
c
      call phase_space_CS_inv(bb,rb,ab,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,bb,ab))
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,i,a)
      map2=real_labels(dsec,i,a)
      call get_wcbar_nlo(map1,map2,rb)
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ira_bra = ans(0)*wcbar_nlo
c
c     mapping 2: [irb,arb]
      call phase_space_CS_inv(i,r,b,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,b))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      ab = real_labels(a,i,b)
      bb = real_labels(b,i,b)
      jb = real_labels(j,i,b)
      kb = real_labels(k,i,b)
      rb = real_labels(r,i,b)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
      Pb_jkr_irb = CF*(2d0*zbk/zbj+zbj)/sbjk
c
      call phase_space_CS_inv(ab,rb,bb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,ab,bb))
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,i,b)
      map2=real_labels(dsec,i,b)
      call get_wcbar_nlo(map1,map2,rb)
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irb_arb = ans(0)*wcbar_nlo
c
      M2_S_SC_CC_ggq = CF*(r1*Ei_ar*Pb_jkr_ira*BLO_ira_bra+r2*Ei_br*Pb_jkr_irb*BLO_irb_arb)
c
c     soft double-hardcollinear kernel, eq. (C.30) TODO: Q_ij contribution is zero for ee->jj
      M2_S_HCC_ggq = M2_S_CC_ggq - M2_S_SC_CC_ggq
      M2_S_HCC_ggq = M2_S_HCC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj
      M2_S_HCC_ggq = M2_S_HCC_ggq*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
      if(test_sector_function) M2_S_HCC_ggq = ws12_nlo*wcbar_nlo
c
      call ct_log('KS_CC                  ',M2_S_CC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den))
      call ct_log('KS_SC_CC               ',M2_S_SC_CC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den))
c
c     plot
      wgtpl=+M2_S_HCC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_S_HCC_ggq).ge.huge(1d0).or.isnan(M2_S_HCC_ggq))then
         write(77,*)'Exception caught in M2_S_HCC_ggq',M2_S_HCC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
