
      double precision function M2_S_HCC_ggq(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
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
      integer ib,jb,kb,ab,bb,rb,parent
      double precision pref,M2tmp,M2_S_CC_ggq,M2_S_SC_CC_ggq
      double precision wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO_ira_jkr,BLO_irb_jkr,BLO_iba_jkr,BLO_ira_bra,BLO_irb_arb
      double precision sia,sab,sib,sir,sbr,sar,sij,sik,sjk,sjr,skr
      double precision sbjk,sbjr,sbkr,zq
      double precision Ei_ar,Ei_br,Ei_jk,Ei_ab,rho_ab,rho_ba
      double precision Pjk_ira,Pjk_irb,Pjk_iba
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
      flavourmatch_ijjk = leg_PDGs(i).eq.leg_PDGs(j).and.abs(leg_PDGs(k)).le.5.and.leg_PDGs(i).ne.leg_PDGs(k)
      flavourmatch_ijkj = leg_PDGs(i).eq.leg_PDGs(k).and.abs(leg_PDGs(j)).le.5.and.leg_PDGs(i).ne.leg_PDGs(j)
      if(leg_pdgs(i).eq.0 .or. .not.(flavourmatch_ijjk.or.flavourmatch_ijkj)) then
        write(*,*) 'Flavour mismatch in M2_SC_ggq',leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
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
      if(leg_PDGs(a).eq.21.and.abs(leg_PDGs(b)).le.6)then
        rho_ab = rho_gq
        rho_ba = rho_qg
      elseif(leg_PDGs(b).eq.21.and.abs(leg_PDGs(a)).le.6)then
        rho_ab = rho_qg
        rho_ba = rho_gq
      endif
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
c     building S(i)C(i,j,k) according to eq.(15) on dropbox
c
c     mapping 1: [ira,jkr]
c
c     fill mapped labels
      call fill_real_mapped_labels(i,a,leg_pdgs,real_leg_pdgs)
      ab = real_labels(a,i,a)
      bb = real_labels(b,i,a)
      jb = real_labels(j,i,a)
      kb = real_labels(k,i,a)
      rb = real_labels(r,i,a)
      call fill_born_mapped_labels(jb,kb,real_leg_pdgs,born_leg_pdgs)
c
      call phase_space_CS_inv(i,r,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,a))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr+sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      if(leg_PDGs(j).eq.21.and.abs(leg_PDGs(k)).le.6)then
        zq = sbkr/(sbjr+sbkr)
      elseif(leg_PDGs(k).eq.21.and.abs(leg_PDGs(j)).le.6)then
        zq = sbjr/(sbjr+sbkr)
      endif
      Pjk_ira = CF*(2d0*zq/(1d0-zq)+(1d0-zq))/sbjk
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      call get_wcbar_nlo(ab,bb,rb)
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
      ab = real_labels(a,i,b)
      bb = real_labels(b,i,b)
      jb = real_labels(j,i,b)
      kb = real_labels(k,i,b)
      rb = real_labels(r,i,b)
      call fill_born_mapped_labels(jb,kb,real_leg_pdgs,born_leg_pdgs)
c
      call phase_space_CS_inv(i,r,b,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,b))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr+sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      if(leg_PDGs(j).eq.21.and.abs(leg_PDGs(k)).le.6)then
        zq = sbkr/(sbjr+sbkr)
      elseif(leg_PDGs(k).eq.21.and.abs(leg_PDGs(j)).le.6)then
        zq = sbjr/(sbjr+sbkr)
      endif
      Pjk_irb = CF*(2d0*zq/(1d0-zq)+(1d0-zq))/sbjk
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      call get_wcbar_nlo(ab,bb,rb)
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
      ab = real_labels(a,i,a)
      bb = real_labels(b,i,a)
      jb = real_labels(j,i,a)
      kb = real_labels(k,i,a)
      rb = real_labels(r,i,a)
c
      call phase_space_CS_inv(i,b,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,a))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr+sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      if(leg_PDGs(j).eq.21.and.abs(leg_PDGs(k)).le.6)then
        zq = sbkr/(sbjr+sbkr)
      elseif(leg_PDGs(k).eq.21.and.abs(leg_PDGs(j)).le.6)then
        zq = sbjr/(sbjr+sbkr)
      endif
      Pjk_iba = CF*(2d0*zq/(1d0-zq)+(1d0-zq))/sbjk
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      call get_wcbar_nlo(ab,bb,rb)
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
      M2_S_CC_ggq = CF*(rho_ab*Ei_ar*Pjk_ira*BLO_ira_jkr+rho_ba*Ei_br*Pjk_irb*BLO_irb_jkr-rhob_gq*Ei_jk*Pjk_iba*BLO_iba_jkr)
c
c     building S(i)SC(i,a,b)C(i,j,k) according to eq.(18) on dropbox
c
c     mapping 1: [ira,bra]
      ab = real_labels(a,i,a)
      bb = real_labels(b,i,a)
      jb = real_labels(j,i,a)
      kb = real_labels(k,i,a)
      rb = real_labels(r,i,a)
c
      call phase_space_CS_inv(i,r,a,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,a))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr+sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      if(leg_PDGs(j).eq.21.and.abs(leg_PDGs(k)).le.6)then
        zq = sbkr/(sbjr+sbkr)
      elseif(leg_PDGs(k).eq.21.and.abs(leg_PDGs(j)).le.6)then
        zq = sbjr/(sbjr+sbkr)
      endif
      Pjk_ira = CF*(2d0*zq/(1d0-zq)+(1d0-zq))/sbjk
c
      call phase_space_CS_inv(bb,rb,ab,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,bb,ab))
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      call get_wcbar_nlo(ab,bb,rb)
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_ira_bra = ans(0)*wcbar_nlo
c
c     mapping 2: [irb,arb]
      ab = real_labels(a,i,b)
      bb = real_labels(b,i,b)
      jb = real_labels(j,i,b)
      kb = real_labels(k,i,b)
      rb = real_labels(r,i,b)
c
      call phase_space_CS_inv(i,r,b,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,b))
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr+sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_S_HCC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      if(leg_PDGs(j).eq.21.and.abs(leg_PDGs(k)).le.6)then
        zq = sbkr/(sbjr+sbkr)
      elseif(leg_PDGs(k).eq.21.and.abs(leg_PDGs(j)).le.6)then
        zq = sbjr/(sbjr+sbkr)
      endif
      Pjk_irb = CF*(2d0*zq/(1d0-zq)+(1d0-zq))/sbjk
c
      call phase_space_CS_inv(ab,rb,bb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,born_labels(:,ab,bb))
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      call get_wcbar_nlo(ab,bb,rb)
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO_irb_arb = ans(0)*wcbar_nlo
c
      M2_S_SC_CC_ggq = CF*(rho_ab*Ei_ar*Pjk_ira*BLO_ira_bra+rho_ba*Ei_br*Pjk_irb*BLO_irb_arb)
c
c     soft double-hardcollinear kernel, eq. (C.30) TODO: Q_ij contribution is zero for ee->jj
      M2_S_CC_ggq = M2_S_CC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_S_SC_CC_ggq = M2_S_SC_CC_ggq*pref*ws12_nlo*extra*%(proc_prefix_rr)s_fl_factor*xj*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_S_HCC_ggq = M2_S_CC_ggq - M2_S_SC_CC_ggq
c
      if(test_sector_function)then
        M2_S_HCC_ggq = ws12_nlo*wcbar_nlo
      endif
c
      call ct_log('KS_CC                  ',M2_S_CC_ggq)
      call ct_log('KS_SC_CC               ',M2_S_SC_CC_ggq)
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
