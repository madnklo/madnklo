
      double precision function M2_HC_CC_GGQ(itmp,jtmp,ktmp,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j) C(i,j,k) * Wcollinear-doublecollinear - S(i) C(i,j) C(i,j,k) * Wsoftcollinear-doublecollinear
c     i, j are a g-g pair while k is a q (or qb) with any flavour
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
      integer i,j,k,r,ierr,nit,itmp,jtmp,ktmp
      integer jb,kb,rb
      double precision pref,M2_C_CC_ggq,M2_SC_CC_ggq,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,KKBLO
      double precision dot
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ktb(0:3),ktb2,kt(0:3),kt2
      double precision x,y,xinit
      double precision ans(0:NSQSO_BORN)
      double precision sij,sir,sjr,sbjk,sbjr,sbkr
      double precision zi,zj,zbj,zbk,zbki
      double precision Pij,Qij,Pb_jk,Ei_jr,Eb_jkr,Eb_kjr
      double precision alphas,alpha_qcd
      integer, parameter :: hel = - 1
      logical flavourmatch_ijjk,flavourmatch_ijkj
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      double precision %(proc_prefix_Born)s_GET_KKBLO
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
      integer map1,map2,parent_leg
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
      M2_HC_CC_ggq=0d0
      M2_C_CC_ggq=0d0
      M2_SC_CC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec .and. bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_HC_CC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch_ijjk = leg_PDGs(itmp).eq.leg_PDGs(jtmp).and.abs(leg_PDGs(ktmp)).le.5.and.leg_PDGs(itmp).ne.leg_PDGs(ktmp)
      flavourmatch_ijkj = leg_PDGs(itmp).eq.leg_PDGs(ktmp).and.abs(leg_PDGs(jtmp)).le.5.and.leg_PDGs(itmp).ne.leg_PDGs(jtmp)
      if(.not.(flavourmatch_ijjk.or.flavourmatch_ijkj)) then
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
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=(8d0*pi*alphas)**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
      zi   = sir/(sir+sjr)
      zj   = 1d0-zi
c
c     fill mapped labels
      call fill_real_mapped_labels(i,j,leg_pdgs,real_leg_pdgs)
      call fill_born_mapped_labels(real_labels(j,i,j),real_labels(k,i,j),real_leg_pdgs,born_leg_pdgs)
      jb = real_labels(j,i,j)
      kb = real_labels(k,i,j)
      rb = real_labels(r,i,j)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        write(77,*)'Inaccuracy 2 in M2_HC_CC_ggq',sbjk,sbjr,sbkr
        goto 999
      endif
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
      parent_leg = born_labels(kb,jb,kb)
c
      kt(:) = zj*xp(:,i) - zi*xp(:,j) -(zj-zi)*sij/(sir+sjr)*xp(:,r)
      kt2 = dot(kt,kt)
      ktb(:) = zbk*xpb(:,jb) - zbj*xpb(:,kb) + (zbk-zbj)*sbjk/(sbjr+sbkr)*xpb(:,rb)
      ktb2 = dot(ktb,ktb)
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ans(0)
c
      KKBLO = %(proc_prefix_Born)s_GET_KKBLO(parent_leg,xpbb,kt)
c
c     collinear double-collinear kernel, eq. (C.39) of 2212.11190v2
c     (same as dropbox eq. 23) (term 2 is zero for k=q)
      Pij = 2d0*CA*(zi/zj+zj/zi+zi*zj)
      Qij = -2d0*CA*zi*zj
      Pb_jk = CF*(2d0*zbk/zbj+zbj)
      Eb_jkr = sbkr/sbjk/sbjr
      Eb_kjr = sbjr/sbjk/sbkr
      M2tmp = Pij*Pb_jk/sbjk*BLO - 2d0*CF*Eb_jkr*Qij*(-1d0+2d0*dot(kt,ktb)**2/kt2/ktb2)*BLO
      M2_C_CC_ggq = M2tmp/sij
c
c     compute collinear double-collinear sector function eq. (C.82) of 2212.11190v2
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(i,j,k,r)
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,i,j)
      map2=real_labels(dsec,i,j)
      call get_wcbar_nlo(map1,map2,rb)
      M2_C_CC_ggq=M2_C_CC_ggq*wc_nlo*wcbar_nlo
c
c     soft-collinear double-collinear kernel, eq. (C.40) of 2212.11190v2
c     (same as dropbox eq. 31)
      Ei_jr = sjr/sij/sir
      M2_SC_CC_ggq = 2d0*CA*Ei_jr*Pb_jk/sbjk*BLO
c
c     compute soft-collinear double-collinear sector function eq. (C.83) of 2212.11190v2
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,i,j)
      map2=real_labels(dsec,i,j)
      call get_wcbar_nlo(map1,map2,rb)
      M2_SC_CC_ggq=M2_SC_CC_ggq*wcbar_nlo
c
      M2_HC_CC_ggq = M2_C_CC_ggq - M2_SC_CC_ggq
c
c     include correct multiplicity and flavour factors
      M2_HC_CC_ggq = M2_HC_CC_ggq*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_HC_CC_ggq = M2_HC_CC_ggq*%(proc_prefix_rr)s_fl_factor
      M2_HC_CC_ggq = M2_HC_CC_ggq*pref*xj*extra
c
      if(test_sector_function) M2_HC_CC_ggq = wc_nlo*wcbar_nlo - wcbar_nlo
c
      call ct_log('KC_CC                  ',M2_C_CC_ggq*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)*%(proc_prefix_rr)s_fl_factor*pref*xj*extra)
      call ct_log('KS_C_CC                ',M2_SC_CC_ggq*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)*%(proc_prefix_rr)s_fl_factor*pref*xj*extra)
c
c     plot
      wgtpl=-M2_HC_CC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_HC_CC_ggq).ge.huge(1d0).or.isnan(M2_HC_CC_ggq))then
         write(77,*)'Exception caught in M2_HC_CC_ggq',M2_HC_CC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
