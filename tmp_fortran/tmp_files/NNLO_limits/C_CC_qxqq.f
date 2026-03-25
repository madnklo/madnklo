      double precision function M2_C_CC_qxqq(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j) C(i,j,k) kernel times WC_CC: i, j are a q-qb pair with same flavour
c     while k is a q (or qb) with any flavour
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer i,j,k,r,ierr,nit
      integer jb,kb,rb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision dot
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ktb(0:3),ktb2,kt(0:3),kt2
      double precision x,y,xinit
      double precision ans(0:NSQSO_BORN)
      double precision sij,sir,sjr,sbjk,sbjr,sbkr
      double precision zi,zj,zbj,zbk,zbki
      double precision Pij,Qij,Pbjk,Ebjkr
      double precision alphas,alpha_qcd
      integer, parameter :: hel = - 1
      logical flavourmatch
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
      integer map1,map2,parent_leg
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
c
c     initialise
      M2_C_CC_qxqq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec .and. bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_C_CC_qxqqp',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_C_CC_qxqq', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
c
c     safety checks
      IF(sij.lt.0d0.or.sir.lt.0d0.or.sjr.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_C_CC_qxqqp',SIJ,SIR,SJR
        GOTO 999
      ENDIF
c
c
      zi   = sir/(sir+sjr)
      zj   = 1d0-zi
c
c     check reshuffled real flavour -> not needed anymore?
c      if(real_leg_pdgs(j).ne.21)then
c         write(*,*) 'Wrong parent particle label 1 in M2_C_CC_qxqqp', j, real_leg_pdgs(j)
c         stop
c      endif
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      jb = real_mapped_labels(j)
      kb = real_mapped_labels(k)
      rb = real_mapped_labels(r)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      IF(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        WRITE(77,*)'Inaccuracy 2 in M2_C_CC_qxqq',SBJK,SBJR,SBKR
        GOTO 999
      ENDIF
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
      parent_leg = real_mapped_labels(jb)
c
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     calculate kt between i and j, as well as ktb between jb and kb
c     TODO: check if labels are fine after reshufflings
      kt(:) = zj*xp(:,i) - zi*xp(:,j) -(zj-zi)*sij/(sir+sjr)*xp(:,r)
      kt2 = -zi*zj*sij
      ktb(:) = zbk*xpb(:,jb) - zbj*xpb(:,kb) + (zbk-zbj)*sbjk/(sbjr+sbkr)*xpb(:,rb)
      ktb2 = -zbj*zbk*sbjk
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ANS(0)
c
c     collinear double-collinear kernel, eq. (C.39) of 2212.11190v2
      Pij = TR*(1d0-2d0*zi*zj)
      Qij = TR*2d0*zi*zj
      Pbjk = CF*(1d0+zbk**2)/zbj
      Ebjkr = sbkr/sbjk/sbjr
      M2TMP = Pij*Pbjk/sbjk-2d0*CF*Ebjkr*Qij*(-1d0+2d0*dot(kt,ktb)**2/kt2/ktb2)
      M2TMP = M2TMP/sij*BLO
c
c     compute collinear triple-collinear sector function eq. (C.82) of 2212.11190v2
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(i,j,ksec,r)
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wcbar_nlo(map1,map2,rb)
      M2tmp=M2tmp*wc_nlo*wcbar_nlo
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_C_CC_qxqq = M2tmp*pref*xj*extra
c
c     plot
      wgtpl=-M2_C_CC_qxqq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_C_CC_qxqq).ge.huge(1d0).or.isnan(M2_C_CC_qxqq))then
         write(77,*)'Exception caught in M2_C_CC_qxqq',M2_C_CC_qxqq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_C_SS_qqx_CC_qxqq(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C_(i,j) S(i,j) C(i,j,k) kernel times WSS_C_CC: i, j are a q-qb pair with same flavour
c     while k is a q (or qb) with any flavour
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer i,j,k,r,ierr,nit
      integer jb,kb,rb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision dot
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ktb(0:3),ktb2,kt(0:3),kt2
      double precision x,y,xinit
      double precision ans(0:NSQSO_BORN)
      double precision sij,sir,sjr,sbjk,sbjr,sbkr
      double precision zi,zj,zbj,zbk
      double precision Pij,Qij,Pbjk,Ebjkr
      integer, parameter :: hel = - 1
      logical flavourmatch
c     set logical doplot
      double precision alphas,alpha_qcd
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
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
c
c     initialise
      M2_C_SS_qqx_CC_qxqq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology(only appears in ijjk)
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_C_SS_qqx_CC_qxqq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_C_SS_qqx_CC_qxqq', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
c
c     safety checks
      if(sij.lt.0d0.or.sir.lt.0d0.or.sjr.lt.0d0)then
        write(77,*)'Inaccuracy 1 in m2_c_ss_qqx_cc_qxqq',sij,sir,sjr
        goto 999
      endif
      zi = sir/(sir+sjr)
      zj = 1d0-zi
      jb = real_mapped_labels(j)
      kb = real_mapped_labels(k)
      rb = real_mapped_labels(r)
      sbjk = xsb(jb,kb)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      if(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
         write(77,*)'Inaccuracy 2 in m2_c_ss_qqx_cc_qxqq',sbjk,sbjr,sbkr
         goto 999
      endif
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
c
c     calculate kt between i and j, as well as ktb between jb and kb
c     TODO: check if labels are fine after reshufflings
      kt(:) = zj*xp(:,i) - zi*xp(:,j) -(zj-zi)*sij/(sir+sjr)*xp(:,r)
      kt2 = -zi*zj*sij
      ktb(:) = zbk*xpb(:,jb) - zbj*xpb(:,kb) + (zbk-zbj)*sbjk/(sbjr+sbkr)*xpb(:,rb)
      ktb2 = -zbj*zbk*sbjk
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ANS(0)
c
c     collinear double-soft double-collinear kernel, eq. (C.41) of 2212.11190v2
      Pij = TR*(1d0-2d0*zi*zj)
      Qij = TR*2d0*zi*zj
      Pbjk = CF*(1d0+zbk**2)/zbj
      Ebjkr = sbkr/sbjk/sbjr
      M2TMP = 2d0*CF*Ebjkr*(Pij-Qij*(-1d0+2d0*dot(kt,ktb)**2/kt2/ktb2))
      M2TMP = M2TMP/sij*BLO
c
c     compute soft-collinear triple-collinear sector function eq. (C.84) of 2212.11190v2
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(i,j,ksec,r)
      M2TMP=M2TMP*wc_nlo
c
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_C_SS_qqx_CC_qxqq=M2tmp*pref*xj*extra
c
c     plot
      wgtpl=-M2_C_SS_qqx_CC_qxqq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_C_SS_qqx_CC_qxqq).ge.huge(1d0).or.isnan(M2_C_SS_qqx_CC_qxqq))then
         write(77,*)'Exception caught in M2_C_SS_qqx_CC_qxqq', M2_C_SS_qqx_CC_qxqq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


