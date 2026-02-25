      double precision function M2_C_CC_qxqq(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j) C(i,j,k) kernel times WC_CC: i, j are a q-qb pair with same flavour
c     while k is a q (or qb) with any flavour
      use sectors2_module
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision dot
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision xpbsave(0:3,nexternal-1),xpbbsave(0:3,nexternal-2)
      double precision xsbsave(nexternal-1,nexternal-1),xsbbsave(nexternal-2,nexternal-2)
      double precision ktb(0:3),ktb2,kt(0:3),kt2,WCCC_NNLO
      double precision x,y,xinit
      double precision ans(0:NSQSO_BORN)
      double precision sij,sir,sjr,sbjk,sbjr,sbkr
      double precision zi,zj,zbj,zbk,zbki
      double precision Pij,Qij,Pbjk,Ebjkr
      double precision alphas,alpha_qcd
      double precision alphaz
      parameter(alphaz=2d0)
      integer nlo_mapped_labels(nexternal), nlo_mapped_flavours(nexternal)
      integer lo_mapped_labels(nexternal), lo_mapped_flavours(nexternal)
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
      double precision alphas,alpha_qcd
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
      M2_C_CC_qxqq=0d0
      M2tmp=0d0
      ierr=0
      xpbsave  = xpb
      xpbbsave = xpbb
      xsbsave  = xsb
      xsbbsave = xsbb
c
c     check sector topology
      if(dsec.ne.bsec .and. csec.ne.bsec) then
        write (*,*) 'Wrong topology in M2_C_CC_qxqqp',isec,jsec,ksec,lsec,asec,bsec,csec,dsec
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
c     reshuffle NLO momenta and labels according to real_leg_pdgs and check
      call reshuffle_momenta(nexternal,real_leg_pdgs,NLO_mapped_flavours,NLO_mapped_labels,xpbsave)
      call invariants_from_p(xpbsave,nexternal-1,xsbsave,ierr)
      if(ierr.eq.1)goto 999
c
c     reshuffle LO momenta and labels according to Born_leg_pdgs and check
      call reshuffle_momenta(nexternal-1,Born_leg_pdgs,LO_mapped_flavours,LO_mapped_labels,xpbbsave)
      call invariants_from_p(xpbbsave,nexternal-2,xsbbsave,ierr)
      if(ierr.eq.1)goto 999
c
c     possible cuts
      call get_underlying_pdgs(asec,bsec,csec,dsec,nexternal-1,real_leg_pdgs)
      call get_underlying_pdgs(asec,bsec,csec,dsec,nexternal-2,Born_leg_pdgs)
      if(docut(xpbb,nexternal-2,born_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
      zi   = sir/(sir+sjr)
      zj   = 1d0-zi
      jb = NLO_mapped_labels(j)
      kb = NLO_mapped_labels(k)
      rb = NLO_mapped_labels(r)
      sbjk = xsbsave(jb,kb)
      sbjr = xsbsave(jb,rb)
      sbkr = xsbsave(kb,rb)
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
c
c     calculate kt between i and j, as well as ktb between jb and kb
c     TODO: check if labels are fine after reshufflings
      kt(:) = zj*xp(:,i) - zi*xp(:,j) -(zj-zi)*sij/(sir+sjr)*xp(:,r)
      kt2 = -zi*zj*sij
      ktb(:) = zbk*xpbsave(:,jb) - zbj*xpbsave(:,kb) + (zbk-zbj)*sbjk/(sbjr+sbkr)*xpbsave(:,rb)
      ktb2 = -zbj*zbk*sbjk
c
c     safety checks
      IF(sij.lt.0d0.or.sir.lt.0d0.or.sjr.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_C_CC_qxqq',SIJ,SIR,SJR,ZI,ZJ
        GOTO 999
      ENDIF
      IF(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0.or.zbj.lt.0d0.or.zbk.lt.0d0)then
        WRITE(77,*)'Inaccuracy 2 in M2_C_CC_qxqq',SBJK,SBJR,SBKR,ZBJ,ZBK
        GOTO 999
      ENDIF
C
C     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbbsave,hel,alphas,ans)
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
      call get_wc_nlo(xs,i,j,r,alphaz,nexternal)
      M2tmp=M2tmp*wc_nlo
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wc_nlo(xsbsave,csec,dsec,rb,1d0,nexternal-1)
      M2tmp=M2tmp*wc_nlo
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_C_CC_qxqq = M2tmp*pref*xj*extra
c
c     plot
      wgtpl=-M2_C_CC_qxqq*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbbsave,xsbbsave,nexternal-2,born_leg_pdgs,wgtpl)
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
      use sectors2_module
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision dot
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision xpbsave(0:3,nexternal-1),xpbbsave(0:3,nexternal-2)
      double precision xsbsave(nexternal-1,nexternal-1),xsbbsave(nexternal-2,nexternal-2)
      double precision ktb(0:3),ktb2,kt(0:3),kt2,WSSCCC_NNLO
      double precision x,y,xinit
      double precision ans(0:NSQSO_BORN)
      double precision sij,sir,sjr,sbjk,sbjr,sbkr
      double precision zi,zj,zbj,zbk
      double precision Pij,Qij,Pbjk,Ebjkr
      integer nlo_mapped_labels(nexternal), nlo_mapped_flavours(nexternal)
      integer lo_mapped_labels(nexternal), lo_mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      logical flavourmatch
c     set logical doplot
      double precision alphas,alpha_qcd
      double precision alphaz
      parameter(alphaz=2D0)
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
      xpbsave  = xpb
      xpbbsave = xpbb
      xsbsave  = xsb
      xsbbsave = xsbb
c
c     check sector topology
      if(dsec.ne.bsec .and. bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_C_SS_qqx_CC_qxqqp',isec,jsec,ksec,lsec,asec,bsec,csec,dsec
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
c     reshuffle NLO momenta and labels according to real_leg_pdgs and check
      call reshuffle_momenta(nexternal,real_leg_pdgs,NLO_mapped_flavours,NLO_mapped_labels,xpbsave)
      call invariants_from_p(xpbsave,nexternal-1,xsbsave,ierr)
      if(ierr.eq.1)goto 999
c
c     reshuffle LO momenta and labels according to Born_leg_pdgs and check
      call reshuffle_momenta(nexternal-1,Born_leg_pdgs,LO_mapped_flavours,LO_mapped_labels,xpbbsave)
      call invariants_from_p(xpbbsave,nexternal-2,xsbbsave,ierr)
      if(ierr.eq.1)goto 999
c
c     possible cuts
      call get_underlying_pdgs(asec,bsec,csec,dsec,nexternal-1,real_leg_pdgs)
      call get_underlying_pdgs(asec,bsec,csec,dsec,nexternal-2,Born_leg_pdgs)
      if(docut(xpbb,nexternal-2,born_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
      zi   = sir/(sir+sjr)
      zj   = 1d0-zi
      jb = NLO_mapped_labels(j)
      kb = NLO_mapped_labels(k)
      rb = NLO_mapped_labels(r)
      sbjk = xsbsave(jb,kb)
      sbjr = xsbsave(jb,rb)
      sbkr = xsbsave(kb,rb)
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
c
c     calculate kt between i and j, as well as ktb between jb and kb
c     TODO: check if labels are fine after reshufflings
      kt(:) = zj*xp(:,i) - zi*xp(:,j) -(zj-zi)*sij/(sir+sjr)*xp(:,r)
      kt2 = -zi*zj*sij
      ktb(:) = zbk*xpbsave(:,jb) - zbj*xpbsave(:,kb) + (zbk-zbj)*sbjk/(sbjr+sbkr)*xpbsave(:,rb)
      ktb2 = -zbj*zbk*sbjk
c
c     safety checks
      IF(sij.lt.0d0.or.sir.lt.0d0.or.sjr.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_C_CC_qxqq',SIJ,SIR,SJR,ZI,ZJ
        GOTO 999
      ENDIF
      IF(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0.or.zbj.lt.0d0.or.zbk.lt.0d0)then
        WRITE(77,*)'Inaccuracy 2 in M2_C_CC_qxqq',SBJK,SBJR,SBKR,ZBJ,ZBK
        GOTO 999
      ENDIF
C
C     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbbsave,hel,alphas,ans)
      BLO = ANS(0)
c
c     collinear double-soft double-collinear kernel, eq. (C.41) of 2212.11190v2
      Pij = TR*(1d0-2d0*zi*zj)
      Qij = TR*2d0*zi*zj
      Pbjk = CF*(1d0+zbk**2)/zbj
      Ebjkr = sbkr/sbjk/sbjr
      M2TMP = 2d0*CF*Ebjkr*(Pij*-Qij*(-1d0+2d0*dot(kt,ktb)**2/kt2/ktb2))
      M2TMP = M2TMP/sij*BLO
c
c     compute soft-collinear triple-collinear sector function eq. (C.83) of 2212.11190v2
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wc_nlo(xsbsave,csec,dsec,rb,1d0,nexternal-1)
      M2TMP=M2TMP*wc_nlo
c
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_C_SS_qqx_CC_qxqq=M2tmp*pref*xj*extra
c
c     plot
      wgtpl=-M2_C_SS_qqx_CC_qxqq*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbbsave,xsbbsave,nexternal-2,born_leg_pdgs,wgtpl)
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


