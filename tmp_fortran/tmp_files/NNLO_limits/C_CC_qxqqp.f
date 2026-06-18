      double precision function M2_C_CC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
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
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ktb(0:3),ktb2,kt(0:3),kt2
      double precision x,y,xinit
      double precision ans(0:NSQSO_BORN)
      double precision sij,sir,sjr,sbjk,sbjr,sbkr
      double precision zi,zj,zbj,zbk
      double precision Pij,Qij,Pbjk,Ebjkr
      double precision dot
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
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
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_C_CC_qxqqp=0d0
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
        write(*,*) 'Flavour mismatch in M2_C_CC_qxqqp', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
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
      zi   = sir/(sir+sjr)
      zj   = 1d0-zi
c
      jb = real_mapped_labels(j)
      kb = real_mapped_labels(k)
      rb = real_mapped_labels(r)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      IF(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0) then
        WRITE(77,*)'Inaccuracy 2 in M2_C_CC_qxqqp',SBJK,SBJR,SBKR
        GOTO 999
      ENDIF
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
c
c     calculate kt between i and j, as well as ktb between jb and kb
      kt(:)  = zj*xp(:,i) - zi*xp(:,j) - (zj-zi)*sij/(sir+sjr)*xp(:,r)
      kt2    = -zi*zj*sij
      ktb(:) = zbk*xpb(:,jb) - zbj*xpb(:,kb) - (zbk-zbj)*sbjk/(sbjr+sbkr)*xpb(:,rb)
      ktb2   = -zbj*zbk*sbjk
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ANS(0)
c
c     collinear double-collinear kernel, eq. (C.39) of 2212.11190v2
      Pij = TR*(1d0-2d0*zi*zj)
      Qij = TR*2d0*zi*zj
      Pbjk = CF*(1d0+zbk**2)/(1d0-zbk)
      Ebjkr = sbkr/sbjk/sbjr
      M2tmp = Pij*Pbjk/sbjk-2d0*CF*Ebjkr*Qij*(-1d0+2d0*dot(kt,ktb)**2/kt2/ktb2)
      M2tmp = M2tmp/sij*BLO
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
      M2_C_CC_qxqqp = M2tmp*pref*xj*extra
      if(test_sector_function) M2_C_CC_qxqqp = WC_NLO*WCBAR_NLO
c
c     plot
      wgtpl=-M2_C_CC_qxqqp*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_C_CC_qxqqp).ge.huge(1d0).or.isnan(M2_C_CC_qxqqp))then
         write(77,*)'Exception caught in M2_C_CC_qxqqp',M2_C_CC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
