
      double precision function M2_HC_SC_ggq_CC_ggq(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j) SC(i,j,k) C(k,i,j) * W1 - S(i) C(i,j) SC(i,j,k) C(k,i,j) * W2
c     i, k are a g-g while j is a q (or qb) with any flavour
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,r,ir,ierr,nit
      integer jb,kb,rb
      double precision pref,M2_C_SC_ggq_CC_ggq,M2_SC_SC_ggq_CC_ggq,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra
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
      double precision Pij,Ei_jr,Eb_kjr
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
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_HC_SC_ggq_CC_ggq=0d0
      M2_C_SC_ggq_CC_ggq=0d0
      M2_SC_SC_ggq_CC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology (only appears in ijkj)
      if(bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_HC_SC_ggq_CC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.leg_PDGs(k).and.leg_PDGs(i).eq.21.and.abs(leg_PDGs(j)).le.5
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_HC_SC_ggq_CC_ggq', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=(8d0*pi*alphas)**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
      zi = sir/(sir+sjr)
      zj = 1d0-zi
c
c     safety checks
      if(sij.lt.0d0.or.sir.lt.0d0.or.sjr.lt.0d0)then
        write(77,*)'Inaccuracy 1 in M2_HC_SC_ggq_CC_ggq',sij,sir,sjr
        goto 999
      endif
c
c     getting PDG's
      jb = real_mapped_labels(j)
      kb = real_mapped_labels(k)
      rb = real_mapped_labels(r)
      sbjk = xsb(jb,kb)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
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
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
c
c     collinear soft-collinear double-collinear kernel, eq. (C.43) of 2212.11190v2
      Pij = CF*(2d0*zj/zi+zi)
      Eb_kjr = sbjr/sbjk/sbkr
      M2tmp = 2d0*CF*Eb_kjr*Pij
      M2_C_SC_ggq_CC_ggq = M2tmp/sij*BLO
c
c     compute collinear soft-collinear double-collinear sector function eq. (C.86) of 2212.11190v2
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(i,j,ksec,r)
      M2_C_SC_ggq_CC_ggq=M2_C_SC_ggq_CC_ggq*wc_nlo
c
c     soft-collinear soft-collinear double-collinear kernel, eq. (C.44) of 2212.11190v2
      Ei_jr = sjr/sij/sir
      M2_SC_SC_ggq_CC_ggq = 4d0*CA*CF*Ei_jr*Eb_kjr*BLO
c
c     soft-collinear soft-collinear double-collinear sector function eq. (C.87) of 2212.11190v2 is 1
c
      M2_HC_SC_ggq_CC_ggq = M2_C_SC_ggq_CC_ggq - M2_SC_SC_ggq_CC_ggq
c
c     Including correct multiplicity factor
      M2_HC_SC_ggq_CC_ggq = M2_HC_SC_ggq_CC_ggq*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_HC_SC_ggq_CC_ggq = M2_HC_SC_ggq_CC_ggq*%(proc_prefix_rr)s_fl_factor
      M2_HC_SC_ggq_CC_ggq = M2_HC_SC_ggq_CC_ggq*pref*xj*extra
c
      if(test_sector_function) M2_HC_SC_ggq_CC_ggq = wc_nlo-1d0
c
c     plot
      wgtpl=-M2_HC_SC_ggq_CC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_HC_SC_ggq_CC_ggq).ge.huge(1d0).or.isnan(M2_HC_SC_ggq_CC_ggq))then
         write(77,*)'Exception caught in M2_HC_SC_ggq_CC_ggq', M2_HC_SC_ggq_CC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
