      double precision function M2_CC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j,k) kernel times WCC: i, j are a q-qb pair with same flavour
c     while k is a q (or qb) with different flavour
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
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ans(0:NSQSO_BORN)
      double precision sijk,sij,sik,sjk,sir,sjr,skr
      double precision zi,zj,zk,zij
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

c
c     initialise
      M2_CC_qxqqp=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec .and. bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_CC_qxqqp',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(i))
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_CC_qxqqp', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
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
      sjk  = xs(j,k)
      sik  = xs(i,k)
      sijk = sij+sik+sjk
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
      zi   = sir/(sir+sjr+skr)
      zj   = sjr/(sir+sjr+skr)
      zk   = skr/(sir+sjr+skr)
      zij  = zi+zj
c
c     safety check
      IF(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0.or.zk.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_CC_qxqqp',SIJ,SIK,SJK,ZI,ZJ,ZK
        GOTO 999
      ENDIF
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ANS(0)
c
c     double-collinear kernel, eq. (B.16) of 2212.11190
      M2tmp = CF*TR*(-sijk**2/(2D0*sij**2)*(sjk/sijk-sik/sijk+(zi-zj)/zij)**2+sijk/sij*(2D0*(zk-zi*zj)/zij+zij)-1D0/2D0)
      M2tmp = M2tmp*BLO
c
c     include double-collinear sector function
      call get_hatsignnlo(r,xs,nexternal)
      call get_wcc_nnlo(asec,bsec,csec,dsec)
      M2tmp=M2tmp*wcc_nnlo
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_CC_qxqqp = M2tmp*pref/sijk**2*xj*extra ! eq (C.15)
c
c     plot
      wgtpl=-M2_CC_qxqqp*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_CC_qxqqp).ge.huge(1d0).or.isnan(M2_CC_qxqqp))then
         write(77,*)'Exception caught in M2_CC_qxqqp',M2_CC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_SS_qqx_CC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i,j) C(i,j,k) kernel times WSSCC: i, j are a q-qb pair with same flavour
c     while k is a q (or qb) with different flavour
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
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ans(0:NSQSO_BORN)
      double precision sijk,sij,sik,sjk,sir,sjr,skr
      double precision zi,zj,zk,zij,eijkr
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
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
c
c     initialise
      M2_SS_qqx_CC_qxqqp=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec .and. bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_SS_qqx_CC_qxqqp',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(i))
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_SS_qqx_CC_qxqqp', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
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
      sjk  = xs(j,k)
      sik  = xs(i,k)
      sijk = sij+sik+sjk
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
      zi   = sir/(sir+sjr+skr)
      zj   = sjr/(sir+sjr+skr)
      zk   = skr/(sir+sjr+skr)
      zij  = zi+zj
c
c     safety check
      IF(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0.or.zk.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_SS_qqx_CC_qxqqp',SIJ,SIK,SJK,ZI,ZJ,ZK
        GOTO 999
      ENDIF
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ANS(0)
c
c     double-soft double-collinear kernel, eq. (C.16) of 2212.11190
      Eijkr = (1/sij**2)*((sik*zj+zi*sjk)/((sik+sjk)*zij)-sik*sjk/(sik+sjk)**2-zi*zj/zij**2)-zk/sij/(sik+sjk)/zij
      M2tmp = CF*(-2d0*TR*Eijkr)
      M2tmp = M2tmp*BLO
c
c     include double-soft double-collinear sector function
      call get_hatsignnlo(r,xs,nexternal)
      call get_wss_cc_nnlo(asec,bsec,csec,dsec)
      M2tmp=M2tmp*wss_cc_nnlo
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_SS_qqx_CC_qxqqp=M2tmp*pref*xj*extra ! eq (C.16)
c
c     plot
      wgtpl=+M2_SS_qqx_CC_qxqqp*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_SS_qqx_CC_qxqqp).ge.huge(1d0).or.isnan(M2_SS_qqx_CC_qxqqp))then
         write(77,*)'Exception caught in M2_SS_qqx_CC_qxqqp',M2_SS_qqx_CC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
