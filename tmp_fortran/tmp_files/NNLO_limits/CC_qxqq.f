
      double precision function M2_CC_qxqq(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j,k) kernel times WCC: i, j are a q-qb pair with same flavour
c     while k is a q (or qb) with same flavour
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
      double precision zi,zj,zk,zij,zik,zjk
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
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_CC_qxqq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec .and. bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_CC_qxqq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).eq.abs(leg_PDGs(i))
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_CC_qxqq', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
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
      zik  = zi+zk
      zjk  = zj+zk
      zij  = zi+zj
c
c     safety check
      if(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0.or.zk.lt.0d0)then
        write(77,*)'Inaccuracy 1 in M2_CC_qxqq',sij,sik,sjk,zi,zj,zk
        goto 999
      endif
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ANS(0)
c
c     double-collinear kernel, using eq. (B.17) in eq. (B.14) of 2212.11190
c     the choice follows the order
c     1. IJK = QX Q Q
c     2. IJK = QX Q QX
c
      if (leg_pdgs(i).EQ.-leg_pdgs(k).and.leg_pdgs(i).eq.-leg_pdgs(j)) then
        M2tmp = CF*TR*(-sijk**2/(2d0*sik**2)*(sjk/sijk-sij/sijk+(zi-zk)/zik)**2+sijk/sik*(2d0*(zj-zi*zk)/zik+zik)-1d0/2d0)
        M2tmp = M2tmp + CF*TR*(-sijk**2/(2d0*sij**2)*(sjk/sijk-sik/sijk+(zi-zj) /zij)**2+sijk/sij*(2d0*(zk-zi*zj)/zij+zij)-1d0/2d0)
        M2tmp = M2tmp + (2d0*CF**2-CA*CF)*(-sijk**2*zi/(2d0*sik*sij)*(1d0+zi**2)/(zik*zij)+(sjk/sik+sjk/sij)+sijk/(2d0*sik)*((1d0+zi**2)/zij-2d0*zk/zik)+ sijk/(2d0*sij)*((1d0+zi**2)/zik-2d0*zj/zij ))
      else if (leg_pdgs(i).eq.leg_pdgs(k).and.leg_pdgs(i).eq.-leg_pdgs(j)) then
        M2tmp = CF*TR*(-sijk**2/(2d0*sij**2)*(sjk/sijk-sik/sijk+(zi-zj)/zij)**2+sijk/sij*(2d0*(zk-zi*zj)/zij+zij)-1d0/2d0)
        M2tmp = M2tmp +CF*TR*(-sijk**2/(2d0*sjk**2)*(sij/sijk-sik/sijk+(zk-zj)/zjk)**2+sijk/sjk*(2d0*(zi-zk*zj)/zjk+zjk)-1d0/2d0)
        M2tmp = M2tmp +CF*(2d0*CF-CA)*(-sijk**2*zj/(2d0*sjk*sij)*(1d0+zj**2)/zjk/zij+sik/sjk+sik/sij+sijk/2d0/sjk*((1d0+zj**2)/zij-2d0*zk/zjk)+sijk/2d0/sij*((1d0+zj**2)/zjk-2d0*zi/zij))
      endif
c
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
      M2_CC_qxqq = M2tmp*pref/sijk**2*xj*extra ! eq.(C.15)
      if(test_sector_function) M2_CC_qxqq = wcc_nnlo
c
c     plot
      wgtpl=-M2_CC_qxqq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_CC_qxqq).ge.huge(1d0).or.isnan(M2_CC_qxqq))then
         write(77,*)'Exception caught in M2_CC_qxqq',M2_CC_qxqq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
