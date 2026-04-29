
      double precision function M2_SC_GGQ(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     SC(i,j,k) kernel times WSC
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,r,l,m,ierr,nit
      integer jb,lb,mb
      integer jbb,lbb,mbb
      logical isNLOmappedQCDparton(nexternal-1)
      logical isLOmappedQCDparton(nexternal-2)
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO,QUADBLO_mlml,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sil,sim,slm,sij,sjl,sjm,ml2,mm2,y,z,x,damp
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_rr)s_fl_factor
      common/%(proc_prefix_rr)s_flavour_factor/%(proc_prefix_rr)s_fl_factor
c     external
      integer get_color_dipole_index
      external get_color_dipole_index
      integer, parameter :: HEL = - 1
      double precision   %(proc_prefix_Born)s_GET_CCBLO
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
c      integer (proc_prefix_S_g)s_den
c      common/(proc_prefix_S_g)s_iden/(proc_prefix_S_g)s_den
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
      M2_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_SC_ggq', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     TODO: check the PDGs
c     get PDGs
      jb = real_mapped_labels(j)
      do l=1,nexternal
         if(l.eq.isec) cycle
          if(abs(leg_pdgs(l)).le.6.or.leg_pdgs(l).eq.21) isNLOmappedQCDparton(real_mapped_labels(l)) = .true.
      enddo
      do l=1,nexternal-1
         if(l.eq.jb) cycle
          if(abs(real_leg_pdgs(l)).le.6.or.real_leg_pdgs(l).eq.21) isLOmappedQCDparton(Born_mapped_labels(l)) = .true.
      enddo
c
c     call W soft-collinear, eq. (C.55) of 2212.11190
c     a small detail is that sig2 is always called with alpha=1 in the limit
c     the necessary sig2's are raised to the respective alpha in the soft-collinear sector functions
      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      call get_wsc_nnlo(asec,bsec,csec,dsec,iref)
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=32d0*pi**2*alphas**2
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(m))cycle
         if(m.eq.i.or.m.eq.j)cycle
         do l=1,nexternal
            if(.not.ISNNLOQCDPARTON(l))cycle
            if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c
            lb = real_mapped_labels(l)
            mb = real_mapped_labels(m)
            lbb = Born_mapped_labels(lb)
            mbb = Born_mapped_labels(mb)
c
c         check labels and pdgs
            if(.not.(isnlomappedqcdparton(lb).and.isnlomappedqcdparton(mb)))then
               write(*,*)'Wrong indices 1 in M2_SC_ggq',lb,mb
               stop
            endif
            if(.not.(islomappedqcdparton(lbb).and.islomappedqcdparton(mbb)))then
               write(*,*)'Wrong indices 2 in M2_SC_ggq',lbb,mbb
               stop
            endif
            if(leg_pdgs(l).ne.born_leg_pdgs(lbb).or.leg_pdgs(m).ne.born_leg_pdgs(mbb))then
               write(*,*)'Wrong indices 3 in M2_SC_ggq',l,m,lbb,mbb
               stop
            endif
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the double-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            call phase_space_CS_inv(jb,lb,mb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c     invariant quantities: (c,d) in the paper --> (m,l)
            sij = xs(i,j)
            sil = xs(i,l)
            sim = xs(i,m)
            sjl = xs(j,l)
            sjm = xs(j,m)
            slm = xs(l,m)
c
c     safety check
            if(sij.le.0d0.or.(sil+sjl).le.0d0.or.(sim+sjm).le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_SC_ggq',sij, sil+sjl, sim+sjm
               goto 999
            endif
c
c     call colour-connected Born
c     TODO: fix strings for the associated underlying Born
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_Born)s_GET_CCBLO(lbb,mbb)
c
c
c     eikonal
c     TODO: write the kernel C.13
c     See file K2_I2_G_v2.pdf in the DropBox directory
c     (c,d) -> (m,l) (verified)
            M2tmp = -2d0*CA*CCBLO
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_SC_ggq=M2_SC_ggq+pref*M2tmp*wsc_nnlo*extra
c
c     plot
            wgtpl=-pref*M2tmp*wsc_nnlo*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo
c
c     apply flavour factor
      M2_SC_ggq = M2_SC_ggq * %(proc_prefix_rr)s_fl_factor
c
c     sanity check
      if(abs(M2_SC_ggq).ge.huge(1d0).or.isnan(M2_SC_ggq))then
         write(77,*)'Exception caught in M2_SC_ggq',M2_SC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_SC_GGQ_CC_GGQ(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     SC(i,j,k)SC(i,j,k) kernel times WSC_CC: i, j is a g-g pair
c     while k is a q (or qb)
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
c
c     initialise
      M2_SC_ggq_CC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.leg_PDGs(j).and.abs(leg_PDGs(k)).le.5.and.leg_PDGs(i).ne.leg_PDGs(k)
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_SC_ggq_CC_ggq',leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
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
        write(77,*)'Inaccuracy 1 in M2_SC_ggq_CC_ggq',sij,sik,sjk,zi,zj,zk
        goto 999
      endif
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLO = ANS(0)
c
c     double-colinear soft-collinear kernel, eq. (C.18) of 2212.11190
c     TODO: write the kernel
      M2tmp = CF
      M2tmp = M2tmp/sjk
c
c     include triple-collinear soft-collinear sector function, eq. (C.62-C.64) of 2212.11190
c     a small detail is that sig2 is always called with alpha=1 in the limit
c     the necessary sig2's are raised to the respective alpha in the soft-collinear sector functions
      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      call get_wsc_cc_nnlo(asec,bsec,csec,dsec,iref)
      M2tmp=M2tmp*wsc_cc_nnlo
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_SC_ggq_CC_ggq = M2tmp*pref/sijk**2*xj*extra ! eq.(C.15)
c
c     plot
      wgtpl=-M2_SC_ggq_CC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_SC_ggq_CC_ggq).ge.huge(1d0).or.isnan(M2_SC_ggq_CC_ggq))then
         write(77,*)'Exception caught in M2_SC_ggq_CC_ggq',M2_SC_ggq_CC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
