
      double precision function M2_SS_GG_SC_GGQ(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i,j)SC(i,j,k) kernel times WSS_SC
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
      integer i,j,k,r,m,ierr,nit
      integer ib,jb,mb,rb,kb
      integer jbb,mbb,rbb,kbb
      logical isNLOmappedQCDparton(nexternal-1)
      logical isLOmappedQCDparton(nexternal-2)
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLOlrkimk,ccBLOkrliml,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sij,sjk,sjm,sik,sim,sjr,skr,skm,xk,xl,ml2,mm2,y,z,x,damp
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
      M2_SS_gg_SC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_SS_gg_SC_ggq', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     check sector topology
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_SS_gg_SC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     call W double-soft soft-collinear, eq. (C.57-C.58) of 2212.11190
c     a small detail is that sig2 is always called with alpha=1 in the limit
c     the necessary sig2's are raised to the respective alpha in the soft-collinear sector functions
      call get_sig2(xs,1d0,nexternal)
      call get_wss_sc_nnlo(asec,bsec,csec,dsec)
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=32d0*pi**2*alphas**2
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(m))cycle
         if(m.eq.i.or.m.eq.k.or.m.eq.j)cycle
c
            mb = real_mapped_labels(m)
            mbb = Born_mapped_labels(mb)
c
c           Mapping 1 for B[lrk,imk]
c
c           get PDGs
            ib = real_mapped_labels(i)
            rb = real_mapped_labels(r)
            kb = real_mapped_labels(k)
            rbb = Born_mapped_labels(rb)
            kbb = Born_mapped_labels(kb)
c
c           underlying Born configuration is remapped
            call phase_space_CS_inv(j,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            call phase_space_CS_inv(ib,mb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c           call colour-connected Born matrix element with the mapping [lrk,imk]
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLOlrkimk = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)

c
c           Mapping 2 for B[krl,icl]
c
c           get PDGs
            ib = real_mapped_labels(i)
            jb = real_mapped_labels(j)
            kb = real_mapped_labels(k)
            rbb = Born_mapped_labels(rb)
            kbb = Born_mapped_labels(kb)
c
c           underlying Born configuration is remapped
            call phase_space_CS_inv(k,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            call phase_space_CS_inv(ib,mb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
            if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c           call colour-connected Born matrix element with the mapping [krl,iml]
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
            ccBLOkrliml = %(proc_prefix_Born)s_GET_CCBLO(kbb,mbb)
c
c           possible cuts
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c           invariant quantities: (c --> m)
            sij = xs(i,j)
            sjk = xs(j,k)
            sjm = xs(j,m)
            sik = xs(i,k)
            sik = xs(i,k)
            sim = xs(i,m)
            skr = xs(k,r)
            sjr = xs(j,r)
            skm = xs(k,m)
c
c           safety check
            if(sij.le.0d0.or.skr.le.0d0.or.sjk.le.0d0.or.sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_SS_gg_SC_ggq',sij, skr, sjk, sim
               goto 999
            endif
c
c           Double-soft soft-collinear kernel according to the eq.(C.14)
            M2tmp = -2d0*sjr/skr/sjk*(CA*skm/sik/sim*ccBLOlrkimk+(2d0*CF-CA)*sjm/sim/sij*ccBLOkrliml)
c           Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_SS_gg_SC_ggq=M2_SS_gg_SC_ggq+pref*M2tmp*wss_sc_nnlo*extra
c
c           plot
            wgtpl=-pref*M2tmp*wss_sc_nnlo*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
c     apply flavour factor
      M2_SS_gg_SC_ggq = M2_SS_gg_SC_ggq * %(proc_prefix_rr)s_fl_factor
c
c     sanity check
      if(abs(M2_SS_gg_SC_ggq).ge.huge(1d0).or.isnan(M2_SS_gg_SC_ggq))then
         write(77,*)'Exception caught in M2_SS_SC_ggq',M2_SS_gg_SC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_SS_GG_SC_GGQ_CC_GGQ(i,j,k,r,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     S(i,j)C(i,j,k)SC(i,j,k) kernel times WSS_CC_SC: i, j is a g-g pair
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
      integer ib,jb,kb,rb
      integer jbb,kbb,rbb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjb,extra,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLOkrjirj,BLOjrkirk
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
      M2_SS_gg_SC_ggq_CC_ggq=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(bsec.ne.csec) then
        write (*,*) 'Wrong topology in M2_SS_gg_SC_ggq_CC_ggq',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.leg_PDGs(j).and.abs(leg_PDGs(k)).le.5.and.leg_PDGs(i).ne.leg_PDGs(k)
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_SS_gg_SC_ggq_CC_ggq', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     possible cuts
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
c
c     Mapping 1 for B[krj,irj]
c
c     get PDGs
      ib = real_mapped_labels(i)
      rb = real_mapped_labels(r)
      jb = real_mapped_labels(j)
      rbb = Born_mapped_labels(rb)
      jbb = Born_mapped_labels(jb)
c
c     underlying Born configuration is remapped
      call phase_space_CS_inv(k,r,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
      call phase_space_CS_inv(ib,rb,jb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element with the mapping [krj,irj]
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLOkrjirj = ANS(0)
c
c     Mapping 2 for B[jrk,irk]
c
c     get PDGs
      ib = real_mapped_labels(i)
      rb = real_mapped_labels(r)
      kb = real_mapped_labels(k)
      rbb = Born_mapped_labels(rb)
      kbb = Born_mapped_labels(kb)
c
c     underlying Born configuration is remapped
      call phase_space_CS_inv(j,r,k,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
      call phase_space_CS_inv(ib,rb,kb,xpb,xpbb,nexternal-1,real_leg_PDGs,xjCS2,Born_mapped_labels)
      if(xjCS1.eq.0d0.or.xjCS2.eq.0d0)goto 999
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     call Born matrix element with the mapping [jrk,irk]
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
      BLOjrkirk = ANS(0)
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
        write(77,*)'Inaccuracy 1 in M2_SS_gg_SC_ggq_CC_ggq',sij,sik,sjk,zi,zj,zk
        goto 999
      endif
c
c     double-soft double-colinear soft-collinear kernel, eq. (C.19) of 2212.11190
      M2tmp = 2d0*CF*skr/sjk/sjr*(CA*sjr/sij/sir*BLOkrjirj+(2*CF-CA)*skr/sik/sir*BLOjrkirk)
c
c     include double-soft double-collinear soft-collinear sector function, eq.(C.65-C.67) of 2212.11190
c     a small detail is that sig2 is always called with alpha=1 in the limit
c     the necessary sig2's are raised to the respective alpha in the soft-collinear sector functions
      call get_w(xs,nexternal)
      call get_sig2(xs,1d0,nexternal)
      call get_wss_sc_cc_nnlo(asec,bsec,csec,dsec,iref)
      M2tmp=M2tmp*wss_sc_cc_nnlo
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_SS_gg_SC_ggq_CC_ggq = M2tmp*pref*xj*extra
c
c     plot
      wgtpl=-M2_SS_gg_SC_ggq_CC_ggq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_SS_gg_SC_ggq_CC_ggq).ge.huge(1d0).or.isnan(M2_SS_gg_SC_ggq_CC_ggq))then
         write(77,*)'Exception caught in M2_SS_gg_SC_ggq_CC_ggq',M2_SS_gg_SC_ggq_CC_ggq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

