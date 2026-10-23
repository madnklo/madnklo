
      double precision function M2_S_SS_gg(i,j,r,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     S(i) S(i,j) kernel times sector function
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'mapped_labels_common.inc'
      include 'colored_partons.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,b,r,l,m,ierr,nit,k
      integer ib,jb,bb,lb,mb,kb,rb
      integer jbb,bbb,lbb,mbb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2,damp,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1),xsbb(nexternal-2,nexternal-2)
      double precision sjm,sij,sim,sml,sil,sbm,sib
      double precision sbml,sbjm,sbjl,sbbm,sbjb
      double precision Ei_jm,Ei_ml,Ei_bm,Ebj_ml,Ebj_bm
      double precision BLO,ccBLO,ccBLO_imj_mjl,ccBLO_imj_ljm
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),xpbb(0:3,nexternal-2)
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
      logical test_sector_function
      common/ctestsecfun/test_sector_function
      logical consistency_check
      common/cconscheck/consistency_check
c
c     initialise
      M2_S_SS_gg=0d0
      M2tmp=0d0
      ierr=0
      damp=1d0
c
c     check sector topology
      if(bsec.ne.csec.and.bsec.ne.dsec) then
        write (*,*) 'Wrong topology in M2_S_SS_gg',asec,bsec,csec,dsec
        stop 1
      endif
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_S_SS_gg', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=(8d0*pi*alphas)**2/2d0
c
c     compute soft double-soft sector function eq. (C.72)
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)

      b = dsec
c
c     Common label filling
c
c     fill mapped labels
      do m=1,nexternal
         if(.not.isNNLOQCDparton(m)) cycle
         if(m.eq.i.or.m.eq.j) cycle
         ! this covers [i c dsec],[i c j] -> [i m x]
         call fill_real_mapped_labels(i,m,leg_pdgs,real_leg_pdgs)
         ! born-mappings:
         ! [i c dsec, j dsec c] ->
         call fill_born_mapped_labels(real_labels(j,i,m),real_labels(dsec,i,m),real_leg_pdgs,born_leg_pdgs)
         ! [i c j, c j d]       ->
         ! [i c j, d j c]       -> Use [jbar, cbar, xbar] -> [jbar, mbar, xbar]
         ! Note: jbb never appears as Born index
         call fill_born_mapped_labels(real_labels(j,i,m),real_labels(m,i,m),real_leg_pdgs,born_leg_pdgs)
      enddo

c
c     3rd line_ sum_[c!=i,a,b]_[d!=i,a,c] -> 1 mapping [icd,acd] = [iml,jml]
c

      do m=1,nexternal
        if(.not.isNNLOqcdparton(m))cycle
        if(m.eq.i.or.m.eq.j.or.m.eq.b)cycle
        do l=1,nexternal
          if(.not.isNNLOqcdparton(l))cycle
          if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle

c         invariant quantities:
          sml = xs(m,l)
          sim = xs(i,m)
          sil = xs(i,l)
          Ei_ml = sml/sim/sil
c         safety check
          if(sim.le.0d0.or.sil.le.0d0)then
             write(77,*)'Inaccuracy 1 in M2_S_SS_gg',sim,sil
             goto 999
          endif
c
          lb = real_labels(l,i,m)
          mb = real_labels(m,i,m)
          jb = real_labels(j,i,m)
          lbb = born_labels(lb,jb,mb)
          mbb = born_labels(mb,jb,mb)
c         RR->R [iml]
          call phase_space_cs_inv(i,m,l,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,i,m))
          call invariants_from_p(xpb,nexternal-1,xsb,ierr)
          if(ierr.eq.1)goto 999
          sbml = xsb(mb,lb)
          sbjm = xsb(jb,mb)
          sbjl = xsb(jb,lb)
          Ebj_ml = sbml/sbjm/sbjl
c         call wsbar
          call get_sig2(xsb,1d0,nexternal-1)
          map1=real_labels(csec,i,m)
          map2=real_labels(dsec,i,m)
          call get_wsbar_nlo(map1,map2)

c     R->B [iml,jml] -> [jb,mb,lb]
          call phase_space_cs_inv(jb,mb,lb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,mb))
          if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c     possible cuts
          if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c     call colour-connected born
          call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
          ccBLO = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
          M2tmp = 2d0*Ei_ml*Ebj_ml*(2d0*CF**2*ans(0)+CA*ccBLO)*wsbar_nlo
          M2_S_SS_gg = M2_S_SS_gg + M2tmp
c
c     plot
          wgtpl = -M2tmp*ws_nlo
          wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
          wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
          wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
          wgts=wgtpl
c     if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,born_leg_pdgs,wgtpl)
          if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
        enddo
      enddo
c
c     4th line_ sum_[c!=i,a,b] -> 1 mapping [icb,abc] = [imb,jbm]
c

      do m=1,nexternal
        if(.not.isNNLOqcdparton(m))cycle
        if(m.eq.i.or.m.eq.j.or.m.eq.b)cycle

c       invariant quantities:
        sbm = xs(b,m)
        sim = xs(i,m)
        sib = xs(i,b)
        Ei_bm = sbm/sim/sib
c       safety check
        if(sim.le.0d0.or.sib.le.0d0)then
           write(77,*)'Inaccuracy 1 in M2_S_SS_gg',sim,sib
           goto 999
        endif
c
        bb = real_labels(b,i,m)
        mb = real_labels(m,i,m)
        jb = real_labels(j,i,m)
        bbb = born_labels(bb,jb,bb)
        mbb = born_labels(mb,jb,bb)
c       RR->R [imb]
        call phase_space_cs_inv(i,m,b,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,i,m))
        call invariants_from_p(xpb,nexternal-1,xsb,ierr)
        if(ierr.eq.1)goto 999
        sbbm = xsb(bb,mb)
        sbjm = xsb(jb,mb)
        sbjb = xsb(jb,bb)
        Ebj_bm = sbbm/sbjm/sbjb
c       call wsbar
        call get_sig2(xsb,1d0,nexternal-1)
        map1=real_labels(csec,i,m)
        map2=real_labels(dsec,i,m)
        call get_wsbar_nlo(map1,map2)

c       R->B [imb,jbm] -> [jb,bb,mb]
        call phase_space_cs_inv(jb,bb,mb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,bb))
        if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c       possible cuts
        if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c       call colour-connected born
        call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
        ccBLO = %(proc_prefix_Born)s_GET_CCBLO(mbb,bbb)
        M2tmp = 2d0*Ei_bm*Ebj_bm*(2d0*CF**2*ans(0)+CA*ccBLO)*wsbar_nlo
        M2_S_SS_gg = M2_S_SS_gg + M2tmp
c
c       plot
        wgtpl = -M2tmp*ws_nlo
        wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
        wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
        wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
        wgts=wgtpl
c       if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,born_leg_pdgs,wgtpl)
        if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo

c
c     7th line_ sum_[c!=i,a][d!=i,a,c] -> 2 mappings [ica,cad]&[ica,dac] = [imj,mjl]&[imj,ljm]
c

      do m=1,nexternal
         if(.not.ISNNLOQCDPARTON(m))cycle
         if(m.eq.i.or.m.eq.j)cycle
         do l=1,nexternal
            if(.not.ISNNLOQCDPARTON(l))cycle
            if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c
c           invariant quantities: (c,d) in the paper --> (m,l)
            sjm = xs(j,m)
            sij = xs(i,j)
            sim = xs(i,m)
            Ei_jm = sjm/sij/sim
c           safety check
            if(sij.le.0d0.or.sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S_SS_gg',sij,sim
               goto 999
            endif
c
            jb = real_labels(j,i,m)
            lb = real_labels(l,i,m)
            mb = real_labels(m,i,m)

c           RR->R [imj]
            call phase_space_CS_inv(i,m,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_labels(:,i,m))
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
            sbml = xsb(mb,lb)
            sbjm = xsb(jb,mb)
            sbjl = xsb(jb,lb)
            Ebj_ml = sbml/sbjm/sbjl
c
c           call wsbar
            call get_sig2(xsb,1d0,nexternal-1)
            map1=real_labels(csec,i,m)
            map2=real_labels(dsec,i,m)
            call get_wsbar_nlo(map1,map2)

c           R->B mapping (a) : [imj,mjl]
            lbb = born_labels(lb,jb,mb)
            mbb = born_labels(mb,jb,mb)
            call phase_space_cs_inv(mb,jb,lb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,mb))
            if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))goto 77
c           call colour-connected Born matrix element for [ica,cad]=[imj,mjl]
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO_imj_mjl = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
            M2tmp = -2d0*CA*Ei_jm*Ebj_ml*ccBLO_imj_mjl*wsbar_nlo
            M2_S_SS_gg = M2_S_SS_gg + M2tmp

c           R->B mapping (b) : [imj,ljm]
            lbb = born_labels(lb,jb,lb)
            mbb = born_labels(mb,jb,lb)
 77         call phase_space_cs_inv(lb,jb,mb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,jb,lb))
            if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
            if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c           call colour-connected Born matrix element for [ica,dac]=[imj,ljm]
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO_imj_ljm = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
            M2tmp = -2d0*CA*Ei_jm*Ebj_ml*ccBLO_imj_ljm*wsbar_nlo
            M2_S_SS_gg = M2_S_SS_gg + M2tmp
c
c           plot
            wgtpl = -M2tmp*ws_nlo
            wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
            wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c           if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo

c
      M2_S_SS_gg = M2_S_SS_gg*pref*ws_nlo*xj*damp*extra*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_S_SS_gg = M2_S_SS_gg * %(proc_prefix_rr)s_fl_factor
c
      if(test_sector_function) M2_S_SS_gg = wsbar_nlo*ws_nlo
c
      call ct_log('KS_SS                  ',M2_S_SS_gg)
c
c     sanity check
      if(abs(M2_S_SS_gg).ge.huge(1d0).or.isnan(M2_S_SS_gg))then
         write(77,*)'Exception caught in M2_S_SS_gg',M2_S_SS_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
