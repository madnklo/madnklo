
      double precision function M2_S_SS_gg(i,j,r,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     S(i) S(i,j) kernel times sector function
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
      integer i,j,b,r,l,m,ierr,nit
      integer ib,jb,bb,lb,mb
      integer jbb,bbb,lbb,mbb
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2,damp,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1),xsbb(nexternal-2,nexternal-2)
      double precision sjm,sij,sim,sml,sil,sbm,sib
      double precision sbml,sbjm,sbjl,sbbm,sbjb
      double precision Ei_jm,Ei_ml,Ei_bm,Ebj_ml,Ebj_bm
      double precision blo,ccblo,ccblo_imj_mjl,ccblo_imj_ljm
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
      integer Born_2_mapped_labels(nexternal-1)
      common/c_NNLO_2_mapped_labels/Born_2_mapped_labels
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
c     get PDGs
      ib = real_mapped_labels(i)
      jb = real_mapped_labels(j)
c
c     compute soft double-soft sector function eq. (C.72)
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
c
c     mapping 1: [imj,jml]; third block eq. 115 in Dropbox
c     eikonal double sum (c,d) -> (m,l)
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
               write(77,*)'inaccuracy 1 in m2_s_ss_gg',sij,sim
               goto 999
            endif
c
            lb = real_mapped_labels(l)
            mb = real_mapped_labels(m)
            if(m.ne.iref) then
               lbb = born_mapped_labels(lb)
               mbb = born_mapped_labels(mb)
            elseif(m.eq.iref) then
               lbb = born_2_mapped_labels(lb)
               mbb = born_2_mapped_labels(mb)
            endif
c
            call phase_space_CS_inv(i,m,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
            sbml = xsb(mb,lb)
            sbjm = xsb(jb,mb)
            sbjl = xsb(jb,lb)
            Ebj_ml = sbml/sbjm/sbjl
c
c           call wsbar
            call get_sig2(xsb,1d0,nexternal-1)
            map1=real_mapped_labels(csec)
            map2=real_mapped_labels(dsec)
            call get_wsbar_nlo(map1,map2)
c
            if(m.ne.iref) then
               call phase_space_cs_inv(mb,jb,lb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_mapped_labels)
               if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
               if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
            elseif(m.eq.iref) then
               call phase_space_cs_inv(mb,jb,lb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_2_mapped_labels)
               if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
               if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
            endif
c
c           call colour-connected Born matrix element
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO_imj_mjl = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
            M2tmp = -2d0*CA*Ei_jm*Ebj_ml*ccBLO_imj_mjl*wsbar_nlo
            M2_s_ss_gg = M2_s_ss_gg + M2tmp
c
c           plot
            wgtpl = -m2tmp*ws_nlo
            wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
            wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo
c
c     mapping 2: [ijl,ljm]; third block eq. 115 in Dropbox
c     eikonal double sum (c,d) -> (m,l)
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
               write(77,*)'inaccuracy 1 in m2_s_ss_gg',sij,sim
               goto 999
            endif
c
            lb = real_mapped_labels(l)
            mb = real_mapped_labels(m)
            if(m.eq.iref) then
               lbb = born_mapped_labels(lb)
               mbb = born_mapped_labels(mb)
            elseif(m.ne.iref) then
               lbb = born_2_mapped_labels(lb)
               mbb = born_2_mapped_labels(mb)
            endif
c
            call phase_space_CS_inv(i,m,j,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
            sbml = xsb(mb,lb)
            sbjm = xsb(jb,mb)
            sbjl = xsb(jb,lb)
            Ebj_ml = sbml/sbjm/sbjl
c
c           call wsbar
            call get_sig2(xsb,1d0,nexternal-1)
            map1=real_mapped_labels(csec)
            map2=real_mapped_labels(dsec)
            call get_wsbar_nlo(map1,map2)
c
            if(m.eq.iref) then
               call phase_space_cs_inv(lb,jb,mb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_mapped_labels)
               if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
               if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
            elseif(m.ne.iref) then
               call phase_space_cs_inv(lb,jb,mb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_2_mapped_labels)
               if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
               if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
            endif
c
c           call colour-connected Born matrix element
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO_imj_ljm = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
            M2tmp = -2d0*CA*Ei_jm*Ebj_ml*ccBLO_imj_ljm*wsbar_nlo
            M2_s_ss_gg = M2_s_ss_gg + M2tmp
c
c           plot
            wgtpl = -m2tmp*ws_nlo
            wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
            wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo
c
c     mapping 3: [iml,jml]; first block eq (115) in Dropbox
      b = dsec
      do m=1,nexternal
        if(.not.isNNLOqcdparton(m))cycle
        if(m.eq.i.or.m.eq.j.or.m.eq.b)cycle
        do l=1,nexternal
          if(.not.isNNLOqcdparton(l))cycle
          if(l.eq.i.or.l.eq.j.or.l.eq.m)cycle
c
c         invariant quantities: (c,d) in the paper --> (m,l)
          sml = xs(m,l)
          sim = xs(i,m)
          sil = xs(i,l)
          ei_ml = sml/sim/sil
c         safety check
          if(sim.le.0d0.or.sil.le.0d0)then
            write(77,*)'inaccuracy 1 in m2_s_ss_gg',sim,sil
            goto 999
          endif
c
          lb = real_mapped_labels(l)
          mb = real_mapped_labels(m)
          lbb = born_mapped_labels(lb)
          mbb = born_mapped_labels(mb)
c
          call phase_space_cs_inv(i,m,l,xp,xpb,nexternal,leg_pdgs,xjcs1,real_mapped_labels)
          call invariants_from_p(xpb,nexternal-1,xsb,ierr)
          if(ierr.eq.1)goto 999
          sbml = xsb(mb,lb)
          sbjm = xsb(jb,mb)
          sbjl = xsb(jb,lb)
          ebj_ml = sbml/sbjm/sbjl
c         call wsbar
          call get_sig2(xsb,1d0,nexternal-1)
          map1=real_mapped_labels(csec)
          map2=real_mapped_labels(dsec)
          call get_wsbar_nlo(map1,map2)
c
          call phase_space_cs_inv(jb,mb,lb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_mapped_labels)
          if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c         possible cuts
          if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c         call colour-connected born
          call epem_ccx_me_accessor_hook(xpbb,hel,alphas,ans)
          ccblo = %(proc_prefix_Born)s_GET_CCBLO(mbb,lbb)
          m2tmp = 2d0*ei_ml*ebj_ml*(2d0*cf**2*ans(0)+ca*ccblo)*wsbar_nlo
          M2_s_ss_gg = M2_s_ss_gg + M2tmp
c
c         plot
          wgtpl = -m2tmp*ws_nlo
          wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
          wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
          wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
          wgts=wgtpl
c          if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,born_leg_pdgs,wgtpl)
          if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
        enddo
      enddo
c
c     mapping 4: [imb,jbm]; second block eq (115) in Dropbox
      b = dsec
      do m=1,nexternal
        if(.not.isNNLOqcdparton(m))cycle
        if(m.eq.i.or.m.eq.j.or.m.eq.b)cycle
c
c       invariant quantities: (c,d) in the paper --> (m,l)
        sbm = xs(b,m)
        sib = xs(i,b)
        sim = xs(i,m)
        ei_bm = sbm/sib/sim
c       safety check
        if(sib.le.0d0.or.sim.le.0d0)then
          write(77,*)'inaccuracy 1 in m2_s_ss_gg',sib,sim
          goto 999
        endif
c
        bb = real_mapped_labels(b)
        mb = real_mapped_labels(m)
        bbb = born_mapped_labels(bb)
        mbb = born_mapped_labels(mb)
c
        call phase_space_cs_inv(i,m,b,xp,xpb,nexternal,leg_pdgs,xjcs1,real_mapped_labels)
        call invariants_from_p(xpb,nexternal-1,xsb,ierr)
        if(ierr.eq.1)goto 999
        sbbm = xsb(bb,mb)
        sbjb = xsb(jb,bb)
        sbjm = xsb(jb,mb)
        ebj_bm = sbbm/sbjb/sbjm
c       call wsbar
        call get_sig2(xsb,1d0,nexternal-1)
        map1=real_mapped_labels(csec)
        map2=real_mapped_labels(dsec)
        call get_wsbar_nlo(map1,map2)
c
        call phase_space_cs_inv(jb,bb,mb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,Born_mapped_labels)
        if(xjcs1.eq.0d0.or.xjcs2.eq.0d0)goto 999
c       possible cuts
        if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))cycle
c
c       call colour-connected born
        call epem_ccx_me_accessor_hook(xpbb,hel,alphas,ans)
        ccblo = %(proc_prefix_Born)s_GET_CCBLO(bbb,mbb)
        m2tmp = 2d0*ei_bm*ebj_bm*(2d0*cf**2*ans(0)+ca*ccblo)*wsbar_nlo
        M2_s_ss_gg = M2_s_ss_gg + M2tmp
c
c       plot
        wgtpl = -m2tmp*ws_nlo
        wgtpl = wgtpl*pref*extra*damp*xj*wgt/nit*wgt_chan
        wgtpl = wgtpl*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
        wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
        wgts=wgtpl
c        if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,born_leg_pdgs,wgtpl)
        if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      enddo
c
      m2_s_ss_gg = m2_s_ss_gg*pref*ws_nlo*xj*damp*extra*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      m2_s_ss_gg = m2_s_ss_gg * %(proc_prefix_rr)s_fl_factor
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
