      double precision function M2_S_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     single-soft limit S_(i) * Wsoft
c     it returns 0 if i is not a gluon
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
      integer i,j,l,m,ierr,nit
      integer jb,lb,mb
      integer jbb,lbb,mbb
      logical isNLOmappedQCDparton(nexternal-1)
      logical isLOmappedQCDparton(nexternal-2)
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS1,xjCS2
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision sil,sim,slm,sij,sjl,sjm,ml2,mm2,y,z,x,damp
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
      double precision pmass(nexternal)
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
      double precision  %(proc_prefix_S_g)s_GET_CCBLO
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      integer %(proc_prefix_S_g)s_den
      common/%(proc_prefix_S_g)s_iden/%(proc_prefix_S_g)s_den
      integer isec,jsec,ksec,lsec,iref
      common/cpartindices/isec,jsec,ksec,lsec,iref
      integer asec,bsec,csec,dsec
      common/csecindices/asec,bsec,csec,dsec
      integer map1,map2
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer real_mapped_labels(nexternal),Born_mapped_labels(nexternal-1)
      common/c_NNLO_mapped_labels/real_mapped_labels,Born_mapped_labels
      include 'pmass.inc'
c
c     initialise
      M2_S_g=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c     return if not gluon
      if(leg_pdgs(i).ne.21)return
c
c     safety check on PDGs
      if(size(leg_pdgs).ne.nexternal)then
        write(*,*) 'M2_s_g:'
        write(*,*) 'Wrong dimension for leg_pdgs',size(leg_pdgs),nexternal
        stop
      endif
c
c     call soft limit of sector function according to eq. (C.51)
      call get_sig2(xs,alpha_mod,nexternal)
      call get_WS_NLO(asec,bsec)
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wbar_nlo(map1,map2)
c
c     overall kernel prefix
      alphas=alpha_qcd(asmz,nloop,scale)
      pref=-8d0*pi*alphas
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isnnloQCDparton(m))cycle
         if(m.eq.i)cycle
         do l=1,nexternal
            if(.not.isnnloQCDparton(l))cycle
            if(l.eq.i)cycle
            if(l.eq.m)cycle
c
            lb=real_mapped_labels(l)
            mb=real_mapped_labels(m)
c
c         check labels and pdgs
          if(.not.(isnlomappedqcdparton(lb).and.isnlomappedqcdparton(mb)))then
            write(*,*)'wrong indices 1 in M2_S_g',lb,mb
            stop
          endif
          if(leg_pdgs(l).ne.born_leg_pdgs(lb).or.born_leg_pdgs(m).ne.born_leg_pdgs(mb))then
            write(*,*)'wrong indices 2 in M2_S_g',l,m,lb,mb
            stop
          endif
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the double-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS1,real_mapped_labels)
            if(xjCS1.eq.0d0)goto 999
            call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            if(docut(xpb,nexternal-1,born_leg_pdgs,0))cycle
c
c     invariant quantities
            sil=xs(i,l)
            sim=xs(i,m)
            slm=xs(l,m)
            ml2=pmass(l)**2
            mm2=pmass(m)**2
c
c     safety check
            if(sil*sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S_g',sil,sim
               goto 999
            endif
c
c     call colour-connected Born
            call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S_g)s_GET_CCBLO(lb,mb)
c
c     eikonal
            M2tmp = slm/(sil*sim) - ml2/sil**2 - mm2/sim**2
            M2tmp = ccBLO*M2tmp
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_rr)s_den)
c
c     damping factors
            if(m.gt.2.and.l.gt.2)then
               y=sil/(sil+sim+slm)
               z=sim/(sim+slm)
               damp=((1d0-y)*(1d0-z))**alpha
            elseif(m.gt.2.and.l.le.2)then
               z=sim/(sim+slm)
               x=1d0-sil/(sim+slm)
               damp=((1d0-z)*x)**alpha
            elseif(m.le.2.and.l.le.2)then
               x=1d0-(sil+sim)/slm
               damp=x**alpha
            endif
            M2tmp=M2tmp*damp*xj
            M2_S_g=M2_S_g+pref*M2tmp*WS_NLO*Wbar_NLO*extra
c
c     plot
            wgtpl=-pref*M2tmp*extra*wgt/nit*wgt_chan
            wgts = wgtpl
c            if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
            if(doplot)call analysis_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgts)
c
         enddo
      enddo
c
c     apply flavour factor
      M2_S_g = M2_S_g * %(proc_prefix_rr)s_fl_factor
c
c     sanity check
      if(abs(M2_S_g).ge.huge(1d0).or.isnan(M2_S_g))then
         write(77,*)'Exception caught in M2_S_g',M2_S_g
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

