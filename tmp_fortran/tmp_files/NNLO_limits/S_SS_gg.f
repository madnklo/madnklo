
      double precision function M2_S_SS_gg(i,j,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     S(i) S(i,k) kernel times Wsector
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'leg_pdgs.inc'
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
      logical test_sector_function
      common/ctestsecfun/test_sector_function
c
c     initialise
      M2_S_SS_gg=0d0
      M2tmp=0d0
      ierr=0
c
c     check flavour match
      if(leg_pdgs(i).eq.0 .or. leg_pdgs(i).ne.leg_pdgs(j)) then
        write(*,*) 'Flavour mismatch in M2_S_SS_gg', leg_PDGs(i),leg_PDGs(j)
        stop 1
      endif
c
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
c     compute soft double-soft sector function from eq. (C.72)
      call get_sig2(xs,alpha_mod,nexternal)
      call get_ws_nlo(asec,bsec)
      call get_sig2(xsb,1d0,nexternal-1)
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wsbar_nlo(map1,map2)
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
               write(*,*)'Wrong indices 1 in M2_S_SS_gg',lb,mb
               stop
            endif
            if(.not.(islomappedqcdparton(lbb).and.islomappedqcdparton(mbb)))then
               write(*,*)'Wrong indices 2 in M2_S_SS_gg',lbb,mbb
               stop
            endif
            if(leg_pdgs(l).ne.born_leg_pdgs(lbb).or.leg_pdgs(m).ne.born_leg_pdgs(mbb))then
               write(*,*)'Wrong indices 3 in M2_S_SS_gg',l,m,lbb,mbb
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
               write(77,*)'Inaccuracy 1 in M2_S_SS_gg',sij, sil+sjl, sim+sjm
               goto 999
            endif
c
c     call colour-connected Born
c     TODO: fix strings for the associated underlying Born
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_Born)s_GET_CCBLO(lbb,mbb)
c
c     call Born matrix element
            call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ans)
            BLO = ANS(0)
c
c     TODO: call quadruple born B[icd,jcd]_{cdcd}
c            QUADBLO_mlml = 0d0
c
c     eikonal
c     See file K2_I2_G_v2.pdf in the DropBox directory
c     (c,d) -> (m,l) (verified)
            M2tmp = -2d0*CA*CCBLO*(((sim*sjl+sil*sjm)/((sim+sjm)*(sil+sjl))-sim*sjm/(sim+sjm)**2-sil*sjl/(sil+sjl)**2)/sij**2-2d0*slm/(sij*(sim+sjm)*(sil+sjl))+slm*(sim*sjl+sil*sjm-sij*slm)/(sij*sim*sjl*sil*sjm)*(1d0-1d0/2d0*(sim*sjl+sil*sjm)/((sim+sjm)*(sil+sjl))))
            M2tmp = M2tmp + 2d0*(slm/sim/sil)*(slm/sjm/sjl)*2d0*CF**2*BLO
c     Including correct multiplicity factor
            M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
            damp=1d0
            M2tmp=M2tmp*damp*xj
            M2_S_SS_gg=M2_S_SS_gg+pref*M2tmp*extra*ws_nlo*wsbar_nlo
c
c     plot
            wgtpl=-pref*M2tmp*extra*wgt/nit*wgt_chan*wsbar_nlo*ws_nlo
            wgtpl = wgtpl*%(proc_prefix_rr)s_fl_factor
            wgts=wgtpl
c            if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgtpl)
            if(doplot)call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
         enddo
      enddo
c
c     apply flavour factor
      M2_S_SS_gg = M2_S_SS_gg * %(proc_prefix_rr)s_fl_factor
      if(test_sector_function) M2_S_SS_gg = wsbar_nlo*ws_nlo
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
