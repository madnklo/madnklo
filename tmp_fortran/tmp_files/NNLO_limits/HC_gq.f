
      double precision function M2_HC_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     collinear limit C_(ia,ib) * Wcollinear  - S_(ia)C_(ia,ib) * Wsoftcollinear
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'mapped_labels_common.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      include 'input.inc'
      include 'run.inc'
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,M2_SC_gq,M2_C_gq,wgt,wgts(1),wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,damp
      double precision rnlo,ans(0:nsqso_born)
      integer, parameter :: hel = - 1
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
      integer %(proc_prefix_HC_gq)s_den
      common/%(proc_prefix_HC_gq)s_iden/%(proc_prefix_HC_gq)s_den
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
      M2_HC_gq=0d0
      M2_SC_gq=0d0
      M2_C_gq=0d0
      ierr=0
c
c     checks
      if(.not.(abs(leg_pdgs(ib)).le.6.and.leg_pdgs(ia).eq.21))then
         write(*,*)'Wrong pdgs in M2_HC_gq',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.(ia.eq.isec.and.ib.eq.jsec))then
         write(*,*)'Wrong indices in M2_HC_gq',ia,ib,isec,jsec
         stop
      endif
c
c     possible cuts
      if(docut(xpb,nexternal-1,real_leg_pdgs,1))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c
c     fill mapped labels (ijr)
      call fill_real_mapped_labels(ia,ib,leg_pdgs,real_leg_pdgs)
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sbr/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HC_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Real
      call %(proc_prefix_HC_gq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      rnlo = ans(0)
c
      M2_C_gq=CF*((1d0-x)+2d0*x/(1d0-x))*rnlo
      if(ia.eq.isec) M2_SC_gq=CF*2d0*x/(1d0-x)*rnlo
c
c     compute soft-collinear limit of sector function, eq. (C.53)
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_labels(csec,ia,ib)
      map2=real_labels(dsec,ia,ib)
      call get_wbar_nlo(map1,map2)
      M2_SC_gq=M2_SC_gq*wbar_nlo
c
c     compute collinear limit of sector function, eq. (C.52)
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(ia,ib,ksec,ir)
      M2_C_gq=M2_C_gq*wc_nlo*wbar_nlo
c
c     include prefactors
      M2_C_gq = M2_C_gq*dble(%(proc_prefix_HC_gq)s_den)/dble(%(proc_prefix_rr)s_den)*%(proc_prefix_rr)s_fl_factor*pref/sab*xj*extra
      M2_SC_gq = M2_SC_gq*dble(%(proc_prefix_HC_gq)s_den)/dble(%(proc_prefix_rr)s_den)*%(proc_prefix_rr)s_fl_factor*pref/sab*xj*extra
      M2_HC_gq = M2_C_gq-M2_SC_gq
c
      if(test_sector_function) M2_HC_gq = (wc_nlo-1d0)*wbar_nlo
c
      call ct_log('KC                     ',M2_C_gq)
      call ct_log('KS_C                   ',M2_SC_gq)
c
c     plot
      wgtpl=-M2_HC_gq*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_HC_gq).ge.huge(1d0).or.isnan(M2_HC_gq))then
         write(77,*)'Exception caught in M2_HC_gq',M2_HC_gq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
