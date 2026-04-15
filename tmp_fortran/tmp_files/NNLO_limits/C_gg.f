      double precision function M2_C_gg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     collinear limit C_(ia,ib) * Wcollinear  - S_(ia)C_(ia,ib) * Wsoftcollinear
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,M2_SC_gg,M2_HC_gg,wgt,wgts(1),wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision RNLO,KKRNLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr,xa,xb
      double precision ans(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer nlo_mapped_labels(nexternal), nlo_mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision %(proc_prefix_C_gg)s_get_kkblo
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
      integer %(proc_prefix_C_gg)s_den
      common/%(proc_prefix_C_gg)s_iden/%(proc_prefix_C_gg)s_den
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
      M2_C_gg=0d0
      M2_SC_gg=0d0
      M2_HC_gg=0d0
      ierr=0
      damp=0d0
c
c     checks
      if(.not.(leg_pdgs(ia).eq.21.and.leg_pdgs(ib).eq.21))then
         write(*,*)'Wrong pdgs in M2_HC_gg',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.(ia.eq.isec.and.ib.eq.jsec))then
         write(*,*)'Wrong indices in M2_HC_gg',ia,ib,isec,jsec
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
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     coefficients of kt
c     kt = wa pa + wb pb + wr pr
      wa= 1d0 - x
      wb= - x
      wr= - (1d0-2d0*x)*sab/(sar+sbr)
      kt(:) = wa*xp(:,ia) + wb*xp(:,ib) + wr*xp(:,ir)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HC_gg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Real as per eq. (C.9)
      call %(proc_prefix_C_gg)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      RNLO = ANS(0)
c
      parent_leg = real_mapped_labels(ib)
      KKRNLO = %(proc_prefix_C_gg)s_GET_KKBLO(parent_leg,xpb,kt)
c     TODO: improve ktmuktnuRmunu / kt^2
      M2_C_gg=CA*2d0*(2d0/sab*KKRNLO+x/(1d0-x)*(1d0+1d0-x**alpha)*RNLO+(1d0-x)/x*(1d0+1d0-(1d0-x)**alpha)*RNLO)
      M2_S_C_gg=CA*2d0*((1d0-x)/x*(1d0+1d0-(1d0-x)**alpha)*RNLO)
c
c     compute soft-collinear limit of sector function
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wbar_nlo(map1,map2)
      M2_S_C_gg=M2_S_C_gg*wbar_nlo
c
c     compute collinear limit of sector function
      call get_sig2(xs,alpha_mod,nexternal)
      call get_wc_nlo(ia,ib,ksec,ir)
      call get_sig2(xsb,alpha_mod_bar,nexternal-1)
      map1=real_mapped_labels(csec)
      map2=real_mapped_labels(dsec)
      call get_wbar_nlo(map1,map2)
      M2_C_gg=M2_C_gg*wc_nlo*wbar_nlo
c
c     account for different damping factors according to recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
c
      M2_HC_gg = M2_C_gg-M2_SC_gg
      M2_C_gg = M2_HC_gg
c     include prefactors
      M2_C_gg = M2_C_gg*dble(%(proc_prefix_C_gg)s_den)/dble(%(proc_prefix_rr)s_den)*%(proc_prefix_rr)s_fl_factor*damp*pref/sab*xj*extra
c
c     plot
      wgtpl=-M2_C_gg*wgt/nit*wgt_chan
      wgts=wgtpl
c      if(doplot)call histo_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgtpl)
      if(doplot)call analysis_fill(xpb,xsb,nexternal-1,real_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_C_gg).ge.huge(1d0).or.isnan(M2_C_gg))then
         write(77,*)'Exception caught in M2_HC_gg',M2_C_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
