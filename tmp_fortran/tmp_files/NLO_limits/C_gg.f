

      double precision function M2_C_gg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     collinear limit C_(ia,ib)
c     this is meant to represent the full collinear
c     for sector (ia,ib)
      use sectors2_module
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision alphaz
      parameter(alphaz=1d0)
      double precision %(proc_prefix_C_gg)s_GET_KKBLO
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_C_gg)s_den
      common/%(proc_prefix_C_gg)s_iden/%(proc_prefix_C_gg)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      integer mapped_flavours(nexternal-1),mapped_indices_shuff(nexternal-1)
      common/c_mapped_quantities_c/mapped_labels,mapped_flavours,mapped_indices_shuff
      double precision xpb_to_ME(0:3,nexternal-1)
c     
c     initialise
      M2_C_gg=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpb_to_ME=0d0
c     
c     checks
      if(leg_pdgs(ia).ne.21.or.leg_pdgs(ib).ne.21)then
         write(*,*)'Wrong pdgs in M2_C_gg',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.(ia.eq.isec.and.ib.eq.jsec))then
         write(*,*)'Wrong indices in M2_C_gg',ia,ib,isec,jsec
         stop
      endif
c
c     possible cuts
      IF(DOCUT(XPB,NEXTERNAL-1,MAPPED_FLAVOURS,0))RETURN
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
         write(77,*)'Inaccuracy 1 in M2_C_gg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      XPB_TO_ME(0:3,MAPPED_INDICES_SHUFF(:))=XPB(0:3,:)
      call %(proc_prefix_C_gg)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
      BLO = ANS(0)
c
      parent_leg = mapped_indices_shuff(mapped_labels(ib))
      if(underlying_leg_pdgs(parent_leg).ne.21)then
         write(*,*)'wrong parent label in M2_C_gg'
         write(*,*)ib,parent_leg,mapped_labels(ib),underlying_leg_pdgs(parent_leg)
         stop
      endif
c
      KKBLO = %(proc_prefix_C_gg)s_GET_KKBLO(parent_leg,xpb_to_ME,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=CA*2d0*(2d0/sab*KKBLO+x/(1d0-x)*(1d0+1d0-x**alpha)*BLO+(1d0-x)/x*(1d0+1d0-(1d0-x)**alpha)*BLO)
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_C_gg)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir) 
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp

c     compute collinear limit of sector function
      call get_wc_nlo(xs,isec,jsec,iref,alphaz,nexternal)
      
      M2_C_gg=M2tmp*pref/sab*xj*extra*wc_nlo
c     apply flavour factor
      M2_C_gg=M2_C_gg*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_C_gg*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,mapped_flavours,wgtpl)
c
c     sanity check
      if(abs(M2_C_gg).ge.huge(1d0).or.isnan(M2_C_gg))then
         write(77,*)'Exception caught in M2_C_gg',M2_C_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_SC_gg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     soft-collinear limit S_(ia)C_(ia,ib)
c     this is meant to represent the soft-collinear
c     for sector (ia,ib)
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,damp,y,xinit
      double precision ANS(0:NSQSO_BORN)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_SC_gg)s_den
      common/%(proc_prefix_SC_gg)s_iden/%(proc_prefix_SC_gg)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      integer mapped_flavours(nexternal-1),mapped_indices_shuff(nexternal-1)
      common/c_mapped_quantities_c/mapped_labels,mapped_flavours,mapped_indices_shuff
      double precision xpb_to_ME(0:3,nexternal-1)
c     
c     initialise
      M2_SC_gg=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpb_to_ME=0d0
c     
c     checks
      if(leg_pdgs(ia).ne.21.or.leg_pdgs(ib).ne.21)then
         write(*,*)'Wrong pdgs in M2_SC_gg',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.(ia.eq.isec.and.ib.eq.jsec))then
         write(*,*)'Wrong indices in M2_SC_gg',ia,ib,isec,jsec
         stop
      endif
c
c     possible cuts
      IF(DOCUT(XPB,NEXTERNAL-1,MAPPED_FLAVOURS,0))RETURN      
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
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_SC_gg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      XPB_TO_ME(0:3,MAPPED_INDICES_SHUFF(:))=XPB(0:3,:)
      call %(proc_prefix_SC_gg)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
      BLO = ANS(0)
c
      M2tmp=CA*2d0*((1d0-x)/x*(1d0+1d0-(1d0-x)**alpha)*BLO)
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_SC_gg)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir) 
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      
      M2_SC_gg=M2tmp*pref/sab*xj*extra
c     apply flavour factor
      M2_SC_gg=M2_SC_gg*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=+M2_SC_gg*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,mapped_flavours,wgtpl)
c
c     sanity check
      if(abs(M2_SC_gg).ge.huge(1d0).or.isnan(M2_SC_gg))then
         write(77,*)'Exception caught in M2_SC_gg',M2_SC_gg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end  
