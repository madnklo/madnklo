
      double precision function M2_HC_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     collinear limit C_(ia,ib) * Wcollinear - S_(ia)C_(ia,ib)
      use sectors2_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      include 'input.inc'
      include 'run.inc'
      integer ia,ib,ir,ierr,nit
      double precision pref,M2_C_gq,M2_SC_gq,wgt,wgts(1),wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision ans(0:nsqso_born)
      integer,parameter :: hel = - 1
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
      integer %(proc_prefix_HC_gq)s_den
      common/%(proc_prefix_HC_gq)s_iden/%(proc_prefix_HC_gq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
c
c     initialise
      M2_C_gq=0d0
      M2_SC_gq=0d0
      M2_HC_gq=0d0
      ierr=0
      damp=0d0
c
c     checks
      if(.not.(abs(leg_pdgs(ib)).le.6.and.leg_pdgs(ia).eq.21))then
         write(*,*)'Wrong pdgs in M2_HC_gq',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)))then
         write(*,*)'Wrong indices in M2_HC_gq',ia,ib,isec,jsec
         stop
      endif
c
c     possible cuts
      if(docut(xpb,nexternal-1,underlying_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sbr/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HC_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_HC_gq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c     In the following equation the x variable is related to the quark energy
      M2_C_gq  = CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*BLO
      if(ia.eq.isec)M2_SC_gq = CF*(2d0*x/(1d0-x)*(1d0+1d0-x**alpha))*BLO
c     compute collinear limit of sector function
      call get_wc_nlo(isec,jsec,iref)
      M2_C_gq =  M2_C_gq*wc_nlo
c     account for different damping factors according to recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2_HC_gq = M2_C_gq-M2_SC_gq
c     include prefactors
      M2_HC_gq = M2_HC_gq *dble(%(proc_prefix_HC_gq)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor*damp*pref/sab*xj*extra
c
c     plot
      wgtpl=-M2_HC_gq*wgt*wgt_chan
c     if(doplot)call histo_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgtpl)
      wgts=wgtpl
      if(doplot)call analysis_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgts)
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

