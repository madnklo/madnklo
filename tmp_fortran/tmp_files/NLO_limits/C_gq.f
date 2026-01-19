
                  
      double precision function M2_C_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
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
      integer ia,ib,ir,ierr,nit
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision ans(0:nsqso_born)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      double precision alphas,alpha_qcd
      double precision alphaz
      parameter(alphaz=1d0)
      integer,parameter :: HEL = - 1
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_C_gq)s_den
      common/%(proc_prefix_C_gq)s_iden/%(proc_prefix_C_gq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      double precision xpbsave(0:3,nexternal-1)
c
c     initialise
      M2_C_gq=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpbsave=xpb
c
c     possible cuts
      call get_underlying_pdgs(isec,jsec,ksec,lsec,nexternal-1,underlying_leg_pdgs)
      call get_collinear_mapped_labels(ia,ib,nexternal,leg_pdgs,mapped_labels,mapped_flavours)
c     reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal,underlying_leg_pdgs,mapped_flavours,mapped_labels,xpbsave)

      if(docut(xpbsave,nexternal-1,underlying_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      if(leg_pdgs(ia).eq.21) then
         x=sbr/(sar+sbr)
      else
         x=sar/(sar+sbr)
      endif
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_C_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_C_gq)s_ME_ACCESSOR_HOOK(xpbsave,hel,alphas,ANS)
      BLO = ANS(0)
c     In the following equation the x variable is related to the quark energy
      M2tmp=BLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0+1d0-x**alpha))
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_C_gq)s_den)/dble(%(proc_prefix_real)s_den)
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

      M2_C_gq=M2tmp*pref/sab*xj*extra*wc_nlo
c     apply flavour factor
      M2_C_gq=M2_C_gq*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_C_gq*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,UNDERLYING_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_C_gq).ge.huge(1d0).or.isnan(M2_C_gq))then
         write(77,*)'Exception caught in M2_C_gq',M2_C_gq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

                  
      double precision function M2_SC_gq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     soft-collinear limit SC_(ia,ib)
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
      integer ia,ib,ir,ierr,nit
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,damp,y,xinit
      double precision ans(0:nsqso_born)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      double precision alphas,alpha_qcd
      integer,parameter :: HEL = - 1
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_SC_gq)s_den
      common/%(proc_prefix_SC_gq)s_iden/%(proc_prefix_SC_gq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      double precision xpbsave(0:3,nexternal-1)
c
c     initialise
      M2_SC_gq=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpbsave=xpb
c     
c     return if not gluon
      if(leg_pdgs(ia).ne.21) then
        write(77,*)'not a gluon in m2_sc_gq',leg_pdgs(ia)
        goto 999
      endif
c     
c     possible cuts
      call get_underlying_pdgs(isec,jsec,ksec,lsec,nexternal-1,underlying_leg_pdgs)
      call get_collinear_mapped_labels(ia,ib,nexternal,leg_pdgs,mapped_labels,mapped_flavours)
c     reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal,underlying_leg_pdgs,mapped_flavours,mapped_labels,xpbsave)

      if(docut(xpbsave,nexternal-1,underlying_leg_pdgs,0))return
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      if(leg_pdgs(ia).eq.21) then
         x=sbr/(sar+sbr)
      else
         x=sar/(sar+sbr)
      endif
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_SC_gq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_SC_gq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c     In the following equation the x variable is related to the quark energy
      M2tmp=BLO*CF*(2d0*x/(1d0-x)*(1d0+1d0-x**alpha))
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_SC_gq)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp

      M2_SC_gq=M2tmp*pref/sab*xj*extra
c     apply flavour factor
      M2_SC_gq=M2_SC_gq*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=+M2_SC_gq*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,UNDERLYING_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_SC_gq).ge.huge(1d0).or.isnan(M2_SC_gq))then
         write(77,*)'Exception caught in M2_SC_gq',M2_SC_gq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
