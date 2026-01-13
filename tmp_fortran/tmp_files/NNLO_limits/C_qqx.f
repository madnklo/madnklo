

      double precision function M2_C_qqx(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     collinear limit C_(ia,ib)
c     this is meant to represent the collinear
c     for sector (ia,ib)
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
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
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
      double precision %(proc_prefix_C_qqx)s_get_kkblo
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_C_qqx)s_den
      common/%(proc_prefix_C_qqx)s_iden/%(proc_prefix_C_qqx)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      double precision xpbsave(0:3,nexternal-1)
c
c     initialise
      M2_C_qqx=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      xpbsave=xpb
c     Check over (ia,ib)
c     They must be equal to (isec,jsec)

      if(.not.((ia.eq.isec.and.ib.eq.jsec))) then
         write (*,*) 'Wrong indices in M2_C_qqx:'
         write(*,*) 'ia, ib, isec, jsec = ', ia, ib, isec, jsec
         stop
      endif


c
c     possible cuts
c      call GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,UNDERLYING_LEG_PDGS)
      call get_collinear_mapped_labels(ia,ib,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
c     Reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal,underlying_leg_pdgs,mapped_flavours,mapped_labels,xpbsave)
      call invariants_from_p(xpbsave,nexternal-1,xsb,ierr)
      if(ierr.eq.1) goto 999



      IF(DOCUT(XPBSAVE,NEXTERNAL-1,UNDERLYING_LEG_PDGS,0))RETURN
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
         write(77,*)'Inaccuracy 1 in M2_C_qqx',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_C_qqx)s_ME_ACCESSOR_HOOK(xpbsave,hel,alphas,ANS)
      BLO = ANS(0)
c
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      KKBLO = %(proc_prefix_C_qqx)s_GET_KKBLO(parent_leg,xpbsave,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=TR*(BLO-4d0/sab*KKBLO)
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_C_qqx)s_den)/dble(%(proc_prefix_real)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp

c     compute collinear limit of sector function
      call get_wc_nnlo(xs,ia,ib,ir,alphaz,n_ext_in)

      M2_C_qqx=M2tmp*pref/sab*xj*extra*wc_nnlo
c     apply flavour factor
      M2_C_qqx=M2_C_qqx*%(proc_prefix_real)s_fl_factor
c
c     plot
      wgtpl=-M2_C_qqx*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbsave,xsb,nexternal-1,UNDERLYING_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_C_qqx).ge.huge(1d0).or.isnan(M2_C_qqx))then
         write(77,*)'Exception caught in M2_C_qqx',M2_C_qqx
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

