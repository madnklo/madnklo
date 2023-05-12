      double precision function M2_S(i,xs,xp,wgt,Wsoft,xj,nit,extra,ierr)
c     single-soft limit S_(i) * Wsoft
c     it returns 0 if i is not a gluon
      implicit none
      include 'nexternal.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      integer i,l,m,lb,mb,ierr,nit
      double precision pref,M2tmp,wgt,wgtpl,Wsoft,xj,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,y,z,x,damp
      integer mapped_labels(nexternal), mapped_flavours(NEXTERNAL)
      logical isLOQCDparton(nexternal-1)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
c     external
      integer get_color_dipole_index
      external get_color_dipole_index
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
      logical docut
      integer, parameter :: HEL = - 1
      double precision  %(proc_prefix_S)s_GET_CCBLO
c
c     initialise
      M2_S=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(AS,NLOOP,MU_R)
      pref=-8d0*pi*alphas
c
c     eikonal double sum
      do m=1,nexternal-1
         if(.not.isNLOQCDparton(m))cycle
         if(m.eq.i)cycle
         do l=m+1,nexternal
            if(.not.isNLOQCDparton(l))cycle
            if(l.eq.i)cycle
c
c     determine indices in the n-1 body kinematics
c     TODO: check new get_mapped_labels()
            call get_soft_mapped_labels(i,l,m,nexternal,leg_PDGs,mapped_labels,mapped_flavours,isLOQCDparton)
            lb=mapped_labels(l)
            mb=mapped_labels(m)
c     check on LO color labels 
            if(.not.(isLOQCDparton(lb).and.isLOQCDparton(mb)))then
               write(*,*)'Wrong LO indices in soft kernel',lb,mb
               stop
            endif
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the single-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
c     check on leg_PDGs
            if(size(leg_PDGs).ne.nexternal)then
               write(*,*) 'Wrong dimension for leg_PDGs',size(leg_PDGs), nexternal
               stop
            endif
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,'S',xjCS)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
c            if(docut(xpb,nexternal-1))cycle
c
c     invariant quantities
            sil=xs(i,l)
            sim=xs(i,m)
            slm=xs(l,m)
c
c     safety check
            if(sil*sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S',sil,sim
               goto 999
            endif
c
c     call colour-connected Born
            call %(proc_prefix_S)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S)s_GET_CCBLO(lb,mb)
c
c     eikonal
c     TODO: check for DIS
            M2tmp=ccBLO*2d0*slm/(sil*sim)
            if(m.gt.2.and.l.gt.2)then
               y=sil/(sil+sim+slm)
               z=sim/(sim+slm)
               damp=((1d0-y)*(1d0-z))**alpha
            elseif(m.gt.2.and.l.le.2)then
               z=sim/(sim+slm)
               x=1d0 - sil/(sim+slm)
               damp=((1d0-z)*x)**alpha
            elseif(m.le.2.and.l.le.2)then
               x=1d0 - (sil+sim)/slm
               damp=x**alpha
            endif
            M2tmp=M2tmp*damp
            M2_S=M2_S+pref*M2tmp*Wsoft*extra
c
c     plot
            wgtpl=-pref*M2tmp*Wsoft*extra*xj*wgt/nit/2d0/sCM
            if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
         enddo
      enddo
c
c     sanity check
      if(abs(M2_S).ge.huge(1d0).or.isnan(M2_S))then
         write(77,*)'Exception caught in M2_S',M2_S
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_H_C_FgFg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
c     hard-collinear limit C_(ia,ib) - S_(ia)C_(ia,ib) - S_(ib)C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      implicit none
      include 'nexternal.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,M2tmp,wgt,wgtpl,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      double precision alphas,alpha_qcd
      logical docut
      double precision %(proc_prefix_H_C_FgFg)s_GET_KKBLO
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
c
c     initialise
      M2_H_C_FgFg=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c     possible cuts
c      if(docut(xpb,nexternal-1))return
c
c     overall kernel prefix
      alphas=alpha_QCD(as,nloop,mu_R)
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
         write(77,*)'Inaccuracy 1 in M2_H_C_FgFg',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_H_C_FgFg)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(ia,ib,ir,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      KKBLO = %(proc_prefix_H_C_FgFg)s_GET_KKBLO(parent_leg,xpb,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=CA*2d0*(2d0/sab*KKBLO+x/(1d0-x)*(1d0-x**alpha)*BLO+(1d0-x)/x*(1d0-(1d0-x)**alpha)*BLO)
c     account for different damping factors according to
c     recoiler position (ir) 
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_H_C_FgFg=M2tmp*pref/sab*extra
c
c     plot
      wgtpl=-M2_H_C_FgFg*xj*wgt/nit/2d0/sCM
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C_FgFg).ge.huge(1d0).or.isnan(M2_H_C_FgFg))then
         write(77,*)'Exception caught in M2_H_C_FgFg',M2_H_C_FgFg
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
  
                  
      double precision function M2_H_C_FgFq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
c     hard-collinear limit C_(ia,ib) - S_(ia)C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      implicit none
      include 'nexternal.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      integer ia,ib,ir,ierr,nit
      double precision pref,M2tmp,wgt,wgtpl,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision ans(0:nsqso_born)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      double precision alphas,alpha_qcd
      logical docut
      integer,parameter :: HEL = - 1
c
c     initialise
      M2_H_C_FgFq=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c      if(docut(xpb,nexternal-1))return
c
c     overall kernel prefix
      alphas=alpha_QCD(as,nloop,mu_R)
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
         write(77,*)'Inaccuracy 1 in M2_H_C_FgFq',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_H_C_FgFq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c     In the following equation the x variable is related to the quark energy
      M2tmp=BLO*CF*((1d0-x)+2d0*x/(1d0-x)*(1d0-x**alpha))
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_H_C_FgFq=M2tmp*pref/sab*extra
c
c     plot
      wgtpl=-M2_H_C_FgFq*xj*wgt/nit/2d0/sCM
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C_FgFq).ge.huge(1d0).or.isnan(M2_H_C_FgFq))then
         write(77,*)'Exception caught in M2_H_C_FgFq',M2_H_C_FgFq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_H_C_FqFqx(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
c     hard-collinear limit C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      implicit none
      include 'nexternal.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,M2tmp,wgt,wgtpl,xj,extra
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
      double precision alphas,alpha_qcd
      double precision %(proc_prefix_H_C_FqFqx)s_get_kkblo
      logical docut
c
c     initialise
      M2_H_C_FqFqx=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c
c     possible cuts
c      if(docut(xpb,nexternal-1))return
c
c     overall kernel prefix
      alphas=alpha_QCD(as,nloop,mu_R)
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
         write(77,*)'Inaccuracy 1 in M2_H_C_FqFqx',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_H_C_FqFqx)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(ia,ib,ir,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      KKBLO = %(proc_prefix_H_C_FqFqx)s_GET_KKBLO(parent_leg,xpb,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=TR*(BLO-4d0/sab*KKBLO)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_H_C_FqFqx=M2tmp*pref/sab*extra
c
c     plot
      wgtpl=-M2_H_C_FqFqx*xj*wgt/nit/2d0/sCM
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C_FqFqx).ge.huge(1d0).or.isnan(M2_H_C_FqFqx))then
         write(77,*)'Exception caught in M2_H_C_FqFqx',M2_H_C_FqFqx
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      SUBROUTINE DUMMY_ME_ACCESSOR_HOOK(P,HEL,USER_ALPHAS,ANS)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'nsqso_born.inc'
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS(0:NSQSO_BORN)
      INTEGER HEL
      DOUBLE PRECISION USER_ALPHAS

      write(*,*) 'This subroutine should never be called!!!'
      write(*,*) 'This counterterm does not contribute to this sector.'
      write(*,*) 'Exit...'
      STOP
      END


      DOUBLE PRECISION FUNCTION DUMMY_GET_CCBLO(LB,MB)
      IMPLICIT NONE
      INTEGER LB,MB

      WRITE(*,*) 'YOU ARE IN DUMMY_GET_CCBLO'
      WRITE(*,*) 'THIS FUNCTION MUST NEVER BE CALLED!'
      WRITE(*,*) 'EXIT...'

      STOP
      END


      DOUBLE PRECISION FUNCTION DUMMY_GET_KKBLO(PARENT_LEG,XPB,KT)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INTEGER PARENT_LEG
      DOUBLE PRECISION XPB(0:3,NEXTERNAL-1),KT(0:3)

      WRITE(*,*) 'YOU ARE IN DUMMY_GET_KKBLO'
      WRITE(*,*) 'THIS FUNCTION MUST NEVER BE CALLED!'
      WRITE(*,*) 'EXIT...'

      STOP
      END







