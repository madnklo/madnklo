

      double precision function M2_HCC_qxqqp(i,j,k,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     hard-collinear limit C_(ia,ib,ic)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib,ic)+permutations...
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'      
      integer i,j,k,iref,ierr,nit,parent_leg
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2),xsbb(nexternal-2,nexternal-2)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision xpbb(0:3,nexternal-2)
      double precision x,y,xinit,damp
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
      integer %(proc_prefix_rr)s_fl_factor
      common/%(proc_prefix_rr)s_flavour_factor/%(proc_prefix_rr)s_fl_factor
      double precision alphas,alpha_qcd
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      integer %(proc_prefix_HC_qqx)s_den
      common/%(proc_prefix_HC_qqx)s_iden/%(proc_prefix_HC_qqx)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-2)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      double precision sijk, sij, sik, sjk
      double precision zi, zj, zk, zij
      double precision sir, sjr, skr
c
c     initialise
      M2_HC_qxqqp=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      sijk = 0d0
      sij  = 0d0
      sik  = 0d0
      sjk  = 0d0
      sir  = 0d0
      sjr  = 0d0
      skr  = 0d0
      zi   = 0d0
      zj   = 0d0
      zk   = 0d0
      zij  = 0d0

c     Check over flavours

      if(.not.(leg_PDGs(i).eq.(-leg_PDGs(j)).and.abs(leg_PDGs(i)).ne.abs(leg_PDGs(k)))) return
      
      
c
c     possible cuts
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)

      IF(DOCUT(XPBB,NEXTERNAL-2,BORN_LEG_PDGS,0))RETURN


c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     invariant quantities


      sij  = xs(i,j)
      sjk  = xs(j,k)
      sik  = xs(i,k)
      sijk = sij+sik+sjk
      zi   = xs(i,iref)/(xs(i,iref)+xs(j,iref)+xs(k,iref))
      zj   = xs(j,iref)/(xs(i,iref)+xs(j,iref)+xs(k,iref))
      zk   = xs(k,iref)/(xs(i,iref)+xs(j,iref)+xs(k,iref))
      zij  = zi + zj

      sir = xs(i,iref)
      sjr = xs(j,iref)
      skr = xs(k,iref)
c
c     safety check

      IF(SIJ.LE.0D0.OR.SIJK.LE.0d0.or.ZIJ.LE.0D0.OR.ZI.LE.0D0.OR.ZJ.LE.0D0.OR.ZK.LE.0D0)THEN
        WRITE(77,*)'Inaccuracy 1 in M2_HC_qqx',SIJ,SIJK,ZIJ,ZI,ZJ,ZK
        GOTO 999
      ENDIF
C
C     Call Born
      call %(proc_prefix_HC_qqx)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(ia,ib,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
      M2tmp = CF*TR*(-SIJK**2/(2D0*SIJ**2)*(SJK/SIJK-SIK/SIJK+
     $     (ZI-ZJ)/ZIJ)**2+SIJK/SIJ*(2D0*(ZK-ZI*ZJ)/ZIJ+ZIJ)-1D0/2D0)
      M2TMP = M2TMP*BLO
      

c     Subtract double soft limit

      M2TMP = M2tmp - SIJK**2*(-CF*
     $     (2D0*TR*(((SIK*SJR-SIR*SJK)**2-SKR*SIJ*(SIK+SJK)*(SIR
     $     +SJR))/(SIJ**2*(SIK+SJK)**2*(SIR+SJR)**2))))*BLO
C     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_HC_qqx)s_den)/dble(%(proc_prefix_rr)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_HCC_qxqqp=M2tmp*pref/sijk**2*xj*extra
c     apply flavour factor
      M2_HCC_qxqqp=M2_HC_qqx*%(proc_prefix_rr)s_fl_factor
c
c     plot
      wgtpl=-M2_HCC_qxqqp*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_HCC_qxqqp).ge.huge(1d0).or.isnan(M2_HCC_qxqqp))then
         write(77,*)'Exception caught in M2_HCC_qxqqp',M2_HCC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

