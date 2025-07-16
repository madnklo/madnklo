

      double precision function M2_HCC_qxqq(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
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
      integer i,j,k,r,ierr,nit,parent_leg
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
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
      integer %(proc_prefix_Born)s_den
      common/%(proc_prefix_Born)s_iden/%(proc_prefix_Born)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-2)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      double precision sijk, sij, sik, sjk
      double precision zi, zj, zk, zij, zik, zjk
      double precision sir, sjr, skr
c
c     initialise
      M2_HCC_qxqq=0d0
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
      zik  = 0d0
      zjk  = 0d0

c     Check over flavours

      if(.not.(leg_PDGs(i).eq.(-leg_PDGs(j)).and.abs(leg_PDGs(k)).le.6)) return
c
c     possible cuts
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)
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
      zi   = xs(i,r)/(xs(i,r)+xs(j,r)+xs(k,r))
      zj   = xs(j,r)/(xs(i,r)+xs(j,r)+xs(k,r))
      zk   = xs(k,r)/(xs(i,r)+xs(j,r)+xs(k,r))
      zij  = zi + zj
      zik  = zi + zk
      zjk  = zj + zk

      sir = xs(i,r)
      sjr = xs(j,r)
      skr = xs(k,r)
c
c     safety check

      IF(SIJ.LE.0D0.OR.SIJK.LE.0d0.or.ZIJ.LE.0D0.OR.ZI.LE.0D0.OR.ZJ.LE.0D0.OR.ZK.LE.0D0)THEN
        WRITE(77,*)'Inaccuracy 1 in M2_HC_qqx',SIJ,SIJK,ZIJ,ZI,ZJ,ZK
        GOTO 999
      ENDIF
C
C     Call Born
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
      BLO = ANS(0)
c
      call get_collinear_mapped_labels(i,j,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent_leg = mapped_labels(j)
      if(mapped_flavours(j).ne.21)then
         write(*,*) 'Wrong parent particle label!', j, mapped_flavours(j)
         stop
      endif
c
c     If q=qp, the kernel is given by (30) in CG/9810389
c     P^(0g)(I,J,K) + P^(0g)(I,K,J) + P^(0g,id)
c
c     P^(0g) first term in (30) or [B.16 in 2212.11190v2]
      M2TMP = CF*TR*(-SIJK**2/(2D0*SIJ**2)*(SJK/SIJK-SIK/SIJK+(ZI-ZJ)/ZIJ)**2+SIJK/SIJ*(2D0*(ZK-ZI*ZJ)/ZIJ+ZIJ)-1D0/2D0)
c
c     P^(0g) second term in (30) or [B.16 in 2212.11190v2] with J<>K
      M2TMP = M2TMP + CF*TR*(-SIJK**2/(2D0*SIK**2)*(SJK/SIJK-SIJ/SIJK+(ZI-ZK)/ZIK)**2+SIJK/SIK*(2D0*(ZJ-ZI*ZK)/ZIK+ZIK)-1D0/2D0)
c
c     P^(0g,id) third term in (30) or [B.17 in 2212.11190v2]
      if (LEG_PDGS(K).EQ.LEG_PDGS(I)) then
         M2TMP = M2TMP + CF*(2D0*CF-CA)*(-SIJK**2*ZK/(2D0*SJK*SIK)*(1D0+ZK**2)/(ZJK*ZIK)+(SIJ/SJK+SIJ/SIK)+SIJK/(2D0*SJK)*((1D0+ZK**2)/ZIK-2D0*ZJ/ZJK)+SIJK/(2D0*SIK)*((1D0+ZK**2)/ZJK-2D0*ZI/ZIK))
      elseif (LEG_PDGS(K).EQ.-LEG_PDGS(I)) then
         M2TMP = M2TMP + CF*(2D0*CF-CA)*(-SIJK**2*ZI/(2D0*SIJ*SIK)*(1D0+ZI**2)/(ZIJ*ZIK)+(SJK/SIJ+SJK/SIK)+SIJK/(2D0*SIJ)*((1D0+ZI**2)/ZIK-2D0*ZJ/ZIJ)+SIJK/(2D0*SIK)*((1D0+ZI**2)/ZIJ-2D0*ZK/ZIK))
      endif

      M2TMP = M2TMP/SIJK**2*BLO

c     Subtract double soft limit

      M2TMP = M2tmp - (-CF*(2D0*TR*(((SIK*SJR-SIR*SJK)**2-SKR*SIJ*(SIK+SJK)*(SIR+SJR))/(SIJ**2*(SIK+SJK)**2*(SIR+SJR)**2))))*BLO
C     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(r.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_HCC_qxqq=M2tmp*pref*xj*extra
c     apply flavour factor
      M2_HCC_qxqq=M2_HCC_qxqq*%(proc_prefix_rr)s_fl_factor
c
c     plot
      wgtpl=-M2_HCC_qxqq*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,BORN_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_HCC_qxqq).ge.huge(1d0).or.isnan(M2_HCC_qxqq))then
         write(77,*)'Exception caught in M2_HCC_qxqq',M2_HCC_qxqq
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
