      double precision function M2_CC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     triple-collinear limit C_(ia,ib,ic)
c     for sectors (ia,ib,ic)+permutations...
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
      logical flavourmatch_ij, flavourmatch_ik, flavourmatch_jk
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
      double precision zi, zj, zk, zij
      double precision sir, sjr, skr
c
c     initialise
      M2_CC_qxqqp=0d0
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
c
c     Check over flavours
c
      if(.not.(leg_PDGs(i).eq.(-leg_PDGs(j)).and.abs(leg_PDGs(i)).ne.abs(leg_PDGs(k)).and.abs(leg_PDGs(k)).le.6)) return
c
c     possible cuts
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)
c
      IF(DOCUT(XPBB,NEXTERNAL-2,BORN_LEG_PDGS,0))RETURN
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     checking sector topology
      if(lsec.ne.0) then
              write (*,*) 'C_ijk called with four indices'
              stop 1
      endif
c
       if(.not.((ia.eq.isec.and.ib.eq.jsec.and.ic.eq.ksec).or.(ia.eq.isec.and.ib.eq.ksec.and.ic.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec.and.ic.eq.ksec).or.(ia.eq.jsec.and.ib.eq.ksec.and.ic.eq.isec).or.(ia.eq.ksec.and.ib.eq.isec.and.ic.eq.jsec).or.(ia.eq.ksec.and.ib.eq.jsec.and.ic.eq.isec)))  then
              write (*,*) 'Wrong indices in C_ijk', ia, ib, ic, isec, jsec, ksec
        stop 1
      endif
c
c     invariant quantities
      sij  = xs(i,j)
      sjk  = xs(j,k)
      sik  = xs(i,k)
      sijk = sij+sik+sjk
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
      zi   = sir/(sir+sjr+skr)
      zj   = sjr/(sir+sjr+skr)
      zk   = skr/(sir+sjr+skr)
      zik  = zi + zk
      zjk  = zj + zk
      zij  = zi + zj
c
c     safety check

      IF(SIJ.LE.0D0.OR.SIJK.LE.0d0.or.ZIJ.LE.0D0.OR.ZI.LE.0D0.OR.ZJ.LE.0D0.OR.ZK.LE.0D0)THEN
        WRITE(77,*)'Inaccuracy 1 in M2_CC_qxqqp',SIJ,SIJK,ZIJ,ZI,ZJ,ZK
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
      flavourmatch_ij=leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(i))
      flavourmatch_ik=leg_PDGs(k).eq.-leg_PDGs(i).and.abs(leg_PDGs(j)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(j))
      flavourmatch_jk=leg_PDGs(k).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(i))
c
      if(flavourmatch_ij) then
        M2tmp = CF*TR*(-SIJK**2/(2D0*SIJ**2)*(SJK/SIJK-SIK/SIJK+(ZI-ZJ)/ZIJ)**2+SIJK/SIJ*(2D0*(ZK-ZI*ZJ)/ZIJ+ZIJ)-1D0/2D0)
      elseif(flavourmatch_ik)
        M2tmp = CF*TR*(-SIJK**2/(2D0*SIK**2)*(SJK/SIJK-SIJ/SIJK+(ZI-ZK)/ZIK)**2+SIJK/SIK*(2D0*(ZJ-ZI*ZK)/ZIK+ZIK)-1D0/2D0)
      elseif(flavourmatch_jk)
        M2tmp = CF*TR*(-SIJK**2/(2D0*SJK**2)*(SIJ/SIJK-SIK/SIJK+(ZK-ZJ)/ZJK)**2+SIJK/SJK*(2D0*(ZI-ZK*ZJ)/ZJK+ZJK)-1D0/2D0)
      else
        stop 1
      endif
c
c
      M2TMP = M2TMP*BLO
c
c     compute limit of sector function
      call get_wcc_nnlo(xs,ia,ib,ic,ir,alphaz,nexternal)
c
      M2TMP=M2TMP*wcc_nnlo
c
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(r.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_CC_qxqqp=M2tmp*pref/sijk**2*xj*extra ! see [eq.(C.15); is consistent]
c     apply flavour factor
      M2_CC_qxqqp=M2_CC_qxqqp*%(proc_prefix_rr)s_fl_factor
c
c     plot
      wgtpl=-M2_CC_qxqqp*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,BORN_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_CC_qxqqp).ge.huge(1d0).or.isnan(M2_CC_qxqqp))then
         write(77,*)'Exception caught in M2_CC_qxqqp',M2_CC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      double precision function M2_SSCC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     double-soft collinear limit SS_CC_(ia,ib,ic)
c     for sectors (ia,ib,ic)+permutations...
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
      logical flavourmatch_ij, flavourmatch_ik, flavourmatch_jk
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
      double precision zi, zj, zk, zij
      double precision sir, sjr, skr
c
c     initialise
      M2_SSCC_qxqqp=0d0
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
c
c     Check over flavours
c
      if(.not.(leg_PDGs(i).eq.(-leg_PDGs(j)).and.abs(leg_PDGs(i)).ne.abs(leg_PDGs(k)).and.abs(leg_PDGs(k)).le.6)) return
c
c     possible cuts
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)
c
      IF(DOCUT(XPBB,NEXTERNAL-2,BORN_LEG_PDGS,0))RETURN
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     checking sector topology
      if(lsec.ne.0) then
              write (*,*) 'S_ij*C_ijk called with four indices'
              stop 1
      endif
c
       if(.not.((ia.eq.isec.and.ib.eq.jsec.and.ic.eq.ksec).or.(ia.eq.isec.and.ib.eq.ksec.and.ic.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec.and.ic.eq.ksec).or.(ia.eq.jsec.and.ib.eq.ksec.and.ic.eq.isec).or.(ia.eq.ksec.and.ib.eq.isec.and.ic.eq.jsec).or.(ia.eq.ksec.and.ib.eq.jsec.and.ic.eq.isec)))  then
              write (*,*) 'Wrong indices in S_ij*C_ijk', ia, ib, ic, isec, jsec, ksec
        stop 1
      endif
c
c     invariant quantities
      sij  = xs(i,j)
      sjk  = xs(j,k)
      sik  = xs(i,k)
      sijk = sij+sik+sjk
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
      zi   = sir/(sir+sjr+skr)
      zj   = sjr/(sir+sjr+skr)
      zk   = skr/(sir+sjr+skr)
      zik  = zi + zk
      zjk  = zj + zk
      zij  = zi + zj
c
c     safety check

      IF(SIJ.LE.0D0.OR.SIJK.LE.0d0.or.ZIJ.LE.0D0.OR.ZI.LE.0D0.OR.ZJ.LE.0D0.OR.ZK.LE.0D0)THEN
        WRITE(77,*)'Inaccuracy 1 in M2_SSCC_qxqqp',SIJ,SIJK,ZIJ,ZI,ZJ,ZK
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
      flavourmatch_ij=leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(i))
      flavourmatch_ik=leg_PDGs(k).eq.-leg_PDGs(i).and.abs(leg_PDGs(j)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(j))
c
      if(ib.eq.jsec.and.ic.eq.jsec) then     !ijjk
            Eijkr = (1/sij**2)*((sik*sjr+sir*sjk)/((sik+sjk)*(sir+sjr))-sik*sjk/(sik+sjk)**2-sir*sjr/(sir+sjr)**2) - skr/sij/(sik+sjk)/(sir+sjr)
            M2TMP = SIJK**2*(CF*(-2d0*TR*Eijkr))
      elseif(ib.eq.jsec.and.ic.eq.ksec) then !ijkj
            Eikjr = (1/sik**2)*((sij*skr+sir*sjk)*((sij+sjk)*(sir+skr))-sij*sjk/(sij+sjk)**2-sir*skr/(sir+skr)**2) - sjr/sik/(sij+sjk)/(sir+skr)
            M2TMP = SIJK**2*(CF*(-2d0*TR*Eikjr))
      endif
c
      M2TMP = M2TMP*BLO
c
c     compute limit of sector function
      call get_wsscc_nnlo(xs,ia,ib,ic,ir,alphaz,nexternal)
c
      M2TMP=M2TMP*wsscc_nnlo
c
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(r.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
c
      M2tmp=M2tmp*damp
      M2_SSCC_qxqqp=M2tmp*pref*CF/xj*extra ! see [eq.(C.16)]
c     apply flavour factor
      M2_SSCC_qxqqp=M2_SSCC_qxqqp*%(proc_prefix_rr)s_fl_factor
c
c     plot
      wgtpl=-M2_SSCC_qxqqp*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,BORN_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_SSCC_qxqqp).ge.huge(1d0).or.isnan(M2_SSCC_qxqqp))then
         write(77,*)'Exception caught in M2_SSCC_qxqqp',M2_SSCC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

