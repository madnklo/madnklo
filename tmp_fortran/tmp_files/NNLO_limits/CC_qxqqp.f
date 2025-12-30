      double precision function M2_CC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C(i,j,k) limit for sectors (i,j,k) + permutations
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer i,j,k,r,ierr,nit,parent_leg
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      logical flavourmatch_ij,flavourmatch_ik,flavourmatch_jk
      logical match_ijk,match_ikj,match_jik,match_jki,match_kij,match_kji
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
      double precision sijk,sij,sik,sjk,sir,sjr,skr
      double precision zi,zj,zk,zij,zik,zjk
c
c     initialise
      M2_CC_qxqqp=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(lsec.ne.0) then
        write (*,*) 'M2_CC_qxqqp called with four indices',isec,jsec,ksec,lsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch_ij = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(i))
      flavourmatch_ik = leg_PDGs(i).eq.-leg_PDGs(k).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(j)).le.5.and.abs(leg_PDGs(j)).ne.abs(leg_PDGs(i))
      flavourmatch_jk = leg_PDGs(j).eq.-leg_PDGs(k).and.abs(leg_PDGs(j)).le.5.and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(i)).ne.abs(leg_PDGs(j))
      if(.not.(flavourmatch_ij.or.flavourmatch_ik.or.flavourmatch_jk))then
        write(*,*) 'flavour mismatch in M2_CC_qxqqp', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     check match of sector indices with function arguments
      match_ijk = i.eq.isec.and.j.eq.jsec.and.k.eq.ksec
      match_ikj = i.eq.isec.and.k.eq.jsec.and.j.eq.ksec
      match_jik = j.eq.isec.and.i.eq.jsec.and.k.eq.ksec
      match_jki = j.eq.isec.and.k.eq.jsec.and.i.eq.ksec
      match_kij = k.eq.isec.and.i.eq.jsec.and.j.eq.ksec
      match_kji = k.eq.isec.and.j.eq.jsec.and.i.eq.ksec
      if(.not.(match_ijk.or.match_ikj.or.match_jik.or.match_jki.or.match_kij.or.match_kji))then
        write (*,*) 'Wrong indices in M2_CC_qxqqp',i,j,k,isec,jsec,ksec
        stop 1
      endif
c
c     possible cuts
!     TODO: CHECK CORRECT VALUES FOR THE FIRST FOUR ARGUMENTS OF GET_UNDERLYING_PDGS!
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
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
      zi   = sir/(sir+sjr+skr)
      zj   = sjr/(sir+sjr+skr)
      zk   = skr/(sir+sjr+skr)
      zik  = zi+zk
      zjk  = zj+zk
      zij  = zi+zj
c
c     safety check
      IF(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0.or.zk.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_CC_qxqqp',SIJ,SIK,SJK,ZI,ZJ,ZK
        GOTO 999
      ENDIF
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
      BLO = ANS(0)
c
c$$$      call get_collinear_mapped_labels(i,j,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
c$$$      parent_leg = mapped_labels(j)
c$$$      if(mapped_flavours(j).ne.21)then
c$$$         write(*,*) 'Wrong parent particle label!', j, mapped_flavours(j)
c$$$         stop
c$$$      endif
c
c     double-collinear kernel, eq. (B.16) of 2212.11190
      if(flavourmatch_ij) then
         M2tmp = CF*TR*(-SIJK**2/(2D0*SIJ**2)*(SJK/SIJK-SIK/SIJK+(ZI-ZJ)/ZIJ)**2+SIJK/SIJ*(2D0*(ZK-ZI*ZJ)/ZIJ+ZIJ)-1D0/2D0)
      elseif(flavourmatch_ik)
         M2tmp = CF*TR*(-SIJK**2/(2D0*SIK**2)*(SJK/SIJK-SIJ/SIJK+(ZI-ZK)/ZIK)**2+SIJK/SIK*(2D0*(ZJ-ZI*ZK)/ZIK+ZIK)-1D0/2D0)
      elseif(flavourmatch_jk)
         M2tmp = CF*TR*(-SIJK**2/(2D0*SJK**2)*(SIJ/SIJK-SIK/SIJK+(ZK-ZJ)/ZJK)**2+SIJK/SJK*(2D0*(ZI-ZK*ZJ)/ZJK+ZK)-1D0/2D0)
      endif
      M2TMP = M2TMP*BLO
c
c     include double-collinear sector function
      call get_wcc_nnlo(xs,ia,ib,ic,ir,alphaz,nexternal)
      M2TMP=M2TMP*WCC_NNLO
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_CC_qxqqp = M2tmp*pref/sijk**2*xj*extra ! see [eq.(C.15); is consistent]
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
c     S(i,j) C(i,j,k) limit for sectors (i,j,k) + permutations
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer i,j,k,r,ierr,nit,parent_leg
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
      logical flavourmatch_ij, flavourmatch_ik, flavourmatch_jk
      logical match_ijk,match_ikj,match_jik,match_jki,match_kij,match_kji
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
      double precision sijk,sij,sik,sjk,sir,sjr,skr
      double precision zi,zj,zk,zij,zik,zjk
c
c     initialise
      M2_SSCC_qxqqp=0d0
      M2tmp=0d0
      ierr=0
c
c     check sector topology
      if(lsec.ne.0) then
        write (*,*) 'M2_SSCC_qxqqp called with four indices',isec,jsec,ksec,lsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch_ij = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5.and.abs(leg_PDGs(k)).ne.abs(leg_PDGs(i))
      flavourmatch_ik = leg_PDGs(i).eq.-leg_PDGs(k).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(j)).le.5.and.abs(leg_PDGs(j)).ne.abs(leg_PDGs(i))
      flavourmatch_jk = leg_PDGs(j).eq.-leg_PDGs(k).and.abs(leg_PDGs(j)).le.5.and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(i)).ne.abs(leg_PDGs(j))
!     TODO: BELOW WE JUST HAVE FLAVOURMATCH_IJ and FLAVOURMATCH_JK. IF ONLY THOSE ARE NEEDED, MODIFY THIS CHECK
      if(.not.(flavourmatch_ij.or.flavourmatch_ik.or.flavourmatch_jk))then
        write(*,*) 'flavour mismatch in M2_SSCC_qxqqp', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     check match of sector indices with function arguments
      match_ijk = i.eq.isec.and.j.eq.jsec.and.k.eq.ksec
      match_ikj = i.eq.isec.and.k.eq.jsec.and.j.eq.ksec
      match_jik = j.eq.isec.and.i.eq.jsec.and.k.eq.ksec
      match_jki = j.eq.isec.and.k.eq.jsec.and.i.eq.ksec
      match_kij = k.eq.isec.and.i.eq.jsec.and.j.eq.ksec
      match_kji = k.eq.isec.and.j.eq.jsec.and.i.eq.ksec
      if(.not.(match_ijk.or.match_ikj.or.match_jik.or.match_jki.or.match_kij.or.match_kji))then
        write (*,*) 'Wrong indices in M2_CC_qxqqp',i,j,k,isec,jsec,ksec
        stop 1
      endif
c
c     possible cuts
!     TODO: CHECK CORRECT VALUES FOR THE FIRST FOUR ARGUMENTS OF GET_UNDERLYING_PDGS!
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
      sir  = xs(i,r)
      sjr  = xs(j,r)
      skr  = xs(k,r)
      zi   = sir/(sir+sjr+skr)
      zj   = sjr/(sir+sjr+skr)
      zk   = skr/(sir+sjr+skr)
      zik  = zi+zk
      zjk  = zj+zk
      zij  = zi+zj
c
c     safety check
      IF(sij.lt.0d0.or.sik.lt.0d0.or.sjk.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0.or.zk.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_CC_qxqqp',SIJ,SIK,SJK,ZI,ZJ,ZK
        GOTO 999
      ENDIF
c
c     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
      BLO = ANS(0)
c
c$$$      call get_collinear_mapped_labels(i,j,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
c$$$      parent_leg = mapped_labels(j)
c$$$      if(mapped_flavours(j).ne.21)then
c$$$         write(*,*) 'Wrong parent particle label!', j, mapped_flavours(j)
c$$$         stop
c$$$      endif
c
!     TODO: CROSS CHECK THIS
c     double-soft double-collinear kernel, eq. (C.16) of 2212.11190
      if(j.eq.jsec.and.k.eq.jsec) then     !ijjk
        Eijkr = (1/sij**2)*((sik*sjr+sir*sjk)/((sik+sjk)*(sir+sjr))-sik*sjk/(sik+sjk)**2-sir*sjr/(sir+sjr)**2) - skr/sij/(sik+sjk)/(sir+sjr)
        M2TMP = SIJK**2*(CF*(-2d0*TR*Eijkr))
      elseif(j.eq.jsec.and.k.eq.ksec) then !ijkj
        Eikjr = (1/sik**2)*((sij*skr+sir*sjk)*((sij+sjk)*(sir+skr))-sij*sjk/(sij+sjk)**2-sir*skr/(sir+skr)**2) - sjr/sik/(sij+sjk)/(sir+skr)
        M2TMP = SIJK**2*(CF*(-2d0*TR*Eikjr))
      endif
      M2TMP = M2TMP*BLO
c
c     include double-soft double-collinear sector function
      call get_wsscc_nnlo(xs,ia,ib,ic,ir,alphaz,nexternal)
      M2TMP=M2TMP*WSSCC_NNLO
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_SSCC_qxqqp=M2tmp*pref*CF/xj*extra ! see [eq.(C.16)]
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
