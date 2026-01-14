      double precision function M2_C_CC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     C_(i,j) C_(i,j,k) kernel times WC_CC: i, j are a q-qb pair with same flavour
c     while k is a q (or qb) with any flavour
      implicit none
      use sectors2_module
      use sectors4_module
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer i,j,k,r,ierr,nit
      integer jb,kb,rb
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision xpbb(0:3,nexternal-2)
      double precision ktb(0:3),ktb2,kkBLO,WCCC_NNLO
      double precision x,y,xinit
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer nlo_mapped_labels(nexternal), nlo_mapped_flavours(nexternal)
      integer lo_mapped_labels(nexternal), lo_mapped_flavours(nexternal)
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
      INTEGER REAL_LEG_PDGS(NEXTERNAL-1)
      double precision sij,sir,sjr,sbjk,sbjr,sbkr
      double precision zi,zj
      double precision zbj,zbk
      double precision xpbsave(0:3,nexternal-1),xpbbsave(0:3,nexternal-2)
      double precision xsbsave(nexternal-1,nexternal-1),xsbbsave(nexternal-2,nexternal-2)
c
c     initialise
      M2_C_CC_qxqqp=0d0
      M2tmp=0d0
      ierr=0
      xpbsave  = xpb
      xpbbsave = xpbb
      xsbsave  = xsb
      xsbbsave = xsbb
c
c     check sector topology
      if(lsec.ne.jsec .and. lsec.ne.ksec) then
        write (*,*) 'Wrong topology in M2_C_CC_qxqqp',isec,jsec,ksec,lsec
        stop 1
      endif
c
c     check flavour match
      flavourmatch = leg_PDGs(i).eq.-leg_PDGs(j).and.abs(leg_PDGs(i)).le.5.and.abs(leg_PDGs(k)).le.5
      if(.not.(flavourmatch))then
        write(*,*) 'Flavour mismatch in M2_C_CC_qxqqp', leg_PDGs(i),leg_PDGs(j),leg_PDGs(k)
        stop 1
      endif
c
c     reshuffle NLO momenta and labels according to real_leg_pdgs and check
      call reshuffle_momenta(nexternal,real_leg_pdgs,NLO_mapped_flavours,NLO_mapped_labels,xpbsave)
      call get_collinear_mapped_labels(i,j,nexternal,leg_PDGs,NLO_mapped_labels,NLO_mapped_flavours)
      if(NLO_mapped_flavours(j).ne.21)then
         write(*,*) 'Wrong parent particle label 1 in M2_C_CC_qxqqp', j, NLO_mapped_flavours(j)
         stop
      endif
      call invariants_from_p(xpbsave,nexternal-1,xsbsave,ierr)
      if(ierr.eq.1)goto 999
c
c     reshuffle LO momenta and labels according to Born_leg_pdgs and check
      call reshuffle_momenta(nexternal-1,Born_leg_pdgs,LO_mapped_flavours,LO_mapped_labels,xpbbsave)
      call get_collinear_mapped_labels(jb,kb,nexternal-1,real_leg_PDGs,LO_mapped_labels,LO_mapped_flavours)
      if(LO_mapped_flavours(kb).ne.NLO_mapped_flavours(k))then
         write(*,*) 'Wrong parent particle label 2 in M2_C_CC_qxqqp', kb,k,LO_mapped_flavours(kb),NLO_mapped_flavours(k)
         stop
      endif
      call invariants_from_p(xpbbsave,nexternal-2,xsbbsave,ierr)
      if(ierr.eq.1)goto 999
c
c     possible cuts
!     TODO: CHECK CORRECT VALUES FOR THE FIRST FOUR ARGUMENTS OF GET_UNDERLYING_PDGS!
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-1,REAL_LEG_PDGS)
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)
      IF(DOCUT(XPBBSAVE,NEXTERNAL-2,BORN_LEG_PDGS,0))RETURN
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
      zi   = sir/(sir+sjr)
      zj   = 1d0-zi
      jb = NLO_mapped_labels(j)
      kb = NLO_mapped_labels(k)
      rb = NLO_mapped_labels(r)
      sbjk = xsbsave(jb,kb)
      sbjr = xsbsave(jb,rb)
      sbkr = xsbsave(kb,rb)
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0-zbj
c
c     coefficients of ktb = wa pab + wb pbb + wr prb
      wa =  1-zbj
      wb = -zbj
      wr = - (1d0-2d0*zbj)*sbjk/(sbjr+sbkr)
      ktb(:) = wa*xpbsave(:,jb) + wb*xpbsave(:,kb) + wr*xpbsave(:,rb)
      ktb2 = -zbj*zbk*zbjk
c     TODO: check proc_prefix here below, and also arguments in KKBLO
      KKBLO = %(proc_prefix_C_CC_qxqq)s_GET_KKBLO(jb,xpbsave,ktb)
c
c     safety checks
      IF(sij.lt.0d0.or.sir.lt.0d0.or.sjr.lt.0d0.or.zi.lt.0d0.or.zj.lt.0d0)then
        WRITE(77,*)'Inaccuracy 1 in M2_C_CC_qxqqp',SIJ,SIR,SJR,ZI,ZJ
        GOTO 999
      ENDIF
      IF(sbjk.lt.0d0.or.sbjr.lt.0d0.or.sbkr.lt.0d0.or.zbj.lt.0d0.or.zbk.lt.0d0)then
        WRITE(77,*)'Inaccuracy 2 in M2_C_CC_qxqqp',SBJK,SBJR,SBKR,ZBJ,ZBK
        GOTO 999
      ENDIF
C
C     call Born matrix element
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbbsave,hel,alphas,ANS)
      BLO = ANS(0)
c
c     collinear double-collinear kernel, eq. (C.39) of 2212.11190v2
      Pij = TR*(1d0-2d0*zi*zj)
      Qij = TR*2d0*zi*zj
      Pbjk = CF*(1d0+zbk**2)/zbj
      Ebjkr = sbkr/sbjk/sbjr
      M2TMP = Pij/sij*Pbjk/sbjk - 2d0*CF*Ebjkr*Qij/sij/ktb2
      M2TMP = M2TMP*BLO
c
c     include collinear double-collinear sector function
c     TODO: check syntax and arguments here!
      call get_wccc_nnlo(xs,isec,jsec,ksec,lsec,r,alphaz,nexternal)
      M2TMP=M2TMP*WCCC_NNLO
c
c     include correct multiplicity and flavour factors
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2tmp = M2tmp*%(proc_prefix_rr)s_fl_factor
      M2_C_CC_qxqqp = M2tmp*pref*xj*extra
c
c     plot
      wgtpl=-M2_C_CC_qxqqp*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbbsave,xsbb,nexternal-2,BORN_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_C_CC_qxqqp).ge.huge(1d0).or.isnan(M2_C_CC_qxqqp))then
         write(77,*)'Exception caught in M2_C_CC_qxqqp',M2_C_CC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


c     TODO: adapt the following limit
      double precision function M2_C_SS_qqx_CC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     soft-collinear triple-collinear limit SC_(ia,ib)C_(ia,ib,ic)
c     this is meant to represent the collinear and appears only in topology ijjk
c     for sectors (ia,ib,ic)+permutations...
      implicit none
      use sectors4_module
      use sectors2_module
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      integer i,j,k,r,ierr,nit,parent_leg
      integer jb,kb,rb
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,xjb,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision xsbb(nexternal-2,nexternal-2)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision xpbb(0:3,nexternal-2)
      double precision x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer nlo_mapped_labels(nexternal), nlo_mapped_flavours(nexternal)
      integer lo_mapped_labels(nexternal), lo_mapped_flavours(nexternal)
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
      INTEGER REAL_LEG_PDGS(NEXTERNAL-1)
      double precision sij, sir, sjr, sbjk, sbjr, sbkr
      double precision zi, zj
      double precision zbj, zbk
      integer ic,id
      double precision xpbsave(0:3,nexternal-1), xpbbsave(0:3,nexternal-2)
c
c     initialise
      M2_C_SS_qqx_CC_qxqqp=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      sij  = 0d0
      sir  = 0d0
      sjr  = 0d0
      sbjk = 0d0
      sbjr = 0d0
      sbkr = 0d0
      zi   = 0d0
      zj   = 0d0
      zbj  = 0d0
      zbk  = 0d0
      ic = 0
      id = 0
      xpbsave  = xpb
      xpbbsave = xpbb

c     Check over flavours

      if(.not.(leg_PDGs(i).eq.(-leg_PDGs(j)).and.abs(leg_PDGs(i)).ne.abs(leg_PDGs(k)).and.abs(leg_PDGs(k)).le.6)) return

c     initial checks and label assignment
      if(lsec.eq.0)then
         if((i.eq.isec.and.j.eq.jsec).or.(i.eq.jsec.and.j.eq.isec)) then
            ic = ksec
         elseif((i.eq.isec.and.j.eq.ksec).or.(i.eq.ksec.and.j.eq.isec)) then
            ic = jsec
         elseif((i.eq.jsec.and.j.eq.ksec).or.(i.eq.ksec.and.j.eq.jsec)) then
            ic = isec
         else
            write(*,*)'Wrong indices 1 in M2_C_SS_qqx_CC_qxqqp'
            write(*,*)i,j,isec,jsec,ksec
            stop
         endif
      else
         if((i.eq.isec.and.j.eq.jsec).or.(i.eq.jsec.and.j.eq.isec)) then
            ic = ksec
            id = lsec
         elseif((i.eq.ksec.and.j.eq.lsec).or.(i.eq.lsec.and.j.eq.ksec)) then
            ic = isec
            id = jsec
         else
            write(*,*)'Wrong indices 2 in M2_C_SS_qqx_CC_qxqqp'
            write(*,*)i,j,isec,jsec,ksec,lsec
            stop
         endif
      endif

c
c     TODO:check over k (why k and not ksec or viceversa?)
c     possible cuts
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-1,REAL_LEG_PDGS)
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)

      IF(DOCUT(XPBBSAVE,NEXTERNAL-2,BORN_LEG_PDGS,0))RETURN
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
      pref=64d0*pi**2*alphas**2
c
c     invariant quantities
      sij  = xs(i,j)
      sir  = xs(i,r)
      sjr  = xs(j,r)
c
c     safety check
c
      IF(SIJ.LE.0D0.OR.SIR.LE.0d0.OR.SJR.LE.0D0)THEN
        WRITE(77,*)'Inaccuracy 1 in M2_C_SS_qqx_CC_qxqqp',SIJ,SIR,SJR
        GOTO 999
      ENDIF
      zi   = sir/(sir+sjr)
      zj   = 1d0-zi
      call get_collinear_mapped_labels(i,j,nexternal,leg_PDGs,NLO_mapped_labels,NLO_mapped_flavours)
      if(NLO_mapped_flavours(j).ne.21)then
         write(*,*) 'Wrong parent particle label 1!', j, NLO_mapped_flavours(j)
         stop
      endif
c     Reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal,real_leg_pdgs,NLO_mapped_flavours,NLO_mapped_labels,xpbsave)
      call invariants_from_p(xpbsave,nexternal-1,xsb,ierr)
      if(ierr.eq.1)goto 999
      jb = NLO_mapped_labels(j)
      kb = NLO_mapped_labels(k)
      rb = NLO_mapped_labels(r)
      sbjr = xsb(jb,rb)
      sbkr = xsb(kb,rb)
      sbjk = xsb(jb,kb)
      zbj = sbjr/(sbjr+sbkr)
      zbk = 1d0 - zbj
      parent_leg = nlo_mapped_labels(jb)
c
      call get_collinear_mapped_labels(jb,kb,nexternal-1,real_leg_PDGs,LO_mapped_labels,LO_mapped_flavours)
      if(LO_mapped_flavours(kb).ne.NLO_mapped_flavours(k))then
         write(*,*) 'Wrong parent particle label 2!', jb,k,LO_mapped_flavours(jb),NLO_mapped_flavours(k)
         stop
      endif
c     Reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal-1,Born_leg_pdgs,LO_mapped_flavours,LO_mapped_labels,xpbbsave)
      call invariants_from_p(xpbbsave,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
c
c     Call Born
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbbsave,hel,alphas,ANS)
      BLO = ANS(0)
c
c     Formula (C.40) of 2212.11190v2 TODO: complete kernel (_y)
      M2tmp = !Pmunu_jk(r)/sbjk.....
      M2tmp = M2tmp * BLO
c
c     include double-collinear sector function Formula (C.83) TODO:  check indices - notation (_y)
      call get_wc_nlo(xs,ia,ib,ir,alphaz,nexternal)
      WC_SS_CC_NNLO = wc_nlo_bar
      M2TMP=M2TMP*WC_SS_CC_NNLO
c
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      M2_C_SS_qqx_CC_qxqqp=M2tmp*pref*xj*extra
c     apply flavour factor
      M2_C_SS_qqx_CC_qxqqp=M2_C_SS_qqx_CC_qxqqp*%(proc_prefix_rr)s_fl_factor
c
c     plot
      wgtpl=-M2_C_SS_qqx_CC_qxqqp*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbbsave,xsbb,nexternal-2,BORN_LEG_PDGS,wgtpl)
c
c     sanity check
      if(abs(M2_C_SS_qqx_CC_qxqqp).ge.huge(1d0).or.isnan(M2_C_SS_qqx_CC_qxqqp))then
         write(77,*)'Exception caught in M2_C_SS_qqx_CC_qxqqp', M2_C_SS_qqx_CC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


