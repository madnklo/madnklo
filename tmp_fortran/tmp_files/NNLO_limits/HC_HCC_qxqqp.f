

      double precision function M2_HC_HCC_qxqqp(i,j,k,r,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     hard-collinear limit C_(ia,ib)C_(ia,ib,ic)
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
      INTEGER REAL_LEG_PDGS(NEXTERNAL-1)
      double precision sij, sir, sjr, sbjk, sbjr 
      double precision zi, zj
      double precision zbj, zbk
c
c     initialise
      M2_HC_HCC_qxqqp=0d0
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
            write(*,*)'Wrong indices 1 in M2_HC_HCC_qxqqp'
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
            write(*,*)'Wrong indices 2 in M2_HC_HCC_qxqqp'
            write(*,*)i,j,isec,jsec,ksec,lsec
            stop
         endif
      endif

c
c     TODO:check over k (why k and not ksec or viceversa?)
c     possible cuts
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-1,REAL_LEG_PDGS)
      call GET_UNDERLYING_PDGS(I,J,KSEC,LSEC,NEXTERNAL-2,BORN_LEG_PDGS)

      IF(DOCUT(XPBB,NEXTERNAL-2,BORN_LEG_PDGS,0))RETURN
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

      IF(SIJ.LE.0D0.OR.SIR.LE.0d0.OR.SJR.LE.0D0)THEN
        WRITE(77,*)'Inaccuracy 1 in M2_HC_HCC_qxqqp',SIJ,SIR,SJR
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
      call reshuffle_momenta(nexternal,real_leg_pdgs,NLO_mapped_flavours,NLO_mapped_labels,xpb)
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
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
      if(LO_mapped_flavours(jb).ne.NLO_mapped_flavours(k))then
         write(*,*) 'Wrong parent particle label 2!', jb,k,LO_mapped_flavours(jb),NLO_mapped_flavours(k)
         stop
      endif
c     Reshuffle momenta and labels according to underlying_leg_pdgs
      call reshuffle_momenta(nexternal-1,Born_leg_pdgs,LO_mapped_flavours,LO_mapped_labels,xpbb)
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.eq.1)goto 999
C
C     Call Born
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,hel,alphas,ANS)
      BLO = ANS(0)
c
c     Formula (C.45) of 2212.11190v2
      M2tmp = TR*(1d0-2d0*zi*zj)/sij * CF*zbj/sbjk
      M2tmp = M2tmp * BLO
c
C     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
c
      M2_HC_HCC_qxqqp=M2tmp*pref*xj*extra
c     apply flavour factor
      M2_HC_HCC_qxqqp=M2_HC_HCC_qxqqp*%(proc_prefix_rr)s_fl_factor
c
c     plot
      wgtpl=-M2_HC_HCC_qxqqp*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpbb,xsbb,nexternal-2,wgtpl)
c
c     sanity check
      if(abs(M2_HC_HCC_qxqqp).ge.huge(1d0).or.isnan(M2_HC_HCC_qxqqp))then
         write(77,*)'Exception caught in M2_HC_HCC_qxqqp',M2_HC_HCC_qxqqp
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

