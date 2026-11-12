C     arXiv:2212.11190v2, (6.12)-(6.19), Appendix E.
C     All five Laurent coefficients -4:0 in the prescription of
C     Section 7.2: positive NLO virtual orders omitted consistently
C     with VV (see footnote 12). No positive virtual orders required.
C     Same I2 structure: wrapper, core, Laurent helpers, kernel tables.
C     CORE returns no a_s**2 prefactor; wrapper multiplies it once.
C     V=VNLO/a_s, VC=VCNLO/a_s; a_s=alpha_s/(2*pi).
C     BC,B4,B3,VC use ORIGINAL leg labels; S,FLAV,REFS are compact.
C     Process wrapper: e+e- -> u ubar (two massless QCD Born legs).
C     General core also accepts explicit higher-multiplicity input.
C     Included one-loop hook calls V_ML5_1_1_SLOOPMATRIX_THRES:
C       IRV_NLO_ACCESSOR_HOOK(P,HEL,ALPHAS,VNLO,VCNLO,IERR)
C     VNLO(-2:0), VCNLO(-2:0,NEXTERNAL,NEXTERNAL) include a_s.
C     Error 1 invalid data; 2 unsupported process; 3 invalid mapping.
C     Last line (6.12) includes B, as implied by (6.7).
      SUBROUTINE INT_COUNTER_IRV_NNLO(P,SLO,IRVNNLO,IERR)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'double_virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix)s.inc'
      INCLUDE 'colored_partons.inc'
      INTEGER IERR
      DOUBLE PRECISION P(0:3,NEXTERNAL),SLO(NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION IRVNNLO(-4:0),VNLO(-2:0)
      DOUBLE PRECISION VCNLO(-2:0,NEXTERNAL,NEXTERNAL),ALPHAS
      IRVNNLO=0D0
      CALL SET_MUR_FROM_MOMENTA(P,NEXTERNAL,LEG_PDGS_%(proc_prefix)s)
      ALPHAS=ALPHAS_CURRENT
      VNLO=HUGE(1D0)
      VCNLO=HUGE(1D0)
      CALL IRV_NLO_ACCESSOR_HOOK(P,-1,ALPHAS,VNLO,VCNLO,IERR)
      IF (IERR.NE.0) RETURN
      CALL INT_COUNTER_IRV_NNLO_WITH_V(P,SLO,VNLO,VCNLO, IRVNNLO,IERR)
      END

      SUBROUTINE INT_COUNTER_IRV_NNLO_WITH_V(P,SLO,VNLO,VCNLO, IRVNNLO,IERR)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'double_virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix)s.inc'
      INCLUDE 'colored_partons.inc'
      INTEGER IERR,I,J,N,A,B,IP,REFS(NEXTERNAL)
      INTEGER LEGS(NEXTERNAL),FLAV(NEXTERNAL),POS(NEXTERNAL)
      INTEGER,PARAMETER :: HEL=-1
      DOUBLE PRECISION P(0:3,NEXTERNAL),SLO(NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION IRVNNLO(-4:0),ALPHAS,ALPHA_QCD,PREF,BORN
      DOUBLE PRECISION ANS(0:NSQSO_BORN),PMASS(NEXTERNAL)
      DOUBLE PRECISION SS(NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION BC(NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION B4(NEXTERNAL,NEXTERNAL,NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION COEFF(-4:0),MU2
      DOUBLE PRECISION GET_CCBLO,GET_QUADRUPLE_BLO
      EXTERNAL ALPHA_QCD,GET_CCBLO,GET_QUADRUPLE_BLO
      DOUBLE PRECISION VNLO(-2:0),VCNLO(-2:0,NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION V(-2:0),VC(-2:0,NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION B3(NEXTERNAL,NEXTERNAL,NEXTERNAL),AA
      INCLUDE 'pmass.inc'
      IERR=0
      IRVNNLO=0D0
      N=0
      POS=0
      REFS=0
      SS=0D0
      BC=0D0
      B4=0D0
      LEGS=0
      B3=0D0
      DO I=1,NEXTERNAL
        IF (.NOT.ISLOQCDPARTON(I)) CYCLE
        IF (I.LE.NINCOMING.OR.PMASS(I).NE.0D0) THEN
           IERR=2
           RETURN
        ENDIF
        IP=ABS(LEG_PDGS_%(proc_prefix)s(I))
        IF (IP.NE.21.AND.(IP.LT.1.OR.IP.GT.6)) THEN
          IERR=2
          RETURN
        ENDIF
        N=N+1
        LEGS(N)=I
        POS(I)=N
        FLAV(N)=IP
      ENDDO
      
c     It was assuming QCD particles in the initial state and
c     at least one QCD particle in the final state at LO
c     IF (N.LT.3) THEN
      IF (N.LT.2) THEN
        IERR=2
        RETURN
      ENDIF
      DO A=1,LEN_IREF
        I=IREF(1,A)
        J=IREF(2,A)
        IF (I.LT.1.OR.I.GT.NEXTERNAL) THEN
          IERR=3
          RETURN
       ENDIF
        IF (POS(I).EQ.0) CYCLE
        IF (J.LT.1.OR.J.GT.NEXTERNAL) THEN
          IERR=3
          RETURN
       ENDIF
       IF (POS(J).EQ.0.OR.I.EQ.J) THEN
          IERR=3
          RETURN
       ENDIF
        IF (REFS(POS(I)).NE.0) THEN
          IF (REFS(POS(I)).NE.POS(J)) THEN
            IERR=3
            RETURN
          ENDIF
       ENDIF
        REFS(POS(I))=POS(J)
      ENDDO
      DO A=1,N
        IF (REFS(A).EQ.0) THEN
          IERR=3
          RETURN
        ENDIF
        DO B=1,N
          IF (A.EQ.B) CYCLE
          SS(A,B)=SLO(LEGS(A),LEGS(B))
          IF (.NOT.(SS(A,B).GT.0D0.AND. SS(A,B).LT.HUGE(1D0))) THEN
          IERR=1
          RETURN
          ENDIF
        ENDDO
      ENDDO

C     renormalisation scale for the n-body kinematics
      CALL SET_MUR_FROM_MOMENTA(P,NEXTERNAL,LEG_PDGS_%(proc_prefix)s)
      ALPHAS=ALPHAS_CURRENT
      MU2=MUR_CURRENT**2
      AA=ALPHAS/(2D0*PI)
      IF (.NOT.(AA.GT.0D0.AND.AA.LT.HUGE(1D0))) THEN
         IERR=1
         RETURN
      ENDIF
      PREF=AA**2
      CALL ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BORN=ANS(0)

C     Current process has two coloured Born legs. Correlator arrays
C     keep ORIGINAL external labels; non-QCD slots remain zero.
      IF (N.NE.2) THEN
         IERR=2
         RETURN
      ENDIF
      IF (.NOT.(ABS(BORN).LT.HUGE(1D0))) THEN
         IERR=1
         RETURN
      ENDIF
      DO A=1,N
         I=LEGS(A)
         DO B=1,N
            IF (A.EQ.B) CYCLE
            J=LEGS(B)
            BC(I,J)=GET_CCBLO(I,J)
C     On the two-leg singlet: T_i.T_j=-C_i, so no division by B.
            IF (FLAV(A).EQ.21) THEN
               B4(I,J,I,J)=-2D0*CA*BC(I,J)
            ELSE
               B4(I,J,I,J)=-2D0*CF*BC(I,J)
            ENDIF
         ENDDO
      ENDDO
C     The loop provider includes its a_s; CORE, like I2, has none.
      V=VNLO/AA
      VC=0D0
      DO I=1,NEXTERNAL
         IF (.NOT.ISLOQCDPARTON(I)) CYCLE
         DO J=1,NEXTERNAL
            IF (.NOT.ISLOQCDPARTON(J)) CYCLE
            IF (I.EQ.J) CYCLE
            VC(:,I,J)=VCNLO(:,I,J)/AA
         ENDDO
      ENDDO
      CALL IRVNNLO_CORE(NEXTERNAL,N,LEGS,FLAV,REFS,SS,MU2, BORN,BC,B4,B3,V,VC,COEFF,IERR)
      IF (IERR.NE.0) RETURN
      DO IP=-4,0
         IF (.NOT.(ABS(PREF*COEFF(IP)).LT.HUGE(1D0))) THEN
            IERR=1
            IRVNNLO=0D0
            RETURN
         ENDIF
         IRVNNLO(IP)=PREF*COEFF(IP)
      ENDDO
      END

      subroutine IRVNNLO_CORE(nd,n,legs,flav,refs,s,mu2, born,bc,b4,b3,v,vc,res,ierr)
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      integer nd,n,legs(nd),flav(nd),refs(nd),ierr
      integer c,d,e,j,k,r,rp,id
      double precision s(nd,nd),mu2,born
      double precision bc(nd,nd),b4(nd,nd,nd,nd),b3(nd,nd,nd)
      double precision v(-2:0),vc(-2:0,nd,nd),res(-4:0)
      double precision logs(nd,nd),cas(nd),gh(nd),ph(nd)
      double precision t(-4:2,24),hc(-4:2),ht(-4:2)
      double precision da(-4:2),db(-4:2),xa(-4:2),xb(-4:2)
      double precision dc(-4:2),w(-4:2),q(-4:2),u(-4:2)
      double precision aa,beta,l,x,IRVNNLO_LI3NEGLOG
      logical IRVNNLO_VALID
      external IRVNNLO_VALID,IRVNNLO_LI3NEGLOG
      ierr=0
      res=0d0
      if (n.lt.2.or.n.gt.nd) then
         ierr=2
         return
      endif
      do c=1,n
         if (legs(c).lt.1.or.legs(c).gt.nd) then
            ierr=3
            return
         endif
         do d=1,c-1
            if (legs(c).eq.legs(d)) then
               ierr=3
               return
            endif
         enddo
      enddo
      if (.not.IRVNNLO_VALID(mu2).or.mu2.le.0d0.or. .not.IRVNNLO_VALID(ca).or.ca.le.0d0.or. .not.IRVNNLO_VALID(cf).or.cf.le.0d0.or. .not.IRVNNLO_VALID(tr).or.tr.le.0d0.or. .not.IRVNNLO_VALID(nf).or.nf.lt.0d0.or. .not.IRVNNLO_VALID(born)) then
         ierr=1
         return
      endif
      if (.not.all(abs(v).lt.huge(1d0))) then
         ierr=1
         return
      endif
c     Coupling factored out as in I2; V inputs are divided by a_s.
      aa=1d0
      beta=(11d0*ca-4d0*tr*nf)/3d0
      logs=0d0
      do c=1,n
         if (flav(c).eq.21) then
            cas(c)=ca
            gh(c)=beta/2d0-2d0*ca
            ph(c)=4d0*beta/3d0-16d0*ca/3d0
         elseif (abs(flav(c)).ge.1.and.abs(flav(c)).le.6) then
            cas(c)=cf
            gh(c)=-cf/2d0
            ph(c)=-cf
         else
            ierr=2
            return
         endif
         if (refs(c).lt.1.or.refs(c).gt.n.or.refs(c).eq.c) then
            ierr=3
            return
         endif
         do d=1,n
            if (d.eq.c) cycle
            if (.not.IRVNNLO_VALID(s(c,d)).or.s(c,d).le.0d0.or. .not.IRVNNLO_VALID(bc(legs(c),legs(d))).or. .not.IRVNNLO_VALID(b4(legs(c),legs(d),legs(c),legs(d))).or. .not.all(abs(vc(:,legs(c),legs(d))).lt.huge(1d0))) then
               ierr=1
               return
            endif
            logs(c,d)=log(s(c,d))-log(mu2)
            do e=1,n
               if (e.eq.c.or.e.eq.d) cycle
               if (.not.IRVNNLO_VALID(b3(legs(c),legs(d),legs(e)))) then
                  ierr=1
                  return
               endif
            enddo
         enddo
      enddo
      do id=1,24
         call IRVNNLO_TABLE(id,ca,cf,tr,nf,t(-4,id))
      enddo
c     Soft line of (6.12), with (6.13)-(6.15).
      do c=1,n
         do d=1,n
            if (d.eq.c) cycle
            l=logs(c,d)
            u=0d0
            u(-2:0)=vc(:,legs(c),legs(d))
            call IRVNNLO_PRODUCT(t(:,1),u,q)
            call IRVNNLO_ADD(res,q,l,-aa)
c     -JsRV*Bcd: the defining minus in (6.13) reverses the sign.
            call IRVNNLO_ADD(res,t(:,5),2d0*l,aa**2*ca*bc(legs(c), legs(d)))
            w=0d0
            w(-4:0)=t(-3:1,1)
            call IRVNNLO_ADD(res,w,l,aa**2*beta*bc(legs(c),legs(d))/2d0)
            call IRVNNLO_ADD(res,t(:,6),l, aa**2*2d0*cas(d)*bc(legs(c),legs(d)))
            w=t(:,7)-t(:,6)
            call IRVNNLO_ADD(res,w,l,-aa**2*b4(legs(c),legs(d), legs(c),legs(d)))
c     Tripole is already finite AFTER the complete ordered colour sum.
            do e=1,n
               if (e.eq.c.or.e.eq.d) cycle
               x=logs(c,e)-logs(d,e)
               res(0)=res(0)+2d0*pi*aa**2*b3(legs(c),legs(d),legs(e))* (x*l*l/2d0+x**3/6d0+IRVNNLO_LI3NEGLOG(x))
            enddo
         enddo
      enddo
c     Hard-collinear line, with (6.16)-(6.19).
      do j=1,n
         r=refs(j)
         if (flav(j).eq.21) then
            hc=nf*t(:,2)+t(:,4)/2d0
            ht=nf*t(:,8)+t(:,10)/2d0
            dc=nf*t(:,11)+t(:,13)/2d0
            da=nf*t(:,14)+t(:,18)/2d0
            db=nf*t(:,15)+t(:,19)/2d0
            xa=2d0*nf*t(:,20)+t(:,23)
            xb=2d0*nf*t(:,21)+t(:,24)
         else
            hc=t(:,3)
            ht=t(:,9)
            dc=t(:,12)
            da=t(:,16)
            db=t(:,17)
            xa=t(:,22)
c     XB(qg)+XB(gq) = 2CF * universal bracket, Eq. (E.20).
            call IRVNNLO_QXB(cf,xb)
         endif
         l=logs(j,r)
         u=0d0
         u(-2:0)=v
         call IRVNNLO_PRODUCT(hc,u,q)
         call IRVNNLO_ADD(res,q,l,aa)
         call IRVNNLO_ADD(res,ht,2d0*l,aa**2*born)
         w=0d0
         w(-4:0)=hc(-3:1)
         call IRVNNLO_ADD(res,w,l,-aa**2*beta*born/2d0)
         w=da+xa
         call IRVNNLO_ADD(res,w,l,-aa**2*cas(j)*born)
         w=db-xa
         call IRVNNLO_ADD(res,w,l,aa**2*bc(legs(j),legs(r)))
         w=db+xb
         do c=1,n
            if (c.eq.j.or.c.eq.r) cycle
            call IRVNNLO_ADD(res,w,logs(j,c),aa**2*bc(legs(j),legs(c)))
            call IRVNNLO_ADD(res,dc,l,aa**2*bc(legs(c),legs(r)))
         enddo
c     Recoiler difference in last line of (6.12), times B (6.7).
c     For two legs there is a unique NLO map. Its compensating
c     difference is identically zero; do not request a third leg.
         if (n.eq.2) cycle
         do k=1,n
            if (k.eq.j) cycle
            do rp=1,n
               if (rp.ne.j.and.rp.ne.k) exit
            enddo
            w=ph(k)*hc
            w(-4:0)=w(-4:0)+gh(k)*hc(-3:1)
            call IRVNNLO_ADD(res,w,l,aa**2*born)
            call IRVNNLO_ADD(res,w,logs(j,rp),-aa**2*born)
         enddo
      enddo
      if (.not.all(abs(res).lt.huge(1d0))) then
         ierr=1
         res=0d0
      endif
      end

      logical function IRVNNLO_VALID(x)
      implicit none
      double precision x
      IRVNNLO_VALID=abs(x).lt.huge(1d0)
      end

c     Keep positive orders until AFTER convolution with V or 1/eps.
      subroutine IRVNNLO_PRODUCT(a,b,c)
      implicit none
      double precision a(-4:2),b(-4:2),c(-4:2)
      integer i,j
      c=0d0
      do i=-4,2
         do j=-4,2
            if (i+j.lt.-4.or.i+j.gt.0) cycle
            c(i+j)=c(i+j)+a(i)*b(j)
         enddo
      enddo
      end

      subroutine IRVNNLO_ADD(res,c,l,weight)
      implicit none
      double precision res(-4:0),c(-4:2),l,weight,term
      integer i,j
      do i=-4,0
         term=weight*c(i)
         res(i)=res(i)+term
         do j=i+1,0
            term=term*(-l)/dble(j-i)
            res(j)=res(j)+term
         enddo
      enddo
      end

c     Li3(-exp(x)) without exponentiating a possibly huge positive x.
c     Inversion: Li3(-e^x)-Li3(-e^-x)=-(x^3+pi^2*x)/6.
c     For |z|<=1 an Euler transform of the alternating series is
c     used directly, requiring at most 64 forward differences.
      double precision function IRVNNLO_LI3NEGLOG(x)
      implicit none
      double precision x,q,a(64),sum,weight,pi
      integer k,m
      pi=acos(-1d0)
      q=exp(-abs(x))
      do k=1,64
         a(k)=q**k/dble(k)**3
      enddo
      sum=0d0
      weight=0.5d0
      do m=1,64
         sum=sum+weight*a(1)
         do k=1,64-m
            a(k)=a(k)-a(k+1)
         enddo
         weight=weight/2d0
      enddo
      IRVNNLO_LI3NEGLOG=-sum
      if (x.gt.0d0) IRVNNLO_LI3NEGLOG= -sum-(x**3+pi*pi*x)/6d0
      end

      subroutine IRVNNLO_QXB(cf,c)
      implicit none
      double precision cf,c(-4:2),p2,z3
      p2=acos(-1d0)**2
      z3=1.2020569031595942854d0
      c=0d0
      c(-3)=-cf
      c(-2)=-2d0*cf
      c(-1)=2d0*cf*(-3d0+7d0*p2/24d0)
      c(0)=2d0*cf*(-9d0+7d0*p2/12d0+25d0*z3/6d0)
      end
c     IDs 1..4: E.2/E.8; 5..7: E.6; 8..10: E.17;
c     11..13: E.18; 14..19: E.19; 20..24: E.20.
c     ID22 combines Delta_hc,A(qg)+Delta_hc,A(gq).
c     IDs 5,8..10 have exp(-2 eps L); others exp(-eps L).
      subroutine IRVNNLO_TABLE(id,ca,cf,tr,nf,c)
      implicit none
      integer id
      double precision ca,cf,tr,nf,c(-4:2),p2,p4,z3
      p2=acos(-1d0)**2
      p4=p2*p2
      z3=1.2020569031595942854d0
      c=0d0
      if (id.ge.1.and.id.le.4) then
         call IRVNNLO_NLO_TABLE(id,ca,cf,tr,c)
         return
      endif
      select case(id)
      case(5)
         c(-4) = 1d0/4d0
         c(-3) = 1d0
         c(-2) = 4d0-7d0*p2/24d0
         c(-1) = 16d0-7d0*p2/6d0-14d0*z3/3d0
         c(0) = 64d0-14d0*p2/3d0-56d0*z3/3d0-7d0*p4/480d0
      case(6)
         c(-2) = 2d0-p2/3d0
         c(-1) = 16d0-2d0*p2/3d0-12d0*z3
         c(0) = 92d0-7d0*p2/2d0-24d0*z3-7d0*p4/18d0
      case(7)
         c(-2) = 2d0-p2/3d0
         c(-1) = 14d0-2d0*p2/3d0-10d0*z3
         c(0) = 74d0-23d0*p2/6d0-20d0*z3-7d0*p4/36d0
      case(8)
         c(-3) = tr*(2d0*cf/3d0-ca/3d0)
         c(-2) = tr*(4d0*nf*tr/9d0+31d0*cf/9d0-31d0*ca/18d0)
         c(-1) = tr*(64d0*nf*tr/27d0+cf*(431d0/27d0-p2)+ca*(-211d0 /27d0+p2/2d0))
         c(0) = tr*(nf*tr*(284d0/27d0-2d0*p2/3d0)+cf*(5506d0/81d0 -31d0*p2/6d0-124d0*z3/9d0)+ca*(-5281d0/162d0+31d0*p2/12d0+62d0 *z3/9d0))
      case(9)
         c(-3) = cf*(ca/4d0)
         c(-2) = cf*(cf*(-5d0/4d0+p2/3d0)+ca/2d0)
         c(-1) = cf*(cf*(-15d0/2d0+2d0*p2/3d0+10d0*z3)+ca*(1d0-p2 /24d0-4d0*z3))
         c(0) = cf*(cf*(-141d0/4d0+109d0*p2/24d0+20d0*z3-7d0*p4/45d0) +ca*(9d0/4d0+7d0*p2/12d0-67d0*z3/6d0-11d0*p4/90d0))
      case(10)
         c(-3) = ca*(ca/6d0)
         c(-2) = ca*(ca*(-28d0/9d0+2d0*p2/3d0))
         c(-1) = ca*(nf*tr/3d0+ca*(-61d0/3d0+7d0*p2/4d0+12d0*z3))
         c(0) = ca*(25d0*nf*tr/9d0+ca*(-15805d0/162d0+38d0*p2/3d0 +221d0*z3/9d0-5d0*p4/9d0))
      case(11)
         c(-1) = tr*(-4d0/3d0+2d0*p2/9d0)
         c(0) = tr*(-104d0/9d0+16d0*p2/27d0+8d0*z3)
      case(12)
         c(-1) = cf*(-1d0+p2/6d0)
         c(0) = cf*(-8d0+p2/3d0+6d0*z3)
      case(13)
         c(-1) = ca*(-2d0/3d0+p2/9d0)
         c(0) = ca*(-52d0/9d0+8d0*p2/27d0+4d0*z3)
      case(14)
         c(-3) = tr*(-4d0/3d0)
         c(-2) = tr*(-32d0/9d0)
         c(-1) = tr*(-280d0/27d0+7d0*p2/9d0)
         c(0) = tr*(-2504d0/81d0+56d0*p2/27d0+100d0*z3/9d0)
      case(15)
         c(-3) = tr*(4d0/3d0)
         c(-2) = tr*(32d0/9d0)
         c(-1) = tr*(244d0/27d0-5d0*p2/9d0)
         c(0) = tr*(1784d0/81d0-40d0*p2/27d0-52d0*z3/9d0)
      case(16)
         c(-3) = cf*(-1d0)
         c(-2) = cf*(-6d0+2d0*p2/3d0)
         c(-1) = cf*(-30d0+23d0*p2/12d0+16d0*z3)
         c(0) = cf*(-126d0+49d0*p2/6d0+121d0*z3/3d0+p4/9d0)
      case(17)
         c(-3) = cf*(1d0)
         c(-2) = cf*(2d0)
         c(-1) = cf*(5d0-5d0*p2/12d0)
         c(0) = cf*(12d0-5d0*p2/6d0-13d0*z3/3d0)
      case(18)
         c(-3) = ca*(-2d0/3d0)
         c(-2) = ca*(-88d0/9d0+4d0*p2/3d0)
         c(-1) = ca*(-1436d0/27d0+55d0*p2/18d0+32d0*z3)
         c(0) = ca*(-18748d0/81d0+406d0*p2/27d0+626d0*z3/9d0+2d0*p4/9d0)
      case(19)
         c(-3) = ca*(2d0/3d0)
         c(-2) = ca*(16d0/9d0)
         c(-1) = ca*(122d0/27d0-5d0*p2/18d0)
         c(0) = ca*(892d0/81d0-20d0*p2/27d0-26d0*z3/9d0)
      case(20)
         c(-3) = tr*(2d0/3d0)
         c(-2) = tr*(16d0/9d0)
         c(-1) = tr*(122d0/27d0-4d0*p2/9d0)
         c(0) = tr*(1108d0/81d0-44d0*p2/27d0-47d0*z3/9d0)
      case(21)
         c(-3) = tr*(-2d0/3d0)
         c(-2) = tr*(-16d0/9d0)
         c(-1) = tr*(-140d0/27d0+7d0*p2/18d0)
         c(0) = tr*(-1252d0/81d0+28d0*p2/27d0+50d0*z3/9d0)
      case(22)
         c(-3) = cf
         c(-2) = 2d0*cf
         c(-1) = (2d0*cf-ca)*(7d0/2d0-11d0*p2/24d0)+ca*(5d0/2d0-7d0 *p2/24d0)
         c(0) = (2d0*cf-ca)*(23d0/2d0-2d0*p2/3d0-32d0*z3/3d0)+ca *(13d0/2d0-5d0*p2/6d0-5d0*z3/3d0)
      case(23)
         c(-3) = ca*(1d0/3d0)
         c(-2) = ca*(8d0/9d0)
         c(-1) = ca*(88d0/27d0-11d0*p2/36d0)
         c(0) = ca*(716d0/81d0-17d0*p2/54d0-64d0*z3/9d0)
      case(24)
         c(-3) = ca*(-1d0/3d0)
         c(-2) = ca*(-8d0/9d0)
         c(-1) = ca*(-70d0/27d0+7d0*p2/36d0)
         c(0) = ca*(-626d0/81d0+14d0*p2/27d0+25d0*z3/9d0)
      end select
      end

      subroutine IRVNNLO_NLO_TABLE(id,ca,cf,tr,c)
      implicit none
      integer id
      double precision ca,cf,tr,c(-4:2),p2,p4,z3
      p2=acos(-1d0)**2
      p4=p2**2
      z3=1.2020569031595942854d0
      c=0d0
      select case(id)
c     E.2: J_s
      case(1)
         c(-2) = 1d0
         c(-1) = 2d0
         c(0) = 6d0 - 7d0 * p2 / 12d0
         c(1) = 18d0 - 7d0 * p2 / 6d0 - 25d0 * z3 / 3d0
         c(2) = 54d0 - 7d0 * p2 / 2d0 - 50d0 * z3 / 3d0 - 71d0 * p4 /1440d0
c     E.8: J_hc^(0g)
      case(2)
         c(-1) =  - 2d0 * tr / 3d0
         c(0) =  - 16d0 * tr / 9d0
         c(1) =  - tr *  ( 140d0 / 27d0 - 7d0 * p2 / 18d0 )
         c(2) =  - tr *  ( 1252d0 / 81d0 - 28d0 * p2 / 27d0 - 50d0 * z3 /9d0 )
c     E.8: J_hc^(1g)
      case(3)
         c(-1) =  - cf / 2d0
         c(0) =  - cf
         c(1) =  - cf *  ( 3d0 - 7d0 * p2 / 24d0 )
         c(2) =  - cf *  ( 9d0 - 7d0 * p2 / 12d0 - 25d0 * z3 / 6d0 )
c     E.8: J_hc^(2g)
      case(4)
         c(-1) =  - ca / 3d0
         c(0) =  - 8d0 * ca / 9d0
         c(1) =  - ca *  ( 70d0 / 27d0 - 7d0 * p2 / 36d0 )
         c(2) =  - ca *  ( 626d0 / 81d0 - 14d0 * p2 / 27d0 - 25d0 * z3 /9d0 )
      end select
      end

C     Process adapter: same MadLoop call and Laurent ordering as
C     NLO_V_%(proc_prefix)s.f. Physical a_s is already included in VNLO.
C     For the qqbar singlet VC_ij=-CF*V for its two QCD labels.
      SUBROUTINE IRV_NLO_ACCESSOR_HOOK(P,HEL,ALPHAS, VNLO,VCNLO,IERR)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'nsquaredSO.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'double_virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix)s.inc'
      INCLUDE 'colored_partons.inc'
      INTEGER HEL,IERR,I,J,N,LEGS(2),RETURNCODE
      INTEGER MATELEM_ARRAY_DIM,NSQUAREDSO_LOOP
      INTEGER NCOLORCORRELATORS
      PARAMETER (NCOLORCORRELATORS=4)
      INTEGER COLOR_CORRELATORS_TO_CONSIDER(0:NCOLORCORRELATORS)
      DOUBLE PRECISION COLOR_CORRELATED_EVALS(NCOLORCORRELATORS, 0:3,0:NSQUAREDSO)
      COMMON/V_ML5_1_1_COLOR_CORRELATIONS/ COLOR_CORRELATORS_TO_CONSIDER,COLOR_CORRELATED_EVALS
      DOUBLE PRECISION P(0:3,NEXTERNAL),ALPHAS,VNLO(-2:0)
      DOUBLE PRECISION VCNLO(-2:0,NEXTERNAL,NEXTERNAL)
      DOUBLE PRECISION,ALLOCATABLE :: MATELEM(:,:),PREC_FOUND(:)
      LOGICAL LOOP_ARRAYS_INITIALISED
      DATA LOOP_ARRAYS_INITIALISED/.FALSE./
      SAVE LOOP_ARRAYS_INITIALISED,MATELEM,PREC_FOUND
      SAVE MATELEM_ARRAY_DIM,NSQUAREDSO_LOOP
      IERR=0
      VNLO=0D0
      VCNLO=0D0
      IF (HEL.NE.-1) THEN
         IERR=2
         RETURN
      ENDIF
      IF (.NOT.(ALPHAS.GT.0D0.AND.ALPHAS.LT.HUGE(1D0))) THEN
         IERR=1
         RETURN
      ENDIF
      N=0
      DO I=1,NEXTERNAL
         IF (.NOT.ISLOQCDPARTON(I)) CYCLE
         N=N+1
         IF (N.GT.2) THEN
            IERR=2
            RETURN
         ENDIF
         LEGS(N)=I
      ENDDO
      IF (N.NE.2) THEN
         IERR=2
         RETURN
      ENDIF
      I=LEGS(1)
      J=LEGS(2)
      IF (ABS(LEG_PDGS_%(proc_prefix)s(I)).NE.2.OR. LEG_PDGS_%(proc_prefix)s(J).NE.-LEG_PDGS_%(proc_prefix)s(I)) THEN
         IERR=2
         RETURN
      ENDIF
      IF (.NOT.LOOP_ARRAYS_INITIALISED) THEN
         CALL V_ML5_1_1_GET_ANSWER_DIMENSION(MATELEM_ARRAY_DIM)
         CALL V_ML5_1_1_GET_NSQSO_LOOP(NSQUAREDSO_LOOP)
         IF (MATELEM_ARRAY_DIM.LT.0.OR.NSQUAREDSO_LOOP.LT.0) THEN
            IERR=4
            RETURN
         ENDIF
         ALLOCATE(MATELEM(0:3,0:MATELEM_ARRAY_DIM))
         ALLOCATE(PREC_FOUND(0:NSQUAREDSO_LOOP))
         LOOP_ARRAYS_INITIALISED=.TRUE.
      ENDIF
      MATELEM=HUGE(1D0)
      PREC_FOUND=0D0
      COLOR_CORRELATED_EVALS=0D0
      CALL V_ML5_1_1_SLOOPMATRIX_THRES(P,MATELEM,-1D0, PREC_FOUND,RETURNCODE)
      VNLO(-2)=MATELEM(3,0)
      VNLO(-1)=MATELEM(2,0)
      VNLO(0)=MATELEM(1,0)
c     From MadLoop to LASS
      VNLO(0) = VNLO(0) - VNLO(-2)*pi**2/12d0
      IF (.NOT.ALL(ABS(VNLO).LT.HUGE(1D0))) THEN
         IERR=4
         VNLO=0D0
         RETURN
      ENDIF
      VCNLO(:,I,J)=-CF*VNLO
      VCNLO(:,J,I)=-CF*VNLO
      END
