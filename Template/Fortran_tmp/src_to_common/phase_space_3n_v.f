c$$$      subroutine phase_space_CS(xx,iU,iS,iB,iA,p,pbar,npart,leg_PDGs,maptype,xjac)
c$$$c     Build n+1 momenta p from n momenta pbar 
c$$$c     iU is the unresolved parton associated
c$$$c     with the soft singularity
c$$$c     npart = n+1
c$$$      implicit none
c$$$      include 'math.inc'
c$$$      double precision xx(3)
c$$$      integer i,j,iU,iS,iB,iA,npart,ierr
c$$$      double precision p(0:3,npart),pbar(0:3,npart-1),xjac
c$$$      double precision zCS,yCS,phCS
c$$$      double precision cosphi,sinphi,cosphiA,sinphiA
c$$$      double precision cosphiU,sinphiU,costhUB,sinthUB
c$$$      double precision Eu,pmod,pAtmod,G
c$$$      double precision nB(0:3),pbB(0:3),pbBsave(0:3)
c$$$      double precision pbS(0:3),pA(0:3),pboost(0:3)
c$$$      double precision pAr(0:3),pUr(0:3)
c$$$      double precision dot
c$$$      double precision ran2
c$$$      integer leg_PDGs(npart)
c$$$      character*1 maptype
c$$$      integer mapped_labels(npart),mapped_flavours(npart)
c$$$c     TODO: change name to isunderlyingqcdparton()
c$$$      logical isLOQCDparton(npart-1)
c$$$      integer idum
c$$$      common/rand/idum
c$$$c
c$$$c     initialise
c$$$      xjac=0d0
c$$$c
c$$$c     input check
c$$$      if(iU.eq.iS.or.iU.eq.iB.or.iU.eq.iA.or.
c$$$     &   iS.eq.iB.or.iS.eq.iA.or.iB.eq.iA)then
c$$$         write(*,*)'Wrong input indices in phase_space_CS'
c$$$         write(*,*)iU,iS,iB,iA
c$$$         stop
c$$$      endif
c$$$c
c$$$c     iA (reference four-vector for definition of the azimuth)
c$$$c     must be != iB, iU, iS;
c$$$      call get_mapped_labels(maptype,iU,iS,iB,npart,leg_pdgs,mapped_labels,
c$$$     $     mapped_flavours,isLOQCDparton)
c$$$      pbB(:)=pbar(:,mapped_labels(iB))
c$$$      pbS(:)=pbar(:,mapped_labels(iS))
c$$$      pA(:) =pbar(:,mapped_labels(iA))
c$$$      pbBsave=pbB
c$$$c
c$$$c     boost to CM of pbS+pbB
c$$$      pboost=pbS+pbB
c$$$      call boostinv(pbB,pboost,pbB,ierr)
c$$$      if(ierr.eq.1)goto 999
c$$$      call boostinv(pA,pboost,pA,ierr)
c$$$      if(ierr.eq.1)goto 999
c$$$c
c$$$c     Catani-Seymour parametrisation
c$$$      zCS =xx(1)
c$$$      yCS =xx(2)
c$$$      phCS=2d0*pi*xx(3)
c$$$c
c$$$c     rotate pA to the frame where p(iB) and pbar(iB) are along (0,0,-1)
c$$$      pmod=sqrt(pbB(1)**2+pbB(2)**2+pbB(3)**2)
c$$$      if(pmod.le.0d0)then
c$$$         write(77,*)'Wrong pmod in phase_space_CS',pmod
c$$$         goto 999
c$$$      endif
c$$$      nB=pbB/pmod
c$$$      call rotate_to_z(nB(1),pA(1),pAr(1))
c$$$c
c$$$c     build NLO unbarred momenta from LO barred ones and from (zCS,yCS,phCS).
c$$$c
c$$$c     phCS is the azimuth of pUr with respect to pAr, where pUr is p(iU)
c$$$c     rotated in the frame where p(iB) is along (0,0,-1)
c$$$      cosphi=cos(phCS)
c$$$      sinphi=sin(phCS)
c$$$c
c$$$c     transverse part of pAr
c$$$      pAtmod=sqrt(pAr(1)**2+pAr(2)**2)
c$$$      if(pAtmod.le.0d0)then
c$$$         write(77,*)'Wrong pAtmod in phase_space_CS',pAtmod
c$$$         goto 999
c$$$      endif
c$$$c
c$$$c     phiA is the azimuth between pAr and the x-z plane
c$$$      cosphiA=pAr(1)/pAtmod
c$$$      sinphiA=pAr(2)/pAtmod
c$$$c
c$$$c     phiU is the azimuth between pUr and the x-z plane
c$$$      cosphiU=cosphi*cosphiA-sinphi*sinphiA
c$$$      sinphiU=sinphi*cosphiA+cosphi*sinphiA
c$$$c
c$$$c     thUB is the polar angle between pUr and (0,0,1),
c$$$c     namely the angle between p(iU) and p(iB)
c$$$      if(1d0-(1d0-yCS)*(1d0-zCS).le.0d0.or.yCS*zCS*(1d0-zCS).le.0d0)then
c$$$         write(77,*)'Wrong CS variables'
c$$$         write(77,*)1d0-(1d0-yCS)*(1d0-zCS),yCS*zCS*(1d0-zCS),yCS,zCS
c$$$         goto 999
c$$$      endif
c$$$      costhUB=1d0-2d0*zCS/(1d0-(1d0-yCS)*(1d0-zCS))
c$$$      sinthUB=2d0*sqrt(yCS*(1d0-zCS)*zCS)/(1d0-(1d0-yCS)*(1d0-zCS))
c$$$c
c$$$c     Eu is the energy of p(iU)
c$$$      Eu=pmod*(1d0-(1d0-yCS)*(1d0-zCS))
c$$$c
c$$$c     construct rotated prU
c$$$      pUr(1)= Eu*sinthUB*cosphiU
c$$$      pUr(2)= Eu*sinthUB*sinphiU
c$$$      pUr(3)=-Eu*costhUB
c$$$c
c$$$c     rotate p(iU) back
c$$$      p(0,iU)=Eu
c$$$      call rotate_to_z_inv(nB(1),pUr(1),p(1,iU))
c$$$c
c$$$c     Boost back
c$$$      call boost(p(0,iU),pboost,p(0,iU),ierr)
c$$$      if(ierr.eq.1)goto 999
c$$$c
c$$$c     construct p from pbar
c$$$      p(:,iB)=(1d0-yCS)*pbBsave(:)
c$$$      p(:,iS)=yCS*pbBsave(:)+pbS(:)-p(:,iU)
c$$$      do j=1,npart
c$$$         if(j.eq.iU.or.j.eq.iB.or.j.eq.iS)cycle
c$$$         p(:,j)=pbar(:,mapped_labels(j))
c$$$      enddo
c$$$c
c$$$c     consistency check
c$$$      do i=1,npart
c$$$         if(p(0,i).lt.0d0)then
c$$$            write(77,*)'Inaccuracy 1 in phase_space_CS',i,p(0,i),npart
c$$$            goto 999
c$$$         endif
c$$$      enddo
c$$$c
c$$$c     construct xjac
c$$$      G=1d0/16d0/pi**3
c$$$      xjac=G*(4d0*pmod**2)*pi*(1-yCS)
c$$$c
c$$$      return
c$$$ 999  xjac=0d0
c$$$      return
c$$$      end


c$$$      subroutine phase_space_CS_inv(iU,iS,iB,p,pbar,npart,leg_PDGs,
c$$$     & maptype,xjac)
c$$$c     Build n momenta pbar from n+1 momenta p
c$$$c     iU is the unresolved parton associated
c$$$c     with the soft singularity
c$$$c
c$$$c     TODO: add mappings for initial state
c$$$c
c$$$      implicit none
c$$$      include 'math.inc'
c$$$      integer i,j,iU,iS,iB,npart
c$$$      double precision p(0:3,npart),pbar(0:3,npart-1)
c$$$      integer leg_PDGs(npart)
c$$$      double precision yCS,xjac,G,Qsq,dot
c$$$      character*1 maptype
c$$$      integer mapped_labels(npart),mapped_flavours(npart)
c$$$c     TODO: change name to isunderlyingqcdparton()
c$$$      logical isLOQCDparton(npart-1)
c$$$c
c$$$c     initialise
c$$$      xjac=0d0
c$$$      mapped_labels = 0
c$$$c
c$$$c     auxiliary quantities
c$$$      Qsq=2d0*(dot(p(0,iU),p(0,iS))+dot(p(0,iU),p(0,iB))+
c$$$     &         dot(p(0,iS),p(0,iB)))
c$$$      yCS=2d0*dot(p(0,iU),p(0,iS))/Qsq
c$$$c
c$$$c     construct pbar from p
c$$$      call get_mapped_labels(maptype,iU,iS,iB,npart,leg_pdgs,mapped_labels,
c$$$     $           mapped_flavours,isLOQCDparton)
c$$$c
c$$$      pbar(:,mapped_labels(iB))=p(:,iB)/(1-yCS)
c$$$      pbar(:,mapped_labels(iS))=p(:,iU)+p(:,iS)-p(:,iB)*yCS/(1-yCS)
c$$$      do j=1,npart
c$$$         if(j.eq.iU.or.j.eq.iB.or.j.eq.iS)cycle
c$$$         pbar(:,mapped_labels(j))=p(:,j)
c$$$      enddo
c$$$c
c$$$c     consistency check
c$$$      do i=1,npart-1
c$$$         if(pbar(0,i).lt.0d0)then
c$$$            write(77,*)'Inaccuracy 1 in phase_space_CS_inv',i,pbar(0,i),npart-1
c$$$            goto 999
c$$$         endif
c$$$      enddo
c$$$c
c$$$c     construct xjac
c$$$      G=1d0/16d0/pi**3
c$$$      xjac=G*Qsq*pi*(1-yCS)
c$$$c
c$$$      return
c$$$ 999  xjac=0d0
c$$$      end


      subroutine phase_space_n(x,shat,p,npart,xjac)
      implicit none
      include 'math.inc'
      integer npart
      double precision x(4),shat
      double precision p(0:3,npart),xjac
c
c     local variables
      double precision app,pmod,costh,sinth,ph,ran2
      double precision e_sup,e_inf,s1,q1(0:3),e1,p1mod
      integer idum
      common/rand/idum
      double precision xm(7),xmsq(7)
      common/c_masses/xm,xmsq
c
c     initialise
      p=0d0
      xm=0d0
      xmsq=0d0
      xjac=GeVtopb
c
c     s1 = (p1+p2)**2
      e_sup=sqrt(shat)-xm(3)
      e_inf=xm(1)+xm(2)      
      if(e_inf.ge.e_sup)then
        xjac=0d0
        return
      endif     
      e1=(e_sup-e_inf)*x(1)+e_inf
      s1=e1*e1
      xjac=xjac*(e_sup-e_inf)*2d0*e1
      xjac=xjac/twopi
c
c     q1 and p3
      app=((shat-s1-xmsq(3))**2-4d0*s1*xmsq(3))/(4d0*shat)
      if(app.le.0d0)then
         xjac=0d0
         return
      endif
      pmod=sqrt(app)
      costh=max(-1d0,min(x(2)*2d0-1d0,1d0))
      sinth=sqrt(1d0-costh*costh)
      ph=RAN2(idum)*twopi
c
c     final-state momenta
c     p3
      q1(0)=sqrt(pmod*pmod+s1)
      q1(1)=pmod*sinth*cos(ph)
      q1(2)=pmod*sinth*sin(ph)
      q1(3)=pmod*costh
      p(0,npart)=sqrt(pmod*pmod+xmsq(3))
      p(1,npart)=-q1(1)
      p(2,npart)=-q1(2)
      p(3,npart)=-q1(3)
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*fourpi
c
c     p1 and p2      
      app=((s1-xmsq(1)-xmsq(2))**2-4d0*xmsq(1)*xmsq(2))/(4d0*s1)
      if(app.le.0d0)then
        xjac=0d0
        return
      endif
      p1mod=sqrt(app)
      call twobody(q1,s1,p1mod,xmsq(1),xmsq(2),x(3),x(4),p(0,3),p(0,4))
c
c     initial-state and center-of-mass momenta
!      p(0,0)=sqrt(shat)
      p(0,2)=dsqrt(shat)/2d0
      p(3,2)=dsqrt(shat)/2d0
      p(0,1)=dsqrt(shat)/2d0
      p(3,1)=-dsqrt(shat)/2d0
c
c     phase-space weight
      xjac=xjac*p1mod/sqrt(s1)/(fourpi)**2
      xjac=xjac*fourpi
c
      return
      end



c$$$      subroutine phase_space_npo(x,shat,iU,iS,iB,iA,p,pbar,xjac,xjacB)
c$$$c     iU is the unresolved parton associated with the soft singularity
c$$$      implicit none
c$$$      include 'math.inc'
c$$$      include 'nexternal.inc'
c$$$      include 'leg_PDGs.inc'
c$$$      double precision x(7),shat
c$$$      double precision p(0:3,nexternal),pbar(0:3,nexternal-1)
c$$$      double precision xjac,xjacB,xjacCS
c$$$      integer i,j,iU,iS,iB,iA
c$$$c
c$$$c     initialise
c$$$      p=0d0
c$$$      pbar=0d0
c$$$      xjac=0d0
c$$$c
c$$$c     call Born phase space
c$$$      call phase_space_n(x(4),shat,pbar,nexternal-1,xjacB)
c$$$c
c$$$c     call radiation phase space
c$$$      call phase_space_CS(x(1),iU,iS,iB,iA,p,pbar,nexternal,leg_PDGs,'C',xjacCS)
c$$$c
c$$$c     total jacobian
c$$$      xjac=xjacB*xjacCS
c$$$c
c$$$      return
c$$$      end


      logical function dotechcut(s,ndim,tiny)
      implicit none
      integer ndim,i,j
      double precision s(ndim,ndim), tiny, smin

      dotechcut=.false.
      smin=1d40

      do i=1,ndim-1
         do j=i+1,ndim
            if(s(i,j).lt.smin) smin=s(i,j)
         enddo
      enddo

      if(smin/s(1,2).lt.tiny) dotechcut=.true.
      return
      end


c$$$      subroutine phase_space_npt(x,shat,iU1,iS1,iB1,iA1,iA2,p,pbar,ptilde,xjac,iU2,iS2,iB2)
c$$$c     iU1 and iU2 are the unresolved partons associated with the soft singularity
c$$$c     iU2, iS2, and iB2 are outputs
c$$$      implicit none
c$$$      include 'dims.inc'
c$$$      include 'setup.inc'
c$$$      include 'mxdim.inc'
c$$$      double precision x(10),shat
c$$$      double precision p(0:3,-2:5),pbar(0:3,-2:4),ptilde(0:3,-2:3)
c$$$      double precision xjac,xjacB,xjacCS1,xjacCS2
c$$$      integer i,j,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
c$$$      double precision dot
c$$$c
c$$$c     initialise
c$$$      p=0d0
c$$$      pbar=0d0
c$$$      ptilde=0d0
c$$$      xjac=1d0
c$$$      p(0,0)=sqrt(shat)
c$$$      pbar(0,0)=sqrt(shat)
c$$$c
c$$$c     call Born phase space
c$$$      call phase_space_n(x(7),shat,ptilde,xjacB)
c$$$c
c$$$      iU2=imap(iS1,iU1,iS1,0,0,npartNNLO)
c$$$      iS2=imap(iB1,iU1,iS1,0,0,npartNNLO)
c$$$      iB2=imap(iA1,iU1,iS1,0,0,npartNNLO)
c$$$c
c$$$c     call radiation phase space from Born to real
c$$$      call phase_space_CS(x(4),iU2,iS2,iB2,iA2,pbar,ptilde,npartNLO,xjacCS2)
c$$$c
c$$$c     call radiation phase space from real to double real
c$$$      call phase_space_CS(x(1),iU1,iS1,iB1,iA1,p,pbar,npartNNLO,xjacCS1)
c$$$c
c$$$c     total jacobian
c$$$      xjac=xjacB*xjacCS1*xjacCS2
c$$$c
c$$$      return
c$$$      end


c$$$      SUBROUTINE RAMBO(N,ET,XXM,P,WT)
c$$$***********************************************************************
c$$$*                       RAMBO                                         *
c$$$*    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)                  *
c$$$*                                                                     *
c$$$*    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR                *
c$$$*    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING                 *
c$$$*    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS                      *
c$$$*    -- ADJUSTED BY HANS KUIJF, WEIGHTS ARE LOGARITHMIC (20-08-90)    *
c$$$*                                                                     *
c$$$*    N  = NUMBER OF FINAL-STATE PARTICLES                             *
c$$$*    ET = TOTAL CENTRE-OF-MASS ENERGY                                 *
c$$$*    XXM = PARTICLE MASSES ( DIM=npartLO )                            *
c$$$*    P  = PARTICLE MOMENTA ( DIM=(4,npartLO) )                        *
c$$$*    WT = WEIGHT OF THE EVENT                                         *
c$$$***********************************************************************
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$      include 'setup.inc'
c$$$      DIMENSION XXM(npartLO),P(4,npartLO),Q(4,npartLO),Z(npartLO),R(4),
c$$$     .   B(3),P2(npartLO),XXM2(npartLO),E(npartLO),V(npartLO),IWARN(5)
c$$$      SAVE ACC,ITMAX,IBEGIN,IWARN
c$$$      DATA ACC/1.D-14/,ITMAX/6/,IBEGIN/0/,IWARN/5*0/
c$$$*
c$$$* INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
c$$$      IF(IBEGIN.NE.0) GOTO 103
c$$$      IBEGIN=1
c$$$      TWOPI=8.*DATAN(1.D0)
c$$$      PO2LOG=LOG(TWOPI/4.)
c$$$      Z(2)=PO2LOG
c$$$      DO 101 K=3,npartLO
c$$$  101 Z(K)=Z(K-1)+PO2LOG-2.*LOG(DFLOAT(K-2))
c$$$      DO 102 K=3,npartLO
c$$$  102 Z(K)=(Z(K)-LOG(DFLOAT(K-1)))
c$$$*
c$$$* CHECK ON THE NUMBER OF PARTICLES
c$$$  103 IF(N.GT.1.AND.N.LT.101) GOTO 104
c$$$      PRINT 1001,N
c$$$      STOP
c$$$*
c$$$* CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
c$$$  104 XXMT=0.
c$$$      NM=0
c$$$      DO 105 I=1,N
c$$$      IF(XXM(I).NE.0.D0) NM=NM+1
c$$$  105 XXMT=XXMT+ABS(XXM(I))
c$$$      IF(XXMT.LE.ET) GOTO 201
c$$$      PRINT 1002,XXMT,ET
c$$$      STOP
c$$$*
c$$$* THE PARAMETER VALUES ARE NOW ACCEPTED
c$$$*
c$$$* GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
c$$$  201 DO 202 I=1,N
c$$$         r1=rn(1)
c$$$      C=2.*r1-1.
c$$$      S=SQRT(1.-C*C)
c$$$      F=TWOPI*RN(2)
c$$$      r1=rn(3)
c$$$      r2=rn(4)
c$$$      Q(4,I)=-LOG(r1*r2)
c$$$      Q(3,I)=Q(4,I)*C
c$$$      Q(2,I)=Q(4,I)*S*COS(F)
c$$$  202 Q(1,I)=Q(4,I)*S*SIN(F)
c$$$*
c$$$* CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
c$$$      DO 203 I=1,4
c$$$  203 R(I)=0.
c$$$      DO 204 I=1,N
c$$$      DO 204 K=1,4
c$$$  204 R(K)=R(K)+Q(K,I)
c$$$      RMAS=SQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
c$$$      DO 205 K=1,3
c$$$  205 B(K)=-R(K)/RMAS
c$$$      G=R(4)/RMAS
c$$$      A=1./(1.+G)
c$$$      X=ET/RMAS
c$$$*
c$$$* TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
c$$$      DO 207 I=1,N
c$$$      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
c$$$      DO 206 K=1,3
c$$$  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
c$$$  207 P(4,I)=X*(G*Q(4,I)+BQ)
c$$$*
c$$$* CALCULATE WEIGHT AND POSSIBLE WARNINGS
c$$$      WT=PO2LOG
c$$$      IF(N.NE.2) WT=(2.*N-4.)*LOG(ET)+Z(N)
c$$$      IF(WT.GE.-180.D0) GOTO 208
c$$$      IF(IWARN(1).LE.5) PRINT 1004,WT
c$$$      IWARN(1)=IWARN(1)+1
c$$$  208 IF(WT.LE. 174.D0) GOTO 209
c$$$      IF(IWARN(2).LE.5) PRINT 1005,WT
c$$$      IWARN(2)=IWARN(2)+1
c$$$*
c$$$* RETURN FOR WEIGHTED MASSLESS MOMENTA
c$$$  209 IF(NM.NE.0) GOTO 210
c$$$* RETURN LOG OF WEIGHT
c$$$      WT=WT
c$$$      RETURN
c$$$*
c$$$* MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
c$$$  210 XMAX=SQRT(1.-(XXMT/ET)**2)
c$$$      DO 301 I=1,N
c$$$      XXM2(I)=XXM(I)**2
c$$$  301 P2(I)=P(4,I)**2
c$$$      ITER=0
c$$$      X=XMAX
c$$$      ACCU=ET*ACC
c$$$  302 F0=-ET
c$$$      G0=0.
c$$$      X2=X*X
c$$$      DO 303 I=1,N
c$$$      E(I)=SQRT(XXM2(I)+X2*P2(I))
c$$$      F0=F0+E(I)
c$$$  303 G0=G0+P2(I)/E(I)
c$$$      IF(ABS(F0).LE.ACCU) GOTO 305
c$$$      ITER=ITER+1
c$$$      IF(ITER.LE.ITMAX) GOTO 304
c$$$      PRINT 1006,ITMAX
c$$$      GOTO 305
c$$$  304 X=X-F0/(X*G0)
c$$$      GOTO 302
c$$$  305 DO 307 I=1,N
c$$$      V(I)=X*P(4,I)
c$$$      DO 306 K=1,3
c$$$  306 P(K,I)=X*P(K,I)
c$$$  307 P(4,I)=E(I)
c$$$*
c$$$* CALCULATE THE MASS-EFFECT WEIGHT FACTOR
c$$$      WT2=1.
c$$$      WT3=0.
c$$$      DO 308 I=1,N
c$$$      WT2=WT2*V(I)/E(I)
c$$$  308 WT3=WT3+V(I)**2/E(I)
c$$$      WTM=(2.*N-3.)*LOG(X)+LOG(WT2/WT3*ET)
c$$$*
c$$$* RETURN FOR  WEIGHTED MASSIVE MOMENTA
c$$$      WT=WT+WTM
c$$$      IF(WT.GE.-180.D0) GOTO 309
c$$$      IF(IWARN(3).LE.5) PRINT 1004,WT
c$$$      IWARN(3)=IWARN(3)+1
c$$$  309 IF(WT.LE. 174.D0) GOTO 310
c$$$      IF(IWARN(4).LE.5) PRINT 1005,WT
c$$$      IWARN(4)=IWARN(4)+1
c$$$* RETURN LOG OF WEIGHT
c$$$  310 WT=WT
c$$$      RETURN
c$$$*
c$$$ 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
c$$$ 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',
c$$$     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
c$$$ 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
c$$$ 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
c$$$ 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',
c$$$     . ' DESIRED ACCURACY =',D15.6)
c$$$      END
c$$$
c$$$      FUNCTION RN(IDUMMY)
c$$$      REAL*8 RN,RAN
c$$$      SAVE INIT
c$$$      DATA INIT /1/
c$$$      IF (INIT.EQ.1) THEN
c$$$        INIT=0
c$$$        CALL RMARIN(1802,9373)
c$$$      END IF
c$$$*
c$$$  10  CALL RANMAR(RAN)
c$$$      IF (RAN.LT.1D-16) GOTO 10
c$$$      RN=RAN
c$$$*
c$$$      END
c$$$
c$$$
c$$$
c$$$      SUBROUTINE RANMAR(RVEC)
c$$$*     -----------------
c$$$* Universal random number generator proposed by Marsaglia and Zaman
c$$$* in report FSU-SCRI-87-50
c$$$* In this version RVEC is a double precision variable.
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
c$$$      COMMON/ RASET2 / IRANMR,JRANMR
c$$$      SAVE /RASET1/,/RASET2/
c$$$      UNI = RANU(IRANMR) - RANU(JRANMR)
c$$$      IF(UNI .LT. 0D0) UNI = UNI + 1D0
c$$$      RANU(IRANMR) = UNI
c$$$      IRANMR = IRANMR - 1
c$$$      JRANMR = JRANMR - 1
c$$$      IF(IRANMR .EQ. 0) IRANMR = 97
c$$$      IF(JRANMR .EQ. 0) JRANMR = 97
c$$$      RANC = RANC - RANCD
c$$$      IF(RANC .LT. 0D0) RANC = RANC + RANCM
c$$$      UNI = UNI - RANC
c$$$      IF(UNI .LT. 0D0) UNI = UNI + 1D0
c$$$      RVEC = UNI
c$$$      END
c$$$ 
c$$$      SUBROUTINE RMARIN(IJ,KL)
c$$$*     -----------------
c$$$* Initializing routine for RANMAR, must be called before generating
c$$$* any pseudorandom numbers with RANMAR. The input values should be in
c$$$* the ranges 0<=ij<=31328 ; 0<=kl<=30081
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
c$$$      COMMON/ RASET2 / IRANMR,JRANMR
c$$$      SAVE /RASET1/,/RASET2/
c$$$* This shows correspondence between the simplified input seeds IJ, KL
c$$$* and the original Marsaglia-Zaman seeds I,J,K,L.
c$$$* To get the standard values in the Marsaglia-Zaman paper (i=12,j=34
c$$$* k=56,l=78) put ij=1802, kl=9373
c$$$      I = MOD( IJ/177 , 177 ) + 2
c$$$      J = MOD( IJ     , 177 ) + 2
c$$$      K = MOD( KL/169 , 178 ) + 1
c$$$      L = MOD( KL     , 169 )
c$$$      DO 300 II = 1 , 97
c$$$        S =  0D0
c$$$        T = .5D0
c$$$        DO 200 JJ = 1 , 24
c$$$          M = MOD( MOD(I*J,179)*K , 179 )
c$$$          I = J
c$$$          J = K
c$$$          K = M
c$$$          L = MOD( 53*L+1 , 169 )
c$$$          IF(MOD(L*M,64) .GE. 32) S = S + T
c$$$          T = .5D0*T
c$$$  200   CONTINUE
c$$$        RANU(II) = S
c$$$  300 CONTINUE
c$$$      RANC  =   362436D0 / 16777216D0
c$$$      RANCD =  7654321D0 / 16777216D0
c$$$      RANCM = 16777213D0 / 16777216D0
c$$$      IRANMR = 97
c$$$      JRANMR = 33
c$$$      END