      subroutine phase_space_CS(xx,iU,iS,iB,iA,p,pbar,npart,xjac)
c     Build n+1 momenta p from n momenta pbar 
c     iU is the unresolved parton associated
c     with the soft singularity
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      double precision xx(3)
      double precision p(0:3,-2:npart),pbar(0:3,-2:npart-1),xjac
      double precision zCS,yCS,phCS
      double precision cosphi,sinphi,cosphiA,sinphiA
      double precision cosphiU,sinphiU,costhUB,sinthUB
      double precision Eu,pmod,pAtmod,G
      integer i,j,iU,iS,iB,iA,npart,ierr
      double precision nB(0:3),pbB(0:3),pbBsave(0:3)
      double precision pbS(0:3),pA(0:3),pboost(0:3)
      double precision pAr(0:3),pUr(0:3)
      double precision dot
      double precision ran2
      integer idum
      common/rand/idum
c
c     initialise
      xjac=0d0
c
c     input check
      if(iU.eq.iS.or.iU.eq.iB.or.iU.eq.iA.or.
     &   iS.eq.iB.or.iS.eq.iA.or.iB.eq.iA)then
         write(*,*)'Wrong input indices in phase_space_CS'
         write(*,*)iU,iS,iB,iA
         stop
      endif
c
c     iA (reference four-vector for definition of the azimuth)
c     must be != iB, iU, iS;
      do i=0,3
         pbB(i)=pbar(i,imap(iB,iS,iU,0,0,npart))
         pbS(i)=pbar(i,imap(iS,iS,iU,0,0,npart))
         pA(i) =pbar(i,imap(iA,iS,iU,0,0,npart))
      enddo
      pbBsave=pbB
c
c     boost to CM of pbS+pbB
      pboost=pbS+pbB
      call boostinv(pbB,pboost,pbB,ierr)
      if(ierr.eq.1)goto 999
      call boostinv(pA,pboost,pA,ierr)
      if(ierr.eq.1)goto 999
c
c     Catani-Seymour parametrisation
      zCS =xx(1)
      yCS =xx(2)
      phCS=2d0*pi*xx(3)
c
c     rotate pA to the frame where p(iB) and pbar(iB) are along (0,0,-1)
      pmod=sqrt(pbB(1)**2+pbB(2)**2+pbB(3)**2)
      if(pmod.le.0d0)then
         write(77,*)'Wrong pmod in phase_space_CS',pmod
         goto 999
      endif
      nB=pbB/pmod
      call rotate_to_z(nB(1),pA(1),pAr(1))
c
c     build NLO unbarred momenta from LO barred ones and from (zCS,yCS,phCS).
c
c     phCS is the azimuth of pUr with respect to pAr, where pUr is p(iU)
c     rotated in the frame where p(iB) is along (0,0,-1)
      cosphi=cos(phCS)
      sinphi=sin(phCS)
c
c     transverse part of pAr
      pAtmod=sqrt(pAr(1)**2+pAr(2)**2)
      if(pAtmod.le.0d0)then
         write(77,*)'Wrong pAtmod in phase_space_CS',pAtmod
         goto 999
      endif
c
c     phiA is the azimuth between pAr and the x-z plane
      cosphiA=pAr(1)/pAtmod
      sinphiA=pAr(2)/pAtmod
c
c     phiU is the azimuth between pUr and the x-z plane
      cosphiU=cosphi*cosphiA-sinphi*sinphiA
      sinphiU=sinphi*cosphiA+cosphi*sinphiA
c
c     thUB is the polar angle between pUr and (0,0,1),
c     namely the angle between p(iU) and p(iB)
      if(1d0-(1d0-yCS)*(1d0-zCS).le.0d0.or.yCS*zCS*(1d0-zCS).le.0d0)then
         write(77,*)'Wrong CS variables'
         write(77,*)1d0-(1d0-yCS)*(1d0-zCS),yCS*zCS*(1d0-zCS),yCS,zCS
         goto 999
      endif
      costhUB=1d0-2d0*zCS/(1d0-(1d0-yCS)*(1d0-zCS))
      sinthUB=2d0*sqrt(yCS*(1d0-zCS)*zCS)/(1d0-(1d0-yCS)*(1d0-zCS))
c
c     Eu is the energy of p(iU)
      Eu=pmod*(1d0-(1d0-yCS)*(1d0-zCS))
c
c     construct rotated prU
      pUr(1)= Eu*sinthUB*cosphiU
      pUr(2)= Eu*sinthUB*sinphiU
      pUr(3)=-Eu*costhUB
c
c     rotate p(iU) back
      p(0,iU)=Eu
      call rotate_to_z_inv(nB(1),pUr(1),p(1,iU))
c
c     Boost back
      call boost(p(0,iU),pboost,p(0,iU),ierr)
      if(ierr.eq.1)goto 999
c
c     construct p from pbar
      do i=0,3
         p(i,iB)=(1d0-yCS)*pbBsave(i)
         p(i,iS)=yCS*pbBsave(i)+pbS(i)-p(i,iU)
         p(i,0)=pbar(i,0)
         do j=-2,npart
            if(j.eq.iU.or.j.eq.iB.or.j.eq.iS.or.j.eq.0)cycle
            p(i,j)=pbar(i,imap(j,iS,iU,0,0,npart))
         enddo
      enddo
c
c     consistency check
      do i=-2,npart
         if(p(0,i).lt.0d0)then
            write(77,*)'Inaccuracy 1 in phase_space_CS',i,p(0,i),npart
            goto 999
         endif
      enddo
c
c     construct xjac
      G=1d0/16d0/pi**3
      xjac=G*(4d0*pmod**2)*pi*(1-yCS)
c
      return
 999  xjac=0d0
      return
      end


      subroutine phase_space_CS_inv(iU,iS,iB,p,pbar,npart,xjac)
c     Build n momenta pbar from n+1 momenta p
c     iU is the unresolved parton associated
c     with the soft singularity
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      double precision p(0:3,-2:npart),pbar(0:3,-2:npart-1)
      double precision yCS,xjac,G,Qsq,dot
      integer i,j,iU,iS,iB,npart
c
c     initialise
      xjac=0d0
c
c     auxiliary quantities
      Qsq=2d0*(dot(p(0,iU),p(0,iS))+dot(p(0,iU),p(0,iB))+
     &         dot(p(0,iS),p(0,iB)))
      yCS=2d0*dot(p(0,iU),p(0,iS))/Qsq
c
c     construct pbar from p
      do i=0,3
         pbar(i,imap(iB,iS,iU,0,0,npart))=p(i,iB)/(1-yCS)
         pbar(i,imap(iS,iS,iU,0,0,npart))=p(i,iU)+p(i,iS)-p(i,iB)*yCS/(1-yCS)
         pbar(i,0)=p(i,0)
         do j=-2,npart
            if(j.eq.iU.or.j.eq.iB.or.j.eq.iS.or.j.eq.0)cycle
            pbar(i,imap(j,iS,iU,0,0,npart))=p(i,j)
         enddo
      enddo
c
c     consistency check
      do i=-2,npart-1
         if(pbar(0,i).lt.0d0)then
            write(77,*)'Inaccuracy 1 in phase_space_CS_inv',i,pbar(0,i),npart-1
            goto 999
         endif
      enddo
c
c     construct xjac
      G=1d0/16d0/pi**3
      xjac=G*Qsq*pi*(1-yCS)
c
      return
 999  xjac=0d0
      end


      subroutine phase_space_n(x,shat,p,xjac)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      double precision x(1),shat
      double precision p(0:3,-2:2),xjac
c
c     local variables
      double precision app,pmod,costh,sinth,ph,ran2
      integer idum
      common/rand/idum
c
c     initialise
      p=0d0
      xjac=1d0
c
c     app is modulus squared of three-momentum
      app=((shat-xmsq(1)-xmsq(2))**2-4d0*xmsq(1)*xmsq(2))/(4d0*shat)
      if(app.le.0d0)then
        xjac=0d0
        return
      endif
      pmod=sqrt(app)
      costh=max(-1d0,min(x(1)*2d0-1d0,1d0))
      sinth=sqrt(1d0-costh*costh)
      ph=RAN2(idum)*twopi
c
c     final-state momenta
      p(0,1)=sqrt(pmod*pmod+xmsq(1))
      p(1,1)=pmod*sinth*cos(ph)
      p(2,1)=pmod*sinth*sin(ph)
      p(3,1)=pmod*costh
      p(0,2)=sqrt(pmod*pmod+xmsq(2))
      p(1,2)=-p(1,1)
      p(2,2)=-p(2,1)
      p(3,2)=-p(3,1)
c
c     initial-state and center-of-mass momenta
      p(0,0)=sqrt(shat)
      p(0,-1)=p(0,0)/2d0
      p(3,-1)=p(0,0)/2d0
      p(0,-2)=p(0,0)/2d0
      p(3,-2)=-p(0,0)/2d0
c
c     phase-space weight
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*fourpi
c
      return
      end


      subroutine phase_space_npo(x,shat,iU,iS,iB,iA,p,pbar,xjac,xjacB)
c     iU is the unresolved parton associated with the soft singularity
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      double precision x(4),shat
      double precision p(0:3,-2:3),pbar(0:3,-2:2)
      double precision xjac,xjacB,xjacCS
      integer i,j,iU,iS,iB,iA
c
c     initialise
      p=0d0
      pbar=0d0
      xjac=1d0
      p(0,0)=sqrt(shat)
c
c     call Born phase space
      call phase_space_n(x(4),shat,pbar,xjacB)
c
c     labels are assigned so to have (iU1 iS1 iB1)
c     as an n+1 -> n mapping
c
c     call radiation phase space
      call phase_space_CS(x(1),iU,iS,iB,iA,p,pbar,npartNLO,xjacCS)
c
c     total jacobian
      xjac=xjacB*xjacCS
c
      return
      end


      subroutine phase_space_npt(x,shat,iU1,iS1,iB1,iA1,iA2,p,pbar,ptilde,xjac,xjacB,iU2,iS2,iB2)
c     iU1 and iU2 are the unresolved partons associated with the soft singularity
c     iU2, iS2, and iB2 are outputs
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      double precision x(7),shat
      double precision p(0:3,-2:4),pbar(0:3,-2:3),ptilde(0:3,-2:2)
      double precision xjac,xjacB,xjacCS1,xjacCS2
      integer i,j,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2
      double precision dot
c
c     initialise
      p=0d0
      pbar=0d0
      ptilde=0d0
      xjac=1d0
      p(0,0)=sqrt(shat)
      pbar(0,0)=sqrt(shat)
c
c     call Born phase space
      call phase_space_n(x(7),shat,ptilde,xjacB)
c
c     map according to (iU1 iS1 iB1 , iS1 iA1 iB1) = (iU1 iS1 iA1 iB1)
      iU2=imap(iS1,iU1,iS1,0,0,npartNNLO)
      iS2=imap(iA1,iU1,iS1,0,0,npartNNLO)
      iB2=imap(iB1,iU1,iS1,0,0,npartNNLO)
c
c     call radiation phase space from Born to real
      call phase_space_CS(x(4),iU2,iS2,iB2,iA2,pbar,ptilde,npartNLO,xjacCS2)
c
c     call radiation phase space from real to double real
      call phase_space_CS(x(1),iU1,iS1,iB1,iA1,p,pbar,npartNNLO,xjacCS1)
c
c     total jacobian
      xjac=xjacB*xjacCS1*xjacCS2
c
      return
      end


      subroutine check_phsp_consistency(x,npart,xs,xsb,iU1,iS1,iB1,iA1,ierr)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      include 'mxdim.inc'
      integer i,iU1,iS1,iB1,iA1,iU2,iS2,iB2,iA2,npart,ierr
      double precision x(mxdim)
      double precision xs(-2:maxdim,-2:maxdim),xsb(-2:maxdim,-2:maxdim)
      double precision sab,sac,sad,sbc,sbd,scd,sabc,sabcd
      double precision sb_bc,sb_bd,sb_cd,sb_bcd
      double precision y,z,yp,zp,y2,z2
c
c     initialise
      ierr=0
c
      if(npart.eq.npartNLO)then     
c     this routine assumes mapping labels (iU1 iS1 iB1)
c     in passing from n+1 to n
         sab=xs(iU1,iS1)
         sac=xs(iU1,iB1)
         sbc=xs(iS1,iB1)
         sabc=sab+sac+sbc
         y=sab/sabc
         z=sac/(sabc-sab)
c
c     check
         if(sabc.le.0d0.or.sabc-sab.le.0d0)then
            write(77,*)'Inaccuracy in check_phsp_consistency NLO',
     &      sabc,sabc-sab
         endif
c
c     actual test of Catani-Seymour mapping and parametrisation
         if(abs(z-x(1)).gt.tiny0)then
            write(77,*)'Inaccurate z in check_phsp_consistency NLO',abs(z-x(1))
            goto 999
         endif
         if(abs(y-x(2)).gt.tiny0)then
            write(77,*)'Inaccurate y in check_phsp_consistency NLO',abs(y-x(2))
            goto 999
         endif
c
      elseif(npart.eq.npartNNLO)then     
c     this routine assumes mapping labels (iU1 iS1 iB1, iS1 iA1 iB1)
c     in passing from n+2 to n+1 to n
         sab=xs(iU1,iS1)
         sac=xs(iU1,iB1)
         sad=xs(iU1,iA1)
         sbc=xs(iS1,iB1)
         sbd=xs(iS1,iA1)
         scd=xs(iB1,iA1)
         sabc=sab+sac+sbc
         sabcd=sabc+sad+sbd+scd
         iU2=imap(iS1,iU1,iS1,0,0,npartNNLO)
         iS2=imap(iA1,iU1,iS1,0,0,npartNNLO)
         iB2=imap(iB1,iU1,iS1,0,0,npartNNLO)
         sb_bc=xsb(iU2,iB2)
         sb_bd=xsb(iU2,iS2)
         sb_cd=xsb(iS2,iB2)
         sb_bcd=sb_bc+sb_bd+sb_cd
c     y = sbar_bd(abc) / sbar_bdc(abc)
c     z = sbar_bc(abc) / (sbar_bc(abc) + sbar_cd(abc))
         y=(sad+sbd-sab*scd/(sac+sbc))/sabcd
         z=(sac+sbc)/(sac+sbc+scd)
         yp=sab/sabc
         zp=sac/(sabc-sab)
         y2=sb_bd/sb_bcd
         z2=sb_bc/(sb_bc+sb_cd)
c
c     check
         if(sabcd.le.0d0.or.sabcd-sabc.le.0d0.or.
     &   sabc-sab.le.0d0.or.sabc.le.0d0)then
            write(77,*)'Inaccuracy in check_phsp_consistency NNLO',
     &      sabcd,sabcd-sabc,sabc-sab,sabc
         endif
c
c     actual test of Catani-Seymour mapping and parametrisation
         if(abs(zp-x(1)).gt.tiny0)then
            write(77,*)'Inaccurate zp  in check_phsp_consistency NNLO',abs(zp-x(1))
            goto 999
         endif
         if(abs(yp-x(2)).gt.tiny0)then
            write(77,*)'Inaccurate yp  in check_phsp_consistency NNLO',abs(yp-x(2))
            goto 999
         endif
         if(abs(z-x(4)).gt.tiny0)then
            write(77,*)'Inaccurate z   in check_phsp_consistency NNLO',abs(z-x(4))
            goto 999
         endif
         if(abs(y-x(5)).gt.tiny0)then
            write(77,*)'Inaccurate y   in check_phsp_consistency NNLO',abs(y-x(5))
            goto 999
         endif
         if(abs(y-y2).gt.tiny0)then
            write(77,*)'Inaccurate y 2 in check_phsp_consistency NNLO',abs(y-y2)
            goto 999
         endif
         if(abs(z-z2).gt.tiny0)then
            write(77,*)'Inaccurate z 2 in check_phsp_consistency NNLO',abs(z-z2)
            goto 999
         endif
      else
         write(*,*)'Unknown npart in check_phsp_consistency',
     &   npart,npartNLO,npartNNLO
      endif
c
      return
 999  ierr=1
      return
      end





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
