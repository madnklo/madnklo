      subroutine phase_space_CS(xx,iU,iS,iB,iA,p,pbar,npart,leg_PDGs,maptype,xjac)
c     Build n+1 momenta p from n momenta pbar 
c     iU is the unresolved parton associated
c     with the soft singularity
c     npart = n+1
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision xx(3)
      integer i,j,iU,iS,iB,iA,npart,ierr
      double precision p(0:3,npart),pbar(0:3,npart-1),xjac
      double precision zCS,yCS,phCS
      double precision cosphi,sinphi,cosphiA,sinphiA
      double precision cosphiU,sinphiU,costhUB,sinthUB
      double precision Eu,pmod,pAtmod,GG
      double precision nB(0:3),pbB(0:3),pbBsave(0:3)
      double precision pbS(0:3),pA(0:3),pboost(0:3)
      double precision pAr(0:3),pUr(0:3)
      double precision dot
      double precision ran2
      integer leg_PDGs(npart)
      character*1 maptype
      integer mapped_labels(npart),mapped_flavours(npart)
c     TODO: change name to isunderlyingqcdparton()
      logical isLOQCDparton(npart-1)
      integer idum
      common/rand/idum
      double precision yplus,z_minus,z_plus
      double precision Q2,sdip,mb2,mc2,lam1,lam2,lambda
      double precision vel
c
c     initialise
      xjac=1d0
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
      call get_mapped_labels(maptype,iU,iS,iB,npart,leg_pdgs,mapped_labels,
     $     mapped_flavours,isLOQCDparton)
      pbB(:)=pbar(:,mapped_labels(iB))
      pbS(:)=pbar(:,mapped_labels(iS))
      pA(:) =pbar(:,mapped_labels(iA))
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
c     for massive emitter and massive recoiler


      mb2 = dot(pbS(:),pbS(:))
      mc2 = dot(pbB(:),pbB(:))
      
      Q2 = dot(pboost(:),pboost(:))
      sdip = Q2-mb2-mc2
      yplus = 1d0 + 2d0*mc2/sdip-(2d0*dsqrt(mc2*(sdip+mb2+mc2)))/sdip
      yCS = xx(1)*yplus
      xjac = xjac*yplus
      vel = dsqrt((2d0*mc2+sdip*(1d0-yCS))**2-4d0*mc2*Q2)/sdip/(1d0-yCS)
      z_minus = sdip*yCS/2d0/(sdip*yCS+mb2)*(1d0-vel)
      z_plus  = sdip*yCS/2d0/(sdip*yCS+mb2)*(1d0+vel)
      zCS = z_minus + (z_plus - z_minus)*xx(2)
      xjac = xjac*(z_plus - z_minus)
      phCS=2d0*pi*xx(3)
      xjac = xjac*2d0*pi
      
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
      lam1=lambda(Q2,mb2,mc2)

      
c
c     Eu is the energy of p(iU)

      Eu=sdip*(1d0-(1d0-zCS)*(1d0-yCS))/2d0/dsqrt(Q2)
      costhUB=(dsqrt(mc2 + pmod**2) - 
     -    ((Eu*(sdip*(2*mc2 + sdip)*vel*(-1 + yCS) + 
     -            dsqrt(lam1)*(2*mc2 + sdip - sdip*yCS)))
     -         /(dsqrt(Q2)*sdip*(-1 + yCS)) + 
     -       dsqrt(lam1)*zCS)/(2d0*Eu*vel))/pmod
      sinthUB= dsqrt(1d0-costhUB**2)

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
c      p(:,iB)=(1d0-yCS)*pbBsave(:)
c      p(:,iS)=yCS*pbBsave(:)+pbS(:)-p(:,iU)

c      lam2=(2d0*mc2+sdip*(1d0-yCS))**2-4d0*mc2*Q2
      lam2 = lambda(Q2,mb2+sdip*yCS,mc2)
      p(:,iB)=dsqrt(lam2/lam1)*
     $     (pbBsave(:)-(Q2+mc2-mb2)/2d0/Q2*pboost(:))+
     $     (sdip*(1d0-yCS)+2d0*mc2)/2d0/Q2*pboost(:)
      
      p(:,iS)=pboost(:)-p(:,iU)-p(:,iB)

      
      do j=1,npart
         if(j.eq.iU.or.j.eq.iB.or.j.eq.iS)cycle
         p(:,j)=pbar(:,mapped_labels(j))
      enddo
c
c     consistency check
      do i=1,npart
         if(p(0,i).lt.0d0)then
            write(77,*)'Inaccuracy 1 in phase_space_CS',i,p(0,i),npart
            goto 999
         endif
      enddo
c
c     construct xjac
c      GG=1d0/16d0/pi**3
c      xjac=GG*(4d0*pmod**2)*pi*(1-yCS)
      xjac = xjac * (sdip**2*(1d0-yCS))/(16d0*dsqrt(lam1)*Pi**2)
c
      return
 999  xjac=0d0
      return
      end


      subroutine phase_space_CS_inv(iU,iS,iB,p,pbar,npart,leg_PDGs,
     & maptype,xjac)
c     Build n momenta pbar from n+1 momenta p
c     iU is the unresolved parton associated
c     with the soft singularity
c
c     TODO: add mappings for initial state
c
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      integer i,j,iU,iS,iB,npart
      double precision p(0:3,npart),pbar(0:3,npart-1),Q(0:3)
      integer leg_PDGs(npart)
      double precision yCS,xjac,GG,Qsq,dot
      character*1 maptype
      integer mapped_labels(npart),mapped_flavours(npart)
c     TODO: change name to isunderlyingqcdparton()
      logical isLOQCDparton(npart-1)
      double precision lam1,lam2,mS2,mB2,miUiS
      double precision lambda

c
c     initialise
      xjac=0d0
      mapped_labels = 0
c
c     auxiliary quantities
      mS2=dot(p(0,iS),p(0,iS))
      mB2=dot(p(0,iB),p(0,iB))
      Q(:) = p(:,iU) + p(:,iS) + p(:,iB)
      Qsq = dot(Q(:),Q(:))
c      Qsq=2d0*(dot(p(0,iU),p(0,iS))+dot(p(0,iU),p(0,iB))+
c     &         dot(p(0,iS),p(0,iB))) + mS2 + mB2
      yCS=2d0*dot(p(0,iU),p(0,iS))/(Qsq-mS2-mB2)

      lam1 = lambda(Qsq,mB2,mS2)
c     miUiS is the invariant mass of (p(:,iU)+p(:,iS))
      miUiS = dot(p(0,iU),p(0,iU))+dot(p(0,iS),p(0,iS))
     $     +2d0*dot(p(0,iU),p(0,is))
      lam2 = lambda(Qsq,miUiS,mB2)
c
c     construct pbar from p
      call get_mapped_labels(maptype,iU,iS,iB,npart,leg_pdgs,mapped_labels,
     $           mapped_flavours,isLOQCDparton)
c


      pbar(:,mapped_labels(iB))=
     $     dsqrt(lam1/lam2)*(p(:,iB)-dot(Q(0),p(:,iB))*Q(:)/Qsq)+
     $     (Qsq+mB2-mS2)/2d0/Qsq*Q(:)
      pbar(:,mapped_labels(iS))=Q(:)-pbar(:,mapped_labels(iB))
      do j=1,npart
         if(j.eq.iU.or.j.eq.iB.or.j.eq.iS)cycle
         pbar(:,mapped_labels(j))=p(:,j)
      enddo
c
c     consistency check
      do i=1,npart-1
         if(pbar(0,i).lt.0d0)then
            write(77,*)'Inaccuracy 1 in phase_space_CS_inv',i,pbar(0,i),npart-1
            goto 999
         endif
      enddo
c
c     construct xjac
c      GG=1d0/16d0/pi**3
c      xjac=GG*Qsq*pi*(1-yCS)

      xjac = ((Qsq-mB2-mS2)**2*(1d0-yCS))/(16d0*dsqrt(lam1)*Pi**2)
c
      return
 999  xjac=0d0
      end


      
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
