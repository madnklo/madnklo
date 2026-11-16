      subroutine phase_space_CS(xx,iU,iS,iB,iA,p,pbar,npart,leg_PDGs,xjac,mapped_labels)
c     Build n+1 momenta p from n momenta pbar 
c     iU is the unresolved parton associated
c     with the soft singularity
c     npart = n+1
      implicit none
      include 'nexternal.inc'
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
      integer idum
      common/rand/idum
      double precision yplus,z_minus,z_plus
      double precision Q2,sdip,mb2,mc2,lam1,lam2,lambda
      double precision vel
      integer mapped_labels(npart)
      double precision pmass(nexternal)
      double precision CSpow(2)
      data CSpow /1d0,1d0/
      common/cCSpow/CSpow
      include 'pmass.inc'
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
c     recoiler in the initial state: FFI or FII mappings (below)
      if(iB.le.2)then
         call phase_space_CS_initial(xx,iU,iS,iB,iA,p,pbar,npart,xjac,
     &        mapped_labels)
         return
      endif
c
c     iA (reference four-vector for definition of the azimuth)
c     must be != iB, iU, iS;
c      call get_mapped_labels_reduced(npart,iU,iS,mapped_labels)
c      call get_collinear_mapped_labels(npart,iU,iS,leg_pdgs,idum1,mapped_labels,idum2,idum3)
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

c     rotate pA to the frame where p(iB) and pbar(iB) are along (0,0,-1)
      pmod=sqrt(pbB(1)**2+pbB(2)**2+pbB(3)**2)
      if(pmod.le.0d0)then
         write(77,*)'Wrong pmod in phase_space_CS',pmod
         goto 999
      endif
      nB=pbB/pmod
      call rotate_to_z(nB(1),pA(1),pAr(1))
c
      mb2 = pmass(iS)**2
      mc2 = pmass(iB)**2
      if(mb2 .ne. 0d0 .or. mc2 .ne. 0d0) then
c     Catani-Seymour parametrisation 
c     for massive emitter and massive recoiler
         Q2 = dot(pboost(:),pboost(:))
         sdip = Q2-mb2-mc2
         yplus = 1d0 + 2d0*mc2/sdip-(2d0*dsqrt(abs(mc2*(sdip+mb2+mc2))))/sdip
         yCS = xx(2)*yplus
         xjac = xjac*yplus
         if((2d0*mc2+sdip*(1d0-yCS))**2-4d0*mc2*Q2.lt.0d0)then
            write(77,*)'Negative sqrt in phase space'
            stop
         endif
         vel = dsqrt((2d0*mc2+sdip*(1d0-yCS))**2-4d0*mc2*Q2)/sdip/(1d0-yCS)
         z_minus = sdip*yCS/2d0/(sdip*yCS+mb2)*(1d0-vel)
         z_plus  = sdip*yCS/2d0/(sdip*yCS+mb2)*(1d0+vel)
         zCS = z_minus + (z_plus - z_minus)*xx(1)
         xjac = xjac*(z_plus - z_minus)
         phCS=2d0*pi*xx(3)
         xjac = xjac*2d0*pi
         lam1=lambda(Q2,mb2,mc2)
         lam2=lambda(Q2,mb2+sdip*yCS,mc2)
c        Eu is the energy of p(iU)
         Eu=sdip*(1d0-(1d0-zCS)*(1d0-yCS))/2d0/dsqrt(Q2)
         costhUB=(dsqrt(mc2 + pmod**2) - 
     -        ((Eu*(sdip*(2*mc2 + sdip)*vel*(-1d0 + yCS) + 
     -        dsqrt(lam1)*(2*mc2 + sdip - sdip*yCS)))
     -        /(dsqrt(Q2)*sdip*(-1d0 + yCS)) + 
     -        dsqrt(lam1)*zCS)/(2d0*Eu*vel))/pmod
         sinthUB= dsqrt(abs(1d0-costhUB**2))
      else
c     Catani-Seymour parametrisation for massless case
         zCS =xx(1)**CSpow(1)
         xjac=xjac* CSpow(1)*xx(1)**(CSpow(1)-1d0)
         yCS =xx(2)**CSpow(2)
         xjac=xjac* CSpow(2)*xx(2)**(CSpow(2)-1d0)
         if(1d0-yCS.eq.1d0 .or. 1d0-zCS.eq.1d0) return
         phCS=2d0*pi*xx(3)
         costhUB=1d0-2d0*zCS/(1d0-(1d0-yCS)*(1d0-zCS))
         sinthUB=2d0*sqrt(yCS*(1d0-zCS)*zCS)/(1d0-(1d0-yCS)*(1d0-zCS))
c
c     Eu is the energy of p(iU)
         Eu=pmod*(1d0-(1d0-yCS)*(1d0-zCS))
      endif
c
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

      if(mb2 .ne.0d0 .or. mc2 .ne. 0d0) then
c     CS massive mapping
         p(:,iB)=dsqrt(lam2/lam1)*
     $        (pbBsave(:)-(Q2+mc2-mb2)/2d0/Q2*pboost(:))+
     $        (sdip*(1d0-yCS)+2d0*mc2)/2d0/Q2*pboost(:)      
         p(:,iS)=pboost(:)-p(:,iU)-p(:,iB)
c     construct xjac
         xjac = xjac * (sdip**2*(1d0-yCS))/(4d0*dsqrt(lam1)*(2d0*Pi)**3)
      else
c     CS massless mapping         
         p(:,iB)=(1d0-yCS)*pbBsave(:)
         p(:,iS)=yCS*pbBsave(:)+pbS(:)-p(:,iU)
c     construct xjac
         GG=1d0/16d0/pi**3
         xjac=xjac*GG*(4d0*pmod**2)*pi*(1-yCS)
      endif
      
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
c
      return
 999  xjac=0d0
      return
      end


      subroutine phase_space_CS_inv(iU,iS,iB,p,pbar,npart,leg_PDGs,xjac,mapped_labels)
c     Build n momenta pbar from n+1 momenta p
c     iU is the unresolved parton associated
c     with the soft singularity
c     (recoiler iB in the initial state: phase_space_CS_initial_inv)
c
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      integer i,j,iU,iS,iB,npart
      double precision p(0:3,npart),pbar(0:3,npart-1),Q(0:3)
      integer leg_PDGs(npart)
      double precision yCS,zCS,xjac,GG,Qsq,dot
      double precision lam1,lam2,mS2,mB2,miUiS2
      double precision lambda,vel,sdip,z_minus,z_plus
      integer mapped_labels(npart)
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      xjac=0d0
c
c     recoiler in the initial state: FFI or FII mappings
      if(iB.le.2)then
         call phase_space_CS_initial_inv(iU,iS,iB,p,pbar,npart,xjac,
     &        mapped_labels)
         return
      endif
c     
c     auxiliary quantities
      mS2=pmass(iS)**2
      mB2=pmass(iB)**2

      Q(:) = p(:,iU) + p(:,iS) + p(:,iB)
      Qsq = dot(Q(:),Q(:))
c
c     construct pbar from p
c     since this routine is to be called with iU = isec = gluon, then
c     mapped labels should be the ones assigned at the beginning of the
c     run for this sector, and stored in c_mapped_labels
c      call get_mapped_labels_reduced(npart,iU,iS,mapped_labels)
c      call get_collinear_mapped_labels(npart,iU,iS,leg_pdgs,idum1,mapped_labels,idum2,idum3)
c
      if(mB2.ne.0d0 .or. mS2.ne.0d0) then
c     CS massive case
         sdip=Qsq-mS2-mB2
         yCS=2d0*dot(p(0,iU),p(0,iS))/sdip
         zCS=dot(p(0,iU),p(0,iB))/(dot(p(0,iU),p(0,iB))+dot(p(0,iS),p(0,iB)))
         vel = dsqrt((2d0*mB2+sdip*(1d0-yCS))**2-4d0*mB2*Qsq)/sdip/(1d0-yCS)
         z_minus = sdip*yCS/2d0/(sdip*yCS+mS2)*(1d0-vel)
         z_plus  = sdip*yCS/2d0/(sdip*yCS+mS2)*(1d0+vel)
         if(zCS.gt.z_plus .or. zCS.lt.z_minus)return         
c     miUiS2 is the invariant mass squared of (p(:,iU)+p(:,iS))
         miUiS2 = mS2+2d0*dot(p(0,iU),p(0,iS))
         lam1 = lambda(Qsq,mS2,mB2)
         lam2 = lambda(Qsq,miUiS2,mB2)
         if(lam1.le.0d0.or.lam2.le.0d0)then
            write(*,*)'negative lambdas in phase_space_CS_inv',lam1,lam2
            stop
         endif

         pbar(:,mapped_labels(iB))=
     $        dsqrt(lam1/lam2)*(p(:,iB)-dot(Q(0),p(0,iB))*Q(:)/Qsq)+
     $        (Qsq+mB2-mS2)/2d0/Qsq*Q(:)
         pbar(:,mapped_labels(iS))=Q(:)-pbar(:,mapped_labels(iB))
c     construct xjac
         xjac = ((Qsq-mB2-mS2)**2*(1d0-yCS))/(16d0*dsqrt(lam1)*Pi**2)
      else
c     CS massless case
         yCS=2d0*dot(p(0,iU),p(0,iS))/Qsq
         pbar(:,mapped_labels(iB))=p(:,iB)/(1-yCS)
         pbar(:,mapped_labels(iS))=p(:,iU)+p(:,iS)-p(:,iB)*yCS/(1-yCS)

c     construct xjac
         GG=1d0/16d0/pi**3
         xjac=GG*Qsq*pi*(1-yCS)
      endif
         
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


c     subroutine get_mapped_labels_reduced(npart,a,b,mapped_labels)
c     implicit none
c     integer npart,a,b,mapped_labels(npart)
c     integer i,i1,i2,ii
c
c     preliminary checks
c      if(a.lt.3.or.b.lt.3)then
c         write(*,*)'wrong a, b in reduced_mapped_labels',a,b
c         stop
c      endif
cc
c      i1=min(a,b)
c      i2=max(a,b)
c      ii=0
c      do i=1,npart
c         if(i.eq.i1)cycle
c         ii=ii+1
c         mapped_labels(i)=ii
c      enddo
c      mapped_labels(i1)=mapped_labels(i2)
c
c      return
c      end


c=======================================================================
c     Phase-space mappings with initial-state partons
c     (2209.09123, App. B; notation (abc) = (iU,iS,iB) = (emitted,
c     emitter, recoiler), all massless):
c
c     FFI: iS final, iB initial
c          kbar_S = k_U + k_S - (1-x) k_B ,   kbar_B = x k_B
c          x = (s_UB + s_SB - s_US)/(s_UB + s_SB)
c          z = s_UB/(s_UB + s_SB)
c          radiative measure  sbar/(16 pi^2) dx dz dphi/(2 pi),
c          sbar = 2 kbar_S.k_B = s_UB + s_SB
c
c     FII: iS initial, iB initial
c          kbar_S = x k_S ,   kbar_B = k_B ,
c          final states (other than iU) Lorentz-transformed from
c          K = k_S + k_B - k_U  to  Kbar = x k_S + k_B
c          x = (s_SB - s_US - s_UB)/s_SB
c          v = s_US/(s_US + s_UB)
c          radiative measure  sbar (1-x)/(16 pi^2) dx dv dphi/(2 pi),
c          sbar = s_SB
c
c     The initial-state parton rescaled by the mapping (iB for FFI, iS
c     for FII) has momentum fraction xbk_npo = xbk_bar/x in the real
c     kinematics, so that, at fixed Born, x ranges in [xbk_bar,1]. The
c     hadronic weight is then (Born weight of gen_mom) x (radiative
c     measure), with the PDF of the rescaled beam evaluated at the real
c     momentum fraction, and no further factor 1/x.
c
c     The momentum fractions of the real kinematics (xbk_npo) and of
c     the mapped Born kinematics (xbk_bar) are passed in the common
c     block c_xbk_mapped:
c       phase_space_CS     sets xbk_npo from the Born fractions xbk
c                          (run.inc) of gen_mom;
c       phase_space_CS_inv sets xbk_bar from xbk_npo.
c=======================================================================
      subroutine phase_space_CS_initial(xx,iU,iS,iB,iA,p,pbar,npart,
     &     xjac,mapped_labels)
c     Build n+1 momenta p from n momenta pbar (lab frame), for a
c     recoiler iB in the initial state: FFI (iS final) or FII (iS
c     initial). Integration variables:
c       xx(1) = z (FFI; z/z+ with a massive emitter) or v (FII),
c       xx(2) = (1-x)/(1-xbk_bar), i.e. xx(2) -> 0 is the x -> 1 limit,
c       xx(3) = phi/(2 pi), azimuth around the emitter-recoiler axis,
c               measured from the reference momentum pbar(iA).
c     The momentum fractions xbk_npo(a) of the initial-state partons
c     a = 1,2 (labels of the real emission) are set by phase_space_npo
c     to the Born ones; here that of the rescaled parton is divided by x.
c     The Born labels of the two initial-state partons may be exchanged
c     (mapped_labels(1) = 2): see get_mapped_labels_initial (imap.f).
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'run.inc'
      integer npart,iU,iS,iB,iA,j,ibeam,ierr
      integer mapped_labels(npart)
      double precision xx(3),p(0:3,npart),pbar(0:3,npart-1),xjac
      double precision xCS,zCS,vCS,phCS,xmin,sbar,kt2
      double precision pbS(0:3),pbB(0:3),pA(0:3),kB(0:3),kS(0:3)
      double precision kt(0:3),K(0:3),Kbar(0:3)
      double precision mS2,zp,P2,bet,PP(0:3),nS(0:3)
      double precision dot
      double precision xbk_npo(2),xbk_bar(2)
      common/c_xbk_mapped/xbk_npo,xbk_bar
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      xjac=0d0
      p=0d0
c
c     massless emitted parton and initial-state partons; the final-state
c     emitter of FFI can be massive (MadNkLO paper with massive
c     particles, App. A.3: FmI mapping)
      if(pmass(iU).ne.0d0.or.(iS.le.2.and.pmass(iS).ne.0d0))then
         write(*,*)'phase_space_CS_initial: massive emitted parton ',
     &        'or initial-state emitter not implemented',iU,iS,iB
         stop
      endif
c
c     Born momenta of emitter, recoiler and reference
      pbS(:)=pbar(:,mapped_labels(iS))
      pbB(:)=pbar(:,mapped_labels(iB))
      pA(:) =pbar(:,mapped_labels(iA))
c
c     initial-state parton rescaled by the mapping, and the minimum x
      if(iS.le.2)then
         ibeam=iS
      else
         ibeam=iB
      endif
      if(mapped_labels(ibeam).gt.2)then
         write(*,*)'phase_space_CS_initial: initial-state label ',ibeam,
     &        ' mapped to ',mapped_labels(ibeam)
         stop
      endif
      xmin=xbk_npo(ibeam)
      if(xmin.le.0d0.or.xmin.ge.1d0)goto 999
c
c     x from xx(2): flat in [xmin,1]
      xCS=1d0-(1d0-xmin)*xx(2)
      xjac=1d0-xmin
      if(xCS.ge.1d0.or.xCS.le.0d0)goto 999
      phCS=2d0*pi*xx(3)
c
      if(iS.gt.2)then
c
c     FFI: emitter iS final, recoiler iB initial
         mS2=pmass(iS)**2
         kB(:)=pbB(:)/xCS
         sbar=2d0*dot(pbS,kB)
         if(mS2.eq.0d0)then
c        massless emitter: z = s_Uk/sbar in [0,1]
         zCS=xx(1)
         if(zCS.le.0d0.or.zCS.ge.1d0)goto 999
         kt2=zCS*(1d0-zCS)*(1d0-xCS)*sbar
         call get_transverse_vector(pbS,kB,pA,kt2,phCS,kt,ierr)
         if(ierr.ne.0)goto 999
         p(:,iU)=zCS*pbS(:)+(1d0-zCS)*(1d0-xCS)*kB(:)+kt(:)
         p(:,iS)=(1d0-zCS)*pbS(:)+zCS*(1d0-xCS)*kB(:)-kt(:)
         else
c        massive emitter (App. A.3): z = s_Uk/sbar in [0,z+], with
c        z+ = sbar(1-x)/(sbar(1-x)+m^2); P = k_U + k_S = pbS + (1-x) kB,
c        light-like n = P - P^2/sbar kB, and
c        k_U = z n + beta kB + kt, beta = (1-x) - z P^2/sbar,
c        so that s_Uk = z sbar, s_US = (1-x) sbar, k_S^2 = m^2
         zp=sbar*(1d0-xCS)/(sbar*(1d0-xCS)+mS2)
         zCS=zp*xx(1)
         xjac=xjac*zp
         if(zCS.le.0d0.or.zCS.ge.zp)goto 999
         PP(:)=pbS(:)+(1d0-xCS)*kB(:)
         P2=mS2+(1d0-xCS)*sbar
         nS(:)=PP(:)-P2/sbar*kB(:)
         bet=(1d0-xCS)-zCS*P2/sbar
         kt2=zCS*bet*sbar
         call get_transverse_vector(nS,kB,pA,kt2,phCS,kt,ierr)
         if(ierr.ne.0)goto 999
         p(:,iU)=zCS*nS(:)+bet*kB(:)+kt(:)
         p(:,iS)=PP(:)-p(:,iU)
         endif
         p(:,iB)=kB(:)
         do j=1,npart
            if(j.eq.iU.or.j.eq.iS.or.j.eq.iB)cycle
            p(:,j)=pbar(:,mapped_labels(j))
         enddo
c        radiative measure
         xjac=xjac*sbar/(16d0*pi**2)
      else
c
c     FII: emitter iS initial, recoiler iB initial
         vCS=xx(1)
         if(vCS.le.0d0.or.vCS.ge.1d0)goto 999
         kS(:)=pbS(:)/xCS
         sbar=2d0*dot(kS,pbB)
         kt2=(1d0-xCS)**2*vCS*(1d0-vCS)*sbar
         call get_transverse_vector(kS,pbB,pA,kt2,phCS,kt,ierr)
         if(ierr.ne.0)goto 999
         p(:,iU)=(1d0-xCS)*(1d0-vCS)*kS(:)+(1d0-xCS)*vCS*pbB(:)+kt(:)
         p(:,iS)=kS(:)
         p(:,iB)=pbB(:)
c        final-state momenta: Lorentz transformation Kbar -> K
         Kbar(:)=pbS(:)+pbB(:)
         K(:)=kS(:)+pbB(:)-p(:,iU)
         do j=1,npart
            if(j.eq.iU.or.j.eq.iS.or.j.eq.iB)cycle
            call lorentz_KbartoK(K,Kbar,pbar(0,mapped_labels(j)),
     &           p(0,j))
         enddo
c        radiative measure
         xjac=xjac*sbar*(1d0-xCS)/(16d0*pi**2)
      endif
c
c     momentum fraction of the rescaled beam in the real kinematics
      xbk_npo(ibeam)=xmin/xCS
c
c     consistency check
      do j=1,npart
         if(p(0,j).lt.0d0)then
            write(77,*)'Inaccuracy 1 in phase_space_CS_initial',j,
     &           p(0,j),npart
            goto 999
         endif
      enddo
c
      return
 999  xjac=0d0
      return
      end


      subroutine phase_space_CS_initial_inv(iU,iS,iB,p,pbar,npart,xjac,
     &     mapped_labels)
c     Build n momenta pbar from n+1 momenta p, for a recoiler iB in
c     the initial state: FFI (iS final) or FII (iS initial). Returns
c     in xjac the radiative measure factor (sbar/(16 pi^2) for FFI,
c     sbar (1-x)/(16 pi^2) for FII), and sets xbk_bar, the momentum
c     fractions of the mapped Born kinematics.
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      integer npart,iU,iS,iB,j,ibeam
      integer mapped_labels(npart)
      double precision p(0:3,npart),pbar(0:3,npart-1),xjac
      double precision sUS,sUB,sSB,xCS,K(0:3),Kbar(0:3)
      double precision dot
      double precision xbk_npo(2),xbk_bar(2)
      common/c_xbk_mapped/xbk_npo,xbk_bar
c
c     initialise
      xjac=0d0
      xbk_bar(1)=xbk_npo(1)
      xbk_bar(2)=xbk_npo(2)
c
      sUS=2d0*dot(p(0,iU),p(0,iS))
      sUB=2d0*dot(p(0,iU),p(0,iB))
      sSB=2d0*dot(p(0,iS),p(0,iB))
c
      if(iS.gt.2)then
c
c     FFI: emitter iS final, recoiler iB initial
         ibeam=iB
         if(sUB+sSB.le.0d0)goto 999
         xCS=(sUB+sSB-sUS)/(sUB+sSB)
         if(xCS.le.0d0.or.xCS.gt.1d0)goto 999
         pbar(:,mapped_labels(iS))=p(:,iU)+p(:,iS)-(1d0-xCS)*p(:,iB)
         pbar(:,mapped_labels(iB))=xCS*p(:,iB)
         do j=1,npart
            if(j.eq.iU.or.j.eq.iS.or.j.eq.iB)cycle
            pbar(:,mapped_labels(j))=p(:,j)
         enddo
         xjac=(sUB+sSB)/(16d0*pi**2)
      else
c
c     FII: emitter iS initial, recoiler iB initial
         ibeam=iS
         if(sSB.le.0d0)goto 999
         xCS=(sSB-sUS-sUB)/sSB
         if(xCS.le.0d0.or.xCS.gt.1d0)goto 999
         pbar(:,mapped_labels(iS))=xCS*p(:,iS)
         pbar(:,mapped_labels(iB))=p(:,iB)
c        final-state momenta: Lorentz transformation K -> Kbar
         K(:)=p(:,iS)+p(:,iB)-p(:,iU)
         Kbar(:)=xCS*p(:,iS)+p(:,iB)
         do j=1,npart
            if(j.eq.iU.or.j.eq.iS.or.j.eq.iB)cycle
            call lorentz_KbartoK(Kbar,K,p(0,j),pbar(0,mapped_labels(j)))
         enddo
         xjac=sSB*(1d0-xCS)/(16d0*pi**2)
      endif
c
c     momentum fraction of the rescaled beam in the mapped Born
      xbk_bar(ibeam)=xbk_npo(ibeam)*xCS
c
      return
 999  xjac=0d0
      return
      end


      subroutine get_transverse_vector(pb,pc,pref,kt2,phi,kt,ierr)
c     Space-like vector kt orthogonal to the massless momenta pb and
c     pc, with kt.kt = -kt2 and azimuth phi around the pb-pc axis,
c     measured from the component of pref orthogonal to pb and pc
c     (e1); e2 completes the basis (e2^mu = eps^{mu nu rho sigma}
c     pb_nu pc_rho e1_sigma, normalised).
      implicit none
      integer ierr,mu
      double precision pb(0:3),pc(0:3),pref(0:3),kt2,phi,kt(0:3)
      double precision e1(0:3),e2(0:3),e2up(0:3),pbpc,n1,n2
      double precision dot
c
      ierr=0
      pbpc=dot(pb,pc)
      if(pbpc.le.0d0.or.kt2.lt.0d0)goto 999
c
c     e1: pref minus its components along pb and pc
      e1(:)=pref(:)-dot(pref,pc)/pbpc*pb(:)-dot(pref,pb)/pbpc*pc(:)
      n1=-dot(e1,e1)
      if(n1.le.0d0)goto 999
      e1(:)=e1(:)/sqrt(n1)
c
c     e2^mu = eps^{mu nu rho sigma} pb_nu pc_rho e1_sigma (eps^{0123}=1)
      call levi_civita_contract(pb,pc,e1,e2up)
      do mu=0,3
         e2(mu)=e2up(mu)
      enddo
      n2=-dot(e2,e2)
      if(n2.le.0d0)goto 999
      e2(:)=e2(:)/sqrt(n2)
c
      kt(:)=sqrt(kt2)*(cos(phi)*e1(:)+sin(phi)*e2(:))
      return
 999  ierr=1
      return
      end


      subroutine levi_civita_contract(a,b,c,r)
c     r^mu = eps^{mu nu rho sigma} a_nu b_rho c_sigma, with eps^{0123}=1
c     and lower indices a_nu = g_{nu nu} a^nu (metric +,-,-,-)
      implicit none
      double precision a(0:3),b(0:3),c(0:3),r(0:3)
      double precision al(0:3),bl(0:3),cl(0:3)
      integer i
      double precision g(0:3)
      data g/1d0,-1d0,-1d0,-1d0/
      do i=0,3
         al(i)=g(i)*a(i)
         bl(i)=g(i)*b(i)
         cl(i)=g(i)*c(i)
      enddo
      r(0)= al(1)*(bl(2)*cl(3)-bl(3)*cl(2))
     &     -al(2)*(bl(1)*cl(3)-bl(3)*cl(1))
     &     +al(3)*(bl(1)*cl(2)-bl(2)*cl(1))
      r(1)=-al(0)*(bl(2)*cl(3)-bl(3)*cl(2))
     &     +al(2)*(bl(0)*cl(3)-bl(3)*cl(0))
     &     -al(3)*(bl(0)*cl(2)-bl(2)*cl(0))
      r(2)= al(0)*(bl(1)*cl(3)-bl(3)*cl(1))
     &     -al(1)*(bl(0)*cl(3)-bl(3)*cl(0))
     &     +al(3)*(bl(0)*cl(1)-bl(1)*cl(0))
      r(3)=-al(0)*(bl(1)*cl(2)-bl(2)*cl(1))
     &     +al(1)*(bl(0)*cl(2)-bl(2)*cl(0))
     &     -al(2)*(bl(0)*cl(1)-bl(1)*cl(0))
      return
      end


      subroutine lorentz_KbartoK(K,Kbar,pin,pout)
c     Lorentz transformation that maps Kbar onto K (K.K = Kbar.Kbar):
c     pout = pin - 2 pin.(K+Kbar)/(K+Kbar)^2 (K+Kbar)
c                + 2 pin.Kbar/Kbar^2 K
c     (initial-initial mapping of 2209.09123, App. B; with K and Kbar
c     exchanged it is the inverse transformation)
      implicit none
      double precision K(0:3),Kbar(0:3),pin(0:3),pout(0:3)
      double precision KpKb(0:3),tmp(0:3)
      double precision dot
c
      KpKb(:)=K(:)+Kbar(:)
      tmp(:)=pin(:)-2d0*dot(pin,KpKb)/dot(KpKb,KpKb)*KpKb(:)
     &     +2d0*dot(pin,Kbar)/dot(Kbar,Kbar)*K(:)
      pout(:)=tmp(:)
      return
      end
