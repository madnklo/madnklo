      subroutine get_W_NLO(xsNLO,alpha,beta,iI,iJ,W_NLO,WS_NLO,WC_NLO,ierr)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer i,k,l,iI,iJ,ierr
      double precision xsNLO(-2:maxdim,-2:maxdim)
      double precision W_NLO(maxdim,maxdim)
      double precision WS_NLO(maxdim,maxdim),WC_NLO(maxdim,maxdim)
      double precision sigma(maxdim,maxdim),sigma2
      double precision sigmas(maxdim,maxdim),sigma2s(maxdim)
      double precision EE(maxdim),omega(maxdim,maxdim)
      double precision alpha,beta,sumW
      double precision sumWsoft(maxdim),sumWcoll(maxdim,maxdim)
      logical :: is_first = .true.
c
c     initialise
      ierr=0
      EE=0d0
      omega=0d0
      sigma=0d0
      sigmas=0d0
      sigma2=0d0
      sigma2s=0d0
      W_NLO=0d0
      WS_NLO=0d0
c
c     input checks
      if(is_first)then
         if(iI.eq.iJ.or.
     &      iI.le.0.or.iI.gt.npartNLO.or.
     &      iJ.le.0.or.iJ.gt.npartNLO)then
            write(*,*)'Input error 1 in get_W_NLO',iI,iJ
            stop
         endif
         if(alpha.le.0d0.or.beta.le.0d0)then
            write(*,*)'Input error 3 in get_W_NLO',alpha,beta
            stop
         endif
         is_first=.false.
      endif
c
      do k=1,npartNLO
         EE(k)=xsNLO(0,k)/sCM
         if(EE(k).le.0d0)then
            write(77,*)'Error 1 in get_W_NLO',k,EE(k)
            goto 999
         endif
         do l=1,npartNLO
            if(l.eq.k)cycle
            omega(k,l)=sCM*xsNLO(k,l)/xsNLO(0,k)/xsNLO(0,l)
            if(omega(k,l).lt.0d0)then
               write(77,*)'Error 2 in get_W_NLO',k,l,omega(k,l)
               goto 999
            endif
c     build the full sigma denominator
            sigma(k,l)=1d0/(EE(k)**alpha*omega(k,l)**beta)
            sigma2=sigma2+sigma(k,l)
c     build the soft sigma denominator
            sigmas(k,l)=1d0/omega(k,l)**beta
            sigma2s(k)=sigma2s(k)+sigmas(k,l)
         enddo
         if(sigma2s(k).le.0d0)then
            write(77,*)'Error 3 in get_W_NLO',k,sigma2s(k)
            goto 999
         endif
      enddo
      if(sigma2.le.0d0)then
         write(77,*)'Error 4 in get_W_NLO',sigma2
         goto 999
      endif
c
c     full NLO W function
      W_NLO(iI,iJ)=sigma(iI,iJ)/sigma2
      W_NLO(iJ,iI)=sigma(iJ,iI)/sigma2
c     single-soft NLO W function
      WS_NLO(iI,iJ)=sigmas(iI,iJ)/sigma2s(iI)
      WS_NLO(iJ,iI)=sigmas(iJ,iI)/sigma2s(iJ)
c     single-collinear NLO W function
      WC_NLO(iI,iJ)=EE(iJ)**alpha/(EE(iI)**alpha+EE(iJ)**alpha)
      WC_NLO(iJ,iI)=EE(iI)**alpha/(EE(iJ)**alpha+EE(iI)**alpha)
c
c     output checks commented out since they would need this routine
c     to compute the whole matrix Wab, while in this version it computes
c     only elements ab = IJ, JI
c$$$c
c$$$c     output checks
c$$$      sumW=0d0
c$$$      sumWsoft=0d0
c$$$      sumWcoll=0d0
c$$$      do k=1,npartNLO
c$$$         do l=1,npartNLO
c$$$            if(l.eq.k)cycle
c$$$            sumW=sumW+W_NLO(k,l)
c$$$            sumWsoft(k)=sumWsoft(k)+WS_NLO(k,l)
c$$$            sumWcoll(k,l)=WC_NLO(k,l)+WC_NLO(l,k)
c$$$         enddo
c$$$      enddo
c$$$      if(abs(sumW-1d0).gt.tiny1)then
c$$$         write(77,*)'Output error 1 in get_W_NLO',sumW
c$$$         goto 999
c$$$      endif 
c$$$      do k=1,npartNLO
c$$$         if(abs(sumWsoft(k)-1d0).gt.tiny1)then
c$$$            write(77,*)'Output error 2 in get_W_NLO',k,sumWsoft(k)
c$$$            goto 999
c$$$         endif
c$$$         do l=1,npartNLO
c$$$            if(l.eq.k)cycle
c$$$            if(abs(sumWcoll(k,l)-1d0).gt.tiny1)then
c$$$               write(77,*)'Output error 3 in get_W_NLO',k,l,sumWcoll(k,l)
c$$$               goto 999
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$c
      return
 999  ierr=1
      return
      end


      subroutine get_W_NNLO(xsNNLO,xp,xsb,xpb,alpha,beta,alpha2,beta2,iI,iJ,iK,iL,
     &iR,do1unres,Wsum,WsumS,WsumC,WsumSS,WsumSSC,ierr)
c     indices iI, iJ, iK, iL, and iR are relevant only for the single-unresolved
c     limits WS_NNLO and WC_NNLO
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer i,j,k,l,iii,khat,lhat,kmap,kb,lb,ierr,itopo
      integer iI,iJ,iK,iL,iR,iA,iB
      double precision xsNNLO(-2:maxdim,-2:maxdim),xsb(-2:maxdim,-2:maxdim)
      double precision xp(0:3,-2:npartNNLO),xpb(0:3,-2:npartNLO)
      double precision Wsum,WsumS,WsumC,WsumSS,WsumSSC
      double precision W_NLO(maxdim,maxdim)
      double precision WS_NLO(maxdim,maxdim),WC_NLO(maxdim,maxdim)
      double precision sigma(maxdim,maxdim,maxdim,maxdim),sigma4
      double precision sigmass(maxdim,maxdim),sigmass_tmp(maxdim,maxdim)
      double precision sigmacc(maxdim,maxdim,maxdim),sigmasscc(maxdim,maxdim,maxdim)
      double precision sigmacc2(maxdim,maxdim,maxdim,maxdim),sigmasscc2(maxdim,maxdim,maxdim,maxdim)
      double precision EE(maxdim),EEb(maxdim),omega(maxdim,maxdim)
      double precision alpha,beta,alpha2,beta2
      double precision sumW,sumWdsoft(maxdim,maxdim)
      double precision EEalpha(maxdim),omegabeta(maxdim,maxdim),sum_inv_omegabeta(maxdim)
      logical do1unres
      logical :: is_first = .true.
      logical :: is_first_1unres = .true.
c
c     initialise
      ierr=0
      EE=0d0
      EEalpha=0d0
      omega=0d0
      omegabeta=0d0
      sum_inv_omegabeta=0d0
      sigma=0d0
      sigma4=0d0
      sigmass=0d0
      sigmacc=0d0
      sigmasscc=0d0
      sigmass_tmp=0d0
      Wsum=0d0
      WsumS=0d0
      WsumC=0d0
      WsumSS=0d0
      WsumSSC=0d0
      itopo=topology(iI,iJ,iK,iL)
c
c     input checks
      if(is_first)then
         if(iI.eq.iJ.or.iI.eq.iK.or.iI.eq.iL.or.iK.eq.iL.or.
     &      iI.le.0.or.iI.gt.npartNNLO.or.
     &      iJ.le.0.or.iJ.gt.npartNNLO.or.
     &      iK.le.0.or.iK.gt.npartNNLO.or.
     &      iL.le.0.or.iL.gt.npartNNLO)then
            write(*,*)'Input error 1 in get_W_NNLO',iI,iJ,iK,iL
            stop
         endif
         if(alpha.le.beta.or.beta.le.1d0)then
            write(*,*)'Input error 2 in get_W_NNLO',alpha,beta
            stop
         endif
         is_first=.false.
      endif
c
      do k=1,npartNNLO
         EE(k)=xsNNLO(0,k)/sCM
         EEalpha(k)=EE(k)**alpha
         if(EE(k).le.0d0)then
            write(77,*)'Error 1 in get_W_NNLO',k,EE(k)
            goto 999
         endif
         do l=1,npartNNLO
            if(l.eq.k)cycle
            omega(k,l)=sCM*xsNNLO(k,l)/xsNNLO(0,k)/xsNNLO(0,l)
            omegabeta(k,l)=omega(k,l)**beta
            sum_inv_omegabeta(k)=sum_inv_omegabeta(k)+1d0/omegabeta(k,l)
            if(omega(k,l).lt.0d0)then
               write(77,*)'Error 2 in get_W_NNLO',k,l,omega(k,l)
               goto 999
            endif
         enddo
      enddo
c
      do i=1,npartNNLO
         do j=1,npartNNLO
            if(j.eq.i)cycle
            do k=1,npartNNLO
               if(k.eq.i)cycle
               do l=1,npartNNLO
                  if(l.eq.i.or.l.eq.k)cycle
                  if(j.eq.k)then
                     sigma(i,j,k,l)=(EE(k)+EE(i))*omega(k,l)
                  else
                     sigma(i,j,k,l)=EE(k)*omega(k,l)
                  endif
c     build the full sigma denominator
                  sigma(i,j,k,l)=1d0/(EEalpha(i)*omegabeta(i,j)*sigma(i,j,k,l))
                  sigma4=sigma4+sigma(i,j,k,l)
c     build half of the double-soft sigma denominator
                  sigmass_tmp(i,k)=sigmass_tmp(i,k)+sigma(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      if(sigma4.le.0d0)then
         write(77,*)'Error 3 in get_W_NNLO',sigma4
         goto 999
      endif
c
      do i=1,npartNNLO-1
         do k=i+1,npartNNLO
c     build the double-soft sigma denominator
            sigmass(i,k)=sigmass_tmp(i,k)+sigmass_tmp(k,i)
            sigmass(k,i)=sigmass(i,k)
         enddo
      enddo
c
      do i=1,npartNNLO
         do k=1,npartNNLO
            if(k.eq.i)cycle
            do j=1,npartNNLO
               if(j.eq.i.or.j.eq.k)cycle
c     build the double-collinear sigma denominator for topologies 1 and 2
               sigmacc(i,j,k)=sigma(i,j,j,k)+sigma(i,j,k,j)+sigma(i,k,k,j)
     &                       +sigma(i,k,j,k)+sigma(j,i,i,k)+sigma(j,i,k,i)
     &                       +sigma(j,k,k,i)+sigma(j,k,i,k)+sigma(k,i,i,j)
     &                       +sigma(k,i,j,i)+sigma(k,j,j,i)+sigma(k,j,i,j)
c     build the double-soft double-collinear sigma denominator for topologies 1 and 2
               sigmasscc(i,j,k)=sigma(i,j,j,k)+sigma(i,k,j,k)
     &                         +sigma(j,i,i,k)+sigma(j,k,i,k)
               sigmasscc(i,k,j)=sigma(i,k,k,j)+sigma(i,j,k,j)
     &                         +sigma(k,i,i,j)+sigma(k,j,i,j)
               sigmasscc(k,j,i)=sigma(k,j,j,i)+sigma(k,i,j,i)
     &                         +sigma(j,k,k,i)+sigma(j,i,k,i)
            enddo
         enddo
      enddo
c
      do i=1,npartNNLO
         do j=1,npartNNLO
            if(j.eq.i)cycle
            do k=1,npartNNLO
               if(k.eq.i)cycle
               do l=1,npartNNLO
                  if(l.eq.i.or.l.eq.k)cycle
c     build the double-collinear sigma denominator for topology 3
                  sigmacc2(i,j,k,l)=sigma(i,j,k,l)+sigma(i,j,l,k)+sigma(j,i,k,l)+sigma(j,i,l,k)
     &                             +sigma(k,l,i,j)+sigma(k,l,j,i)+sigma(l,k,i,j)+sigma(l,k,j,i)
c     build the double-soft double-collinear sigma denominator for topology 3
                  sigmasscc2(i,j,k,l)=sigma(i,j,k,l)+sigma(k,l,i,j)
               enddo
            enddo
         enddo
      enddo
c
c     full NNLO W function, together with its double-unresolved limits
         if(itopo.eq.1)then
            Wsum=(sigma(iI,iJ,iJ,iL)+sigma(iI,iL,iL,iJ)+sigma(iJ,iI,iI,iL)+
     &            sigma(iJ,iL,iL,iI)+sigma(iL,iI,iI,iJ)+sigma(iL,iJ,iJ,iI)+
     &            sigma(iI,iJ,iL,iJ)+sigma(iI,iL,iJ,iL)+sigma(iJ,iI,iL,iI)+
     &            sigma(iJ,iL,iI,iL)+sigma(iL,iI,iJ,iI)+sigma(iL,iJ,iI,iJ))/sigma4
c     This WsumSS assumes sector labels are such that, in itopo 1, iI and iJ are the ones that become soft,
c     thus keep only sectors with iI and iJ in first and third position. If more than one pair could go soft,
c     one would need a different WsumSS per pair
            WsumSS=(sigma(iI,iJ,iJ,iL)+sigma(iJ,iI,iI,iL)+
     &              sigma(iI,iL,iJ,iL)+sigma(iJ,iL,iI,iL))/sigmass(iI,iJ)
         elseif(itopo.eq.3)then
            Wsum=(sigma(iI,iJ,iK,iL)+sigma(iI,iJ,iL,iK)+sigma(iI,iK,iJ,iL)+
     &            sigma(iI,iK,iL,iJ)+sigma(iI,iL,iJ,iK)+sigma(iI,iL,iK,iJ)+
     &            sigma(iJ,iI,iK,iL)+sigma(iJ,iI,iL,iK)+sigma(iJ,iK,iI,iL)+
     &            sigma(iJ,iK,iL,iI)+sigma(iJ,iL,iI,iK)+sigma(iJ,iL,iK,iI)+
     &            sigma(iK,iI,iJ,iL)+sigma(iK,iI,iL,iJ)+sigma(iK,iJ,iI,iL)+
     &            sigma(iK,iJ,iL,iI)+sigma(iK,iL,iI,iJ)+sigma(iK,iL,iJ,iI)+ 
     &            sigma(iL,iI,iJ,iK)+sigma(iL,iI,iK,iJ)+sigma(iL,iJ,iI,iK)+ 
     &            sigma(iL,iJ,iK,iI)+sigma(iL,iK,iI,iJ)+sigma(iL,iK,iJ,iI))/sigma4
c     This WsumSS assumes sector labels are such that, in itopo 3, iI and iJ are the ones that become soft,
c     thus keep only sectors with iI and iJ in first and third position. If more than one pair could go soft,
c     one would need a different WsumSS per pair
            WsumSS=(sigma(iI,iK,iJ,iL)+sigma(iJ,iK,iI,iL)+
     &              sigma(iI,iL,iJ,iK)+sigma(iJ,iL,iI,iK))/sigmass(iI,iJ)
         else
            write(*,*)'Only topology 1, 3 expected in W functions'
            stop
         endif
c
c     single-unresolved limits of the NNLO W function, if needed
      if(.not.do1unres)goto 888
c
      if(itopo.eq.1)then
         iA=imap(iI,iI,iJ,0,0,npartNNLO)
         iB=imap(iL,iI,iJ,0,0,npartNNLO)
         call get_W_NLO(xsb,alpha2,beta2,iA,iB,W_NLO,WS_NLO,WC_NLO,ierr)
         if(ierr.eq.1)goto 999
c     This WsumC assumes sector labels are such that, in itopo 1, iI and iJ are the ones that become collinear,
c     thus keep only sectors with iI and iJ in first and second position.  If more than one pair could go
c     collinear, one would need a different WsumC per pair
         WsumC=W_NLO(iA,iB)+W_NLO(iB,iA)
c     Sectors with iI and iJ in first, second and third position
         WsumSSC=WS_NLO(iA,iB) 
      elseif(itopo.eq.3)then
         iA=imap(iK,iI,iJ,0,0,npartNNLO)
         iB=imap(iL,iI,iJ,0,0,npartNNLO)
         call get_W_NLO(xsb,alpha2,beta2,iA,iB,W_NLO,WS_NLO,WC_NLO,ierr)
         if(ierr.eq.1)goto 999
c     This WsumC assumes sector labels are such that, in itopo 3, iI and iJ are the ones that become collinear,
c     thus keep only sectors with iI and iJ in first and second position.  If more than one pair could go
c     collinear, one would need a different WsumC per pair
         WsumC=W_NLO(iA,iB)+W_NLO(iB,iA)
      endif
 888  continue
c
c     output checks commented out since they would need this routine
c     to compute the whole matrix Wab, while in this version it computes
c     only elements abcd = iI,iJ,iK,iL
c$$$      sumW=0d0
c$$$      sumWdsoft=0d0
c$$$      do i=1,npartNNLO
c$$$         do k=1,npartNNLO
c$$$            if(k.eq.i)cycle
c$$$            do l=1,npartNNLO
c$$$               if(l.eq.i.or.l.eq.k)cycle
c$$$               do j=1,npartNNLO
c$$$                  if(j.eq.i)cycle
c$$$                  sumW=sumW+W_NNLO(i,j,k,l)
c$$$                  sumWdsoft(i,k)=sumWdsoft(i,k)+WSS_NNLO(i,j,k,l)
c$$$               enddo
c$$$               do j=1,npartNNLO
c$$$                  if(j.eq.k)cycle
c$$$                  sumWdsoft(i,k)=sumWdsoft(i,k)+WSS_NNLO(k,j,i,l)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      if(abs(sumW-1d0).gt.tiny1)then
c$$$         write(77,*)'Output error 1 in get_W_NNLO',sumW
c$$$         goto 999
c$$$      endif
c$$$      do i=1,npartNNLO
c$$$         do k=1,npartNNLO
c$$$            if(k.eq.i)cycle
c$$$            if(abs(sumWdsoft(i,k)-1d0).gt.tiny1)then
c$$$               write(77,*)'Output error 2 in get_W_NNLO',i,k,sumWdsoft(i,k)
c$$$               goto 999
c$$$            endif
c$$$         enddo
c$$$      enddo
c
      return
 999  ierr=1
      return
      end
