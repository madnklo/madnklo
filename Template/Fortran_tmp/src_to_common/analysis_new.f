      subroutine analysis_begin(nwgt,weights_info)
      implicit none
      integer nwgt
      character*(*) weights_info(*)
c
      call HwU_inithist(nwgt,weights_info)
c
      call HwU_book(1,'total  ',2,0d0,2d0)
c      call HwU_book(2,'thrust ',17,0.68d0,1.02d0)
      call HwU_book(3,'pt j1  ',50,2d0,502d0)
      call HwU_book(4,'pt j2  ',50,2d0,502d0)
      call HwU_book(5,'pt j3  ',50,2d0,502d0)
      call HwU_book(6,'pt j4  ',50,2d0,502d0)
c
      call HwU_book(7,'abs eta j1  ',30,0d0,6d0)
      call HwU_book(8,'abs eta j2  ',30,0d0,6d0)
      call HwU_book(9,'abs eta j3  ',30,0d0,6d0)
      call HwU_book(10,'abs eta j4 ',30,0d0,6d0)
c
      call HwU_book(11,'abs y j1  ',30,0d0,6d0)
      call HwU_book(12,'abs y j2  ',30,0d0,6d0)
      call HwU_book(13,'abs y j3  ',30,0d0,6d0)
      call HwU_book(14,'abs y j4  ',30,0d0,6d0)
c
      call HwU_book(15,'njet      ',4,0d0,4d0)
c
      return
      end


      subroutine analysis_end(fname,xresc)
      implicit none
      double precision xresc
      character*(*) fname
      open(unit=98,file=trim(fname),status='unknown')
      call HwU_output(98,xresc)
      close(98)
      return
      end


      subroutine analysis_fill(p,xs,nexternal,wgts)
      implicit none
      include 'jets.inc'
      integer nexternal,i,j,nQCD
      double precision xs(nexternal,nexternal)
      double precision p(0:3,nexternal),wgts(*)
      double precision xsec,thrust
      double precision getthrust_3body,getrapidity,getpseudorap
      double precision rfj,sycut,palg,pQCD(0:3,nexternal)
c
c     observables
      xsec=1d0
c$$$      thrust=getthrust_3body(p,nexternal)
c$$$c     jets
c$$$      pjet=0d0
c$$$      jet=0
c$$$      njet=0
c$$$      ptjet=0d0
c$$$      etajet=-100d0
c$$$      yjet=-100d0
c$$$c
c$$$c     cluster partons into jets
c$$$      nQCD=0
c$$$      do j=3,nexternal
c$$$         nQCD=nQCD+1
c$$$         do i=0,3
c$$$            pQCD(i,nQCD)=p(i,j)
c$$$         enddo
c$$$      enddo
c$$$c     
c$$$c     clustering parameters
c$$$      palg=-1d0
c$$$      rfj=0.4d0
c$$$      sycut=10d0
c$$$      call fastjetppgenkt(pQCD,nQCD,rfj,sycut,palg,pjet,njet,jet)
c$$$c
c$$$c     check on jet pt ordering
c$$$      do i=1,njet
c$$$         ptjet(i)=sqrt(pjet(1,i)**2+pjet(2,i)**2)
c$$$         etajet(i)=getpseudorap(pjet(0,i),pjet(1,i),pjet(2,i),pjet(3,i))
c$$$         yjet(i)=getrapidity(pjet(0,i),pjet(3,i))
c$$$         if(i.gt.1)then
c$$$            if (ptjet(i).gt.ptjet(i-1)) then
c$$$               write (*,*) 'Error 1 in docut: jets unordered in pt'
c$$$               stop
c$$$            endif
c$$$         endif
c$$$      enddo
c
c     fill histograms
      call HwU_fill(1,xsec,wgts)
c$$$      call HwU_fill(2,thrust,wgts)
      if(njet.ge.1)then
         call HwU_fill(3,ptjet(1),wgts)
         call HwU_fill(7,abs(etajet(1)),wgts)
         call HwU_fill(11,abs(yjet(1)),wgts)
      endif
      if(njet.ge.2)then
         call HwU_fill(4,ptjet(2),wgts)
         call HwU_fill(8,abs(etajet(2)),wgts)
         call HwU_fill(12,abs(yjet(2)),wgts)
      endif
      if(njet.ge.3)then
         call HwU_fill(5,ptjet(3),wgts)
         call HwU_fill(9,abs(etajet(3)),wgts)
         call HwU_fill(13,abs(yjet(3)),wgts)
      endif
      if(njet.ge.4)then
         call HwU_fill(6,ptjet(4),wgts)
         call HwU_fill(10,abs(etajet(4)),wgts)
         call HwU_fill(14,abs(yjet(4)),wgts)
      endif
      call HwU_fill(15,dble(njet),wgts)
c
      return
      end


      function getthrust_3body(xp,n)
      implicit none
      integer n,i
      double precision xp(0:3,n)
      double precision getthrust_3body
      double precision tiny, sCM
      parameter(tiny=1d-5)
      double precision dot
c
      sCM = 2d0*dot(xp(:,1),xp(:,2))
      getthrust_3body=0d0
      do i=3,n
         if(2*xp(0,i)/sqrt(sCM).gt.getthrust_3body)then
            getthrust_3body=2*xp(0,i)/sqrt(sCM)
         endif
      enddo
c
c     avoid misbinning
      getthrust_3body=max(getthrust_3body,0.5d0+tiny)
      getthrust_3body=min(getthrust_3body,1d0-tiny)
c
      return
      end


      function getrapidity(en,pl)
      implicit none
      double precision getrapidity,en,pl
      double precision xplus,xminus,y
      double precision tiny
      parameter (tiny=1d-8)
c
      xplus=en+pl
      xminus=en-pl
      if(xplus.gt.tiny.and.xminus.gt.tiny)then
         if( (xplus/xminus).gt.tiny.and.(xminus/xplus).gt.tiny)then
            y=0.5d0*log(xplus/xminus)
         else
            y=sign(1d0,pl)*1d8
         endif
      else
         y=sign(1d0,pl)*1d8
      endif
      getrapidity=y
c
      return
      end


      function getpseudorap(en,ptx,pty,pl)
      implicit none
      real*8 getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th
      parameter (tiny=1.d-5)
c
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      getpseudorap=eta
      return
      end


      function getdelphi(ptx1,pty1,ptx2,pty2)
      implicit none
      real*8 getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
      parameter (tiny=1.d-5)
c
      pt1=sqrt(ptx1**2+pty1**2)
      pt2=sqrt(ptx2**2+pty2**2)
      if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
        tmp=ptx1*ptx2+pty1*pty2
        tmp=tmp/(pt1*pt2)
        if(abs(tmp).gt.1.d0+tiny)then
          write(*,*)'Cosine larger than 1'
          stop
        elseif(abs(tmp).ge.1.d0)then
          tmp=sign(1.d0,tmp)
        endif
        tmp=acos(tmp)
      else
        tmp=1.d8
      endif
      getdelphi=tmp
      return
      end


      function getdr(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
      implicit none
      real*8 getdr,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi,
     & getpseudorap,getdelphi
c
      deta=getpseudorap(en1,ptx1,pty1,pl1)-
     &     getpseudorap(en2,ptx2,pty2,pl2)
      dphi=getdelphi(ptx1,pty1,ptx2,pty2)
      getdr=sqrt(dphi**2+deta**2)
c
      return
      end


      function getdry(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
      implicit none
      real*8 getdry,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi,
     & getrapidity,getdelphi
c
      deta=getrapidity(en1,pl1)-
     &     getrapidity(en2,pl2)
      dphi=getdelphi(ptx1,pty1,ptx2,pty2)
      getdry=sqrt(dphi**2+deta**2)
c
      return
      end


      function getdrysq(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
      implicit none
      real*8 getdrysq,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi,
     & getrapidity,getdelphi
c
      deta=getrapidity(en1,pl1)-
     &     getrapidity(en2,pl2)
      dphi=getdelphi(ptx1,pty1,ptx2,pty2)
      getdrysq=dphi**2+deta**2
c
      return
      end
