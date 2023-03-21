      function docut(p,npart)
      implicit none
      integer npart,i,j,nQCD
      double precision p(0:3,npart)
      double precision rfj,sycut,palg,pQCD(0:3,npart)
      logical dojets
c      parameter(dojets=.true.)
      parameter(dojets=.false.)
c     
c     initialise
      docut=.true.
      pjet=0d0
      jet=0
      njet=0
      ptjet=0d0
c
c     cluster partons into jets
      if(dojets)then
         nQCD=0
         do j=1,npart
            nQCD=nQCD+1
            do i=0,3
               pQCD(i,nQCD)=p(i,j)
            enddo
         enddo
c     
c     clustering parameters
         palg=-1d0
         rfj=0.4d0
         sycut=10d0
c******************************************************************************
c     call FASTJET to get all the jets
c     
c     INPUT:
c     input momenta:               pQCD(0:3,n), energy is 0th component
c     number of input momenta:     nQCD
c     radius parameter:            rfj
c     minumum jet pt:              sycut
c     jet algorithm:               palg, 1.0=kt, 0.0=C/A, -1.0 = anti-kt
c     
c     OUTPUT:
c     jet momenta:                           pjet(0:3,n), E is 0th cmpnt
c     the number of jets (with pt > SYCUT):  njet
c     the jet for a given particle 'i':      jet(i),   note that this is the
c     particle in pQCD, which doesn't
c     necessarily correspond to the particle
c     label in the process
c     
         call fastjetppgenkt(pQCD,nQCD,rfj,sycut,palg,pjet,njet,jet)
c     
c******************************************************************************
         do i=1,njet
            ptjet(i)=sqrt(pjet(1,i)**2+pjet(2,i)**2)
            if(i.gt.1)then
               if (ptjet(i).gt.ptjet(i-1)) then
                  write (*,*) 'Error 1 in docut: jets unordered in pt'
                  stop
               endif
            endif
         enddo
      endif
c     
c     user defined cuts
      docut=.false.
c
      return
      end


      function dotechcut(x,npart,tiny)
      implicit none
      include 'setup.inc'
      integer mxdim
      parameter(mxdim=30)
      integer npart,i
      double precision x(mxdim),tiny
      logical dotechcut
c
c     initialise
      dotechcut=.false.
c
      if(npart.eq.npartNLO)then
         do i=1,2
            if(x(i).le.tiny.or.x(i).ge.1d0-tiny)dotechcut=.true.
         enddo
      elseif(npart.eq.npartNNLO)then
         do i=1,5
            if(i.eq.3)cycle
            if(x(i).le.tiny.or.x(i).ge.1d0-tiny)dotechcut=.true.
         enddo
      else
         write(*,*)'Wrong npart in dotechcut',npart
         stop
      endif
c
      return
      end
