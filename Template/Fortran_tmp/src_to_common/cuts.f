      logical function docut(p,npart)
      implicit none

      include 'jets.inc'

      integer npart,i,j,nQCD
      double precision p(0:3,npart)
      double precision rfj,sycut,palg,pQCD(0:3,maxdim),etamax
      logical dojets
      parameter(dojets=.true.)
      double precision dot
      double precision, parameter :: tiny=1d-8
c
c     local variables
      double precision x1,x2
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
         do j=3,npart
            nQCD=nQCD+1
            do i=0,3
               pQCD(i,nQCD)=p(i,j)
            enddo
         enddo
c     
c     clustering parameters
         palg=-1d0
         rfj=0.4d0
         sycut=20d0
         etamax = 5d0

c******************************************************************************
c     call FASTJET to get all the jets
c     
c     INPUT:
c     input momenta:               pQCD(0:3,n), energy is 0th component
c     number of input momenta:     nQCD
c     radius parameter:            rfj
c     minumum jet pt:              sycut
c     jet algorithm:               palg, 1.0=kt, 0.0=C/A, -1.0 = anti-kt
c     maximum jet rapidity:        etamax    
c     
c     OUTPUT:
c     jet momenta:                           pjet(0:3,n), E is 0th cmpnt
c     the number of jets (with pt > SYCUT):  njet
c     the jet for a given particle 'i':      jet(i),   note that this is the
c     particle in pQCD, which doesn't
c     necessarily correspond to the particle
c     label in the process
c
         call fastjetppgenkt_etamax(pQCD,nQCD,rfj,sycut,etamax,palg,pjet,njet,jet)

         
C     call fastjetppgenkt(pQCD,nQCD,rfj,sycut,palg,pjet,njet,jet)

c     
c******************************************************************************
         do i=1,njet
            ptjet(i)=sqrt(pjet(1,i)**2+pjet(2,i)**2)
            if(i.gt.1)then
               if (ptjet(i)-ptjet(i-1).gt.tiny) then
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
      if(njet.le.2)docut=.true.
c
      return
      end


