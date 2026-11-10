c=======================================================================
c     MadNkLO, hadron beams (lpp != 0): parton luminosities from LHAPDF
c
c     The processes of a P1_* directory share one matrix element
c     (evaluated for the first of them), and possibly have mirror
c     processes, with the two beams exchanged (the matrix element is
c     then evaluated with the momenta of the two initial partons
c     exchanged). initial_states.inc, written at generation time, lists
c     their initial-state flavours (beam 1, beam 2):
c       n_init_states, init_states(2,n_init_states),
c       init_mirror(n_init_states) (.true. for the mirror processes),
c       has_mirror
c     The PDFs are those of the LHAPDF set lhaid of the run_card
c     (pdfwrap, called by setrun), always through LHAPDF.
c=======================================================================
      double precision function pdf_lumi(x1,x2,muf,imirror)
c     Sum over the processes of this directory (imirror = 1) or over
c     their mirrors (imirror = 2) of f_a(x1,muf) f_b(x2,muf), where
c     (a,b) are the flavours from beam 1 and beam 2 (antiprotons,
c     lpp = -1: charge-conjugated proton PDFs).
      implicit none
      include 'run.inc'
      include 'initial_states.inc'
      double precision x1,x2,muf
      integer imirror,k
      logical mirror
      double precision pdg2pdf
      external pdg2pdf
      logical firsttime
      data firsttime/.true./
      save firsttime
c
      if(firsttime)then
         call check_beams()
         firsttime=.false.
      endif
      mirror=(imirror.eq.2)
      pdf_lumi=0d0
      do k=1,n_init_states
         if(init_mirror(k).eqv.mirror)then
            pdf_lumi=pdf_lumi
     &           +pdg2pdf(abs(lpp(1)),sign(1,lpp(1))*init_states(1,k),
     &                    x1,muf)
     &           *pdg2pdf(abs(lpp(2)),sign(1,lpp(2))*init_states(2,k),
     &                    x2,muf)
         endif
      enddo
c
      return
      end


      subroutine check_beams()
c     Only lepton beams without PDF (lpp = 0) and (anti)protons with
c     PDFs (lpp = 1, -1) are supported
      implicit none
      include 'run.inc'
      integer i
      do i=1,2
         if(abs(lpp(i)).gt.1)then
            write(*,*)'run_card.dat: lpp',i,' =',lpp(i),
     &           ' is not supported by MadNkLO: use 0 (no PDF), ',
     &           '1 (proton) or -1 (antiproton)'
            stop 1
         endif
      enddo
      return
      end


      double precision function pdf_swap_initial(p,pm,npart)
c     pm = p with the momenta of the two initial partons exchanged
c     (matrix element of the mirror processes); returns 1
      implicit none
      integer npart,i,j
      double precision p(0:3,npart),pm(0:3,npart)
      do i=1,npart
         do j=0,3
            pm(j,i)=p(j,i)
         enddo
      enddo
      do j=0,3
         pm(j,1)=p(j,2)
         pm(j,2)=p(j,1)
      enddo
      pdf_swap_initial=1d0
      return
      end


      subroutine set_tau_min(npart,pdgs)
c     Lower bound tau_min on tau = x1 x2 = shat/stot for the sampling
c     of x1, x2 (genps.f, GENCMS_LOGTAU), from the final-state masses
c     and the generation cuts that every phase-space point of this
c     (Born) process must pass (cuts.f, docut with nUnres = 0):
c       sqrt(shat) >= sum_i sqrt(m_i^2 + ptmin_i^2), with ptmin = ptl
c       for charged leptons, ptj for jet partons (each parton is a jet
c       at Born level), ptgmin for photons;
c       shat >= mll^2 (opposite-sign lepton pair), mll_sf^2 (same
c       flavour).
c     If the bound is 0 (no masses nor cuts), tau_min = 0 and genps.f
c     keeps the flat sampling in tau (GENCMS), with a warning.
      implicit none
      include 'nexternal.inc'
      include 'run.inc'
      include 'cuts.inc'
      integer npart,pdgs(npart)
      integer i,j,id
      double precision rs,pt0,shmin,stot,mb(2),pb(2)
      double precision tau_min
      common/to_tau_min/tau_min
      logical lep_m,lep_p
      double precision pmass(nexternal)
      double precision, parameter :: ZERO=0d0
      include 'coupl.inc'
c
      include 'pmass.inc'
c
      if(npart.ne.nexternal)then
         write(*,*)'set_tau_min: npart =',npart,' != nexternal'
         stop 1
      endif
c     hadronic centre-of-mass energy squared, as in gen_mom (genps.f)
      do i=1,2
         mb(i)=0d0
         if(abs(lpp(i)).eq.1.or.abs(lpp(i)).eq.2)mb(i)=0.938d0
         if(abs(lpp(i)).eq.3)mb(i)=0.000511d0
         pb(i)=sqrt(max(max(ebeam(i),mb(i))**2-mb(i)**2,0d0))
      enddo
      stot=mb(1)**2+mb(2)**2+2d0*(max(ebeam(1),mb(1))*
     &     max(ebeam(2),mb(2))+pb(1)*pb(2))
c
      rs=0d0
      do i=3,npart
         id=pdgs(i)
         pt0=0d0
         if(abs(id).eq.11.or.abs(id).eq.13.or.abs(id).eq.15)then
            pt0=max(ptl,0d0)
         elseif((abs(id).le.maxjetflavor.and.id.ne.0).or.id.eq.21)then
            pt0=max(ptj,0d0)
         elseif(id.eq.22)then
            pt0=max(ptgmin,0d0)
         endif
         rs=rs+sqrt(pmass(i)**2+pt0**2)
      enddo
      shmin=rs**2
      do i=3,npart
         lep_m=pdgs(i).eq.11.or.pdgs(i).eq.13.or.pdgs(i).eq.15
         if(.not.lep_m)cycle
         do j=3,npart
            lep_p=pdgs(j).eq.-11.or.pdgs(j).eq.-13.or.pdgs(j).eq.-15
            if(.not.lep_p)cycle
            if(mll.gt.0d0)shmin=max(shmin,mll**2)
            if(pdgs(j).eq.-pdgs(i).and.mll_sf.gt.0d0)
     &           shmin=max(shmin,mll_sf**2)
         enddo
      enddo
      tau_min=min(shmin/stot,1d0)
      if(tau_min.ge.1d0)then
         write(*,*)'set_tau_min: the cuts and masses require ',
     &        'sqrt(shat) >=',sqrt(shmin),' > sqrt(s) =',sqrt(stot)
         stop 1
      endif
      if(tau_min.gt.0d0)then
         write(*,*)'Hadron beams: x1 x2 sampled in log, sqrt(shat) >=',
     &        sqrt(shmin),' GeV (masses and cuts)'
      else
         write(*,*)'WARNING: hadron beams, but no lower bound on ',
     &        'shat from masses and cuts: x1 x2 sampled uniformly'
      endif
c
      return
      end
