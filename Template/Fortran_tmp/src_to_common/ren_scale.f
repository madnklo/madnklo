c=======================================================================
c     Renormalisation scale and alpha_s
c
c     Every piece of the calculation sets the scale right before
c     evaluating a (Born-like) matrix element, using the momenta and
c     flavours passed to that matrix element:
c       RR            -> n+2 body kinematics
c       R, RV         -> n+1 body kinematics
c       B, V, VV, I   -> n   body kinematics
c       counterterms  -> their own remapped kinematics (e.g. one setting
c                        per dipole inside the soft eikonal sums)
c
c     After a call to set_mur(_from_momenta), the common block in
c     math.inc holds
c       muR_current    : the renormalisation scale just set
c       alphas_current : alpha_s(muR_current)
c     and the model parameters (G, AS, MU_R and the MadLoop couplings,
c     including quadruple precision) are updated consistently.
c=======================================================================


      double precision function get_mur(p,npart,pdgs)
c     Renormalisation scale for momenta p(0:3,npart) with flavours pdgs.
c     run_card switches:
c       fixed_ren_scale = T : muR = muR_over_ref * muR_ref_fixed
c       fixed_ren_scale = F : muR = muR_over_ref * mu_dyn(p,pdgs), with
c                             mu_dyn selected by dynamical_scale_choice
c     Any dynamical choice must be infrared safe, i.e. mu(p) -> mu(pbar)
c     in all soft and collinear limits, so that each local counterterm
c     (evaluated at the scale of its own mapped kinematics) still
c     cancels the singularities of the matrix element it subtracts.
      implicit none
      include 'run.inc'
      integer npart
      integer pdgs(npart)
      double precision p(0:3,npart)
      double precision mu
      integer i
c
      if(fixed_ren_scale)then
         mu=muR_ref_fixed
      else
         select case(dynamical_scale_choice)
c        Add further dynamical scale choices here (MadGraph numbering).
c        Each must be infrared safe: mu(p) -> mu(pbar) in all soft and
c        collinear limits.
         case(3)
c           HT/2: half the sum of the transverse masses of the
c           final-state QCD partons (massless: |pT|; soft partons drop
c           out and collinear massless partons add up to their parent)
            mu=0d0
            do i=3,npart
               if(abs(pdgs(i)).le.6.or.pdgs(i).eq.21)
     &            mu=mu+sqrt(max(0d0,p(0,i)**2-p(3,i)**2))
            enddo
            mu=mu/2d0
         case(4)
c           partonic centre-of-mass energy sqrt(shat) = sqrt((p1+p2)^2)
c           (constant for e+e- without ISR: useful as a check against
c           the fixed-scale result)
            mu=sqrt(max(0d0,(p(0,1)+p(0,2))**2-(p(1,1)+p(1,2))**2
     &              -(p(2,1)+p(2,2))**2-(p(3,1)+p(3,2))**2))
         case default
            write(*,*)'get_mur: dynamical_scale_choice =',
     &           dynamical_scale_choice,' is not implemented.'
            write(*,*)'Set fixed_ren_scale = True in run_card.dat'
            stop
         end select
      endif
c
      get_mur=muR_over_ref*mu
      if(get_mur.le.0d0)then
         write(*,*)'get_mur: non-positive scale',get_mur
         stop
      endif
c
      return
      end


      subroutine set_mur(mu)
c     Set the renormalisation scale to mu: fills muR_current and
c     alphas_current (math.inc) and updates the model couplings.
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'run.inc'
      double precision mu
      double precision alpha_qcd
      external alpha_qcd
      logical firsttime
      data firsttime/.true./
      save firsttime
c
      if(mu.le.0d0)then
         write(*,*)'set_mur: non-positive scale',mu
         stop
      endif
c
c     On the first call MU_R still holds the param_card value (it is
c     only changed below, by update_as_param2): warn once that it is
c     not the argument of alpha_s and of the couplings.
      if(firsttime)then
         call warn_param_card_mur(mu,MU_R)
         firsttime=.false.
      endif
c
c     update only if mu changed, or if G was modified elsewhere
c     (relative tolerance: the matrix-element hooks recompute G with a
c     different but equivalent formula, which may differ in the last
c     bit and must not trigger a full coupling update at every call)
      if(mu.ne.muR_current.or.
     &   abs(G-dsqrt(4d0*pi*alphas_current)).gt.1d-10*abs(G))then
         muR_current=mu
         alphas_current=alpha_qcd(asmz,nloop,mu)
c        keep the legacy run.inc variable in sync
         scale=mu
c        model: sets MU_R, G, AS and recomputes the couplings
         call update_as_param2(mu,alphas_current)
         call mp_update_as_param()
      endif
c
      return
      end


      subroutine set_mur_from_momenta(p,npart,pdgs)
c     Set the renormalisation scale for momenta p(0:3,npart) with
c     flavours pdgs (same arguments as docut).
      implicit none
      integer npart
      integer pdgs(npart)
      double precision p(0:3,npart)
      double precision get_mur
      external get_mur
c
      call set_mur(get_mur(p,npart,pdgs))
c
      return
      end


      subroutine warn_param_card_mur(mu,mu_card)
c     One-time warnings about the param_card MU_R, which is never used
c     as argument of alpha_s and of the couplings: the renormalisation
c     scale is always taken from the run_card (see get_mur).
      implicit none
      include 'run.inc'
      double precision mu,mu_card
c
      if(fixed_ren_scale)then
         if(abs(mu-mu_card).gt.1d-8*abs(mu))then
            write(*,*)
            write(*,*)'WARNING: fixed renormalisation scale from ',
     &           'run_card.dat'
            write(*,*)'  muR = muR_over_ref*muR_ref_fixed =',mu
            write(*,*)'  differs from MU_R in param_card.dat =',mu_card
            write(*,*)'  The param_card MU_R is NOT used as argument ',
     &           'of alpha_s and of the couplings:'
            write(*,*)'  the run_card value muR =',mu,' is used.'
            write(*,*)
         endif
      else
         write(*,*)
         write(*,*)'WARNING: dynamical renormalisation scale from ',
     &        'run_card.dat'
         write(*,*)'  (dynamical_scale_choice =',
     &        dynamical_scale_choice,').'
         write(*,*)'  MU_R in param_card.dat =',mu_card,' is NOT used ',
     &        'as argument of alpha_s and of the couplings:'
         write(*,*)'  muR is computed point by point from the ',
     &        'kinematics (times muR_over_ref).'
         write(*,*)
      endif
c
      return
      end
