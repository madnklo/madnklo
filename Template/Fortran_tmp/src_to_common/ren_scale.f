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


c=======================================================================
c     Renormalisation-scale variations by reweighting
c
c     The scale factors k = sv_fac(1..sv_n) (muR -> k*muR) are read by
c     sv_init from rw_rscale in run_card.dat (reweight_scale = True).
c     Every contribution is a pure power of alpha_s, so that for each
c     phase-space point its value at k*muR is obtained from the one at
c     muR without new matrix-element evaluations:
c       B                : r_k^b
c       R, local CTs     : r_k^(b+1), each at the scale of its own
c                          (mapped) kinematics
c       V+I              : r_k^(b+1) * [ (V+I)(muR)
c                          + b beta0/2 alphas/(2 pi) B log(k^2) ]
c     with r_k = alpha_s(k*muR)/alpha_s(muR) and b = born_qcd_power
c     (qcd_power.inc, written at generation time).
c     The integrands add their pieces with sv_add (vegas integral) and
c     build the histogram weights with sv_wgts. The vegas routine
c     accumulates the integrals and their Monte Carlo errors for each k
c     (sv_vegas_*), with the same statistical treatment as the central
c     value; the grid is adapted on the central value only. Index 0 of
c     the accumulators is the central value rebuilt from the sv_add
c     contributions, which is checked against the vegas result.
c=======================================================================
      module scale_var_module
      implicit none
      integer sv_max
      parameter (sv_max=10)
c     number of muR factors (0: no scale variations) and factors
      integer sv_n
      double precision sv_fac(sv_max)
c     muF factors and dynamical-scale label, only used to lay out the
c     histogram weights as MG5_aMC does (muF does not enter: no PDFs)
      integer sv_nf,sv_dyn
      double precision sv_facf(sv_max)
c     alpha_s(sv_fac(k)*muR_current), set by set_mur
      double precision sv_as(sv_max)
c     current phase-space point
      double precision sv_pt(0:sv_max)
c     vegas accumulators: supercell, iteration, vegas call
      double precision sv_fb(0:sv_max),sv_f2b(0:sv_max)
      double precision sv_ti(0:sv_max),sv_tsi(0:sv_max)
      double precision sv_si(0:sv_max),sv_vs(0:sv_max),sv_sw
      integer sv_nmis
c     sums over channels of the main vegas runs (value, error^2)
      double precision sv_sum(0:sv_max),sv_sum2(0:sv_max)
      integer sv_nmis_tot
      data sv_n/0/
      save
      end module scale_var_module


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
c     With scale variations, also fills sv_as(k) = alpha_s(sv_fac(k)*mu).
      use scale_var_module
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'run.inc'
      double precision mu
      integer k
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
         do k=1,sv_n
            sv_as(k)=alpha_qcd(asmz,nloop,sv_fac(k)*mu)
         enddo
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


      subroutine sv_init()
c     Read the factors of the renormalisation-scale variations from
c     run_card.dat (filled by setrun, which must be called first):
c     variations are switched on if reweight_scale is true (as in
c     MG5_aMC). Only muR is varied (no PDFs, so muF does not enter); the
c     muF factors (rw_fscale) are only used to lay out the histogram
c     weights exactly as MG5_aMC does.
      use scale_var_module
      implicit none
      include 'run.inc'
      include '../SubProcesses/reweight0.inc'
      integer k,n,nf
c
      sv_n=0
      sv_nf=0
c     label of the histogram weights: as in MG5_aMC 3.x, whose run card
c     sets dynamical_scale_choice = -2 when both scales are fixed
      sv_dyn=dyn_scale(1)
      if(fixed_ren_scale.and.fixed_fac_scale)sv_dyn=-2
      if(do_rwgt_scale.and.lscalevar(1))then
         n=nint(scalevarR(0))
         nf=nint(scalevarF(0))
         if(n.gt.sv_max.or.nf.gt.sv_max)then
            write(*,*)'sv_init: too many factors in rw_rscale or ',
     &           'rw_fscale',n,nf,' (maximum',sv_max,')'
            stop
         endif
         do k=1,n
            if(scalevarR(k).le.0d0)then
               write(*,*)'sv_init: non-positive factor in rw_rscale',
     &              scalevarR(k)
               stop
            endif
            sv_fac(k)=scalevarR(k)
         enddo
         do k=1,nf
            sv_facf(k)=scalevarF(k)
         enddo
         sv_n=n
         sv_nf=max(nf,1)
         if(nf.lt.1)sv_facf(1)=1d0
      endif
      if(sv_n.gt.0)then
         write(*,*)'Renormalisation-scale variations, muR/muR0 =',
     &        (sv_fac(k),k=1,sv_n)
      endif
c
      sv_sum=0d0
      sv_sum2=0d0
      sv_nmis_tot=0
      call sv_vegas_reset()
      call sv_point_reset()
c     force set_mur to fill sv_as at its next call
      call sv_force_mur_update()
c
      return
      end


      subroutine sv_force_mur_update()
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      muR_current=-1d0
      return
      end


      subroutine sv_weights_info(nwgt,info,nmax)
c     Number and names of the histogram weights, laid out exactly as in
c     MG5_aMC fixed order (madfks_plot.f, initplot): the central value,
c     then one weight per (muF,muR) pair, muF outer and muR inner loop,
c     named e.g. 'dyn=  -1 muR= 2.000 muF= 1.000'. Without PDFs muF does
c     not enter, so the weights with the same muR coincide.
      use scale_var_module
      implicit none
      integer nwgt,nmax,ii,jj
      character*(*) info(nmax)
c
      nwgt=1+sv_n*sv_nf
      if(nwgt.gt.nmax)then
         write(*,*)'sv_weights_info: too many weights',nwgt,nmax
         stop
      endif
      info(1)='central value               '
      do ii=1,sv_nf
         do jj=1,sv_n
            write(info(1+(ii-1)*sv_n+jj),'(a4,i4,x,a4,f6.3,x,a4,f6.3)')
     &           'dyn=',sv_dyn,'muR=',sv_fac(jj),'muF=',sv_facf(ii)
         enddo
      enddo
c
      return
      end


      subroutine sv_check_power(npow)
      implicit none
      integer npow
      if(npow.lt.0)then
         write(*,*)'The diagrams of the underlying Born do not all ',
     &        'have the same QCD order (born_qcd_power < 0 in ',
     &        'qcd_power.inc):'
         write(*,*)'renormalisation-scale variations by reweighting ',
     &        'are not possible. Set reweight_scale = False in ',
     &        'run_card.dat'
         stop
      endif
      return
      end


      subroutine sv_point_reset()
c     reset the scale-variation contributions of the current point
      use scale_var_module
      implicit none
      sv_pt=0d0
      return
      end


      subroutine sv_add(v,dv,npow)
c     Add to the current phase-space point the contribution v, in the
c     normalisation of the vegas integrand, proportional to alpha_s^npow
c     at the scale last set by set_mur (call it after the set_mur for
c     the kinematics of v). dv is the coefficient of log(k^2) in the
c     explicit muR dependence of v at fixed alpha_s: 0, except for V+I.
      use scale_var_module
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision v,dv
      integer npow,k
c
      if(sv_n.eq.0)return
      call sv_check_power(npow)
      sv_pt(0)=sv_pt(0)+v
      do k=1,sv_n
         sv_pt(k)=sv_pt(k)+(sv_as(k)/alphas_current)**npow*
     &        (v+dv*2d0*log(sv_fac(k)))
      enddo
c
      return
      end


      subroutine sv_wgts(w,dw,npow,wgts)
c     Histogram weights of a contribution w (same conventions as
c     sv_add), in the layout of sv_weights_info: wgts(1) = w, then the
c     value at the scale sv_fac(jj)*muR for each (muF_ii, muR_jj) pair
      use scale_var_module
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision w,dw,wgts(*),wk
      integer npow,ii,jj
c
      wgts(1)=w
      if(sv_n.eq.0)return
      call sv_check_power(npow)
      do jj=1,sv_n
         wk=(sv_as(jj)/alphas_current)**npow*(w+dw*2d0*log(sv_fac(jj)))
         do ii=1,sv_nf
            wgts(1+(ii-1)*sv_n+jj)=wk
         enddo
      enddo
c
      return
      end


      subroutine sv_vegas_reset()
c     called by vegas when its results are reset (init <= 1)
      use scale_var_module
      implicit none
      sv_fb=0d0
      sv_f2b=0d0
      sv_ti=0d0
      sv_tsi=0d0
      sv_si=0d0
      sv_vs=0d0
      sv_sw=0d0
      sv_nmis=0
      return
      end


      subroutine sv_vegas_point(f,wgt)
c     called by vegas after each point: f = wgt*integrand. A point with
c     f = 0 (cut, or failed: the integrand returns 0) gives no
c     contribution for any scale. Points whose central value rebuilt
c     from sv_add differs from f are counted (sv_nmis) and reported.
      use scale_var_module
      implicit none
      double precision f,wgt,fk
      integer k
c
      if(sv_n.eq.0)return
      if(f.eq.0d0)return
      if(abs(wgt*sv_pt(0)-f).gt.1d-8*max(abs(f),abs(wgt*sv_pt(0))))
     &     sv_nmis=sv_nmis+1
      do k=0,sv_n
         fk=wgt*sv_pt(k)
         sv_fb(k)=sv_fb(k)+fk
         sv_f2b(k)=sv_f2b(k)+fk**2
      enddo
c
      return
      end


      subroutine sv_vegas_end_cell(npg)
c     called by vegas at the end of each supercell (npg points): same
c     variance estimate as for the central value
      use scale_var_module
      implicit none
      integer npg,k
      double precision t
c
      if(sv_n.eq.0)return
      do k=0,sv_n
         t=sqrt(sv_f2b(k)*npg)
         t=(t-sv_fb(k))*(t+sv_fb(k))
         if(t.le.0d0)t=1d-30
         sv_ti(k)=sv_ti(k)+sv_fb(k)
         sv_tsi(k)=sv_tsi(k)+t
         sv_fb(k)=0d0
         sv_f2b(k)=0d0
      enddo
c
      return
      end


      subroutine sv_vegas_end_iter(dv2g,w)
c     called by vegas at the end of each iteration. Iterations are
c     combined with the weights w = 1/sigma^2 of the central value
c     (as vegas does for the central value), so that the variance of
c     the combination is sum(w^2 sigma_k^2)/(sum w)^2.
      use scale_var_module
      implicit none
      double precision dv2g,w
      integer k
c
      if(sv_n.eq.0)return
      do k=0,sv_n
         sv_si(k)=sv_si(k)+w*sv_ti(k)
         sv_vs(k)=sv_vs(k)+(w*sv_tsi(k)*dv2g)*w
         sv_ti(k)=0d0
         sv_tsi(k)=0d0
      enddo
      sv_sw=sv_sw+w
c
      return
      end


      subroutine sv_add_channel()
c     add the results of the last vegas run to the sums over channels
c     (call it after each main vegas run)
      use scale_var_module
      implicit none
      integer k
c
      if(sv_n.eq.0)return
      if(sv_sw.le.0d0)return
      do k=0,sv_n
         sv_sum(k)=sv_sum(k)+sv_si(k)/sv_sw
         sv_sum2(k)=sv_sum2(k)+(sv_vs(k)/sv_sw)/sv_sw
      enddo
      sv_nmis_tot=sv_nmis_tot+sv_nmis
c
      return
      end


      subroutine sv_write_results(iu,tag,central,central_err)
c     write the results for each scale factor in the results_*.log file,
c     after the central line ' sigma <tag> [pb]  = central +- error'
      use scale_var_module
      implicit none
      integer iu,k
      character*(*) tag
      double precision central,central_err
c
      if(sv_n.eq.0)return
      do k=1,sv_n
         write(iu,'(a,a,a,f9.5,a,es23.15e3,a,es23.15e3)')
     &        ' sigma ',trim(tag),' [pb] muR/muR0 = ',sv_fac(k),' : ',
     &        sv_sum(k),' +- ',sqrt(sv_sum2(k))
      enddo
c     cross-checks
      if(abs(sv_sum(0)-central).gt.1d-8*abs(central))then
         write(iu,*)'WARNING: the scale-variation reweighting does ',
     &        'not reproduce the central value:',sv_sum(0),' vs ',
     &        central,' +-',central_err
      endif
      if(sv_nmis_tot.gt.0)then
         write(iu,*)'WARNING:',sv_nmis_tot,' phase-space points with ',
     &        'inconsistent scale-variation contributions'
      endif
c
      return
      end
