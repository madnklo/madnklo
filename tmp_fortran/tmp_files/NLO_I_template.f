      subroutine int_counter_NLO(p,sLO,INLO,BLOI,ierr)
c     MSbar integrated counterterm
c     FINITE_PART = INLO(0)
c     SINGLE_POLE = INLO(-1)
c     DOUBLE_POLE = INLO(-2)
c     BLOI = Born matrix element (for the renormalisation-scale variations)
      implicit none
      INCLUDE 'nexternal.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix)s.inc'
      INCLUDE 'colored_partons.inc'
      integer i,j,r
      integer ierr
      double precision p(0:3,nexternal)
      double precision sLO(nexternal,nexternal)
      double precision INLO(-2:0),pref
      double precision BLO,ccBLO,BLOI
      double precision A20a,A21a,A20b,A20,A21,A10,A10g
c     incoming QCD partons (2209.09123, Secs. 4.2, 4.3)
      integer ninit,ia,ib
      double precision betaF,gammaI,CC(nexternal)
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      DOUBLE PRECISION  GET_CCBLO
      integer iref1(nexternal)
      double precision vv,ypl,Q2,ddilog
      double precision pmass(nexternal)
      DOUBLE PRECISION SS,MK2,ML2
      DOUBLE PRECISION FF1,FF2,FF3
      PARAMETER(FF1=1D0,FF2=1D0,FF3=0D0)
      double precision res
c     massive-massless (2503.14629, final-state radiation)
      double precision LKL,LSM,rhom
      double precision IS_M,IHC_0G,IHC_1G
      double precision IHC_2G
      logical massive_coloured
      include 'pmass.inc'
c
c     initialise
      ierr = 0
c     renormalisation scale for the n-body kinematics
      call set_mur_from_momenta(p,nexternal,leg_pdgs_%(proc_prefix)s)
      alphas=alphas_current
      pref=alphas/(2d0*pi)
      INLO = 0d0
      iref1 = 0
      CCBLO = 0d0
      BLO = 0d0
      BLOI = 0d0
      res = 0d0

      CALL ME_ACCESSOR_HOOK(P,HEL,ALPHAS,ANS)
      BLO = ANS(0)
      BLOI = BLO
c
c     TODO: add check 
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
c
c     incoming QCD partons, ninit = 0, 1 (a = ia) or 2 (a,b = 1,2), and
c     damping parameters (damping_factors.inc) of the final-state
c     collinear counterterms (beta of 2209.09123: the recoiler is in
c     the initial state whenever there is an incoming QCD parton) and
c     of the initial-state ones (gamma: beta_IF with a final-state
c     recoiler, ninit = 1; beta_II with an initial-state one)
      ninit=0
      ia=0
      do i=1,2
         if(ISLOQCDPARTON(i))then
            ninit=ninit+1
            ia=i
         endif
      enddo
      betaF=beta_FF
      if(ninit.ge.1)betaF=beta_FI
      gammaI=0d0
      if(ninit.eq.1)gammaI=beta_IF
      if(ninit.eq.2)gammaI=beta_II
c     Massive coloured partons with incoming QCD partons: the soft
c     dipoles (incoming a, massive k) need no term of their own. The
c     pair branches below (massless-massive and massive-massless) plus
c     the colour-conservation constant of the incoming parton (I^I_fin,
c     I^II_fin) give exactly -2 I_s^FmI of the MadNkLO massive paper
c     (eq. IsMI): note massive_isr_phase0_map.md, Sec. 4.2
c     Casimir of each parton (0 for colourless ones)
      do i=1,nexternal
         CC(i)=0d0
         if(leg_pdgs_%(proc_prefix)s(i).eq.21)CC(i)=CA
         if(leg_pdgs_%(proc_prefix)s(i).ne.0.and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6)CC(i)=CF
      enddo
c
c     Born contribution
c     (2209.09123, eq. intFF_fin: phi_k - gamma^hc_k L_kr for a
c     final-state parton k, only - gamma^hc_a L_ar for an incoming one)
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(i.le.2)then
            if(leg_pdgs_%(proc_prefix)s(i).eq.21) then
               INLO(0) = INLO(0) - gamma_hc_g*log(sLO(i,iref1(i))/muR_current**2)
               INLO(0) = INLO(0) + pi**2/12d0 * CA
               INLO(-1) = INLO(-1) + gamma_g
               INLO(-2) = INLO(-2) + CA
            elseif(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6) then
               INLO(0) = INLO(0) - gamma_hc_q*log(sLO(i,iref1(i))/muR_current**2)
               INLO(0) = INLO(0) + pi**2/12d0 * CF
               INLO(-1) = INLO(-1) + gamma_q
               INLO(-2) = INLO(-2) + CF
            endif
            cycle
         endif
         if(leg_pdgs_%(proc_prefix)s(i).eq.21) then
            INLO(0) = INLO(0) + (CA/6d0+2*TR*Nf/3d0)*(log(sLO(i,iref1(i))/muR_current**2)-8d0/3d0)+CA*(6d0-7d0/2d0*zeta2)
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])      
            INLO(0) = INLO(0) + pi**2/12d0 * CA
            INLO(-1) = INLO(-1) + gamma_g
            INLO(-2) = INLO(-2) + CA
         elseif(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6) then
            INLO(0) = INLO(0) + (CF/2d0)*(10d0-7d0*zeta2+log(sLO(i,iref1(i))/muR_current**2))
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
            INLO(0) = INLO(0) + pi**2/12d0 * CF
            INLO(-1) = INLO(-1) + gamma_q
            INLO(-2) = INLO(-2) + CF
         endif
c     hard-collinear counterterm with a massive recoiler r=iref1(i):
c     extra finite term -sum_n gamma^hc_(i,n) I^(ng)_hc,M(s_ir,m_r)
c     (2503.14629, Sec. 2.5 and eq. IhcM)
         if(iref1(i).ne.0)then
         if(pmass(iref1(i)).ne.0d0)then
            rhom=pmass(iref1(i))**2/sLO(i,iref1(i))
            if(leg_pdgs_%(proc_prefix)s(i).eq.21)then
               INLO(0) = INLO(0) + 2d0/3d0*TR*Nf*IHC_0G(rhom) + CA/6d0*IHC_2G(rhom)
            elseif(abs(leg_pdgs_%(proc_prefix)s(i)).le.6)then
               INLO(0) = INLO(0) + CF/2d0*IHC_1G(rhom)
            endif
         endif
         endif
      enddo
c
c     The massive integrated counterterms are computed without damping
c     factors: stop if damping is switched on with massive partons
      massive_coloured=.false.
      do i=1,nexternal
         if(ISLOQCDPARTON(i).and.pmass(i).ne.0d0)massive_coloured=.true.
      enddo
      if(massive_coloured.and.(alpha.ne.0d0.or.beta_FF.ne.0d0))then
         write(*,*)'int_counter_NLO: damping factors (alpha, beta_FF)'
         write(*,*)'are not available with massive coloured partons.'
         write(*,*)'Set them to 0 in Cards/damping_factors.inc'
         stop
      endif
c
c     Include damping factors (2209.09123, eq. intFF_fin: the terms
c     gamma^hc_k A2(beta) for final-state partons k only, and
c     2 A2(alpha) C_j L_jr; the corresponding 2 A2(alpha) L_cd B_cd is
c     in the colour-linked loop below)
      A20a=A20(alpha)
      A21a=A21(alpha)
      A20b=A20(betaF)
      do i=1,nexternal
         if(pmass(i).ne.0d0)cycle
         if(leg_pdgs_%(proc_prefix)s(i).eq.21)INLO(0) = INLO(0) + CA*(A20a*(A20a-2d0*A20b)-A21a)
         if(leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6)INLO(0) = INLO(0) + CF*(A20a*(A20a-2d0*A20b)-A21a)
         if(i.gt.2.and.leg_pdgs_%(proc_prefix)s(i).eq.21)INLO(0) = INLO(0) + (gamma_g-2d0*CA)*A20b
         if(i.gt.2.and.leg_pdgs_%(proc_prefix)s(i).ne.0 .and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6)INLO(0) = INLO(0) + (gamma_q-2d0*CF)*A20b
         if(alpha.ne.0d0.and.CC(i).ne.0d0)INLO(0) = INLO(0) + 2d0*A20a*CC(i)*log(sLO(i,iref1(i))/muR_current**2)
      enddo
c
c     finite contributions of the initial-state counterterms that do not
c     depend on x (2209.09123, I^I_fin and I^II_fin)
      A10g=A10(gammaI)
      if(ninit.eq.1)then
         INLO(0) = INLO(0) + 2d0*CC(ia)*(1d0+zeta2/4d0-A20a*(A10g-A20b-1d0)+A21a)
      elseif(ninit.eq.2)then
         INLO(0) = INLO(0) + (2d0+zeta2/2d0+3d0*A21a-A20a*(2d0*A10g-2d0*A20b+A20a))*(CC(1)+CC(2))
      endif
c
      INLO=INLO*BLO
      if(ninit.eq.2)INLO(0) = INLO(0) + 4d0*(zeta2-1d0+A21a)*GET_CCBLO(1,2)
c
c     Colour-linked-Born contribution
      do i=1,nexternal
         if(.not.ISLOQCDPARTON(i))cycle
         do j=1,nexternal
            if(.not.ISLOQCDPARTON(j))cycle
            if(j.eq.i)cycle
            CCBLO = GET_CCBLO(i,j)
            if(pmass(i).eq.0d0.and.pmass(j).eq.0d0)then
               INLO(0) = INLO(0) + ccBLO*log(sLO(i,j)/muR_current**2)*(2d0-log(sLO(i,j)/muR_current**2)/2d0)
               if(alpha.ne.0d0)INLO(0) = INLO(0) + 2d0*A20a*ccBLO*log(sLO(i,j)/muR_current**2)
               INLO(-1) = INLO(-1) + ccBLO*log(sLO(i,j)/muR_current**2)
            elseif(pmass(i).eq.0d0.and.pmass(j).ne.0d0)then
c     massless i, massive j: integrated massive-massless soft
c     counterterm (2503.14629, eq. IsMF), per ordered pair. The sign of
c     the ln(s_ij/m_j^2) terms follows the paper's integrals (the
c     assembled formula for I in the paper has them with opposite sign)
               LKL=log(sLO(i,j)/muR_current**2)
               LSM=log(sLO(i,j)/pmass(j)**2)
               INLO(0) = INLO(0) + ccBLO*(2d0*LKL-LKL**2/2d0-LKL*LSM -2d0*IS_M(pmass(j)**2/sLO(i,j)))
               INLO(-1) = INLO(-1) + ccBLO*(LKL+LSM)
            elseif(pmass(i).ne.0d0.and.pmass(j).eq.0d0)then
c     massive i, massless j: L_ij term, plus the share
c     -(1/eps+4(1-zeta2)) B_ij of the massive-leg term
c     C_i (1/eps+4(1-zeta2)) B, written through colour conservation
c     as in the massive-massive case below
               LKL=log(sLO(i,j)/muR_current**2)
               INLO(0) = INLO(0) + ccBLO*(LKL-4d0*(1d0-zeta2))
               INLO(-1) = INLO(-1) - ccBLO
            elseif(pmass(i).ne.0d0.and.pmass(j).ne.0d0)then
               ML2=PMASS(I)**2
               MK2=PMASS(J)**2
               SS=SLO(I,J)
               VV=DSQRT(SS**2-4D0*ML2*MK2)/SS
               Q2=SS+ML2+MK2
               YPL=1D0+(DSQRT(ML2)-DSQRT(Q2))*2D0*DSQRT(ML2)/SS
               CALL NLO_I_MASS(ss,vv,mk2,ml2,muR_current,ccBLO,res)
               INLO(0) = INLO(0) + res
               INLO(-1) = INLO(-1) + CCBLO*(-1D0/2D0)*(2D0 - 1D0/VV*DLOG((1D0+VV)/(1D0-VV)) )
            endif
         enddo
      enddo
      INLO = INLO*pref
c
      if(abs(INLO(0)).ge.huge(1d0).or.isnan(INLO(0)))then
         write(77,*)'Exception caught in int_counter_NLO',INLO(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end




      subroutine int_convolution_NLO(pp,sB,imirror,u,nmuf,facf,conv,ierr)
c     x-dependent contribution J(x) + C(x) of the integrated initial-state
c     counterterms and of the MS-bar PDF counterterm (2209.09123, eqs.
c     sigNLOsub, PDFct and J^I, J^II_i of Secs. 4.2, 4.3), in d = 4:
c       conv = sum over the incoming QCD partons i, and over the
c              processes of this directory (imirror = 1) or their
c              mirrors (imirror = 2), of
c              f_o(etab_o) int_0^1 dx [J_i(x) + C_i(x)] f(etab_i/x)/x B,
c     at the Born momenta pp (the swapped ones for the mirrors), with
c     invariants sB. Here x = u (one vegas variable, flat in [0,1]),
c     etab are the Born momentum fractions xbk, o is the other beam, and
c     T(x) = f(etab_i/x)/x for x >= etab_i, 0 below, is the test
c     function of the distributions:
c       (F(x))_+         -> F(x) [T(x) - T(1)]   (+ on the whole F),
c       x/(1-x)_+ in P_gg -> [x T(x) - T(1)]/(1-x),
c       delta(1-x)        -> T(1).
c     The logarithms of J are those of the un-rescaled invariants,
c     ln(s_Born/(x mu^2)) (2209.09123, App. B): -ln x multiplies the
c     regular function, since it vanishes at x = 1. Details and the
c     formula coded: note nlo_isr_phase0_map.md, Sec. 13.5.
c     conv(0) is at the factorisation scale muF; conv(ii), ii = 1..nmuf,
c     at facf(ii)*muF (factorisation-scale variations: the PDFs and
c     log(muF) of the PDF counterterm change, the Born does not).
      implicit none
      INCLUDE 'nexternal.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix)s.inc'
      INCLUDE 'colored_partons.inc'
      INCLUDE 'initial_states.inc'
      integer imirror,nmuf,ierr
      double precision pp(0:3,nexternal),sB(nexternal,nexternal),u
      double precision facf(*),conv(0:nmuf)
      integer i,k,lab,ib,io,ir,ninit,iq,ifl,flb,flo,iv
      integer iref1(nexternal)
      double precision x,lx,l1x,lam,muf,muf0,alphas,BLO,Bab,CC(nexternal)
      double precision betaF,gammaI,A10g,A20a,xa,xb
      double precision eta,etao,Fo,Td,Td1,Tg,Tg1,Tq,sumF,sumD
      double precision LBir,LFir,SD1,ccb,LBik,res,rhoB
      double precision pmass(nexternal)
      double precision ANS(0:NSQSO_BORN)
      INTEGER, PARAMETER :: HEL = - 1
      DOUBLE PRECISION GET_CCBLO,A10,A20,get_muf,pdf_beam
      external get_muf,pdf_beam
      include 'pmass.inc'
c
      ierr=0
      conv=0d0
      x=u
      if(x.le.0d0.or.x.ge.1d0)return
      lx=log(x)
      l1x=log(1d0-x)
c
c     recoilers, incoming QCD partons, Casimirs, damping parameters
c     (as in int_counter_NLO)
      iref1=0
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
      ninit=0
      do i=1,2
         if(ISLOQCDPARTON(i))ninit=ninit+1
      enddo
      if(ninit.eq.0)return
      lam=dble(ninit)
      do i=1,nexternal
         CC(i)=0d0
         if(leg_pdgs_%(proc_prefix)s(i).eq.21)CC(i)=CA
         if(leg_pdgs_%(proc_prefix)s(i).ne.0.and.abs(leg_pdgs_%(proc_prefix)s(i)).le.6)CC(i)=CF
      enddo
      betaF=beta_FI
      gammaI=beta_IF
      if(ninit.eq.2)gammaI=beta_II
      A10g=A10(gammaI)
      A20a=A20(alpha)
      xa=x**(1d0+alpha)
      xb=x**(1d0+betaF)
c
c     scales, alpha_s and Born
      call set_mur_from_momenta(pp,nexternal,leg_pdgs_%(proc_prefix)s)
      alphas=alphas_current
      muf0=get_muf(pp,nexternal,leg_pdgs_%(proc_prefix)s)
      CALL ME_ACCESSOR_HOOK(pp,HEL,ALPHAS,ANS)
      BLO=ANS(0)
      Bab=0d0
      if(ninit.eq.2)Bab=GET_CCBLO(1,2)
c
c     central muF (iv = 0) and its variations (iv = 1..nmuf)
      do iv=0,nmuf
         muf=muf0
         if(iv.gt.0)muf=facf(iv)*muf0
         do lab=1,2
            if(.not.ISLOQCDPARTON(lab))cycle
c
c           final-state term of D1: sum of gamma^hc_k - 2 C_k A2(alpha)
c           over the final-state partons k whose hard-collinear counterterm
c           recoils against lab (iref1(k) = lab): its FI mapping rescales
c           only that incoming parton
            SD1=0d0
            do k=3,nexternal
               if(.not.ISLOQCDPARTON(k))cycle
               if(iref1(k).ne.lab)cycle
c              massive partons: no hard-collinear counterterm
               if(pmass(k).ne.0d0)cycle
               if(leg_pdgs_%(proc_prefix)s(k).eq.21)SD1=SD1+gamma_hc_g-2d0*CA*A20a
               if(abs(leg_pdgs_%(proc_prefix)s(k)).le.6)SD1=SD1+gamma_hc_q-2d0*CF*A20a
            enddo
c           physical beam of the incoming parton lab, and the other one
            ib=lab
            if(imirror.eq.2)ib=3-lab
            io=3-ib
            if(lpp(ib).eq.0)cycle
            eta=xbk(ib)
            etao=xbk(io)
            ir=iref1(lab)
            LBir=log(sB(lab,ir)/muR_current**2)
            LFir=log(sB(lab,ir)/muf**2)
c
            do i=1,n_init_states
               if(init_mirror(i).neqv.(imirror.eq.2))cycle
               flb=init_states(ib,i)
               flo=init_states(io,i)
               Fo=1d0
               if(lpp(io).ne.0)Fo=pdf_beam(io,flo,etao,muf)
               if(Fo.eq.0d0)cycle
c
c              test function of the diagonal terms (a' = flb)
               Td=0d0
               if(x.ge.eta)Td=pdf_beam(ib,flb,eta/x,muf)/x
               Td1=pdf_beam(ib,flb,eta,muf)
c
c              F1 + F2: P-bar and P^(lambda)_fin, summed over the PDF
c              flavours a' that give the Born flavour b
               sumF=0d0
               Tg=0d0
               if(x.ge.eta)Tg=pdf_beam(ib,21,eta/x,muf)/x
               Tg1=pdf_beam(ib,21,eta,muf)
               if(flb.eq.21)then
c                 a' = g: P_gg
                  res=2d0*CA*(x*Tg-Tg1)/(1d0-x)+2d0*CA*((1d0-x)/x+x*(1d0-x))*Tg+beta0/2d0*Tg1
                  sumF=sumF+LFir*res-lx*2d0*CA*(x/(1d0-x)+(1d0-x)/x+x*(1d0-x))*Tg
                  sumF=sumF+2d0*CA*((1d0-x)/x+x*(1d0-x))*(lam*l1x-A10g)*Tg
c                 a' = q, qbar: q -> g, C_F (1+(1-x)^2)/x
                  do iq=1,Nf
                     do ifl=-1,1,2
                        Tq=0d0
                        if(x.ge.eta)Tq=pdf_beam(ib,ifl*iq,eta/x,muf)/x
                        res=CF*(1d0+(1d0-x)**2)/x
                        sumF=sumF+(LFir-lx)*res*Tq
                        sumF=sumF+(res*(lam*l1x-A10g)+CF*x)*Tq
                     enddo
                  enddo
               else
c                 a' = b: q -> q, C_F ((1+x^2)/(1-x))_+
                  res=CF*(1d0+x**2)/(1d0-x)
                  sumF=sumF+LFir*res*(Td-Td1)-lx*res*Td
                  sumF=sumF+CF*(1d0-x)*(lam*l1x+1d0-A10g)*Td
c                 a' = g: g -> q, T_R (x^2+(1-x)^2)
                  res=TR*(x**2+(1d0-x)**2)
                  sumF=sumF+(LFir-lx)*res*Tg
                  sumF=sumF+(res*(lam*l1x-A10g)+TR*2d0*x*(1d0-x))*Tg
               endif
c
c              D1, D2: diagonal terms proportional to B
               sumD=-SD1*xb/(1d0-x)*(Td-Td1)
               sumD=sumD+2d0*CC(lab)*(lam*x*l1x/(1d0-x)*(Td-Td1)-(lam-1d0)*xa*l1x/(1d0-x)*(Td-Td1)-A10g*x/(1d0-x)*(Td-Td1)+xa/(1d0-x)*(Td-Td1)*(A10g-A20a-1d0-LBir)+lx*xa/(1d0-x)*Td)
c
               res=(sumF+sumD)*BLO
c              D3: B_ab (two incoming QCD partons), L_ab = LBir - ln x
               if(ninit.eq.2)res=res-2d0*(xa*l1x/(1d0-x)*(Td-Td1)+xa/(1d0-x)*(Td-Td1)*(A20a+1d0+LBir)-lx*xa/(1d0-x)*Td)*Bab
c              D4: final-state partons k, L_ik = L_ik,Born - ln x
               do k=3,nexternal
                  if(.not.ISLOQCDPARTON(k))cycle
                  ccb=GET_CCBLO(lab,k)
                  LBik=log(sB(lab,k)/muR_current**2)
                  res=res-2d0*(xa/(1d0-x)*(Td-Td1)*LBik-lx*xa/(1d0-x)*Td)*ccb
c                 massive k (alpha = 0): the soft integral J_s^FmI differs
c                 from the massless J_s,FI by -[x/(1-x) ln(1+rho/(1-x))]_+,
c                 rho = m_k^2/s with s = s_Born/x (MadNkLO massive paper,
c                 eq. IsMI; note massive_isr_phase0_map.md, Sec. 4.3)
                  if(pmass(k).ne.0d0)then
                     rhoB=pmass(k)**2/sB(lab,k)
                     res=res+2d0*x/(1d0-x)*(log(1d0+x*rhoB/(1d0-x))*Td-log(1d0+rhoB/(1d0-x))*Td1)*ccb
                  endif
               enddo
               conv(iv)=conv(iv)+Fo*res
            enddo
         enddo
         conv(iv)=conv(iv)*alphas/(2d0*pi)
c
         if(abs(conv(iv)).ge.huge(1d0).or.isnan(conv(iv)))then
            write(77,*)'Exception caught in int_convolution_NLO',conv(iv)
            ierr=1
            conv=0d0
            return
         endif
      enddo
c
      return
      end


      FUNCTION A10(W)
C     A10(w) = Psi0(w+1) + eulergamma
      IMPLICIT NONE
      DOUBLE PRECISION A10,W
C
      IF(W.NE.0D0.AND.W.NE.1D0.AND.W.NE.2D0.AND.W.NE.3D0.AND.W.NE.4D0.AND.W.NE.5D0)THEN                            
        WRITE(*,*)'Value not coded in A10',W
        STOP
      ENDIF
C
      IF(W.EQ.0D0)A10=0D0
      IF(W.EQ.1D0)A10=1D0
      IF(W.EQ.2D0)A10=3D0/2D0
      IF(W.EQ.3D0)A10=11D0/6D0
      IF(W.EQ.4D0)A10=25D0/12D0
      IF(W.EQ.5D0)A10=137/60D0
C
      RETURN
      END


      FUNCTION A20(W)
C     A20(w) = Psi0(w+2) - 1 + eulergamma
      IMPLICIT NONE
      DOUBLE PRECISION A20,W
C
      IF(W.NE.0D0.AND.W.NE.1D0.AND.W.NE.2D0.AND.W.NE.3D0.AND.W.NE.4D0.AND.W.NE.5D0)THEN                            
         WRITE(*,*)'Value not coded in A20',W
        STOP
      ENDIF
C
      IF(W.EQ.0D0)A20=0D0
      IF(W.EQ.1D0)A20=1D0/2D0
      IF(W.EQ.2D0)A20=5D0/6D0
      IF(W.EQ.3D0)A20=13D0/12D0
      IF(W.EQ.4D0)A20=77D0/60D0
      IF(W.EQ.5D0)A20=29D0/20D0
C
      RETURN
      END


      FUNCTION A21(W)
C     A21(w) = Psi1(w+2) + 1 - Zeta2
      IMPLICIT NONE
      DOUBLE PRECISION A21,W
C
      IF(W.NE.0D0.AND.W.NE.1D0.AND.W.NE.2D0.AND.W.NE.3D0.AND.W.NE.4D0.AND.W.NE.5D0)THEN                            
        WRITE(*,*)'Value not coded in A21',W
        STOP
      ENDIF
C
      IF(W.EQ.0D0)A21=0D0
      IF(W.EQ.1D0)A21=-1D0/4D0
      IF(W.EQ.2D0)A21=-13D0/36D0
      IF(W.EQ.3D0)A21=-61D0/144D0
      IF(W.EQ.4D0)A21=-1669D0/3600D0
      IF(W.EQ.5D0)A21=-1769D0/3600D0
C
      RETURN
      END



      subroutine NLO_I_MASS(s,v,mk2,ml2,mu,ccBLO,INLO_MASS)
c     Fully massive integrated soft counterterm (paper eq. IsMM).
c     New compact form using eta, eta_k, eta_l (SciPost 2503.14629).
c
c     Convention: INLO_MASS = -ccBLO * (mc_I_s^MM - pole_coeff*ln(s/mu^2))
c       pole_coeff = 1 + s/(2*sqrt(lam)) * ln(eta)  [= -(1/eps) coefficient]
c       mc_I_s^MM  = eq. (IsMM) in appendix, mu-independent finite part
c
c     The 1/eps pole itself is assembled outside this routine in the main
c     int_counter_NLO loop:
c       INLO(-1) += ccBLO*(-1/2)*(2 - (1/v)*ln((1+v)/(1-v)))
c                 = -ccBLO*(1 + s/(2*sqrtlam)*ln(eta))
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision s,v,mk2,ml2,mu,ccBLO,INLO_MASS
      double precision ddilog
      external ddilog
      double precision lam,sqrtlam,Q2
      double precision eta,eta_k,eta_l
      double precision ln_eta,ln_etak,ln_etal
      double precision pole_coeff,mc_finite,ln_mu2
c
c     -- kinematic variables -------------------------------------------------
      lam     = (s*v)**2d0            ! lambda = s^2 - 4 mk^2 ml^2
      sqrtlam = s*v                   ! sqrt(lambda)
      Q2      = s + mk2 + ml2
c
      eta   = (1d0-v)/(1d0+v)
      eta_k = (s+2d0*mk2-sqrtlam)/(s+2d0*mk2+sqrtlam)
      eta_l = (s+2d0*ml2-sqrtlam)/(s+2d0*ml2+sqrtlam)
c
      ln_eta  = dlog(eta)
      ln_etak = dlog(eta_k)
      ln_etal = dlog(eta_l)
      ln_mu2  = dlog(s/mu**2)
c
c     -- mu-independent finite part  mc_I_s^MM  (eq. IsMM) ------------------
c     (no continuation lines here: the template goes through MG5's
c     FortranWriter, which re-indents every line and splits long ones)
      mc_finite = 2d0*ddilog(eta) + (3d0/8d0)*ln_eta**2 + 0.5d0*dlog(Q2*s/(mk2*ml2))*ln_eta - 2d0*zeta2 - ddilog(1d0-eta_k) - ddilog(1d0-eta_l) - (1d0/8d0)*(ln_etak-ln_etal)**2
      mc_finite = (s/sqrtlam)*mc_finite - (mk2-ml2)/(2d0*sqrtlam)*(ln_etak-ln_etal) - Q2/(2d0*sqrtlam)*ln_eta + dlog(Q2*s/lam) + dlog(mk2*ml2/lam) + 4d0
c
c     -- coefficient of the 1/eps pole (used for the mu term) ----------------
      pole_coeff = 1d0 + s/(2d0*sqrtlam)*ln_eta
c
c     -- assemble ------------------------------------------------------------
      INLO_MASS = -ccBLO * (mc_finite - pole_coeff*ln_mu2)
c
      end


      double precision function IS_M(rho)
c     finite part of the massive-massless integrated soft counterterm,
c     rho=m^2/s (2503.14629, eq. IsMF)
      implicit none
      double precision rho,ddilog
      external ddilog
      IS_M = -ddilog(-rho) - log(rho)**2/4d0 + (1d0+rho)*log(1d0+rho) + (0.5d0-rho)*log(rho)
      return
      end


      double precision function IHC_0G(rho)
c     massive-recoiler hard-collinear finite part, g -> q qbar,
c     rho=m_r^2/s (2503.14629, eq. IhcM)
      implicit none
      double precision rho
      IHC_0G = rho + 1.5d0*rho*log(rho) - (1d0+1.5d0*rho)*log(1d0+rho) - 2d0*rho*sqrt(rho)*atan(sqrt(1d0+rho)-sqrt(rho))
      return
      end


      double precision function IHC_1G(rho)
c     massive-recoiler hard-collinear finite part, q -> q g
      implicit none
      double precision rho
      IHC_1G = rho*log(rho) - (1d0+rho)*log(1d0+rho)
      return
      end


      double precision function IHC_2G(rho)
c     massive-recoiler hard-collinear finite part, g -> g g
      implicit none
      double precision rho
      IHC_2G = -2d0*rho - log(1d0+rho) + 4d0*rho*sqrt(rho)*atan(sqrt(1d0+rho)-sqrt(rho))
      return
      end
