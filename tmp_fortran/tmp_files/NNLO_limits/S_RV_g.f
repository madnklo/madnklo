

      function M2_S_RV_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     single-soft limit S_(i) * Wsoft for RV
c     it returns 0 if i is not a gluon
      use sectors2_module
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      double precision m2_s_rv_g(-2:0)
      integer i,l,m,q,lb,mb,qb,ierr,nit,idum
      double precision pref,M2tmp(-2:0),wgt,wgtpl,wgt_chan,xj,xjB,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO,triBLO,VLO(-2:0),ccVLO(-2:0),extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,ml2,mm2,siq,smq,y,z,x,damp
      double precision eik0,eik1(-2:0),eik2(-2:0)
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
c     external
      integer get_color_dipole_index
      external get_color_dipole_index
      double precision alphas,ans(0:NSQSO_BORN)
      double precision alpha_qcd
      double precision alphaZ
      parameter(alphaZ=1d0)
      integer, parameter :: HEL = - 1
      double precision  %(proc_prefix_S_RV_g)s_GET_CCBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_CCVLO
      double precision  %(proc_prefix_S_RV_g)s_GET_TRIBLO
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_S_RV_g)s_den
      common/%(proc_prefix_S_RV_g)s_iden/%(proc_prefix_S_RV_g)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      integer mapped_flavours(nexternal-1),mapped_indices_shuff(nexternal-1)
      common/c_mapped_quantities_s/mapped_labels,mapped_flavours,mapped_indices_shuff
      double precision xpb_to_ME(0:3,nexternal-1)
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INCLUDE 'pmass.inc'

c
c     initialise
      M2_S_RV_g=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      idum=0
      xpb_to_ME=0d0
c
c     checks
      if(leg_pdgs(i).ne.21)then
         write(*,*)'Wrong pdgs in M2_S_RV_g',leg_pdgs(i)
         stop
      endif
      if(.not.(i.eq.isec))then
         write(*,*)'Wrong indices in M2_S_RV_g',i,isec
         stop
      endif
c
c     call W soft
      CALL GET_SIG2(XS,ALPHAZ,NEXTERNAL)
      CALL GET_WS_NLO(ISEC,JSEC)
c
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref=-8d0*pi*alphas
c
c     eikonal double sum
      do m=1,nexternal
         if(.not.isNLOQCDparton(m))cycle
         if(m.eq.i)cycle
         do l=1,nexternal
            if(.not.isNLOQCDparton(l))cycle
            if(l.eq.i)cycle
            if(l.eq.m)cycle
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the single-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            IF(DOCUT(XPB,NEXTERNAL-1,MAPPED_FLAVOURS,0))CYCLE
c
c     invariant quantities
            sil=xs(i,l)
            sim=xs(i,m)
            slm=xs(l,m)
            ml2=pmass(l)**2
            mm2=pmass(m)**2
c
c     safety check
            if(sil*sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S_RV_g',sil,sim
               goto 999
            endif
c
c     call colour-connected Born and Virtual
            lb=mapped_indices_shuff(mapped_labels(l))
            mb=mapped_indices_shuff(mapped_labels(m))
            XPB_TO_ME(0:3,MAPPED_INDICES_SHUFF(:))=XPB(0:3,:)
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S_RV_g)s_GET_CCBLO(lb,mb)
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
            ccVLO = %(proc_prefix_S_RV_g)s_GET_CCVLO(lb,mb)
c
c     eikonals
            EIK0     =  SLM/(SIL*SIM) - ML2/SIL**2 - MM2/SIM**2
            EIK1(-2) =  CA*EIK0
            EIK1(-1) = -CA*EIK0*log(sil*sim/slm/scale**2)
            EIK1( 0) =  CA*EIK0/2d0*(log(sil*sim/slm/scale**2)**2-5d0*zeta2)

            M2TMP(-2:0) = M2TMP(-2:0) + CCVLO(-2:0)*EIK0
            M2TMP(-2:0) = M2TMP(-2:0) - alphas/2d0/pi*CCBLO*EIK1(-2:0)
            M2TMP(-1)   = M2TMP(-1) - alphas*beta0/4d0/pi*CCBLO*EIK0

            do q=1,nexternal
               if(.not.isNLOQCDparton(q))cycle
               if(q.eq.i.or.q.eq.l.or.q.eq.m)cycle
c
c     invariant quantities
               siq=xs(i,q)
               smq=xs(m,q)
c
c     call triple-colour-connected Born
               qb=mapped_indices_shuff(mapped_labels(q))
               call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_to_ME,hel,alphas,ANS)
               TRIBLO = %(proc_prefix_S_RV_g)s_GET_TRIBLO(lb,mb,qb)
c
c     eikonals
               EIK2(-2) = 0d0
               EIK2(-1) = EIK0
               EIK2( 0) = -EIK0*log(sim*siq/smq/scale**2)

               M2TMP(-2:0) = M2TMP(-2:0) + alphas*TRIBLO*EIK2(-2:0)
            enddo


c     Including correct multiplicity factor
            M2tmp(-2:0) = M2tmp(-2:0)*dble(%(proc_prefix_S_RV_g)s_den)/dble(%(proc_prefix_real)s_den)
c
c     damping factors
c     TODO: adapt damping factors
            if(m.gt.2.and.l.gt.2)then
               y=sil/(sil+sim+slm)
               z=sim/(sim+slm)
               damp=((1d0-y)*(1d0-z))**alpha
            elseif(m.gt.2.and.l.le.2)then
               z=sim/(sim+slm)
               x=1d0 - sil/(sim+slm)
               damp=((1d0-z)*x)**alpha
            elseif(m.le.2.and.l.le.2)then
               x=1d0 - (sil+sim)/slm
               damp=x**alpha
            endif
            M2tmp(-2:0)=M2tmp(-2:0)*damp*xj
            M2_S_RV_g(-2:0)=M2_S_RV_g(-2:0)+pref*M2tmp(-2:0)*WS_NLO*extra
c
c     plot
            wgtpl=-pref*M2tmp(0)*WS_NLO*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*%(proc_prefix_real)s_fl_factor
            if(doplot)call histo_fill(xpb,xsb,nexternal-1,mapped_flavours,wgtpl)
c
         enddo
      enddo
c
c     apply flavour factor
      M2_S_RV_g(-2:0) = M2_S_RV_g(-2:0) * %(proc_prefix_real)s_fl_factor
c
c     sanity check
      if(abs(M2_S_RV_g(0)).ge.huge(1d0).or.isnan(M2_S_RV_g(0)))then
         write(77,*)'Exception caught in finite part of M2_S_RV_g',M2_S_RV_g(0)
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
