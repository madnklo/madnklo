

      SUBROUTINE SUB_M2_S_RV_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr,ret)
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
      double precision ret(-2:0)
      integer i,l,m,q,lb,mb,qb,ierr,nit
      double precision pref,M2tmp(-2:0),wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO,triBLO,VLO(-2:0),ccVLO(-2:0),extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,ml2,mm2,siq,smq,y,z,x,damp
      double precision eik0,eik1(-2:0),eik2(-2:0)
      double precision res_delta
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
      common/c_mapped_labels/mapped_labels
      DOUBLE PRECISION PMASS(NEXTERNAL)
      INCLUDE 'pmass.inc'
c
c     initialise
      ret=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
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
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,leg_PDGs,xjCS,mapped_labels)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            if(docut(xpb,nexternal-1,underlying_leg_pdgs,0))cycle
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
            lb=mapped_labels(l)
            mb=mapped_labels(m)
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            ccBLO = %(proc_prefix_S_RV_g)s_GET_CCBLO(lb,mb)
            ANS=0d0
c            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            ccVLO = (/0d0,0d0,0d0/) !%(proc_prefix_S_RV_g)s_GET_CCVLO(lb,mb)
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
               qb=mapped_labels(q)
c               call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
               TRIBLO = 0d0 !%(proc_prefix_S_RV_g)s_GET_TRIBLO(lb,mb,qb)
c
c     eikonals
               EIK2(-2) = 0d0
               EIK2(-1) = EIK0
               EIK2( 0) = -EIK0*log(sim*siq/smq/scale**2)

               M2TMP(-2:0) = M2TMP(-2:0) + alphas*TRIBLO*EIK2(-2:0)
            enddo
c
c     damping factors; TODO: adapt
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
            ret(-2:0)=ret(-2:0)+pref*M2tmp(-2:0)*WS_NLO*extra
c
c     plot
            wgtpl=-pref*M2tmp(0)*WS_NLO*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*dble(%(proc_prefix_S_RV_g)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor
c     if(doplot)call histo_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgtpl)
            wgts=wgtpl
            if(doplot)call analysis_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgts)
c
         enddo
      enddo
c
c     apply flavour factor
      ret(-2:0) = ret(-2:0) * %(proc_prefix_real)s_fl_factor
c
c     sanity check
      if(abs(ret(0)).ge.huge(1d0).or.isnan(ret(0)))then
         write(77,*)'Exception caught in finite part of M2_S_RV_g',ret(0)
         goto 999
      endif
c
c     add delta_s_rv_g (all prefactors included)
      call DELTA_S_RV_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr,res_delta)
      ret = ret + res_delta
c
      return
 999  ierr=1
      return
      end


      subroutine DELTA_S_RV_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr,res_delta)
c     Delta single-soft limit S_(i) * Wsoft for RV
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
      double precision res_delta(-2:0)
      integer i,k,ierr,nit
      double precision pref,M2tmp(-2:0),wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO,triBLO,quadBLO(-2:0),extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,ml2,mm2,y,z,x,damp
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
      integer, parameter :: HEL = - 1
      double precision cl
      double precision  %(proc_prefix_S_RV_g)s_GET_CCBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_TRIBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_QUADBLO
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_S_RV_g)s_den
      common/%(proc_prefix_S_RV_g)s_iden/%(proc_prefix_S_RV_g)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      double precision delta_s(-2:0)
c     Label conventions according to Eq.(5.19) in 2212.11190
c     with (c,d) ---> (l,m)
      integer l,m,p,q,t,lb,mb,pb,qb,kb,rb,tb
      double precision ccBLO_lm,ccBLO_ml
      double precision ccBLO_rl,ccBLO_lr
      double precision ccBLO_kr,ccBLO_rk
      double precision quadblo_pmlm,quadblo_pqlm,quadblo_tmlm
      double precision xpb_to_ME_lm(0:3,nexternal-1),xpb_to_ME_ml(0:3,nexternal-1)
      double precision xpb_to_ME_kr(0:3,nexternal-1),xpb_to_ME_rk(0:3,nexternal-1)
      double precision xpb_lm(0:3,nexternal-1),xpb_ml(0:3,nexternal-1)
      double precision xpb_rk(0:3,nexternal-1),xpb_kr(0:3,nexternal-1)
      double precision xsb_lm(nexternal-1,nexternal-1),xsb_ml(nexternal-1,nexternal-1)
      double precision xsb_kr(nexternal-1,nexternal-1),xsb_rk(nexternal-1,nexternal-1)
      double precision siq,smq,spm,spq,stm,sbpm,sbpq,sbtm,sik,sir,skr
      double precision mk2,mr2,gamma_l
      DOUBLE PRECISION PMASS(NEXTERNAL)
      double precision M2TMP_KR,EIK_KR
      integer mapped_labels_ilm(nexternal), mapped_labels_iml(nexternal)
      INCLUDE 'pmass.inc'
c
c     initialise
      res_delta = 0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      delta_s=0d0
      ccBLO_lm = 0d0
      ccBLO_ml = 0d0
      ccBLO_rl = 0d0
      ccBLO_lr = 0d0
      xpb_lm = 0d0
      xpb_ml = 0d0
      xpb_kr = 0d0
      xpb_rk = 0d0
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
c     overall kernel prefix
      ALPHAS=ALPHA_QCD(ASMZ,NLOOP,SCALE)
      pref=-8d0*pi*alphas
c
c     eikonal double sum
      do l=1,nexternal
         if(.not.isNLOQCDparton(l))cycle
         if(l.eq.i)cycle
         do m=1,nexternal
            if(.not.isNLOQCDparton(m))cycle
            if(m.eq.i)cycle
            if(m.eq.l)cycle
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the single-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
C     Build B_lm^{(ilm)} and B_lm^{(iml)}
c     The structure is such that we have
c     Eik(xs)*(B_lm^{(ilm)}*theta(ilm)-B_lm^{(iml)}*theta(iml))
c     Take care of this in applying kinematical cuts over Born kinematics
c     Build  B_lm^{(ilm)}
            call phase_space_CS_inv(i,l,m,xp,xpb_lm,nexternal,leg_PDGs,xjCS,mapped_labels)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb_lm,nexternal-1,xsb_lm,ierr)
            if(ierr.eq.1)goto 999
c     Build  B_lm^{(iml)}
 777        call phase_space_CS_inv(i,m,l,xp,xpb_ml,nexternal,leg_PDGs,xjCS,mapped_labels)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb_ml,nexternal-1,xsb_ml,ierr)
            if(ierr.eq.1)goto 999

c     invariant quantities
            sil=xs(i,l)
            sim=xs(i,m)
            slm=xs(l,m)
            ml2=pmass(l)**2
            mm2=pmass(m)**2
c     eikonal
            EIK0 =  SLM/(SIL*SIM) - ML2/SIL**2 - MM2/SIM**2
c
c     safety check
            if(sil*sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S_RV_g',sil,sim
               goto 999
            endif
            lb=mapped_labels(l)
            mb=mapped_labels(m)
            if(abs(leg_pdgs(l)).le.6) then
               gamma_l = gamma_q
               Cl = CF
            elseif(leg_pdgs(l).eq.21) then
               gamma_l = gamma_g
               Cl = CA
            else
               write(*,*) 'delta_S_RV_g:'
               write(*,*) 'Error in evaluating gamma_c,C_c'
               write(*,*) 'c, leg_pdgs(c) = ', l, leg_pdgs(l)
               write(*,*) 'Exit...'
               stop
            endif
c
c     call colour-connected B^{(ilm)} and B^{(iml)}
            ANS = 0d0
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_lm,hel,alphas,ANS)
            ccBLO_lm =%(proc_prefix_S_RV_g)s_GET_CCBLO(lb,mb)
            if(docut(xpb_lm,nexternal-1,underlying_leg_pdgs,0)) goto 778
c
            delta_s(-2) = delta_s(-2) + EIK0*2d0*Cl*ccBLO_lm
            delta_s(-1) = delta_s(-1) + EIK0*ccBLO_lm*(4d0*Cl+gamma_l)
c     Sum over e
            do t=1,nexternal
               if(.not.(isNLOQCDparton(t))) cycle
               if(t.eq.l) cycle
               if(t.eq.m) cycle
               tb = mapped_labels(t)
               stm = xs(t,m)
               sbtm = xsb(tb,mb)
               ANS = 0d0
               call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_lm,hel,alphas,ANS)
               QUADBLO_tmlm= 0d0 !%(proc_prefix_S_RV_g)s_GET_QUADBLO(tb,mb,lb,mb)
               delta_s(-1) = delta_s(-1)-EIK0*dlog(stm/sbtm)*QUADBLO_tmlm
               delta_s(0) = delta_s(0) + eik0*1d0/2d0*dlog(stm/sbtm)**2*QUADBLO_tmlm
            enddo
c
c
c     (c d e f) ---> (l m p q)
            do p=1,nexternal
               if(.not.isNLOQCDparton(p))cycle
               if(p.eq.i) cycle
               if(p.eq.l) cycle
               do q=1,nexternal
                  if(.not.isNLOQCDparton(q))cycle
                  if(q.eq.i) cycle
                  if(q.eq.l) cycle
                  if(q.eq.p) cycle
                  pb = mapped_labels(p)
                  qb = mapped_labels(q)
c     call invariants
                  spm = xs(p,m)
                  spq = xs(p,q)
                  sbpm = xsb_lm(pb,mb)
                  sbpq = xsb_lm(pb,qb)
c     call quadruple-colour-connected Born
                  call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_lm,hel,alphas,ANS)
                  QUADBLO_pqlm = 0d0 !%(proc_prefix_S_RV_g)s_GET_QUADBLO(pb,qb,lb,mb)
c
                  delta_s(-1) = delta_s(-1)-eik0*1d0/2d0*dlog(spq/sbpq)*QUADBLO_pqlm
                  delta_s(0) = delta_s(0) + eik0*1d0/4d0*dlog(spq/sbpq)**2*QUADBLO_pqlm
c
c     damping factors; TODO: adapt
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
            M2TMP(-2:0) = M2TMP(-2:0) + alphas/2d0/pi*delta_s(-2:0)
            M2tmp(-2:0)=M2tmp(-2:0)*damp*xj
            res_delta(-2:0)=res_delta(-2:0)+pref*M2tmp(-2:0)*WS_NLO*extra
c
c     plot
            wgtpl=-pref*M2tmp(0)*WS_NLO*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*dble(%(proc_prefix_S_RV_g)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor
c     if(doplot)call histo_fill(xpb_lm,xsb_lm,nexternal-1,underlying_leg_pdgs,wgtpl)
            wgts=wgtpl
            if(doplot)call analysis_fill(xpb_lm,xsb_lm,nexternal-1,underlying_leg_pdgs,wgts)
c
c     close q
               enddo
c     close p
            enddo

 778        ANS = 0d0
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_ml,hel,alphas,ANS)
            ccBLO_ml = %(proc_prefix_S_RV_g)s_GET_CCBLO(mb,lb)
            if(docut(xpb_lm,nexternal-1,underlying_leg_pdgs,0)) cycle
c
            delta_s(-2) = delta_s(-2) - EIK0*2d0*Cl*ccBLO_ml
            delta_s(-1) = delta_s(-1) - EIK0*ccBLO_ml*(4d0*Cl+gamma_l)
c
            M2TMP(-2:0) = M2TMP(-2:0) + alphas/2d0/pi*delta_s(-2:0)
            M2tmp(-2:0)=M2tmp(-2:0)*damp*xj
            res_delta(-2:0)=res_delta(-2:0)+pref*M2tmp(-2:0)*WS_NLO*extra
c
c     plot
            wgtpl=-pref*M2tmp(0)*WS_NLO*extra*wgt/nit*wgt_chan
            wgtpl = wgtpl*dble(%(proc_prefix_S_RV_g)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor
c     if(doplot)call histo_fill(xpb_lm,xsb_lm,nexternal-1,underlying_leg_pdgs,wgtpl)
            wgts=wgtpl
            if(doplot)call analysis_fill(xpb_lm,xsb_lm,nexternal-1,underlying_leg_pdgs,wgts)
c     close m
         enddo
c     Sum over (k,c) ---> (l,k)
         do k=1,nexternal
            if(.not.isNLOQCDparton(k))cycle
            if(k.eq.i) cycle
            if(k.eq.l) cycle
            if(k.eq.iref) cycle
c     Build  B_lm^{(irk)}
            call phase_space_CS_inv(i,iref,k,xp,xpb_rk,nexternal,leg_PDGs,xjCS,mapped_labels)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb_rk,nexternal-1,xsb_rk,ierr)
            if(ierr.eq.1)goto 999
c     Build  B_lm^{(ikr)}
            call phase_space_CS_inv(i,k,iref,xp,xpb_kr,nexternal,leg_PDGs,xjCS,mapped_labels)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb_kr,nexternal-1,xsb_kr,ierr)
            if(ierr.eq.1)goto 999
c     possible cuts
c     invariant quantities
            sik=xs(i,k)
            sir=xs(i,iref)
            skr=xs(k,iref)
            mk2=pmass(k)**2
            mr2=pmass(iref)**2
c     safety check
            if(sik*sir.le.0d0)then
               write(77,*)'Inaccuracy 2 in M2_S_RV_g',sik,sir
               goto 999
            endif
            kb=mapped_labels(k)
            rb=mapped_labels(iref)
c     call colour-connected B^{(icr)} and B^{(irc)}
            ANS = 0d0
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_kr,hel,alphas,ANS)
            ccBLO_kr =%(proc_prefix_S_RV_g)s_GET_CCBLO(kb,rb)
c
            ANS = 0d0
            call %(proc_prefix_S_RV_g)s_ME_ACCESSOR_HOOK(xpb_rk,hel,alphas,ANS)
            ccBLO_rk = %(proc_prefix_S_RV_g)s_GET_CCBLO(rb,kb)
c     eikonals
            EIK_KR =  SKR/(SIK*SIR) - MK2/SIK**2 - MR2/SIR**2
            if(docut(xpb_kr,nexternal-1,underlying_leg_pdgs,0)) goto 779
            M2TMP_KR = M2TMP_KR+ EIK_KR*gamma_l*ccBLO_kr
 779        if(docut(xpb_rK,nexternal-1,underlying_leg_pdgs,0)) cycle
            M2TMP_KR = M2TMP_KR - EIK_KR*gamma_l*ccBLO_rk

            M2tmp(-1) = M2tmp(-1) + alphas/2d0/pi*pref*WS_NLO*extra*M2TMP_kr
            M2TMP_KR = M2TMP_KR*damp*xj
            res_delta(-1) = res_delta(-1) + M2tmp(-1)
c     close k
         enddo
c     close l
      enddo

c     apply flavour factor
      res_delta(-2:0) = res_delta(-2:0) * %(proc_prefix_real)s_fl_factor
c
c     sanity check
      if(abs(res_delta(0)).ge.huge(1d0).or.isnan(res_delta(0)))then
         write(77,*)'Exception caught in finite part of DELTA_S_RV_g',res_delta(0)
         goto 999
      endif
c
 999  ierr=1
      return
      end
