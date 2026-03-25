
      
      double precision function M2_S_g(i,xs,xp,wgt,xj,xjB,nit,extra,wgt_chan,ierr)
c     single-soft limit S_(i) * Wsoft
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
      integer i,l,m,lb,mb,ierr,nit
      double precision pref,M2tmp,wgt,wgts(1),wgtpl,wgt_chan,xj,xjB,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO,extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,ml2,mm2,y,z,x,damp
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
      double precision  %(proc_prefix_S_g)s_GET_CCBLO
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_S_g)s_den
      common/%(proc_prefix_S_g)s_iden/%(proc_prefix_S_g)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
      double precision pmass(nexternal)
      INCLUDE 'pmass.inc'      
c
c     initialise
      M2_S_g=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
c     
c     checks
      if(leg_pdgs(i).ne.21)then
         write(*,*)'Wrong pdgs in M2_S_g',leg_pdgs(i)
         stop
      endif
      if(.not.(i.eq.isec))then
         write(*,*)'Wrong indices in M2_S_g',i,isec
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
               write(77,*)'Inaccuracy 1 in M2_S_g',sil,sim
               goto 999
            endif
c
c     call colour-connected Born
            call %(proc_prefix_S_g)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
            lb=mapped_labels(l)
            mb=mapped_labels(m)
            ccBLO = %(proc_prefix_S_g)s_GET_CCBLO(lb,mb)
c
c     eikonal
            M2TMP = SLM/(SIL*SIM) - ML2/SIL**2 - MM2/SIM**2
            M2TMP = CCBLO*M2TMP
c
c     damping factors
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
            M2tmp=M2tmp*damp*xj
            M2_S_g=M2_S_g+pref*M2tmp*WS_NLO*extra
c
c     plot
            wgtpl=-pref*M2tmp*WS_NLO*extra*wgt*wgt_chan
            wgtpl = wgtpl*dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_real)s_den)*%(proc_prefix_real)s_fl_factor
c     if(doplot)call histo_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgtpl)
            wgts=wgtpl
            if(doplot)call analysis_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgts)
c
         enddo 
      enddo
c
c     apply flavour factor
      M2_S_g = M2_S_g *dble(%(proc_prefix_S_g)s_den)/dble(%(proc_prefix_real)s_den)* %(proc_prefix_real)s_fl_factor
c
c     sanity check
      if(abs(M2_S_g).ge.huge(1d0).or.isnan(M2_S_g))then
         write(77,*)'Exception caught in M2_S_g',M2_S_g
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

