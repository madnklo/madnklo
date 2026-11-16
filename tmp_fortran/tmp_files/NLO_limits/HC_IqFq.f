      double precision function M2_HC_IqFq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     Hard-collinear counterterm of the initial-state splitting
c     q(ib) -> q(ia) + g, the gluon entering the Born,
c     with ia final and ib initial (2209.09123, eq. (eq:HCijRdamped),
c     kernel P^hc of eq. (eq:hcPI), azimuthal kernel Q of eq. (eq:APkernels2b)):
c        HC = (1 - S_ia) C_(ia,ib) R,
c     with the sector function W_c = 1 (ib initial). The Born is the
c     one of the sector mapping: (ia,ib,ir) FII for ir initial,
c     (ia,ir,ib) FFI for ir final.
      use sectors2_module
      implicit none
      double precision fl_factor_npo
      external fl_factor_npo
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      include 'input.inc'
      include 'run.inc'
      include 'qcd_power.inc'
      integer ia,ib,ir,ierr,nit,parent_leg
      double precision pref,wgt,wgts(sv_maxwgt),wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,KKBLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sij,sir,sjr,x,v,z,xi,kt2,damp,eik,Phc,avg
      double precision ans(0:nsqso_born)
      integer,parameter :: hel = - 1
      double precision alphas,alpha_qcd
      double precision lumi_npo
      external lumi_npo
      integer nspincol
      external nspincol
      double precision %(proc_prefix_HC_IqFq)s_GET_KKBLO
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      logical docut
c     channel decomposition: chan_nocut = .true. while the pieces are
c     evaluated for the channel fractions, without the cuts on the
c     mapped Born (local_counter_chan)
      logical chan_nocut
      common/c_chan_nocut/chan_nocut
      integer %(proc_prefix_real)s_fl_factor
      common/%(proc_prefix_real)s_flavour_factor/%(proc_prefix_real)s_fl_factor
      integer %(proc_prefix_real)s_den
      common/%(proc_prefix_real)s_iden/%(proc_prefix_real)s_den
      integer %(proc_prefix_HC_IqFq)s_den
      common/%(proc_prefix_HC_IqFq)s_iden/%(proc_prefix_HC_IqFq)s_den
      integer isec,jsec,ksec,lsec,iref
      common/csecindices/isec,jsec,ksec,lsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/UNDERLYING_LEG_PDGS
      integer mapped_labels(nexternal)
      common/c_mapped_labels/mapped_labels
c
c     initialise
      M2_HC_IqFq=0d0
      ierr=0
c
c     checks
      if(.not.(abs(leg_pdgs(ib)).le.6.and.leg_pdgs(ia).eq.leg_pdgs(ib)))then
         write(*,*)'Wrong pdgs in M2_HC_IqFq',leg_pdgs(ia),leg_pdgs(ib)
         stop
      endif
      if(.not.(ia.eq.isec.and.ib.eq.jsec.and.ir.eq.iref.and.ia.gt.2.and.ib.le.2))then
         write(*,*)'Wrong indices in M2_HC_IqFq',ia,ib,ir,isec,jsec,iref
         stop
      endif
c
c     possible cuts
      if(.not.chan_nocut.and.docut(xpb,nexternal-1,underlying_leg_pdgs,0))return
c
c     overall kernel prefix
c     renormalisation scale for the mapped n-body kinematics
      call set_mur_from_momenta(xpb,nexternal-1,underlying_leg_pdgs)
      alphas=alphas_current
      pref=8d0*pi*alphas
c
c     invariant quantities (i = ia, j = ib, r = ir)
      sij=xs(ia,ib)
      sir=xs(ia,ir)
      sjr=xs(ib,ir)
c
c     variables of the sector mapping (2209.09123, App. B), and the
c     damping factor of the initial-state collinear counterterm
c     (gamma of 2209.09123: beta_II for ir initial, beta_IF for ir
c     final). eik is the eikonal factor I_jr = s_jr/(s_ij s_ir) times
c     (1-v) for ir initial, 1/(1-z) for ir final.
      if(ir.le.2)then
c        FII (ijr)
         if(sjr.le.0d0.or.sij+sir.le.0d0)goto 998
         x=(sjr-sij-sir)/sjr
         v=sij/(sij+sir)
         damp=(1d0-v)**beta_II
         eik=sjr/(sij*sir)*(1d0-v)
      else
c        FFI (irj)
         if(sij+sjr.le.0d0)goto 998
         x=(sij+sjr-sir)/(sij+sjr)
         z=sij/(sij+sjr)
         damp=(1d0-z)**beta_IF
         eik=sjr/(sij*sir)/(1d0-z)
      endif
c
c     safety check
      if(sij.le.0d0.or.sir.le.0d0.or.x.le.0d0.or.x.ge.1d0)goto 998
c
c     call Born
      call %(proc_prefix_HC_IqFq)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      BLO = ANS(0)
c
c     hard-collinear kernel P^hc(x)/x
c     with the azimuthal term Q(x) [ B + (d-2) kt.B.kt/kt^2 ], d = 4,
c     with the transverse vector of 2209.09123 (eq. (eq:remapped_ki))
c        kt = k_i - xi k_j - (s_ij/s_jr) k_r,  xi = s_ir/s_jr,
c     orthogonal to k_j and k_r, kt^2 = - xi s_ij
      Phc=CF*(2d0*(1d0-x)/x+x)
      M2_HC_IqFq=Phc/x*BLO/sij
      xi=sir/sjr
      kt(:)=xp(:,ia)-xi*xp(:,ib)-sij/sjr*xp(:,ir)
      kt2=-xi*sij
      parent_leg=mapped_labels(ib)
      KKBLO = %(proc_prefix_HC_IqFq)s_GET_KKBLO(parent_leg,xpb,kt)
      M2_HC_IqFq=M2_HC_IqFq+(-CF*2d0*(1d0-x)/x)/x*(BLO+2d0*KKBLO/kt2)/sij
c
c     include prefactors: the matrix elements are averaged over the
c     spins and colours of the initial states, and the splitting
c     kernels (2209.09123) relate averaged ones; the ratio of the
c     symmetry factors of the Born and of the real is the one of
c     their IDEN without these averages, which differ for the two
c     initial-state partons ib and its Born parent
      avg=dble(nspincol(leg_pdgs(ib)))/dble(nspincol(underlying_leg_pdgs(mapped_labels(ib))))
      M2_HC_IqFq = M2_HC_IqFq*damp*dble(%(proc_prefix_HC_IqFq)s_den)/dble(%(proc_prefix_real)s_den)*avg*fl_factor_npo(%(proc_prefix_real)s_fl_factor)*pref*xj*extra
c
c     parton luminosity (hadron beams), at the scale of the mapped
c     kinematics
      M2_HC_IqFq = M2_HC_IqFq*lumi_npo(xpb,nexternal-1,underlying_leg_pdgs)
c
c     plot
      wgtpl=-M2_HC_IqFq*wgt*wgt_chan
      call sv_wgts(wgtpl,0d0,born_qcd_power+1,wgts)
c     renormalisation-scale variations: contribution to the integrand
c     (-K), ~ alpha_s^(b+1) at the scale of the mapped kinematics
      call sv_add(-M2_HC_IqFq*wgt_chan,0d0,born_qcd_power+1)
      if(doplot.and.wgtpl.ne.0d0)call analysis_fill(xpb,xsb,nexternal-1,underlying_leg_pdgs,wgts)
c
c     sanity check
      if(abs(M2_HC_IqFq).ge.huge(1d0).or.isnan(M2_HC_IqFq))then
         write(77,*)'Exception caught in M2_HC_IqFq',M2_HC_IqFq
         goto 999
      endif
c
      return
 998  write(77,*)'Inaccuracy 1 in M2_HC_IqFq',sij,sir,sjr
 999  ierr=1
      return
      end
