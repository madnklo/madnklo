      double precision function M2_CC_gluon_parents(i,j,k,l,kind1,kind2,ordered,xs,xp,wgt,xj,xjb,nit,extra,wgt_chan,ierr)
c     Two independent gluon parents: C_ijkl or (C_ij-S_i C_ij) C_ijkl.
c     Eq. (C.22) and the strongly ordered two-pair limits of 2212.11190v2.
      use sectors4_module
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'colored_partons.inc'
      include 'mapped_labels_common.inc'
      include 'leg_PDGs.inc'
      include 'nsqso_born.inc'
      include 'input.inc'
      include 'run.inc'
      integer i,j,k,l,kind1,kind2,nit,ierr,r,a,b,ib,jb,kb,lb,rb,parent1,parent2,legs(5)
      integer isec,jsec,ksec,lsec,iref,asec,bsec,csec,dsec
      common/cpartindices/isec,jsec,ksec,lsec,iref
      common/csecindices/asec,bsec,csec,dsec
      integer real_leg_pdgs(nexternal-1),Born_leg_pdgs(nexternal-2)
      common/c_NNLO_U_PDGs/real_leg_pdgs,Born_leg_pdgs
      integer %(proc_prefix_rr)s_fl_factor,%(proc_prefix_rr)s_den,%(proc_prefix_Born)s_den
      common/%(proc_prefix_rr)s_flavour_factor/%(proc_prefix_rr)s_fl_factor
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      common/%(proc_prefix_Born)s_iden/%(proc_prefix_Born)s_den
      logical ordered,docut,doplot,test_sector_function
      common/cdoplot/doplot
      common/ctestsecfun/test_sector_function
      double precision xs(nexternal,nexternal),xp(0:3,nexternal),pmass(nexternal)
      double precision xpb(0:3,nexternal-1),xsb(nexternal-1,nexternal-1)
      double precision xpbb(0:3,nexternal-2),xsbb(nexternal-2,nexternal-2)
      double precision wgt,xj,xjb,extra,wgt_chan,xjcs1,xjcs2,weight,wc,soft,norm
      double precision z1,z2,t1,t2,p1,p2,q1,q2,v1(0:3),v2(0:3)
      double precision alphas,ans(0:NSQSO_BORN),BLO,KKBLO1,KKBLO2,KKKKBLO,term,wgts(1)
## if (pyamplicol) {
      double precision %(proc_prefix_Born)s_GET_KKBLO,%(proc_prefix_Born)s_GET_KKKKBLO
## }
      include 'pmass.inc'
      M2_CC_gluon_parents=0d0
      ierr=0
## if (pyamplicol) {
      r=iref
      if(min(i,j,k,l,r).le.nincoming.or.max(i,j,k,l,r).gt.nexternal)then
        write(*,*) 'CC gluon parents: a final-state recoiler outside both pairs is required'
        stop 1
      endif
      legs=(/i,j,k,l,r/)
      do a=1,5
        do b=a+1,5
          if(legs(a).eq.legs(b))then
            write(*,*) 'CC gluon parents: repeated mapping legs'
            stop 1
          endif
        enddo
        if(pmass(legs(a)).ne.0d0)then
          write(*,*) 'CC gluon parents: massive mapping legs are unsupported'
          stop 1
        endif
      enddo
      call check_pair(i,j,kind1)
      call check_pair(k,l,kind2)
      do a=1,nexternal
        if(isnnloqcdparton(a).and.(a.le.nincoming.or.pmass(a).ne.0d0))then
          write(*,*) 'CC gluon parents: only massless final-state QCD is implemented'
          stop 1
        endif
      enddo
      call splitting(i,j,kind1,xs,xp,nexternal,p1,q1,z1,t1,v1)
      if(ierr.ne.0)goto 999
      call fill_real_mapped_labels(i,j,leg_pdgs,real_leg_pdgs)
      jb=real_labels(j,i,j)
      kb=real_labels(k,i,j)
      lb=real_labels(l,i,j)
      rb=real_labels(r,i,j)
      call phase_space_CS_inv(i,j,r,xp,xpb,nexternal,leg_pdgs,xjcs1,real_labels(:,i,j))
      if(xjcs1.eq.0d0)goto 999
      call invariants_from_p(xpb,nexternal-1,xsb,ierr)
      if(ierr.ne.0)goto 999
c     The common massless recoiler only rescales: pair-2 invariants and
c     momentum fractions agree before/after map 1. Use mapped inputs.
      ib=r
      r=rb
      call splitting(kb,lb,kind2,xsb,xpb,nexternal-1,p2,q2,z2,t2,v2)
      r=ib
      if(ierr.ne.0)goto 999
      call fill_born_mapped_labels(kb,lb,real_leg_pdgs,Born_leg_pdgs)
      call phase_space_CS_inv(kb,lb,rb,xpb,xpbb,nexternal-1,real_leg_pdgs,xjcs2,born_labels(:,kb,lb))
      if(xjcs2.eq.0d0)goto 999
      if(docut(xpbb,nexternal-2,Born_leg_pdgs,0))return
      call invariants_from_p(xpbb,nexternal-2,xsbb,ierr)
      if(ierr.ne.0)goto 999
      parent1=born_labels(jb,kb,lb)
      parent2=born_labels(lb,kb,lb)
      if(ordered)then
        call get_sig2(xs,alpha_mod,nexternal)
        call get_wc_nlo(i,j,k,r)
        wc=wc_nlo
        soft=0d0
        if(kind1.eq.1)soft=2d0*CA*(1d0-z1)/z1
        p1=p1*wc-soft
        q1=q1*wc
        call get_sig2(xsb,alpha_mod_bar,nexternal-1)
        call get_wcbar_nlo(real_labels(csec,i,j),real_labels(dsec,i,j),rb)
        weight=wcbar_nlo
      else
c       Eight orientations of the same unordered pair of collinear pairs.
c       The existing get_WCC_NNLO assumes a three-leg topology.
        norm=sector(i,j,k,l)+sector(j,i,k,l)+sector(i,j,l,k)+sector(j,i,l,k)
        norm=norm+sector(k,l,i,j)+sector(l,k,i,j)+sector(k,l,j,i)+sector(l,k,j,i)
        if(norm.le.0d0.or.norm.ge.huge(1d0))goto 999
        weight=sector(asec,bsec,csec,dsec)/norm
      endif
      if(test_sector_function)then
        M2_CC_gluon_parents=weight
        if(ordered)M2_CC_gluon_parents=(wc-merge(1d0,0d0,kind1.eq.1))*weight
        return
      endif
      call set_mur_from_momenta(xpbb,nexternal-2,Born_leg_pdgs)
      alphas=alphas_current
      call %(proc_prefix_Born)s_RESET_SPIN_CORRELATION_VECTORS()
      call %(proc_prefix_Born)s_ME_ACCESSOR_HOOK(xpbb,-1,alphas,ans)
      BLO=ans(0)
      KKBLO1=%(proc_prefix_Born)s_GET_KKBLO(xpbb,alphas,parent1,v1)
      KKBLO2=%(proc_prefix_Born)s_GET_KKBLO(xpbb,alphas,parent2,v2)
      KKKKBLO=%(proc_prefix_Born)s_GET_KKKKBLO(xpbb,alphas,parent1,v1,parent2,v2)
      term=(p1+q1)*(p2+q2)*BLO+2d0*q1*(p2+q2)*KKBLO1/t1
      term=term+2d0*q2*(p1+q1)*KKBLO2/t2+4d0*q1*q2*KKKKBLO/t1/t2
      term=term/xs(i,j)/xsb(kb,lb)*(8d0*pi*alphas)**2*weight*extra*xj
      term=term*%(proc_prefix_rr)s_fl_factor*dble(%(proc_prefix_Born)s_den)/dble(%(proc_prefix_rr)s_den)
      if(abs(term).ge.huge(1d0).or.isnan(term))goto 999
      M2_CC_gluon_parents=term
      if(doplot)then
        wgts(1)=-term*wgt/nit*wgt_chan
        call analysis_fill(xpbb,xsbb,nexternal-2,Born_leg_pdgs,wgts)
      endif
      return
 999  ierr=1
      M2_CC_gluon_parents=0d0
      return

      contains

      double precision function sector(a,b,c,d)
      integer a,b,c,d
      double precision ea,ec,wa,wb,wir,wkr
      ea=(xs(a,1)+xs(a,2))/xs(1,2)
      ec=(xs(c,1)+xs(c,2))/xs(1,2)
      wa=xs(a,b)/ea/(xs(b,1)+xs(b,2))
      wb=xs(c,d)/ec/(xs(d,1)+xs(d,2))
      wir=xs(a,r)/ea/(xs(r,1)+xs(r,2))
      wkr=xs(c,r)/ec/(xs(r,1)+xs(r,2))
      sector=(1d0/ea/wa)**alpha_mod/ec/wb/wir/wkr
      end function sector

      subroutine check_pair(a,b,kind)
      integer a,b,kind
      logical valid
      valid=.false.
      if(kind.eq.1)valid=leg_pdgs(a).eq.21.and.leg_pdgs(b).eq.21
      if(kind.eq.2)valid=abs(leg_pdgs(a)).ge.1.and.abs(leg_pdgs(a)).le.5.and.leg_pdgs(a).eq.-leg_pdgs(b)
      if(.not.valid)then
        write(*,*) 'CC gluon parents: invalid splitting flavours',leg_pdgs(a),leg_pdgs(b)
        stop 1
      endif
      end subroutine check_pair

      subroutine splitting(a,b,kind,s,p,n,Pab,Qab,z,t,v)
      integer a,b,kind,n
      double precision s(n,n),p(0:3,n),Pab,Qab,z,t,v(0:3),sab,sar,sbr
c     Check original flavours separately: stage-2 labels belong to pbar.
      sab=s(a,b)
      sar=s(a,r)
      sbr=s(b,r)
      if(min(sab,sar,sbr).le.0d0)then
        ierr=1
        return
      endif
      z=sar/(sar+sbr)
      t=-z*(1d0-z)*sab
      if(t.eq.0d0)then
        ierr=1
        return
      endif
      v=z*p(:,b)-(1d0-z)*p(:,a)-(2d0*z-1d0)*sab/(sar+sbr)*p(:,r)
      if(kind.eq.1)then
        Pab=2d0*CA*(z/(1d0-z)+(1d0-z)/z+z*(1d0-z))
        Qab=-2d0*CA*z*(1d0-z)
      else
        Pab=TR*(1d0-2d0*z*(1d0-z))
        Qab=2d0*TR*z*(1d0-z)
      endif
      end subroutine splitting
## } else {
      write(*,*) 'CC gluon parents require the pyAmpliCol double-spin provider'
      stop 1
## }
      end function M2_CC_gluon_parents
