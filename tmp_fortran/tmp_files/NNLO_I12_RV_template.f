      subroutine int_counter_I12_NNLO_%(isec)d_%(jsec)d(p,sNLO,sLO,I12NNLO,ierr)
c     MSbar integrated counterterm
c     FINITE_PART = I12NNLO(0)
c     SINGLE_POLE = I12NNLO(-1)
c     DOUBLE_POLE = I12NNLO(-2)
      use sectors2_module
      implicit none
      INCLUDE 'nexternal.inc'
      INCLUDE 'mapped_labels_common.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'input.inc'
      INCLUDE 'virtual_recoilers.inc'
      INCLUDE 'leg_PDGs_%(proc_prefix_real)s.inc'
      INCLUDE 'colored_partons.inc'
      integer i,c,d,e,f,irefp,irefpp
      integer cb,db,rb,kb
      integer ierr
      double precision p(0:3,nexternal),pb(0:3,nexternal-1)
      double precision sNLO(nexternal,nexternal),sLO(nexternal-1,nexternal-1)
      double precision I12NNLO(-2:0),pref
      double precision Is(-2:0),Ic(-2:0),Isc(-2:0)
      double precision BLO,ccBLO,triBLO,quadBLO
      DOUBLE PRECISION ALPHAS,ANS(0:NSQSO_BORN)
      DOUBLE PRECISION ALPHA_QCD
      INTEGER, PARAMETER :: HEL = - 1
      double precision  %(proc_prefix_S_RV_g)s_GET_CCBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_TRIBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_QUADBLO
      integer iref1(nexternal)
      double precision res
      integer isec,jsec,iref
      common/csecindices/isec,jsec,iref
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/underlying_leg_PDGs
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      alphas=alpha_qcd(as,nloop,mu_r)
      pref=8d0*pi*alphas
      Is = 0d0
      Ic = 0d0
      Isc = 0d0
      I12NNLO = 0d0

      i = isec
      iref1 = 0
      BLO = 0d0
      CCBLO = 0d0
      TRIBLO = 0d0
      QUADBLO = 0d0
      res = 0d0
c
c     TODO: check
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
c
c     SOFT PART OF I^(12), eq. (4.47)
c
      do c=1,nexternal
        if(.not.isnloqcdparton(c)) cycle
        if(c.eq.i.or.m.eq.j) cycle
        do d=1,nexternal
          if(.not.isnloqcdparton(d)) cycle
          if(d.eq.i.or.d.eq.c) cycle
          call fill_born_mapped_labels(i,c,leg_pdgs,underlying_leg_pdgs)
          call fill_born_mapped_labels(i,d,leg_pdgs,underlying_leg_pdgs)
        enddo
      enddo
c
      do c=1,nexternal
        if(.not.isnloqcdparton(c)) cycle
        if(c.eq.i) cycle
        do d=1,nexternal
          if(.not.isnloqcdparton(d)) cycle
          if(d.eq.i.or.d.eq.c) cycle
c
          scd = snlo(c,d)
          sic = snlo(i,c)
          sid = snlo(i,d)
          if(sic*sid.le.0d0)then
            write(77,*)'Inaccuracy 1 in I12_RV', sic, sid
            goto 999
          endif
          Ei_cd = scd / sic / sid
c
          cb = born_labels(c,i,c)
          db = born_labels(d,i,c)
c
          call phase_space_CS_inv(i,c,d,p,pb,nexternal,leg_PDGs,xjCS1,born_labels(:,i,c))
          if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 10
          call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
          blo = ans(0)
          ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,db)
c
          Is = Is - Ei_cd * Ca*(Js(snlo(i,c))+Js(snlo(i,d))-Js(snlo(c,d)))*ccblo
          Is = Is - Ei_cd * Jhc(i,snlo(i,irefp))*ccblo
c
          do e=1,nexternal
            if(.not.isnloqcdparton(e)) cycle
            if(e.eq.i) cycle
            do f=1,nexternal
               if(.not.isnloqcdparton(f)) cycle
               if(f.eq.i.or.f.eq.e) cycle
c
               if(len(leg_pdgs_%(proc_prefix_real)s).le.5) then
                 quadblo = 2d0 * Cf**2 * blo
               elseif
                 write(77,*)'ee -> >2j is not implemented yet'
                 goto 999
                 !quadblo = %(proc_prefix_S_RV_g)s_get_quadblo(c,d,e,f)
               endif
               Is = Is + Ei_cd * Js(snlo(e,f))*quadblo/2d0
             enddo
c
             if(e.eq.i.or.e.eq.d) cycle
c
             if(len(leg_pdgs_%(proc_prefix_real)s).le.5) then
               quadblo = 2d0 * Cf**2 * blo
             elseif
               write(77,*)'ee -> >2j is not implemented yet'
               goto 999
               !quadblo = %(proc_prefix_S_RV_g)s_get_quadblo(c,d,e,d)
             endif
             Is = Is + Ei_cd * Js(snlo(d,e))*quadblo
c
 10          cb = born_labels(c,i,d)
             db = born_labels(d,i,d)
c
             call phase_space_CS_inv(i,d,c,p,pb,nexternal,leg_PDGs,xjCS1,born_labels(:,i,d))
             if(docut(pb,nexternal-1,underlying_leg_pdgs,0))cycle
             call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
             blo = ans(0)
             ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,db)
c
             if(len(leg_pdgs_%(proc_prefix_real)s).le.5) then
               quadblo = 2d0 * Cf**2 * blo
             elseif
               write(77,*)'ee -> >2j is not implemented yet'
               goto 999
               !quadblo = %(proc_prefix_S_RV_g)s_get_quadblo(c,d,e,d)
             endif
             Is = Is - Ei_cd * Js(snlo(d,e))*quadblo
          enddo
        enddo
      enddo
c
      do k=1,nexternal
        if(.not.isnloqcdparton(k)) cycle
        if(k.eq.i) cycle
          do c=1,nexternal
            if(.not.isnloqcdparton(c)) cycle
            if(c.eq.i.or.c.eq.k.or.c.eq.iref) cycle
            do d=1,nexternal
              if(.not.isnloqcdparton(d)) cycle
              if(d.eq.i.or.d.eq.c.or.d.eq.k.or.d.eq.iref) cycle
c
              scd = snlo(c,d)
              scr = snlo(c,iref)
              sck = snlo(c,k)
              sic = snlo(i,c)
              sid = snlo(i,d)
              sir = snlo(i,iref)
              sik = snlo(i,k)
              if(sic*sid*sir*sik.le.0d0)then
                 write(77,*)'Inaccuracy 1 in I12_RV', sic,sid,sir,sik
                 goto 999
              endif
              Ei_cd = scd / sic / sid
              Ei_cr = scr / sic / sir
              Ei_ck = sck / sic / sik
c
              cb = born_labels(c,i,c)
              db = born_labels(d,i,c)
              rb = born_labels(iref,i,c)
              kb = born_labels(k,i,c)
c
              call phase_space_CS_inv(i,c,d,p,pb,nexternal,leg_PDGs,xjCS1,born_labels(:,i,c))
              if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 11
              call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
              ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,db)
              Is = Is - Jhc(k,snlo(k,irefpp))*Ei_cd*ccblo
c
 11           call phase_space_CS_inv(i,c,iref,p,pb,nexternal,leg_PDGs,xjCS1,born_labels(:,i,c))
              if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 12
              call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
              ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,rb)
              Is = Is - Jhc(k,snlo(k,irefpp))*2d0*Ei_cr*ccblo
c
 12           call phase_space_CS_inv(i,c,k,p,pb,nexternal,leg_PDGs,xjCS1,born_labels(:,i,c))
              if(docut(pb,nexternal-1,underlying_leg_pdgs,0))cycle
              call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
              blo = ans(0)
              ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,kb)
              Is = Is - Jhc(k,snlo(k,irefpp))*2d0*Ei_ck*ccblo               
            enddo
         enddo
      enddo
c
      Ws_ij = 1d0 ! placeholder
      Is = pref * Is * Ws_ij
c
c     Torino to ML conversion factor (gamma[1-eps] -> exp[ eps eulergamma])
      I12NNLO(0)  = I12NNLO(0) + pi**2/12d0 * CF
      I12NNLO(-1) = I12NNLO(-1) + gamma_q
      I12NNLO(-2) = I12NNLO(-2) + CF
c
      if(abs(I12NNLO(0)).ge.huge(1d0).or.isnan(I12NNLO(0)))then
         write(77,*)'Exception caught in int_counter_I12_NNLO',I12NNLO(0)
         goto 999
      endif
c
c     COLLINEAR PART OF I^(12), eq. (4.48)
c
      call fill_born_mapped_labels(i,j,leg_pdgs,underlying_leg_pdgs)
      call fill_born_mapped_labels(i,iref,leg_pdgs,underlying_leg_pdgs)
      call fill_born_mapped_labels(j,iref,leg_pdgs,underlying_leg_pdgs)
c
      do c=1,nexternal
        if(.not.isnloqcdparton(c)) cycle
        if(c.eq.i.or.c.eq.j) cycle
        do d=1,nexternal
          if(.not.isnloqcdparton(d)) cycle
          if(d.eq.i.or.d.eq.j.or.d.eq.c) cycle
          write(77,*)'ee -> >2j is not implemented yet'
          goto 999
        enddo
c

      enddo

      return
 999  ierr=1
      return
      end

      FUNCTION Js(s)
C     eq. (E.2) of 2212.11190v2
      implicit none
      include 'math.inc'
      include 'coupl.inc'
      double precision s,lmu,js(-2:0)
      alphas = alpha_qcd(as,nloop,mu_r)
      lmu    = log(s/mu_r**2)
      Js(-2) = 1d0
      Js(-1) = 2d0 - lmu
      Js( 0) = 6d0 - 7d0/12d0*pi**2 - 2d0*lmu
      Js     = alphas/2d0/pi * Js
      return
      END

      FUNCTION Jhc(k, s)
C     eq. (E.9) of 2212.11190v2
      implicit none
      include 'math.inc'
      include 'coupl.inc'
      include 'leg_PDGs_%(proc_prefix_real)s.inc'
      integer k
      double precision s,lmu,js(-2:0)
      alphas    = alpha_qcd(as,nloop,mu_r)
      lmu       = log(s/mu_r**2)
      if(leg_pdgs_%(proc_prefix_real)s(k).eq.21) then
        jhc(-1) = gamma_hc_g
        jhc( 0) = phi_hc_g - gamma_hc_g*lmu
      elseif(leg_pdgs_%(proc_prefix_real)s(k).ne.0 .and.abs(leg_pdgs_%(proc_prefix_real)s(k)).le.6) then
        jhc(-1) = gamma_hc_q
        jhc( 0) = phi_hc_q - gamma_hc_q*lmu
      endif
      Js        = alphas/2d0/pi * Js
      return
      END
