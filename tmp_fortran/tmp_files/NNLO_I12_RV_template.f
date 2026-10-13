      subroutine int_counter_I12_NNLO_%(isec)d_%(jsec)d(p,sNLO,sLO,I12NNLO,ierr)
c     MSbar integrated counterterm
c     FINITE_PART = I12NNLO(0)
c     SINGLE_POLE = I12NNLO(-1)
c     DOUBLE_POLE = I12NNLO(-2)
      use sectors2_module
      implicit none
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'math.inc'
      INCLUDE 'damping_factors.inc'
      INCLUDE 'colored_partons.inc'
      INCLUDE 'leg_PDGs.inc'
      INCLUDE 'nsqso_born.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'
      INCLUDE 'virtual_recoilers.inc'
      INCLUDE 'rv_mapped_labels_common.inc'
      integer i,j,k,r,c,d,e,f,rp,rpp
      integer cb,db,rb,jb,kb
      integer ierr
      double precision p(0:3,nexternal),pb(0:3,nexternal-1),xjcs1
      double precision sNLO(nexternal,nexternal),sLO(nexternal-1,nexternal-1)
      double precision I12NNLO(-2:0),pref
      double precision Is(-2:0),Ic(-2:0),Isc(-2:0)
      double precision blo,ccblo,quadblo
      double precision Pgq,Pgq_sc,Ei_cd,Ei_ck,Ei_cr
      double precision C_l(nexternal),Cb_l
      double precision rho_ij_l,rhob_ij_l,rho_ji_l,rhob_ji_l
      double precision scd,scr,sck,sic,sid,sir,sik,sjr,sij,x
      double precision alphas,alpha_qcd,ans(0:nsqso_born)
      integer, parameter :: hel = - 1
      double precision  %(proc_prefix_S_RV_g)s_GET_CCBLO
      double precision  %(proc_prefix_S_RV_g)s_GET_QUADBLO
      integer iref1(nexternal)
      double precision res
      integer isec,jsec,ksec,lsec,iref0
      common/csecindices/isec,jsec,ksec,lsec,iref0
      integer underlying_leg_pdgs(nexternal-1)
      common/c_U_PDGs/underlying_leg_PDGs
      logical docut
      double precision JsIJ(-2:0),JsIR(-2:0),JsJR(-2:0)
      double precision JsIC(-2:0),JsID(-2:0),JsJC(-2:0),JsCD(-2:0)
      double precision JsDE(-2:0),JsEF(-2:0)
      double precision JhcIIR(-2:0),JhcIIRp(-2:0)
      double precision JhcJJR(-2:0),JhcKKRpp(-2:0)
      double precision pmass(nexternal)
      include 'pmass.inc'
c
c     initialise
      alphas=alpha_qcd(asmz,nloop,mu_r)
      pref=8d0*pi*alphas
c
c     TODO: check
      iref1 = 0
      do i=1,len_iref
         iref1(iref(1,i)) = iref(2,i)
      enddo
c     TODO: generalise r,rp,rpp
      r = jsec
      rp = iref1(isec)
      rpp = 0
c
      i = isec
      j = jsec
c
      Is = 0d0
      Ic = 0d0
      Isc = 0d0
      I12nnlo = 0d0
      blo = 0d0
      ccblo = 0d0
      quadblo = 0d0
c
c     SOFT PART OF I^(12), eq. (4.47)
c
      do c=1,nexternal
        if(.not.isnloqcdparton(c)) cycle
        if(c.eq.i.or.c.eq.j) cycle
        do d=1,nexternal
          if(.not.isnloqcdparton(d)) cycle
          if(d.eq.i.or.d.eq.c) cycle
          call fill_mapped_labels(i,c,leg_pdgs,underlying_leg_pdgs)
          call fill_mapped_labels(i,d,leg_pdgs,underlying_leg_pdgs)
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
          cb = mapped_labels(c,i,c)
          db = mapped_labels(d,i,c)
c
          call phase_space_CS_inv(i,c,d,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,c))
          if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 10
          call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
          blo = ans(0)
          ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,db)
c
          call get_Js(snlo(i,c), JsIC)
          call get_Js(snlo(i,d), JsID)
          call get_Js(snlo(c,d), JsCD)
          call get_Jhc(i, snlo(i,rp), JhcIIRp)
          Is = Is - Ei_cd * Ca * (JsIC+JsID-JsCD) * ccblo
          Is = Is - Ei_cd * JhcIIRp * ccblo
c
          do e=1,nexternal
            if(.not.isnloqcdparton(e)) cycle
            if(e.eq.i) cycle
            do f=1,nexternal
               if(.not.isnloqcdparton(f)) cycle
               if(f.eq.i.or.f.eq.e) cycle
c
               if(size(leg_pdgs).le.5) then
                 quadblo = 2d0 * Cf**2 * blo
               else
                 write(77,*)'ee -> >2j is not implemented yet'
                 goto 999
                 !quadblo = %(proc_prefix_S_RV_g)s_get_quadblo(c,d,e,f)
               endif
               call get_Js(snlo(e,f), JsEF)
               Is = Is + Ei_cd * JsEF * quadblo/2d0
             enddo
c
             if(e.eq.i.or.e.eq.d) cycle
c
             if(size(leg_pdgs).le.5) then
               quadblo = 2d0 * Cf**2 * blo
             else
               write(77,*)'ee -> >2j is not implemented yet'
               goto 999
               !quadblo = %(proc_prefix_S_RV_g)s_get_quadblo(c,d,e,d)
             endif
             call get_Js(snlo(d,e), JsDE)
             Is = Is + Ei_cd * JsDE * quadblo
c
 10          cb = mapped_labels(c,i,d)
             db = mapped_labels(d,i,d)
c
             call phase_space_CS_inv(i,d,c,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,d))
             if(docut(pb,nexternal-1,underlying_leg_pdgs,0))cycle
             call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
             blo = ans(0)
             ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,db)
c
             if(size(leg_pdgs).le.5) then
               quadblo = 2d0 * Cf**2 * blo
             else
               write(77,*)'ee -> >2j is not implemented yet'
               goto 999
               !quadblo = %(proc_prefix_S_RV_g)s_get_quadblo(c,d,e,d)
             endif
             call get_Js(snlo(d,e), JsDE)
             Is = Is - Ei_cd * JsDE * quadblo
          enddo
        enddo
      enddo
c
      sir = snlo(i,r)
      do k=1,nexternal
        if(.not.isnloqcdparton(k)) cycle
        if(k.eq.i) cycle
c
        sik = snlo(i,k)
c
          do c=1,nexternal
            if(.not.isnloqcdparton(c)) cycle
            if(c.eq.i.or.c.eq.k) cycle
c
            cb = mapped_labels(c,i,c)
            rb = mapped_labels(r,i,c)
            kb = mapped_labels(k,i,c)
c
            sck = snlo(c,k)
            scr = snlo(c,r)
            sic = snlo(i,c)
            if(sic*sir*sik.le.0d0)then
              write(77,*)'Inaccuracy 1 in I12_RV', sic,sir,sik
              goto 999
            endif
            Ei_cr = scr / sic / sir
            Ei_ck = sck / sic / sik
c
            do d=1,nexternal
              if(c.eq.r) cycle
              if(.not.isnloqcdparton(d)) cycle
              if(d.eq.i.or.d.eq.c.or.d.eq.k.or.d.eq.r) cycle
c
              db = mapped_labels(d,i,c)
c
              scd = snlo(c,d)
              sid = snlo(i,d)
              if(sic*sid.le.0d0)then
                 write(77,*)'Inaccuracy 1 in I12_RV', sic,sid
                 goto 999
              endif
              Ei_cd = scd / sic / sid
c
              call phase_space_CS_inv(i,c,d,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,c))
              if(docut(pb,nexternal-1,underlying_leg_pdgs,0))cycle
              call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
              ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,db)
              call get_Jhc(k, snlo(k,rpp), JhcKKRpp)
              Is = Is - JhcKKRpp * Ei_cd * ccblo
            enddo
c
c     TODO: needs rpp in initial state for 2jets
            if(size(leg_pdgs).le.5)cycle
            if(c.eq.r)goto 11
            call phase_space_CS_inv(i,c,r,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,c))
            if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 11
            call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
            ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,rb)
            call get_Jhc(k, snlo(k,rpp), JhcKKRpp)
            Is = Is - JhcKKRpp * 2d0 * Ei_cr * ccblo
c
c     TODO: needs rpp in initial state for 2jets
 11         if(size(leg_pdgs).le.5)cycle
            call phase_space_CS_inv(i,c,k,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,c))
            if(docut(pb,nexternal-1,underlying_leg_pdgs,0))cycle
            call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
            blo = ans(0)
            ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(cb,kb)
            call get_Jhc(k, snlo(k,rpp), JhcKKRpp)
            Is = Is - JhcKKRpp * 2d0 * Ei_ck * ccblo
         enddo
      enddo
c
      Is = pref * Is
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
c     COLLINEAR + SOFT-COLLINEAR PART OF I^(12), eq. (4.48-49), with i=g, j=q
c
      r = iref1(isec)
      C_l = 0d0
      do c=1,nexternal
        if(abs(leg_pdgs(c)).le.6)then
          C_l(c) = Cf
        elseif(abs(leg_pdgs(c)).eq.21)then
          C_l(c) = Ca
        endif
      enddo
      call fill_mapped_labels(i,j,leg_pdgs,underlying_leg_pdgs)
c
      jb = mapped_labels(j,i,j)
      if(abs(underlying_leg_pdgs(jb)).le.6)then
        Cb_l = Cf
        rho_ij_l = (Cb_l + C_l(i) - C_l(j)) / Cb_l
        rhob_ij_l = (Cb_l - C_l(i) - C_l(j)) / Cb_l
        rho_ji_l = (Cb_l + C_l(j) - C_l(i)) / Cb_l
        rhob_ji_l = (Cb_l - C_l(j) - C_l(i)) / Cb_l
      elseif(abs(underlying_leg_pdgs(jb)).eq.21)then
        Cb_l = Ca
        rho_ij_l = (Cb_l + C_l(i) - C_l(j)) / Cb_l
        rhob_ij_l = (Cb_l - C_l(i) - C_l(j)) / Cb_l
        rho_ji_l = (Cb_l + C_l(j) - C_l(i)) / Cb_l
        rhob_ji_l = (Cb_l - C_l(j) - C_l(i)) / Cb_l
      endif
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
        cb = mapped_labels(c,i,j)
c
        call phase_space_CS_inv(j,r,i,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
        if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 20
        call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
        blo = ans(0)
        ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(jb,cb)
        call get_Js(snlo(i,c),JsIC)
        Ic = Ic + rho_ij_l * JsIC * ccblo
        Isc = Isc + Ca * JsIC * ccblo
        if(size(leg_pdgs).le.5) then
          write(77,*)'ee -> >2j is not implemented yet'
          write(77,*)'cal B contribution for i=q, j=qx missing'
          goto 999
        endif
c
 20     call phase_space_CS_inv(i,r,j,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
        if(docut(pb,nexternal-1,underlying_leg_pdgs,0))cycle
        call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
        blo = ans(0)
        ccblo = %(proc_prefix_S_RV_g)s_get_ccblo(jb,cb)
        call get_Js(snlo(j,c), JsJC)
        Ic = Ic + rho_ji_l * JsJC * ccblo
        Isc = Isc + (2d0*C_l(j)-Ca) * JsJC * ccblo
        if(size(leg_pdgs).le.5) then
          write(77,*)'ee -> >2j is not implemented yet'
          write(77,*)'cal B contribution for i=q, j=qx missing'
          goto 999
        endif
      enddo
c
      call phase_space_CS_inv(i,j,r,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
      if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 21
      call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
      blo = ans(0)
      call get_Js(snlo(i,j), JsIJ)
      Ic = Ic + Cb_l * rhob_ij_l * JsIJ * blo
      Isc = Isc - C_l(j) * Ca * JsIJ * blo
      call get_Js(snlo(i,r), JsIR)
      Isc = Isc - Ca * C_l(j) * JsIR * blo
      call get_Js(snlo(j,r), JsJR)
      Isc = Isc - (2d0*C_l(j)-Ca) * C_l(j) * JsJR * blo
c
 21   call phase_space_CS_inv(j,r,i,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
      if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 22
      call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
      blo = ans(0)
      call get_Js(snlo(i,r), JsIR)
      Ic = Ic + Cb_l * rho_ij_l * JsIR * blo
      Isc = Isc + Ca * C_l(j) * JsIR * blo
c
 22   call phase_space_CS_inv(i,r,j,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
      if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 23
      call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
      blo = ans(0)
      call get_Js(snlo(j,r), JsJR)
      Ic = Ic + Cb_l * rho_ji_l * JsJR * blo
      Isc = Isc + (2d0*C_l(j)-Ca) * C_l(j) * JsJR * blo
      call get_Jhc(j, snlo(j,r), JhcJJR)
      Isc = Isc - JhcJJR * blo
c
 23   call phase_space_CS_inv(i,j,r,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
      if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 24
      call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
      blo = ans(0)
      call get_Js(snlo(i,r), JsIR)
      Ic = Ic - Cb_l * rho_ij_l * JsIR * blo
c
 24   call phase_space_CS_inv(j,i,r,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
      if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 25
      call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
      blo = ans(0)
      call get_Js(snlo(j,r), JsJR)
      Ic = Ic - Cb_l * rho_ji_l * JsJR * blo
c
 25   call phase_space_CS_inv(i,j,r,p,pb,nexternal,leg_PDGs,xjCS1,mapped_labels(:,i,j))
      if(docut(pb,nexternal-1,underlying_leg_pdgs,0))goto 999
      call %(proc_prefix_S_RV_g)s_me_accessor_hook(pb,hel,alphas,ans)
      blo = ans(0)
      call get_Jhc(i, snlo(i,r), JhcIIR)
      call get_Jhc(j, snlo(j,r), JhcJJR)
      Ic = Ic - (JhcIIR+JhcJJR) * blo
      Isc = Isc - JhcIIR * blo
c
      if(size(leg_pdgs).le.5) then
        write(77,*)'ee -> >2j is not implemented yet'
        write(77,*)'r prime contribution missing'
        goto 999
      endif
c
      sij=snlo(i,j)
      sir=snlo(i,r)
      sjr=snlo(j,r)
      x=sjr/(sir+sjr)
      Pgq = Cf*(2d0*x/(1d0-x) + (1d0-x))
      Pgq_sc = Cf*(2d0*x/(1d0-x))
c
      call get_wc_nlo(isec,jsec,r)
      Ic  = -pref * Pgq/sij * Ic * Wc_nlo
      Isc = -pref * Pgq_sc * Isc * 1d0
c
c     FINAL RESULT: eq. (4.46)
c
      call get_ws_nlo(isec,jsec)
      I12NNLO = Is*Ws_nlo + Ic - Isc

      return
 999  ierr=1
      return
      end

      SUBROUTINE GET_Js(s, Js)
C     eq. (E.2) of 2212.11190v2
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      double precision s,lmu,Js(-2:0)
      double precision alphas,alpha_qcd
      alphas = alpha_qcd(asmz,nloop,mu_r)
      lmu    = log(s/mu_r**2)
      Js(-2) = 1d0
      Js(-1) = 2d0 - lmu
      Js( 0) = 6d0 - 7d0/12d0*pi**2 - 2d0*lmu
      Js     = alphas/2d0/pi * Js
      return
      END

      SUBROUTINE GET_Jhc(k, s, Jhc)
C     eq. (E.9) of 2212.11190v2
      implicit none
      include 'coupl.inc'
      include 'math.inc'
      include 'leg_PDGs.inc'
      integer k
      double precision s,lmu,Jhc(-2:0)
      double precision alphas,alpha_qcd
      alphas    = alpha_qcd(asmz,nloop,mu_r)
      lmu       = log(s/mu_r**2)
      if(leg_pdgs(k).eq.21) then
        Jhc(-1) = gamma_hc_g
        Jhc( 0) = phi_hc_g - gamma_hc_g*lmu
      elseif(leg_pdgs(k).ne.0 .and.abs(leg_pdgs(k)).le.6) then
        Jhc(-1) = gamma_hc_q
        Jhc( 0) = phi_hc_q - gamma_hc_q*lmu
      endif
      Jhc       = alphas/2d0/pi * Jhc
      return
      END
