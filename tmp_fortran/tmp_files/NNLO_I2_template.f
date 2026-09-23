c     Implementation of arXiv:2212.11190v2, Eqs. (4.41)--(4.45).
c     Full MSbar Laurent series through epsilon**0; Appendix E
c     conventions.
c     I2NNLO(1:5) = finite, single, double, triple, quadruple pole.
c     Includes (alpha_s/(2*pi))**2. No extra Gamma conversion or damping.
c     Applicable to massless final-state QCD partons, colourless incoming
c     legs.
c     Requires GET_CCBLO(c,d)=B_cd and GET_QUADRUPLE_BLO(c,d,e,f)
c     = <{T_c.T_d,T_e.T_f}> (NOT half the anticommutator), also for
c     repeated leg labels. ME_ACCESSOR_HOOK must populate these accessors.
c     r_k comes from iref; r_jl is the smallest QCD leg different from
c     j,l.
c     A pair recoiler requires at least three QCD legs in this
c     realization.
c     Error codes: 1 invalid/NaN input or output; 2 unsupported coloured
c     masses/initial state/flavour or fewer than three QCD legs;
c     3 absent/invalid recoiler. Failure leaves all five outputs zero.
      subroutine INT_COUNTER_I2_NNLO(p,sLO,INLO,ierr)
      implicit none
      include 'nexternal.inc'
      include 'nsqso_born.inc'
      include 'coupl.inc'
      include 'math.inc'
      include 'input.inc'
      include 'virtual_recoilers.inc'
      include 'leg_PDGs_%(proc_prefix)s.inc'
      include 'colored_partons.inc'
      integer ierr,i,j,k,l,n,a,b,c,d,ip,refs(nexternal)
      integer legs(nexternal),flav(nexternal),pos(nexternal)
      integer,parameter :: hel=-1
      double precision p(0:3,nexternal),sLO(nexternal,nexternal)
      double precision INLO(5),alphas,alpha_qcd,pref,born
      double precision ans(0:nsqso_born),pmass(nexternal)
      double precision ss(nexternal,nexternal)
      double precision bc(nexternal,nexternal)
      double precision b4(nexternal,nexternal,nexternal,nexternal)
      double precision coeff(-4:0),mu2
      double precision get_ccblo,get_quadruple_blo
      external alpha_qcd,get_ccblo,get_quadruple_blo
      include 'pmass.inc'
      ierr=0
      I2NNLO=0d0
      n=0
      pos=0
      refs=0
      ss=0d0
      bc=0d0
      b4=0d0
      do i=1,nexternal
         if (.not.ISLOQCDPARTON(i)) cycle
         if (i.le.nincoming.or.pmass(i).ne.0d0) then
            ierr=2
            return
         endif
         ip=abs(leg_pdgs_%(proc_prefix)s(i))
         if (ip.ne.21.and.(ip.lt.1.or.ip.gt.6)) then
            ierr=2
            return
         endif
         n=n+1
         legs(n)=i
         pos(i)=n
         flav(n)=ip
      enddo
      if (n.lt.3) then
         ierr=2
         return
      endif
      do a=1,len_iref
         i=iref(1,a)
         j=iref(2,a)
         if (i.lt.1.or.i.gt.nexternal) then
            ierr=3
            return
         endif
         if (pos(i).eq.0) cycle
         if (j.lt.1.or.j.gt.nexternal) then
            ierr=3
            return
         endif
         if (pos(j).eq.0.or.i.eq.j) then
            ierr=3
            return
         endif
         if (refs(pos(i)).ne.0) then
            if (refs(pos(i)).ne.pos(j)) then
               ierr=3
               return
            endif
         endif
         refs(pos(i))=pos(j)
      enddo
      do a=1,n
         if (refs(a).eq.0) then
            ierr=3
            return
         endif
         do b=1,n
            if (a.eq.b) cycle
            ss(a,b)=sLO(legs(a),legs(b))
            if (.not.(ss(a,b).gt.0d0.and.
     &           ss(a,b).lt.huge(1d0))) then
               ierr=1
               return
            endif
         enddo
      enddo
      mu2=MU_R**2
      alphas=alpha_qcd(AS,NLOOP,MU_R)
      pref=(alphas/(2d0*pi))**2
      call ME_ACCESSOR_HOOK(p,hel,alphas,ans)
      born=ans(0)
      do a=1,n
         i=legs(a)
         do b=1,n
            if (b.eq.a) cycle
            j=legs(b)
            bc(a,b)=get_ccblo(i,j)
            do c=1,n
               k=legs(c)
               do d=1,n
                  if (d.eq.c) cycle
                  l=legs(d)
                  b4(a,b,c,d)=get_quadruple_blo(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      call I2NNLO_CORE(nexternal,n,flav,refs,ss,mu2,
     &     CA,CF,TR,dble(Nf),born,bc,b4,coeff,ierr)
      if (ierr.ne.0) return
      do ip=-4,0
         if (.not.(abs(pref*coeff(ip)).le.huge(1d0))) then
            ierr=1
            I2NNLO=0d0
            return
         endif
         I2NNLO(1-ip)=pref*coeff(ip)
      enddo
      end

c     Core returns epsilon-indexed coefficients without alpha_s prefactor.
c     All arrays use leading dimension ND; only the first N legs are used.
c     FLAV=21 for gluons; 1..6 for quarks or antiquarks (same integrals).
      subroutine I2NNLO_CORE(nd,n,flav,refs,s,mu2,
     &     ca,cf,tr,nf,born,bc,b4,res,ierr)
      implicit none
      integer nd,n,flav(nd),refs(nd),ierr
      integer c,d,e,f,k,r,j,l,ip,id
      double precision s(nd,nd),mu2,ca,cf,tr,nf,born
      double precision bc(nd,nd),b4(nd,nd,nd,nd),res(-4:0)
      double precision logs(nd,nd),tab(-4:2,30)
      double precision hc(-4:2),w(-4:2),a(-4:2),b(-4:2)
      double precision v(-4:2),sc
      ierr=0
      res=0d0
      if (n.lt.3.or.n.gt.nd) then
         ierr=2
         return
      endif
      if (.not.(mu2.gt.0d0.and.mu2.lt.huge(1d0))) then
         ierr=1
         return
      endif
      if (.not.(ca.gt.0d0.and.cf.gt.0d0.and.tr.gt.0d0
     &     .and.nf.ge.0d0)) then
         ierr=1
         return
      endif
      logs=0d0
      do c=1,n
         if (flav(c).ne.21.and.
     &        (abs(flav(c)).lt.1.or.abs(flav(c)).gt.6)) then
            ierr=2
            return
         endif
         if (refs(c).lt.1.or.refs(c).gt.n.or.refs(c).eq.c) then
            ierr=3
            return
         endif
         do d=1,n
            if (d.eq.c) cycle
            if (.not.(s(c,d).gt.0d0.and.s(c,d).lt.huge(1d0))) then
               ierr=1
               return
            endif
            logs(c,d)=log(s(c,d))-log(mu2)
         enddo
      enddo
      do id=1,30
         call I2NNLO_TABLE(id,ca,cf,tr,tab(-4,id))
      enddo
c     Double soft: Eq. (4.42). Every sum is ordered as in the paper.
      do c=1,n
         do d=1,n
            if (d.eq.c) cycle
            do e=1,n
               if (e.eq.c.or.e.eq.d) cycle
               do f=1,n
                  if (f.eq.c.or.f.eq.d.or.f.eq.e) cycle
                  sc=logs(c,d)+logs(e,f)
                  call I2NNLO_ADD(res,tab(-4,5),sc,b4(c,d,e,f)/4d0)
               enddo
               sc=logs(c,d)+logs(e,d)
               call I2NNLO_ADD(res,tab(-4,6),sc,b4(c,d,e,d))
            enddo
            sc=2d0*logs(c,d)
            call I2NNLO_ADD(res,tab(-4,7),sc,b4(c,d,c,d)/2d0)
            w=nf*tr*tab(:,8)-ca*tab(:,9)/2d0
            call I2NNLO_ADD(res,w,sc,bc(c,d))
         enddo
      enddo
c     Soft x hard collinear: Eqs. (4.43),(4.44).
      do k=1,n
         r=refs(k)
         if (flav(k).eq.21) then
            hc=nf*tab(:,2)+tab(:,4)/2d0
            w=ca*(2d0*nf*tab(:,27)+tab(:,30))
            a=2d0*nf*(tab(:,27)-tab(:,21))+tab(:,30)-tab(:,23)
            b=2d0*nf*(tab(:,24)-tab(:,21))+tab(:,26)-tab(:,23)
         else
            hc=tab(:,3)
            w=2d0*cf*tab(:,28)+ca*(tab(:,29)-tab(:,28))
            a=2d0*tab(:,28)+ca/cf*(tab(:,29)-tab(:,28))
     &           -2d0*tab(:,22)
            b=2d0*(tab(:,25)-tab(:,22))
         endif
         call I2NNLO_PRODUCT(hc,tab(-4,1),v)
         do c=1,n
            do d=1,n
               if (d.eq.c) cycle
               sc=logs(k,r)+logs(c,d)
               call I2NNLO_ADD(res,v,sc,-bc(c,d))
            enddo
         enddo
         sc=2d0*logs(k,r)
         call I2NNLO_ADD(res,w,sc,-born)
         call I2NNLO_ADD(res,a,sc,-bc(k,r))
         do c=1,n
            if (c.eq.k.or.c.eq.r) cycle
            sc=logs(k,r)+logs(k,c)
            call I2NNLO_ADD(res,b,sc,-bc(k,c))
            sc=logs(k,r)+logs(c,r)
            call I2NNLO_ADD(res,b,sc,-bc(c,r))
         enddo
c     Three-parton hard double collinear, between (4.44) and (4.45).
         if (flav(k).eq.21) then
            w=nf*tab(:,12)+tab(:,14)/6d0
         else
            w=nf*tab(:,10)+(tab(:,11)+tab(:,13))/2d0
         endif
         call I2NNLO_ADD(res,w,2d0*logs(k,r),born)
      enddo
c     Four-parton hard double collinear: Eq. (4.45).
c     Symmetric reference rule R_2(j,l): smallest allowed leg, Eq. (A.14).
      do j=1,n
         do l=1,n
            if (l.eq.j) cycle
            do r=1,n
               if (r.ne.j.and.r.ne.l) exit
            enddo
            if (flav(j).eq.21.and.flav(l).eq.21) then
               w=nf**2*tab(:,15)+nf*tab(:,17)+tab(:,20)/4d0
            elseif (flav(j).eq.21.or.flav(l).eq.21) then
               w=nf*tab(:,16)+tab(:,19)/2d0
            else
               w=tab(:,18)
            endif
            sc=logs(j,r)+logs(l,r)
            call I2NNLO_ADD(res,w,sc,born/2d0)
         enddo
      enddo
      do ip=-4,0
         if (.not.(abs(res(ip)).le.huge(1d0))) then
            ierr=1
            res=0d0
            return
         endif
      enddo
      end

c     Multiply Laurent polynomials before truncating: positive orders of
c     the NLO integrals are necessary in Eq. (4.43), even for the finite
c     part.
      subroutine I2NNLO_PRODUCT(a,b,c)
      implicit none
      double precision a(-4:2),b(-4:2),c(-4:2)
      integer i,j
      c=0d0
      do i=-4,2
         do j=-4,2
            if (i+j.lt.-4.or.i+j.gt.0) cycle
            c(i+j)=c(i+j)+a(i)*b(j)
         enddo
      enddo
      end

c     Add WEIGHT * exp(-epsilon*L) * C, expanded through epsilon**0.
c     L=log(s/mu^2)+log(s'/mu^2), or L=2*log(s/mu^2).
      subroutine I2NNLO_ADD(res,c,l,weight)
      implicit none
      double precision res(-4:0),c(-4:2),l,weight,term
      integer i,j
      do i=-4,0
         term=weight*c(i)
         res(i)=res(i)+term
         do j=i+1,0
            term=term*(-l)/dble(j-i)
            res(j)=res(j)+term
         enddo
      enddo
      end
c     Constituent Laurent coefficients of Appendix E, exact v2 edition.
c     Common alpha_s factors and exp(-epsilon*L) are applied by the
c     caller.
      subroutine I2NNLO_TABLE(id,ca,cf,tr,c)
      implicit none
      integer id
      double precision ca,cf,tr,c(-4:2),p2,p4,z3
      p2=acos(-1d0)**2
      p4=p2**2
      z3=1.2020569031595942854d0
      c=0d0
      select case(id)
c     E.2: J_s
      case(1)
         c(-2) = 1d0
         c(-1) = 2d0
         c(0) = 6d0 - 7d0 * p2 / 12d0
         c(1) = 18d0 - 7d0 * p2 / 6d0 - 25d0 * z3 / 3d0
         c(2) = 54d0 - 7d0 * p2 / 2d0 - 50d0 * z3 / 3d0 - 71d0 * p4 /
     &        1440d0
c     E.8: J_hc^(0g)
      case(2)
         c(-1) =  - 2d0 * tr / 3d0
         c(0) =  - 16d0 * tr / 9d0
         c(1) =  - tr *  ( 140d0 / 27d0 - 7d0 * p2 / 18d0 )
         c(2) =  - tr *  ( 1252d0 / 81d0 - 28d0 * p2 / 27d0 - 50d0 * z3 /
     &        9d0 )
c     E.8: J_hc^(1g)
      case(3)
         c(-1) =  - cf / 2d0
         c(0) =  - cf
         c(1) =  - cf *  ( 3d0 - 7d0 * p2 / 24d0 )
         c(2) =  - cf *  ( 9d0 - 7d0 * p2 / 12d0 - 25d0 * z3 / 6d0 )
c     E.8: J_hc^(2g)
      case(4)
         c(-1) =  - ca / 3d0
         c(0) =  - 8d0 * ca / 9d0
         c(1) =  - ca *  ( 70d0 / 27d0 - 7d0 * p2 / 36d0 )
         c(2) =  - ca *  ( 626d0 / 81d0 - 14d0 * p2 / 27d0 - 25d0 * z3 /
     &        9d0 )
c     E.4: J_sxs^(4)
      case(5)
         c(-4) = 1d0
         c(-3) = 4d0
         c(-2) = 16d0 - 7d0 * p2 / 6d0
         c(-1) = 60d0 - 14d0 * p2 / 3d0 - 50d0 * z3 / 3d0
         c(0) = 216d0 - 56d0 * p2 / 3d0 - 200d0 * z3 / 3d0 + 29d0 * p4 /
     &        120d0
c     E.4: J_sxs^(3)
      case(6)
         c(-4) = 1d0
         c(-3) = 4d0
         c(-2) = 17d0 - 4d0 * p2 / 3d0
         c(-1) = 70d0 - 16d0 * p2 / 3d0 - 68d0 * z3 / 3d0
         c(0) = 284d0 - 68d0 * p2 / 3d0 - 272d0 * z3 / 3d0 + 13d0 * p4 /
     &        90d0
c     E.4: J_sxs^(2)
      case(7)
         c(-4) = 1d0
         c(-3) = 4d0
         c(-2) = 18d0 - 3d0 * p2 / 2d0
         c(-1) = 76d0 - 6d0 * p2 - 74d0 * z3 / 3d0
         c(0) = 312d0 - 27d0 * p2 - 308d0 * z3 / 3d0 + 49d0 * p4 / 120d0
c     E.4: J_ss^(qqbar)
      case(8)
         c(-3) = 1d0 / 6d0
         c(-2) = 17d0 / 18d0
         c(-1) = 116d0 / 27d0 - 7d0 * p2 / 36d0
         c(0) = 1474d0 / 81d0 - 131d0 * p2 / 108d0 - 19d0 * z3 / 9d0
c     E.4: J_ss^(gg)
      case(9)
         c(-4) = 1d0 / 2d0
         c(-3) = 35d0 / 12d0
         c(-2) = 487d0 / 36d0 - 2d0 * p2 / 3d0
         c(-1) = 1562d0 / 27d0 - 269d0 * p2 / 72d0 - 77d0 * z3 / 6d0
         c(0) = 19351d0 / 81d0 - 3829d0 * p2 / 216d0 - 1025d0 * z3 / 18d0 -
     &        23d0 * p4 / 240d0
c     E.11: J_hcc^(0g)
      case(10)
         c(-2) = cf * tr / 6d0
         c(-1) = cf * tr *  ( 13d0 / 36d0 + p2 / 9d0 )
         c(0) = cf * tr *  (  - 55d0 / 216d0 + 17d0 * p2 / 18d0 + 14d0 * z3
     &        / 3d0 )
c     E.11: J_hcc^(0g,id)
      case(11)
         c(-1) = cf *  ( 2d0 * cf - ca )  *  ( 13d0 / 8d0 - p2 / 4d0 + z3 )
         c(0) = cf *  ( 2d0 * cf - ca )  *  (  - 227d0 / 16d0 + p2 + 17d0 *
     &        z3 / 2d0 - 11d0 * p4 / 120d0 )
c     E.11: J_hcc^(1g)
      case(12)
         c(-3) =  - 2d0 * cf * tr / 3d0 - ca * tr
         c(-2) =  - 31d0 * cf * tr / 9d0 - 89d0 * ca * tr / 18d0
         c(-1) =  - cf * tr *  ( 899d0 / 54d0 - p2 )  - ca * tr *  ( 1211d0
     &        / 54d0 - 3d0 * p2 / 2d0 )
         c(0) = cf * tr *  (  - 23833d0 / 324d0 + 31d0 * p2 / 6d0 + 160d0 *
     &        z3 / 9d0 )  + ca * tr *  (  - 2620d0 / 27d0 + 89d0 * p2 / 12d0 +
     &        80d0 * z3 / 9d0 )
c     E.11: J_hcc^(2g)
      case(13)
         c(-3) =  - 2d0 * cf ** 2 - cf * ca / 2d0
         c(-2) =  - 37d0 * cf ** 2 / 4d0 + 23d0 * cf * ca / 12d0
         c(-1) =  - cf ** 2 *  ( 307d0 / 8d0 - 3d0 * p2 + 4d0 * z3 )  - cf
     &        * ca *  ( 241d0 / 36d0 - p2 / 18d0 - 4d0 * z3 )
         c(0) = cf ** 2 *  (  - 2361d0 / 16d0 + 111d0 * p2 / 8d0 + 136d0 *
     &        z3 / 3d0 - 7d0 * p4 / 3d0 )  + cf * ca *  (  - 4609d0 / 216d0 +
     &        53d0 * p2 / 216d0 - 47d0 * z3 / 6d0 + 7d0 * p4 / 20d0 )
c     E.11: J_hcc^(3g)
      case(14)
         c(-3) =  - 5d0 * ca ** 2 / 2d0
         c(-2) =  - 77d0 * ca ** 2 / 6d0
         c(-1) =  - ca ** 2 *  ( 48d0 - 11d0 * p2 / 4d0 + 3d0 * z3 )
         c(0) = ca ** 2 *  (  - 16943d0 / 108d0 + 61d0 * p2 / 4d0 + 56d0 *
     &        z3 / 3d0 - 9d0 * p4 / 40d0 )
c     E.13: J_hcxhc^(qqqq)
      case(15)
         c(-2) = 4d0 * tr ** 2 / 9d0
         c(-1) = 64d0 * tr ** 2 / 27d0
         c(0) = tr ** 2 *  ( 284d0 / 27d0 - 16d0 * p2 / 27d0 )
c     E.13: J_hcxhc^(qqqg)
      case(16)
         c(-2) = tr * cf / 3d0
         c(-1) = 14d0 * tr * cf / 9d0
         c(0) = tr * cf *  ( 181d0 / 27d0 - 4d0 * p2 / 9d0 )
c     E.13: J_hcxhc^(qqgg)
      case(17)
         c(-2) = 2d0 * tr * ca / 9d0
         c(-1) = 32d0 * tr * ca / 27d0
         c(0) = tr * ca *  ( 142d0 / 27d0 - 8d0 * p2 / 27d0 )
c     E.13: J_hcxhc^(qgqg)
      case(18)
         c(-2) = cf ** 2 / 4d0
         c(-1) = cf ** 2
         c(0) = cf ** 2 *  ( 17d0 / 4d0 - p2 / 3d0 )
c     E.13: J_hcxhc^(qggg)
      case(19)
         c(-2) = ca * cf / 6d0
         c(-1) = 7d0 * ca * cf / 9d0
         c(0) = ca * cf *  ( 181d0 / 54d0 - 2d0 * p2 / 9d0 )
c     E.13: J_hcxhc^(gggg)
      case(20)
         c(-2) = ca ** 2 / 9d0
         c(-1) = 16d0 * ca ** 2 / 27d0
         c(0) = ca ** 2 *  ( 71d0 / 27d0 - 4d0 * p2 / 27d0 )
c     E.15: J_sxhc^4(1g)
      case(21)
         c(-3) =  - 2d0 * tr / 3d0
         c(-2) =  - 28d0 * tr / 9d0
         c(-1) =  - tr *  ( 344d0 / 27d0 - 7d0 * p2 / 9d0 )
         c(0) = tr *  (  - 3928d0 / 81d0 + 98d0 * p2 / 27d0 + 100d0 * z3 /
     &        9d0 )
c     E.15: J_sxhc^4(2g)
      case(22)
         c(-3) =  - cf / 2d0
         c(-2) =  - 2d0 * cf
         c(-1) =  - cf *  ( 8d0 - 7d0 * p2 / 12d0 )
         c(0) = cf *  (  - 30d0 + 7d0 * p2 / 3d0 + 25d0 * z3 / 3d0 )
c     E.15: J_sxhc^4(3g)
      case(23)
         c(-3) =  - ca / 3d0
         c(-2) =  - 14d0 * ca / 9d0
         c(-1) =  - ca *  ( 172d0 / 27d0 - 7d0 * p2 / 18d0 )
         c(0) = ca *  (  - 1964d0 / 81d0 + 49d0 * p2 / 27d0 + 50d0 * z3 /
     &        9d0 )
c     E.15: J_sxhc^3(1g)
      case(24)
         c(-3) =  - 2d0 * tr / 3d0
         c(-2) =  - 28d0 * tr / 9d0
         c(-1) =  - tr *  ( 362d0 / 27d0 - 8d0 * p2 / 9d0 )
         c(0) = tr *  (  - 4504d0 / 81d0 + 112d0 * p2 / 27d0 + 136d0 * z3 /
     &        9d0 )
c     E.15: J_sxhc^3(2g)
      case(25)
         c(-3) =  - cf / 2d0
         c(-2) =  - 2d0 * cf
         c(-1) =  - cf *  ( 17d0 / 2d0 - 2d0 * p2 / 3d0 )
         c(0) = cf *  (  - 35d0 + 8d0 * p2 / 3d0 + 34d0 * z3 / 3d0 )
c     E.15: J_sxhc^3(3g)
      case(26)
         c(-3) =  - ca / 3d0
         c(-2) =  - 14d0 * ca / 9d0
         c(-1) =  - ca *  ( 181d0 / 27d0 - 4d0 * p2 / 9d0 )
         c(0) = ca *  (  - 2252d0 / 81d0 + 56d0 * p2 / 27d0 + 68d0 * z3 /
     &        9d0 )
c     E.15: J_sxhc^(gqq)
      case(27)
         c(-3) =  - 2d0 * tr / 3d0
         c(-2) =  - 28d0 * tr / 9d0
         c(-1) =  - tr *  ( 344d0 / 27d0 - 17d0 * p2 / 18d0 )
         c(0) = tr *  (  - 4225d0 / 81d0 + 128d0 * p2 / 27d0 + 139d0 * z3 /
     &        9d0 )
c     E.15: J_sxhc^(gqg)
      case(28)
         c(-3) =  - cf / 2d0
         c(-2) =  - 2d0 * cf
         c(-1) =  - cf *  ( 9d0 - 5d0 * p2 / 6d0 )
         c(0) = cf *  (  - 38d0 + 19d0 * p2 / 6d0 + 101d0 * z3 / 6d0 )
c     E.15: J_sxhc^(ggq)
      case(29)
         c(-3) =  - cf / 2d0
         c(-2) =  - 2d0 * cf
         c(-1) =  - cf *  ( 8d0 - 2d0 * p2 / 3d0 )
         c(0) = cf *  (  - 32d0 + 17d0 * p2 / 6d0 + 59d0 * z3 / 6d0 )
c     E.15: J_sxhc^(ggg)
      case(30)
         c(-3) =  - ca / 3d0
         c(-2) =  - 14d0 * ca / 9d0
         c(-1) =  - ca *  ( 199d0 / 27d0 - 5d0 * p2 / 9d0 )
         c(0) = ca *  (  - 2477d0 / 81d0 + 119d0 * p2 / 54d0 + 101d0 * z3 /
     &        9d0 )
      end select
      end
