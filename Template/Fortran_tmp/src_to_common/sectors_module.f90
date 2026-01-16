module sectors2_module
  implicit none
  integer, public :: n_ext
  double precision, public :: alpha_mod, W_NLO, WS_NLO, WC_NLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:), public :: sig2
  public :: get_sig2, get_W_NLO, get_WS_NLO, get_WC_NLO
  private

contains

  subroutine get_sig2(xs_in,alpha_in,n_ext_in)
    implicit none
    ! global
    integer :: n_ext_in
    double precision :: alpha_in
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    ! local
    integer :: i,j
    double precision :: ei,ej,wij
    ! set global module variables
    n_ext=n_ext_in
    if (.not.allocated(xs_mod)) allocate(xs_mod(n_ext,n_ext))
    if (.not.allocated(sig2)) allocate(sig2(3:n_ext,3:n_ext))
    xs_mod=xs_in
    alpha_mod=alpha_in
    ! calculate 2-index sigma
    sig2=0d0
    do i=3,n_ext
       do j=3,n_ext
          if(i.eq.j)cycle
          if( (xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2))*xs_mod(i,j)*xs_mod(1,2).ne.0d0 )then
             ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
             ej=(xs_mod(j,1)+xs_mod(j,2))/xs_mod(1,2)
             wij=xs_mod(1,2)*xs_mod(i,j)/(xs_mod(i,1)+xs_mod(i,2))/(xs_mod(j,1)+xs_mod(j,2))
             sig2(i,j)=(1d0/ei/wij)**alpha_mod
          endif
       enddo
    enddo
  end subroutine get_sig2


  subroutine get_W_NLO(i1,i2)
    !     NLO sector functions W(i1,i2)
    implicit none
    include 'all_sector_list.inc'
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    call sector2_global_checks(i1,i2)
    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       sigma = sigma + sig2(a,b)
    enddo
    W_NLO = num/sigma
    call sector2_sanity_checks(sigma,W_NLO)
  end subroutine get_W_NLO

  subroutine get_WS_NLO(i1,i2)
    !     NLO soft sector functions WS(i1,i2) = barS_i1 W(i1,i2)
    implicit none
    include 'all_sector_list.inc'
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    call sector2_global_checks(i1,i2)
    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       if(a.eq.i1) sigma = sigma + sig2(a,b)
       !if(b.eq.i1) sigma = sigma + sig2(b,a)
    enddo
    WS_NLO = num/sigma
    call sector2_sanity_checks(sigma,WS_NLO)
  end subroutine get_WS_NLO

  subroutine get_WC_NLO(xs_in,ia,ib,ir,alphaz,n_ext_in)
    !     NLO collinear sector functions WC(ia,ib,ir)
    implicit none
    include 'all_sector_list.inc'
    integer :: ia,ib,ir
    integer :: n_ext_in
    double precision :: ei,ej,wij,wir,wjr,alphaz
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    ei=(xs_in(ia,1)+xs_in(ia,2))/xs_in(1,2)
    ej=(xs_in(ib,1)+xs_in(ib,2))/xs_in(1,2)
    wij=xs_in(1,2)*xs_in(ia,ib)/(xs_in(ia,1)+xs_in(ia,2))/(xs_in(ib,1)+xs_in(ib,2))
    wir=xs_in(1,2)*xs_in(ia,ir)/(xs_in(ia,1)+xs_in(ia,2))/(xs_in(ir,1)+xs_in(ir,2))
    wjr=xs_in(1,2)*xs_in(ib,ir)/(xs_in(ib,1)+xs_in(ib,2))/(xs_in(ir,1)+xs_in(ir,2))
    wc_nlo=(ej*wjr)**alphaz/((ei*wir)**alphaz+(ej*wjr)**alphaz)
  end subroutine get_WC_NLO



  subroutine sector2_global_checks(i1,i2)
    implicit none
    integer :: i1,i2
    if(alpha_mod.lt.1d0)then
       write(77,*)'Wrong alpha_mod in sectors2',alpha_mod
       stop
    endif
    if(i1.le.2.or.i2.le.2) then
       write(77,*) 'sectors2: indices must be in final state',i1,i2
       stop
    endif
  end subroutine sector2_global_checks


  subroutine sector2_sanity_checks(sigma,Z)
    implicit none
    double precision :: Z,sigma
    if(sigma.le.0d0)then
       write(*,*)'Wrong sigma ',sigma
       stop
    endif
    if(abs(Z).ge.huge(1d0).or.isnan(Z))then
       write(77,*)'Exception caught ',Z
       stop
    endif
  end subroutine sector2_sanity_checks


end module sectors2_module






module sectors4_module
  implicit none
  integer, public :: n_ext,num_sec
  double precision, public :: alpha_mod, Z_NNLO, ZSS_NNLO, Z_HC_NNLO, ZS_NNLO, WCC_NNLO, WSS_NNLO, WSS_CC_NNLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:), public :: sig2
  double precision, allocatable, dimension(:,:,:,:), public :: sigNNLO
  public :: get_sigNNLO, get_Z_NNLO, get_ZHC_NNLO, get_ZS_NNLO, get_WCC_NNLO, get_WSS_NNLO, get_WSS_CC_NNLO
  private

contains

  subroutine get_sigNNLO(xs_in,alpha_in,n_ext_in)
    implicit none
    ! global
    include 'nexternal.inc'
    integer :: n_ext_in
    double precision :: alpha_in
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    ! local
    integer :: i,j,k,l,del_jk
    double precision :: ei,ej,ek
    double precision :: wij, wik, wjk, wkl

    ! set global module variables
    n_ext=n_ext_in
!    num_sec=(n_ext-2)*(n_ext-3)/2
    num_sec=(nexternal-2)*(nexternal-3)/2
    if (.not.allocated(xs_mod)) allocate(xs_mod(n_ext,n_ext))
    if (.not.allocated(sig2)) allocate(sig2(3:n_ext,3:n_ext))
    if (.not.allocated(sigNNLO)) allocate(sigNNLO(3:n_ext,3:n_ext,3:n_ext,3:n_ext))
    xs_mod=xs_in
    alpha_mod=alpha_in
    ! calculate 2-index and 4-index sigma
    sig2=0d0
    sigNNLO=0d0
    do i=3,n_ext
       do j=3,n_ext
          if(i.eq.j)cycle
          if( (xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2))*xs_mod(i,j)*xs_mod(1,2).ne.0d0 ) then
             ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
             wij=xs_mod(1,2)*xs_mod(i,j)/&
                  (xs_mod(i,1)+xs_mod(i,2))/(xs_mod(j,1)+xs_mod(j,2))
             sig2(i,j)=(1d0/ei/wij)**alpha_mod
          endif
          do k=3,n_ext
             if(k.eq.i)cycle
             del_jk=0
             if(j.eq.k) del_jk = 1
             do l=3,n_ext
                if(l.eq.i.or.l.eq.k)cycle
                if( (xs_mod(k,1)+xs_mod(k,2))*(xs_mod(l,1)+xs_mod(l,2))*xs_mod(k,l)*xs_mod(1,2).ne.0d0 )then
                   ek=(xs_mod(k,1)+xs_mod(k,2))/xs_mod(1,2)
                   wkl=xs_mod(1,2)*xs_mod(k,l)/&
                        (xs_mod(k,1)+xs_mod(k,2))/(xs_mod(l,1)+xs_mod(l,2))
                   sigNNLO(i,j,k,l) = sig2(i,j)*1d0/(ek+del_jk*ei)/wkl
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine get_sigNNLO

  double precision function hatsigNNLO(i,j,k,l,r,n_ext_in,xs_in,alpha_in)
    implicit none
    ! global
    include 'nexternal.inc'
    integer :: n_ext_in
    double precision :: alpha_in
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    ! local
    integer :: i,j,k,l,r,del_jk
    double precision :: ei,ej,ek
    double precision :: wij, wik, wjk, wkl, wir, wkr, hatsig2

    ! set global module variables
    n_ext=n_ext_in
!    num_sec=(n_ext-2)*(n_ext-3)/2
    num_sec=(nexternal-2)*(nexternal-3)/2
    if (.not.allocated(xs_mod)) allocate(xs_mod(n_ext,n_ext))
    xs_mod=xs_in
    alpha_mod=alpha_in
    ! calculate 2-index and 4-index sigma hat
    hatsig2=0d0
    hatsigNNLO=0d0
    if(r.eq.i.or.i.eq.j.or.r.eq.j.or.r.eq.k.or.i.eq.k.or.r.eq.l.or.l.eq.i.or.l.eq.k) then
      write(77,*) 'Wrong indices called in hatsigNNLO', i,j,k,l,r
      stop
    endif
    del_jk=0
    if(j.eq.k)del_jk = 1
    if((xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2))*&
         xs_mod(i,j)*xs_mod(1,2)*(xs_mod(r,1)+xs_mod(r,2))*&
         (xs_mod(k,1)+xs_mod(k,2))*(xs_mod(l,1)+xs_mod(l,2))*xs_mod(k,l).ne.0d0) then
      ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
      ek=(xs_mod(k,1)+xs_mod(k,2))/xs_mod(1,2)
      wij=xs_mod(1,2)*xs_mod(i,j)/&
              (xs_mod(i,1)+xs_mod(i,2))/(xs_mod(j,1)+xs_mod(j,2))
      wir=xs_mod(1,2)*xs_mod(i,r)/&
              (xs_mod(i,1)+xs_mod(i,2))/(xs_mod(r,1)+xs_mod(r,2))
      wkr=xs_mod(1,2)*xs_mod(k,r)/&
              (xs_mod(k,1)+xs_mod(k,2))/(xs_mod(r,1)+xs_mod(r,2))
      wkl=xs_mod(1,2)*xs_mod(k,l)/&
              (xs_mod(k,1)+xs_mod(k,2))/(xs_mod(l,1)+xs_mod(l,2))
      hatsig2=(1d0/ei/wij/wir)**alpha_mod
      hatsigNNLO = hatsig2*1d0/(ek*wkr+del_jk*ei*wir)/wkl
    endif
  end function hatsigNNLO


   subroutine get_WCC_NNLO(xs_in,IA,IB,C,D,ir,alphaz,n_ext_in)
    !     NNLO collinear sector functions WCC(ia,ib,ic,ir)
    implicit none
    include 'all_sector_list.inc'
    integer :: ia,ib,ic,ir,c,d
    integer :: n_ext_in
    double precision :: alphaz, num, sigma, wcc_nnlo
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    num = hatsigNNLO(IA,IB,C,D,ir,n_ext_in,xs_in,alphaz)
    if (IB .eq. C) then
       ic = D
    elseif (IB .eq. D) then
       ic = C
    end if
    sigma = hatsigNNLO(ia,ib,ib,ic,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ia,ic,ib,ic,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ib,ia,ia,ic,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ib,ic,ia,ic,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ia,ib,ic,ib,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ia,ic,ic,ib,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ic,ia,ia,ib,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ic,ib,ia,ib,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ib,ia,ic,ia,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ib,ic,ic,ia,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ic,ia,ib,ia,ir,n_ext_in,xs_in,alphaz) + &
            hatsigNNLO(ic,ib,ib,ia,ir,n_ext_in,xs_in,alphaz)
    wcc_nnlo=num/sigma
  end subroutine get_WCC_NNLO

  subroutine get_WSS_CC_NNLO(xs_in,ia,ib,C,D,ir,alphaz,n_ext_in)
          !     NNLO double-soft collinear sector functions WSSCC(ia,ib,ic,ir)
    implicit none
    include 'all_sector_list.inc'
    integer :: ia,ib,ic,ir,C,D
    integer :: n_ext_in
    double precision :: alphaz, num, sigma, wsscc_nnlo
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    num = hatsigNNLO(IA,IB,C,D,ir,n_ext_in,xs_in,alphaz)
    if(IB.eq.C) then
      ic = D
      sigma = hatsigNNLO(ia,ib,ib,ic,ir,n_ext_in,xs_in,alphaz) + &
              hatsigNNLO(ia,ic,ib,ic,ir,n_ext_in,xs_in,alphaz) + &
              hatsigNNLO(ib,ia,ia,ic,ir,n_ext_in,xs_in,alphaz) + &
              hatsigNNLO(ib,ic,ia,ic,ir,n_ext_in,xs_in,alphaz)
    elseif(IB.eq.D) then
      ic = C
      sigma = hatsigNNLO(ia,ib,ic,ib,ir,n_ext_in,xs_in,alphaz) + &
              hatsigNNLO(ia,ic,ic,ib,ir,n_ext_in,xs_in,alphaz) + &
              hatsigNNLO(ic,ia,ia,ib,ir,n_ext_in,xs_in,alphaz) + &
              hatsigNNLO(ic,ib,ia,ib,ir,n_ext_in,xs_in,alphaz)
    endif
      wsscc_nnlo=num/sigma
  end subroutine get_WSS_CC_NNLO

  subroutine get_Z_NNLO(i1,i2,i3,i4)
    !     NNLO sector functions Z(i1,i2,i3,i4)
    implicit none
    include 'all_sector_list.inc'
    integer :: i,a,b,c,d,i1,i2,i3,i4
    double precision :: num,sigma
    call sector4_global_checks(i1,i2,i3,i4)
    if(i4.eq.0) then
       num = sigNNLO(i1,i2,i2,i3) + &
             sigNNLO(i1,i2,i3,i2) + &
             sigNNLO(i1,i3,i3,i2) + &
             sigNNLO(i1,i3,i2,i3) + &
             sigNNLO(i2,i1,i1,i3) + &
             sigNNLO(i2,i1,i3,i1) + &
             sigNNLO(i2,i3,i3,i1) + &
             sigNNLO(i2,i3,i1,i3) + &
             sigNNLO(i3,i1,i1,i2) + &
             sigNNLO(i3,i1,i2,i1) + &
             sigNNLO(i3,i2,i2,i1) + &
             sigNNLO(i3,i2,i1,i2)
    elseif(i4.ne.0) then
       num = sigNNLO(i1,i2,i3,i4)  + &
             sigNNLO(i1,i2,i4,i3)  + &
             sigNNLO(i3,i4,i1,i2)  + &
             sigNNLO(i3,i4,i2,i1)  + &
             sigNNLO(i4,i3,i1,i2)  + &
             sigNNLO(i4,i3,i2,i1)
    else
       write(*,*) 'get_Z_NNLO: error in the numerator construction...'
       write(*,*) 'negative value for 4th sector index i4...'
       write(*,*) 'i4 = ', i4
       write(*,*) 'exit...'
       stop
    endif
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       c=all_sector_list(3,i)
       d=all_sector_list(4,i)
       if(d.eq.0) then
          sigma = sigma + &
               sigNNLO(a,b,b,c) + &
               sigNNLO(a,b,c,b) + &
               sigNNLO(a,c,c,b) + &
               sigNNLO(a,c,b,c) + &
               sigNNLO(b,a,a,c) + &
               sigNNLO(b,a,c,a) + &
               sigNNLO(b,c,c,a) + &
               sigNNLO(b,c,a,c) + &
               sigNNLO(c,a,a,b) + &
               sigNNLO(c,a,b,a) + &
               sigNNLO(c,b,b,a) + &
               sigNNLO(c,b,a,b)
       elseif(d.ne.0) then
          sigma = sigma + &
               sigNNLO(a,b,c,d) + &
               sigNNLO(a,b,d,c) + &
               sigNNLO(b,a,c,d) + &
               sigNNLO(b,a,d,c) + &
               sigNNLO(c,d,a,b) + &
               sigNNLO(c,d,b,a) + &
               sigNNLO(d,c,a,b) + &
               sigNNLO(d,c,b,a)
       else
          write(*,*) 'get_Z_NNLO: error in the denominator construction...'
          write(*,*) 'negative value for 4th sector index d...'
          write(*,*) 'd = ', d
          write(*,*) 'exit...'
          stop
       endif
    enddo
    Z_NNLO = num/sigma
    call sector2_sanity_checks(sigma,Z_NNLO)
  end subroutine get_Z_NNLO

  subroutine get_WSS_NNLO(i1,i2,i3,i4)
    !     NNLO double-soft sector functions WSS(i1,i2,i3,i4) = barS_i1i3 W(i1,i2,i3,i4) [eq.(C.54)]
    implicit none
    include 'all_sector_list.inc'
    integer :: i,a,b,c,d,i1,i2,i3,i4
    double precision :: num,sigma
    call sector4_global_checks(i1,i2,i3,i4)
    if(i4.ge.0) then
       num = sigNNLO(i1,i2,i3,i4)
    else
       write(*,*) 'get_WSS_NNLO: error in the construction of numerator'
       write(*,*) 'Negative value for 4th sector index i4...'
       write(*,*) 'i4 = ', i4
       write(*,*) 'exit...'
       stop
    endif
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       c=all_sector_list(3,i)
       d=all_sector_list(4,i)
       if(d.eq.0) then
          if((a.eq.i1.and.c.eq.i3).or.(a.eq.i3.and.c.eq.i1)) sigma = sigma + &
               sigNNLO(a,b,c,b) + sigNNLO(a,c,c,b) +  &
               sigNNLO(c,b,a,b) + sigNNLO(c,a,a,b)
          if((a.eq.i1.and.b.eq.i3).or.(a.eq.i3.and.b.eq.i1)) sigma = sigma + &
               sigNNLO(a,c,b,c) + sigNNLO(a,b,b,c) + &
               sigNNLO(b,c,a,c) + sigNNLO(b,a,a,c)
          if((b.eq.i1.and.c.eq.i3).or.(b.eq.i3.and.c.eq.i1)) sigma = sigma + &
               sigNNLO(b,a,c,a) + sigNNLO(b,c,c,a) + &
               sigNNLO(c,a,b,a) + sigNNLO(c,b,b,a)

       elseif(d.ne.0) then
          if((a.eq.i1.and.c.eq.i3).or.(a.eq.i3.and.c.eq.i1)) sigma = sigma + &
               sigNNLO(a,b,c,d) + sigNNLO(c,d,a,b)
          if((a.eq.i1.and.d.eq.i3).or.(a.eq.i3.and.d.eq.i1)) sigma = sigma + &
               sigNNLO(a,b,d,c) + sigNNLO(d,c,a,b)
          if((b.eq.i1.and.c.eq.i3).or.(b.eq.i3.and.c.eq.i1)) sigma = sigma + &
               sigNNLO(b,a,c,d) + sigNNLO(c,d,b,a)
          if((b.eq.i1.and.d.eq.i3).or.(b.eq.i3.and.d.eq.i1)) sigma = sigma + &
               sigNNLO(b,a,d,c) + sigNNLO(d,c,b,a)
       else
          write(*,*) 'get_WSS_NNLO: error in the construction of denominator'
          write(*,*) 'Negative value for 4th sector index i4...'
          write(*,*) 'i4 = ', i4
          write(*,*) 'exit...'
          stop
       endif
    enddo
    WSS_NNLO = num/sigma
    call sector2_sanity_checks(sigma,WSS_NNLO)
  end subroutine get_WSS_NNLO



  subroutine get_ZS_NNLO(i1,i2,sec_list)
    !     NNLO 2-index mapped sector function relevant to the barHCijbarSij limit
    implicit none
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    ! This list contains the pairs (bar{isec}, bar{jsec})
    ! The second entry runs over the number of final state NLO particles
    integer, dimension (2,num_sec) :: sec_list

    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,num_sec
       a = sec_list(1,i)
       b = sec_list(2,i)
       if(a.eq.0.or.b.eq.0) cycle
       if(a.eq.i1) sigma = sigma + sig2(a,b)
       if(b.eq.i1) sigma = sigma + sig2(b,a)
    enddo
    ZS_NNLO = num/sigma
    call sector2_sanity_checks(sigma,ZS_NNLO)
  end subroutine get_ZS_NNLO


  subroutine get_ZHC_NNLO(i1,i2,sec_list)
    !     NNLO 2-index mapped sector function relevant to the barHCij limit
    implicit none
    integer :: i,a,b,i1,i2
    double precision num,sigma
    ! This list contains the pairs (bar{isec}, bar{jsec})
    ! The second entry runs over the number of final state NLO particles
    integer, dimension (2,num_sec) :: sec_list

    num = sig2(i1,i2) + sig2(i2,i1)
    sigma = 0d0
    do i=1,num_sec
       a = sec_list(1,i)
       b = sec_list(2,i)
       if(a.eq.0.or.b.eq.0) cycle
       sigma = sigma + sig2(a,b) + sig2(b,a)
    enddo
    Z_HC_NNLO = num/sigma
    call sector2_sanity_checks(sigma,Z_HC_NNLO)
  end subroutine get_ZHC_NNLO


  subroutine sector4_global_checks(i1,i2,i3,i4)
    implicit none
    integer :: i1,i2,i3,i4
    if(alpha_mod.lt.1d0)then
       write(77,*)'Wrong alpha_mod in sectors4',alpha_mod
       stop
    endif
    if(i1.le.2.or.i2.le.2.or.i3.le.2) then
       if(i4.ne.0.and.i4.le.2) then
          write(77,*) 'sectors4: indices must be in final state',i1,i2,i3,i4
       elseif(i4.eq.0) then
          write(77,*) 'sectors4: indices must be in final state',i1,i2,i3
       endif
       stop
    endif
  end subroutine sector4_global_checks


  subroutine sector2_sanity_checks(sigma,Z)
    implicit none
    double precision :: Z,sigma
    if(sigma.le.0d0)then
       write(*,*)'Wrong sigma ',sigma
       stop
    endif
    if(abs(Z).ge.huge(1d0).or.isnan(Z))then
       write(77,*)'Exception caught ',Z
       stop
    endif
  end subroutine sector2_sanity_checks


end module sectors4_module
