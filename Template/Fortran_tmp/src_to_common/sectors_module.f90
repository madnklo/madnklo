module sectors2_module
  implicit none
  integer, public :: n_ext
  double precision, public :: alpha_mod, W_NLO, WS_NLO, WC_NLO
  double precision, public :: Wbar_NLO, WSbar_NLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:), public :: sig2
  public :: get_sig2, get_W_NLO, get_WS_NLO, get_WC_NLO
  public :: get_Wbar_NLO, get_WSbar_NLO
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
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    include 'all_sector_list.inc'
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

  subroutine get_Wbar_NLO(i1,i2)
    !     NLO sector functions W(i1,i2)
    implicit none
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    include 'all_sector_list_real.inc'
    call sector2_global_checks(i1,i2)
    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       sigma = sigma + sig2(a,b)
    enddo
    Wbar_NLO = num/sigma
    call sector2_sanity_checks(sigma,Wbar_NLO)
  end subroutine get_Wbar_NLO

  subroutine get_WS_NLO(i1,i2)
    !     NLO soft sector functions WS(i1,i2) = barS_i1 W(i1,i2)
    implicit none
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    include 'all_sector_list.inc'
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

  subroutine get_WSbar_NLO(i1,i2)
    !     NLO sector functions WSbar(i1,i2)
    implicit none
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    include 'all_sector_list_real.inc'
    call sector2_global_checks(i1,i2)
    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       sigma = sigma + sig2(a,b)
    enddo
    WSbar_NLO = num/sigma
    call sector2_sanity_checks(sigma,WSbar_NLO)
  end subroutine get_WSbar_NLO

  subroutine get_WC_NLO(xs_in,ia,ib,ir,alphaz,n_ext_in)
    !     NLO collinear sector functions WC(ia,ib,ir)
    implicit none
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
  integer, public :: n_ext
  double precision, public :: alpha_mod, W_NNLO, WSS_NNLO, WCC_NNLO, WSS_CC_NNLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:), public :: sig2,hatsig2
  double precision, allocatable, dimension(:,:,:,:), public :: sigNNLO,hatsigNNLO
  public :: get_sigNNLO, get_hatsigNNLO, get_W_NNLO, get_WSS_NNLO, get_WCC_NNLO, get_WSS_CC_NNLO
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

  subroutine get_hatsigNNLO(iref_in,xs_in,alpha_in,n_ext_in)
    implicit none
    ! global
    include 'nexternal.inc'
    integer :: n_ext_in, iref_in
    double precision :: alpha_in
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    ! local
    integer :: i,j,k,l,r, del_jk
    double precision :: ei, ek, wij, wkl, wir, wkr

    ! set global module variables
    n_ext=n_ext_in
    if (.not.allocated(xs_mod)) allocate(xs_mod(n_ext,n_ext))
    if (.not.allocated(hatsig2)) allocate(hatsig2(3:n_ext,3:n_ext))
    if (.not.allocated(hatsigNNLO)) allocate(hatsigNNLO(3:n_ext,3:n_ext,3:n_ext,3:n_ext))
    xs_mod=xs_in
    alpha_mod=alpha_in
    ! calculate 2-index and 4-index hatsigma
    r = iref_in
    hatsig2=0d0
    hatsigNNLO=0d0
    do i=3,n_ext
       if(i.eq.r)cycle
       wir=xs_mod(1,2)*xs_mod(i,r)/&
            (xs_mod(i,1)+xs_mod(i,2))/(xs_mod(r,1)+xs_mod(r,2))
       do j=3,n_ext
          if(j.eq.i.or.j.eq.r)cycle
          if( (xs_mod(i,1)+xs_mod(i,2))*&
              (xs_mod(j,1)+xs_mod(j,2))*&
               xs_mod(i,j)*xs_mod(1,2).ne.0d0 ) then
             ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
             wij=xs_mod(1,2)*xs_mod(i,j)/&
                  (xs_mod(i,1)+xs_mod(i,2))/(xs_mod(j,1)+xs_mod(j,2))
             hatsig2(i,j)=(1d0/ei/wij/wir)**alpha_mod
          endif
          do k=3,n_ext
             if(k.eq.i.or.k.eq.r)cycle
             ek=(xs_mod(k,1)+xs_mod(k,2))/xs_mod(1,2)
             wkr=xs_mod(1,2)*xs_mod(k,r)/&
                  (xs_mod(k,1)+xs_mod(k,2))/(xs_mod(r,1)+xs_mod(r,2))
             if(k.eq.j)then
                do l=3,n_ext
                   if(l.eq.i.or.l.eq.k.or.l.eq.j.or.l.eq.r)cycle
                   if( (xs_mod(k,1)+xs_mod(k,2))*&
                        (xs_mod(l,1)+xs_mod(l,2))*&
                        xs_mod(k,l)*xs_mod(1,2).ne.0d0 ) then
                      wkl=xs_mod(1,2)*xs_mod(k,l)/&
                           (xs_mod(k,1)+xs_mod(k,2))/(xs_mod(l,1)+xs_mod(l,2))
                      hatsigNNLO(i,j,k,l) = hatsig2(i,j)/(ek*wkr+ei*wir)/wkl
                   endif
                enddo
             else
                l=j
                if( (xs_mod(k,1)+xs_mod(k,2))*&
                     (xs_mod(l,1)+xs_mod(l,2))*&
                     xs_mod(k,l)*xs_mod(1,2).ne.0d0 ) then
                   wkl=xs_mod(1,2)*xs_mod(k,l)/&
                        (xs_mod(k,1)+xs_mod(k,2))/(xs_mod(l,1)+xs_mod(l,2))
                   hatsigNNLO(i,j,k,l) = hatsig2(i,j)/(ek*wkr)/wkl
                endif
             endif
          enddo
       enddo
    enddo
  end subroutine get_hatsigNNLO

  subroutine get_W_NNLO(a,b,c,d)
    ! NNLO sector functions
    implicit none
    integer :: i,a,b,c,d,i1,i2,i3,i4
    double precision :: num,sigma
    include 'all_sector_list.inc'
    call sector4_global_checks(a,b,c,d)
    num = sigNNLO(a,b,c,d)
    sigma = 0d0
    do i=1,lensectors
       i1=all_sector_list(1,i)
       i2=all_sector_list(2,i)
       i3=all_sector_list(3,i)
       i4=all_sector_list(4,i)
       sigma = sigma + sigNNLO(i1,i2,i3,i4)
    enddo
    W_NNLO = num/sigma
    call sector2_sanity_checks(sigma,W_NNLO)
  end subroutine get_W_NNLO

  subroutine get_WSS_NNLO(a,b,c,d)
    ! NNLO double-soft sector functions WSS
    implicit none
    integer :: i,a,b,c,d,i1,i2,i3,i4
    double precision :: num,sigma
    include 'all_sector_list.inc'
    call sector4_global_checks(a,b,c,d)
    num = sigNNLO(a,b,c,d)
    sigma = 0d0
    do i=1,lensectors
       i1=all_sector_list(1,i)
       i2=all_sector_list(2,i)
       i3=all_sector_list(3,i)
       i4=all_sector_list(4,i)
       sigma = sigma + &
               sigNNLO(i1,i2,i3,i4) + sigNNLO(i3,i2,i1,i4)
    enddo
    WSS_NNLO = num/sigma
    call sector2_sanity_checks(sigma,WSS_NNLO)
  end subroutine get_WSS_NNLO

  subroutine get_WCC_NNLO(ia,ib,c,d)
    ! NNLO triple-collinear sector functions WCC
    implicit none
    integer :: ia,ib,ic,ir,c,d
    double precision :: num, sigma
    num = hatsigNNLO(ia,ib,c,d)
    if (ib .eq. c) then
       ic = d
    elseif (ib .eq. d) then
       ic = c
    end if
    sigma = hatsigNNLO(ia,ib,ib,ic) + hatsigNNLO(ia,ic,ib,ic) + &
            hatsigNNLO(ib,ia,ia,ic) + hatsigNNLO(ib,ic,ia,ic) + &
            hatsigNNLO(ia,ib,ic,ib) + hatsigNNLO(ia,ic,ic,ib) + &
            hatsigNNLO(ic,ia,ia,ib) + hatsigNNLO(ic,ib,ia,ib) + &
            hatsigNNLO(ib,ia,ic,ia) + hatsigNNLO(ib,ic,ic,ia) + &
            hatsigNNLO(ic,ia,ib,ia) + hatsigNNLO(ic,ib,ib,ia)
    wcc_nnlo=num/sigma
  end subroutine get_WCC_NNLO

  subroutine get_WSS_CC_NNLO(ia,ib,c,d)
    ! NNLO double-soft triple-collinear sector functions WSSCC
    implicit none
    integer :: ia,ib,ic,ir,c,d
    double precision :: num, sigma
    num = hatsigNNLO(ia,ib,c,d)
    if(ib.eq.c) then
      ic = d
      sigma = hatsigNNLO(ia,ib,ib,ic) + hatsigNNLO(ia,ic,ib,ic) + &
              hatsigNNLO(ib,ia,ia,ic) + hatsigNNLO(ib,ic,ia,ic)
    elseif(ib.eq.d) then
      ic = c
      sigma = hatsigNNLO(ia,ib,ic,ib) + hatsigNNLO(ia,ic,ic,ib) + &
              hatsigNNLO(ic,ia,ia,ib) + hatsigNNLO(ic,ib,ia,ib)
    endif
      wss_cc_nnlo=num/sigma
  end subroutine get_WSS_CC_NNLO

  subroutine sector4_global_checks(i1,i2,i3,i4)
    implicit none
    integer :: i1,i2,i3,i4
    if(alpha_mod.lt.1d0)then
       write(77,*)'Wrong alpha_mod in sectors4',alpha_mod
       stop
    endif
    if(i1.le.2.or.i2.le.2.or.i3.le.2.or.i4.le.2) then
       write(77,*) 'sectors4: indices must be in final state',i1,i2,i3,i4
       stop
    endif
  end subroutine sector4_global_checks

  subroutine sector2_sanity_checks(sigma,W)
    implicit none
    double precision :: W,sigma
    if(sigma.le.0d0)then
       write(*,*)'Wrong sigma ',sigma
       stop
    endif
    if(abs(W).ge.huge(1d0).or.isnan(W))then
       write(77,*)'Exception caught ',W
       stop
    endif
  end subroutine sector2_sanity_checks


end module sectors4_module
