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

  subroutine get_WS_NLO(i1,i2)
    !     NLO soft sector functions WS(i1,i2) = barS_i1 W(i1,i2)
    implicit none
    integer :: i,i1,i2,sec(2)
    double precision :: num,sigma
    include 'all_K_sector_list.inc'
    call sector2_global_checks(i1,i2)
    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,len
       sec=s_sector_list(i1,i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma + &
            sig2(sec(1),sec(2))
    enddo
    WS_NLO = num/sigma
    call sector2_sanity_checks(sigma,WS_NLO)
  end subroutine get_WS_NLO

  subroutine get_WC_NLO(i1,i2,ir)
    !     NLO collinear sector functions WC(i1,i2) = barC_i1 W(i1,i2)
    implicit none
    integer :: i,i1,i2,ir,sec(2)
    double precision :: num,sigma
    include 'all_K_sector_list.inc'
    call sector2_global_checks(i1,ir)
    num = sig2(i1,ir)
    sigma = 0d0
    do i=1,len
       sec=c_sector_list(i1,i2,i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma + &
            sig2(sec(1),sec(2))
    enddo
    WC_NLO = num/sigma
    call sector2_sanity_checks(sigma,WC_NLO)
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

end module sectors2_module
