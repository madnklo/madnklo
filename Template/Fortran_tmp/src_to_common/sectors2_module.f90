module sectors2_module
  implicit none

  integer, public :: n_ext
  double precision, public, parameter :: alpha_mod=1d0
  double precision, public :: W_NLO, WS_NLO, WC_NLO
  double precision, public :: Wbar_NLO, WSbar_NLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:), public :: sig2

  public :: get_sig2, get_W_NLO, get_WS_NLO, get_WC_NLO
  private

  ! Denominators belonging to the current phase-space point.  They are
  ! invalidated by every call to get_sig2 and evaluated only when needed.
  double precision :: sigma_nlo
  double precision, allocatable :: sigma_soft(:)
  double precision, allocatable :: sigma_coll(:,:)
  logical :: sigma_nlo_valid=.false.
  logical, allocatable :: sigma_soft_valid(:)
  logical, allocatable :: sigma_coll_valid(:,:)

contains

  subroutine get_sig2(xs_in,n_ext_in)
    implicit none
    integer, intent(in) :: n_ext_in
    double precision, intent(in) :: xs_in(n_ext_in,n_ext_in)
    integer :: i,j
    double precision :: ei,ej,wij,sigma_ij

    n_ext=n_ext_in
    call ensure_workspace(n_ext)
    xs_mod=xs_in

    ! Calculate each symmetric pair only once.  Sector labels are ordered,
    ! but sigma_ij itself is symmetric under i <-> j.
    sig2=0d0
    do i=3,n_ext-1
       do j=i+1,n_ext
          if ((xs_mod(i,1)+xs_mod(i,2)) * &
              (xs_mod(j,1)+xs_mod(j,2)) * &
              xs_mod(i,j)*xs_mod(1,2).ne.0d0) then
             ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
             ej=(xs_mod(j,1)+xs_mod(j,2))/xs_mod(1,2)
             wij=xs_mod(1,2)*xs_mod(i,j) / &
                  (xs_mod(i,1)+xs_mod(i,2)) / &
                  (xs_mod(j,1)+xs_mod(j,2))
             sigma_ij=(1d0/ei/wij)**alpha_mod
             sig2(i,j)=sigma_ij
             sig2(j,i)=sigma_ij
          endif
       enddo
    enddo

    ! A new phase-space point requires new normalization denominators.
    sigma_nlo=0d0
    sigma_nlo_valid=.false.
    sigma_soft=0d0
    sigma_coll=0d0
    sigma_soft_valid=.false.
    sigma_coll_valid=.false.
  end subroutine get_sig2


  subroutine ensure_workspace(n_ext_in)
    implicit none
    integer, intent(in) :: n_ext_in
    logical :: resize

    resize=.not.allocated(xs_mod)
    if (.not.resize) then
       resize=size(xs_mod,1).ne.n_ext_in
    endif

    if (resize) then
       if (allocated(xs_mod)) deallocate(xs_mod)
       if (allocated(sig2)) deallocate(sig2)
       if (allocated(sigma_soft)) deallocate(sigma_soft)
       if (allocated(sigma_coll)) deallocate(sigma_coll)
       if (allocated(sigma_soft_valid)) deallocate(sigma_soft_valid)
       if (allocated(sigma_coll_valid)) deallocate(sigma_coll_valid)

       allocate(xs_mod(n_ext_in,n_ext_in))
       allocate(sig2(3:n_ext_in,3:n_ext_in))
       allocate(sigma_soft(1:n_ext_in))
       allocate(sigma_coll(1:n_ext_in,1:n_ext_in))
       allocate(sigma_soft_valid(1:n_ext_in))
       allocate(sigma_coll_valid(1:n_ext_in,1:n_ext_in))
    endif
  end subroutine ensure_workspace


  subroutine get_W_NLO(i1,i2)
    ! NLO sector function W(i1,i2).
    implicit none
    integer, intent(in) :: i1,i2
    integer :: i,a,b
    include 'all_sector_list.inc'

    if (.not.sigma_nlo_valid) then
       sigma_nlo=0d0
       do i=1,lensectors
          a=all_sector_list(1,i)
          b=all_sector_list(2,i)
          sigma_nlo=sigma_nlo+sig2(a,b)
       enddo
       sigma_nlo_valid=.true.
    endif

    W_NLO=sig2(i1,i2)/sigma_nlo
  end subroutine get_W_NLO


  subroutine get_WS_NLO(i1,i2)
    ! Soft limit WS(i1,i2) = barS_i1 W(i1,i2).
    implicit none
    integer, intent(in) :: i1,i2
    integer :: i,sec(2)
    include 'all_K_sector_list.inc'

    if (.not.sigma_soft_valid(i1)) then
       sigma_soft(i1)=0d0
       do i=1,len
          sec=s_sector_list(i1,i,:)
          if (all(sec.eq.0)) cycle
          sigma_soft(i1)=sigma_soft(i1)+sig2(sec(1),sec(2))
       enddo
       sigma_soft_valid(i1)=.true.
    endif

    WS_NLO=sig2(i1,i2)/sigma_soft(i1)
  end subroutine get_WS_NLO


  subroutine get_WC_NLO(i1,i2,ir)
    ! Collinear limit WC(i1,i2) = barC_i1i2 W(i1,i2).
    implicit none
    integer, intent(in) :: i1,i2,ir
    integer :: i,sec(2)
    include 'all_K_sector_list.inc'

    if (.not.sigma_coll_valid(i1,i2)) then
       sigma_coll(i1,i2)=0d0
       do i=1,len
          sec=c_sector_list(i1,i2,i,:)
          if (all(sec.eq.0)) cycle
          sigma_coll(i1,i2)=sigma_coll(i1,i2)+sig2(sec(1),sec(2))
       enddo
       sigma_coll_valid(i1,i2)=.true.
    endif

    WC_NLO=sig2(i1,ir)/sigma_coll(i1,i2)
  end subroutine get_WC_NLO

end module sectors2_module
