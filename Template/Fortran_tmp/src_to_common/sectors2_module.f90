module sectors2_module
  implicit none
  integer, public :: n_ext
  double precision, public :: alpha_mod, Z_NLO, ZS_NLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:), public :: sig2
  public :: get_sig2, get_Z_NLO, get_ZS_NLO
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
          if( (xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2))* &
               xs_mod(i,j)*xs_mod(1,2).ne.0d0 )then
             ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
             ej=(xs_mod(j,1)+xs_mod(j,2))/xs_mod(1,2)
             wij=xs_mod(1,2)*xs_mod(i,j)/(xs_mod(i,1)+xs_mod(i,2))/ &
                  (xs_mod(j,1)+xs_mod(j,2))
             sig2(i,j)=(1d0/ei/wij)**alpha_mod
          endif
       enddo
    enddo
  end subroutine get_sig2


  subroutine get_Z_NLO(i1,i2)
    !     NLO sector functions Z(i1,i2)
    !     This function is meant to be called with (i1,i2) = perm(isec,jsec)
    !     i1 must be in the final state 
    implicit none
    include 'all_sector_list.inc'
    integer :: i,a,b,i1,i2
    double precision :: num,sigma
    call sector2_global_checks(i1,i2)
    num = sig2(i1,i2) + sig2(i2,i1)
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       sigma = sigma + sig2(a,b) + sig2(b,a)
    enddo
    Z_NLO = num/sigma
    call sector2_sanity_checks(sigma,Z_NLO)
  end subroutine get_Z_NLO


  subroutine get_ZS_NLO(i1,i2)
    !     NLO soft sector functions ZS(i1,i2) = S_i1 Z(i1,i2)
    !     This function is meant to be called with (i1,i2) = perm(isec,jsec)
    !     i1 must be in the final state and is the one associated with the
    !     soft singularity
    implicit none
    include 'all_sector_list.inc'
    integer i,a,b,i1,i2
    double precision num,sigma
    call sector2_global_checks(i1,i2)
    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,lensectors
       a=all_sector_list(1,i)
       b=all_sector_list(2,i)
       ! check
       ! 3=g, 4=g, 5=q, sectors ab = 34, 35, 45
       ! sector 34 -> ZS(i1=3,i2=4) = S3.Z34 = sig2(3,4) / (sig2(3,4) + sig2(3,5))
       ! sector 34 -> ZS(i1=4,i2=3) = S4.Z34 = sig2(4,3) / (sig2(4,3) + sig2(4,5))
       ! sector 35 -> ZS(i1=3,i2=5) = S3.Z35 = sig2(3,5) / (sig2(3,4) + sig2(3,5))
       ! sector 45 -> ZS(i1=4,i2=5) = S4.Z45 = sig2(4,5) / (sig2(4,3) + sig2(4,5))
       if(a.eq.i1) sigma = sigma + sig2(a,b)
       if(b.eq.i1) sigma = sigma + sig2(b,a)
    enddo
    ZS_NLO = num/sigma
    call sector2_sanity_checks(sigma,ZS_NLO)
  end subroutine get_ZS_NLO


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

