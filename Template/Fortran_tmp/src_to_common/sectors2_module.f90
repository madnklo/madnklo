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






module sectors3_module
  implicit none
  integer, public :: n_ext
  double precision, public :: alpha_mod, Z_NNLO, ZSS_NNLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:,:,:), public :: sigNNLO
  public :: get_sigNNLO, get_Z_NNLO, get_ZSS_NNLO
  private

contains

  subroutine get_sigNNLO(xs_in,alpha_in,n_ext_in)
    implicit none
    ! global
    integer :: n_ext_in
    double precision :: alpha_in
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    ! local
    integer :: i,j,k,l
    integer del_jk
    double precision :: ei,ej,ek
    double precision :: wij, wik, wjk, wkl
    
    ! set global module variables
    n_ext=n_ext_in
    if (.not.allocated(xs_mod)) allocate(xs_mod(n_ext,n_ext))
    if (.not.allocated(sigNNLO)) allocate(sigNNLO(3:n_ext,3:n_ext,3:n_ext,3:n_ext))
    xs_mod=xs_in
    alpha_mod=alpha_in
    ! calculate 4-index sigma
    sigNNLO=0d0
    do i=3,n_ext
       do j=3,n_ext
          do k=3,n_ext
              del_jk=0
             if(j.eq.k) del_jk = 1
             do l=3,n_ext
                if(i.eq.j.or.i.eq.k.or.k.eq.l.or.i.eq.l)cycle
                if( (xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2))* &
                     xs_mod(i,j)*xs_mod(1,2).ne.0d0 .or. &
                     (xs_mod(k,1)+xs_mod(k,2))*(xs_mod(l,1)+xs_mod(l,2))* &
                     xs_mod(k,l)*xs_mod(1,2).ne.0d0)then

                   ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
                   ek=(xs_mod(k,1)+xs_mod(k,2))/xs_mod(1,2)
                   wij=xs_mod(1,2)*xs_mod(i,j)/(xs_mod(i,1)+xs_mod(i,2))/ &
                        (xs_mod(j,1)+xs_mod(j,2))
                   wkl=xs_mod(1,2)*xs_mod(k,l)/(xs_mod(k,1)+xs_mod(k,2))/ &
                        (xs_mod(l,1)+xs_mod(l,2))
                   
                   sigNNLO(i,j,k,l) = 1d0/(ei*wij)**alpha_mod &
                        *1d0/(ek+del_jk*ei)/wkl
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine get_sigNNLO


  subroutine get_Z_NNLO(i1,i2,i3,i4)
    !     NNLO sector functions Z(i1,i2,i3,i4)
    !     This function is meant to be called with (i1,i2,i3,i4) = perm(isec,jsec,ksec,lsec)
    !     i1 must be in the final state 
    implicit none
    include 'all_sector_list.inc'
    integer :: i,a,b,c,d,i1,i2,i3,i4
    double precision :: num,sigma
    call sector4_global_checks(i1,i2,i3,i4)
    if(i4.eq.0) then
       num = sigNNLO(i1,i2,i2,i3) + sigNNLO(i1,i2,i3,i2)  &
           + sigNNLO(i1,i3,i3,i2) + sigNNLO(i1,i3,i2,i3)  &
           + sigNNLO(i2,i1,i1,i3) + sigNNLO(i2,i1,i3,i1)  &
           + sigNNLO(i2,i3,i3,i1) + sigNNLO(i2,i3,i1,i3)  &
           + sigNNLO(i3,i1,i1,i2) + sigNNLO(i3,i1,i2,i1)  &
           + sigNNLO(i3,i2,i2,i1) + sigNNLO(i3,i2,i1,i2)
    elseif(i4.ne.0) then
       num = sigNNLO(i1,i2,i3,i4) + sigNNLO(i1,i2,i4,i3) &
            +sigNNLO(i2,i1,i3,i4) + sigNNLO(i2,i1,i4,i3) &
            +sigNNLO(i3,i4,i1,i2) + sigNNLO(i3,i4,i1,i2) &
            +sigNNLO(i4,i3,i1,i2) + sigNNLO(i4,i3,i2,i1)
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
          sigma = sigma &
           + sigNNLO(a,b,b,c) + sigNNLO(a,b,c,b)  &    
           + sigNNLO(a,c,c,b) + sigNNLO(a,c,b,c)  &
           + sigNNLO(b,a,a,c) + sigNNLO(b,a,c,a)  &
           + sigNNLO(b,c,c,a) + sigNNLO(b,c,a,c)  &
           + sigNNLO(c,a,a,b) + sigNNLO(c,a,b,a)  &
           + sigNNLO(c,b,b,a) + sigNNLO(c,b,a,b)
       elseif(d.ne.0) then
          sigma = sigma &
            + sigNNLO(a,b,c,d) + sigNNLO(a,b,d,c) &    
            + sigNNLO(b,a,c,d) + sigNNLO(b,a,d,c) &
            + sigNNLO(c,d,a,b) + sigNNLO(c,d,b,a) &
            + sigNNLO(d,c,a,b) + sigNNLO(d,c,b,a)
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


  subroutine get_ZSS_NNLO(i1,i2,i3,i4)
    !     NNLO double-soft sector functions ZSS(i1,i2) = S_i1 Z(i1,i2)
    !     This function is meant to be called with (i1,i2) = perm(isec,jsec)
    !     i1 must be in the final state and is the one associated with the
    !     soft singularity
    implicit none
    include 'all_sector_list.inc'
    integer i,a,b,c,d,i1,i2,i3,i4
    double precision num,sigma
    call sector4_global_checks(i1,i2,i3,i4)
    if(i4.eq.0) then
       num = sigNNLO(i1,i2,i3,i2) &
           + sigNNLO(i1,i3,i3,i2) &
           + sigNNLO(i3,i1,i1,i2) &
           + sigNNLO(i3,i2,i1,i2)
    elseif(i4.ne.0) then
       num = sigNNLO(i1,i2,i3,i4) &
           + sigNNLO(i3,i4,i1,i2)
    else
       write(*,*) 'get_ZSS_NNLO: error in the construction of numerator'
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
          if(a.eq.i1.and.c.eq.i3) sigma = sigma + sigNNLO(a,b,c,b) + sigNNLO(a,c,c,b) + sigNNLO(c,b,a,b) + sigNNLO(c,a,a,b)
!          if(a.eq.i3.and.c.eq.i1) sigma = sigma + sigNNLO(c,b,a,b) + sigNNLO(c,a,a,b)
          if(a.eq.i1.and.b.eq.i3) sigma = sigma + sigNNLO(a,c,b,c) + sigNNLO(a,b,b,c) + sigNNLO(b,c,a,c) + sigNNLO(b,a,a,c)
!          if(a.eq.i3.and.b.eq.i1) sigma = sigma + sigNNLO(b,c,a,c) + sigNNLO(b,a,a,c)
          if(b.eq.i1.and.c.eq.i3) sigma = sigma + sigNNLO(b,a,c,a) + sigNNLO(b,c,c,a) + sigNNLO(c,a,b,a) + sigNNLO(c,b,b,a)
!          if(b.eq.i3.and.c.eq.i1) sigma = sigma + sigNNLO(c,a,b,a) + sigNNLO(c,b,b,a)
       elseif(d.ne.0) then
          if(a.eq.i1.and.c.eq.i3) sigma = sigma + sigNNLO(a,b,c,d)
          if(a.eq.i1.and.d.eq.i3) sigma = sigma + sigNNLO(a,b,d,c)
          if(b.eq.i1.and.c.eq.i3) sigma = sigma + sigNNLO(b,a,c,d)
          if(b.eq.i1.and.d.eq.i3) sigma = sigma + sigNNLO(b,a,d,c)
          if(a.eq.i3.and.c.eq.i1) sigma = sigma + sigNNLO(c,d,a,b)
          if(b.eq.i3.and.c.eq.i1) sigma = sigma + sigNNLO(c,d,b,a)
          if(a.eq.i3.and.d.eq.i1) sigma = sigma + sigNNLO(d,c,a,b)
          if(b.eq.i3.and.d.eq.i1) sigma = sigma + sigNNLO(d,c,b,a)
       else
          write(*,*) 'get_ZSS_NNLO: error in the construction of denominator'
          write(*,*) 'Negative value for 4th sector index i4...'
          write(*,*) 'i4 = ', i4
          write(*,*) 'exit...'
          stop
       endif
    enddo
    ZSS_NNLO = num/sigma
    call sector2_sanity_checks(sigma,ZSS_NNLO)
  end subroutine get_ZSS_NNLO
  
  
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


end module sectors3_module
