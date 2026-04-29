module sectors4_module
  implicit none
  integer, public :: n_ext
  double precision, public, parameter :: alpha_mod=2d0
  double precision, public, parameter :: alpha_mod_bar=1d0
  double precision, public :: W_NNLO, WSS_NNLO, WCC_NNLO, WSS_CC_NNLO, WSC_NNLO, WSS_SC_NNLO, WSS_SC_CC_NNLO, WSC_CC_NNLO
  double precision, public :: Wbar_NLO, WS_NLO, WSbar_NLO, WC_NLO, WCbar_NLO
  double precision, allocatable, dimension(:,:), public :: xs_mod
  double precision, allocatable, dimension(:,:), public :: w,sig2,hatsig2
  double precision, allocatable, dimension(:,:,:,:), public :: sigNNLO,hatsigNNLO
  public :: get_sigNNLO, get_hatsigNNLO, get_W_NNLO, get_WSS_NNLO, get_WCC_NNLO, get_WSS_CC_NNLO
  public :: get_WSC_NNLO, get_WSS_SC_NNLO, get_WSS_SC_CC_NNLO, get_WSC_CC_NNLO
  public :: get_w, get_sig2, get_Wbar_NLO, get_WS_NLO, get_WSbar_NLO, get_WC_NLO, get_WCbar_NLO
  private

contains

  subroutine get_sigNNLO(xs_in,n_ext_in)
    implicit none
    ! global
    include 'nexternal.inc'
    integer :: n_ext_in
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

  subroutine get_hatsigNNLO(iref_in,xs_in,n_ext_in)
    implicit none
    ! global
    include 'nexternal.inc'
    integer :: n_ext_in, iref_in
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
  end subroutine get_W_NNLO

  subroutine get_WSS_NNLO(a,b,c,d)
    ! NNLO double-soft sector functions WSS
    implicit none
    integer :: i,a,b,c,d,sec(4)
    double precision :: num,sigma
    include 'all_K2_sector_list.inc'
    num = sigNNLO(a,b,c,d)
    sigma = 0d0
    do i=1,len
       sec=ss_sector_list(a,c,i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma +  sigNNLO(sec(1),sec(2),sec(3),sec(4))
    enddo
    WSS_NNLO = num/sigma
  end subroutine get_WSS_NNLO

  subroutine get_WCC_NNLO(a,b,ic,id)
    ! NNLO triple-collinear sector functions WCC
    implicit none
    integer :: i,ii(3),a,b,c,d,ic,id,sec(4)
    double precision :: num, sigma
    include 'all_K2_sector_list.inc'
    num = hatsigNNLO(a,b,ic,id)
    sigma = 0d0
    if(b.eq.ic) then
       c = id
    elseif(b.eq.id) then
       c = ic
    endif

    ! list-reading requires sorted a,b,c
    ii = [a, b, c]
    if (ii(1)>ii(2)) call swap(ii(1),ii(2))
    if (ii(1)>ii(3)) call swap(ii(1),ii(3))
    if (ii(2)>ii(3)) call swap(ii(2),ii(3))
    do i=1,len
       sec=cc_sector_list(ii(1),ii(2),ii(3),i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma + hatsigNNLO(sec(1),sec(2),sec(3),sec(4))
    enddo
    wcc_nnlo=num/sigma
  end subroutine get_WCC_NNLO

  subroutine get_WSS_CC_NNLO(a,b,c,d)
    ! NNLO double-soft triple-collinear sector functions WSSCC
    implicit none
    integer :: i,a,b,c,d,ib,sec(4)
    double precision :: num, sigma
    include 'all_K2_sector_list.inc'
    num = hatsigNNLO(a,b,c,d)
    sigma = 0d0
    if(b.eq.c) then
       ib = b
    elseif(b.eq.d) then
       ib = c
    endif
    do i=1,len
       sec=ss_cc_sector_list(a,ib,d,i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma + &
            hatsigNNLO(sec(1),sec(2),sec(3),sec(4))
    enddo
    wss_cc_nnlo=num/sigma
  end subroutine get_WSS_CC_NNLO

  subroutine get_w(xs_in,n_ext_in)
    implicit none
    ! global
    integer :: n_ext_in
    double precision, dimension (n_ext_in,n_ext_in) :: xs_in
    ! local
    integer :: i,j
    ! set global module variables
    n_ext=n_ext_in
    if (.not.allocated(xs_mod)) allocate(xs_mod(n_ext,n_ext))
    if (.not.allocated(w)) allocate(w(3:n_ext,3:n_ext))
    xs_mod=xs_in
    ! calculate 2-index w
    w=0d0
    do i=3,n_ext
       do j=3,n_ext
          if(i.eq.j)cycle
          if((xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2)).ne.0d0)then
             w(i,j)=xs_mod(1,2)*xs_mod(i,j)/(xs_mod(i,1)+xs_mod(i,2))/(xs_mod(j,1)+xs_mod(j,2))
          endif
       enddo
    enddo
  end subroutine get_w

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
    ! calculate 2-index sigma
    sig2=0d0
    do i=3,n_ext
       do j=3,n_ext
          if(i.eq.j)cycle
          if((xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2))*xs_mod(1,2).ne.0d0)then
             ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
             ej=(xs_mod(j,1)+xs_mod(j,2))/xs_mod(1,2)
             wij=xs_mod(1,2)*xs_mod(i,j)/(xs_mod(i,1)+xs_mod(i,2))/(xs_mod(j,1)+xs_mod(j,2))
             sig2(i,j)=(1d0/ei/wij)**alpha_in
          endif
       enddo
    enddo
  end subroutine get_sig2

!   double precision function hatsig4(i,j,k,l,r,xs_in,alpha_in,n_ext_in)
!     implicit none
!     ! global
!     integer :: n_ext_in
!     double precision :: alpha_in
!     double precision, dimension (n_ext_in,n_ext_in) :: xs_in
!     ! local
!     integer :: i,j,k,l,r,n_ext
!     double precision :: ei,ek,wij,wkl,wir,wkr
!     ! set global module variables
!     n_ext=n_ext_in
!     if (.not.allocated(xs_mod)) allocate(xs_mod(n_ext,n_ext))
!     xs_mod=xs_in
!     ! calculate (sig[i,j]/wir)^alpha*sig[k,l]/wkr where l=i,j,k for the NNLO SC sector functions
!     hatsig4=0d0
!     if((xs_mod(i,1)+xs_mod(i,2))*(xs_mod(j,1)+xs_mod(j,2))*(xs_mod(k,1)+xs_mod(k,2))* &
!        (xs_mod(r,1)+xs_mod(r,2))*xs_mod(1,2).ne.0d0)then
!         ei=(xs_mod(i,1)+xs_mod(i,2))/xs_mod(1,2)
!         ek=(xs_mod(j,1)+xs_mod(j,2))/xs_mod(1,2)
!         wij=xs_mod(1,2)*xs_mod(i,j)/(xs_mod(i,1)+xs_mod(i,2))/(xs_mod(j,1)+xs_mod(j,2))
!         wkl=xs_mod(1,2)*xs_mod(k,l)/(xs_mod(l,1)+xs_mod(l,2))/(xs_mod(k,1)+xs_mod(k,2))
!         wir=xs_mod(1,2)*xs_mod(i,r)/(xs_mod(i,1)+xs_mod(i,2))/(xs_mod(r,1)+xs_mod(r,2))
!         wkr=xs_mod(1,2)*xs_mod(k,r)/(xs_mod(k,1)+xs_mod(k,2))/(xs_mod(r,1)+xs_mod(r,2))
!         hatsig4=((1d0/ei/wij)/wir)**alpha_in*(1d0/ek/wkl)/wkr
!     endif
!   end

  subroutine get_Wbar_NLO(i1,i2)
    !     NLO sector functions W(i1,i2)
    implicit none
    integer :: i,j,a,b,i1,i2,sec1,sec2
    double precision :: num,sigma
    include 'all_sector_list.inc'
    num = sig2(i1,i2)
    sigma = 0d0
    j=1
    do i=1,lensectors
       sec1=bar_indices(j)
       sec2=bar_indices(j+1)
       j=j+2
       if(sec1.eq.0.and.sec2.eq.0) cycle
       sigma = sigma + sig2(sec1,sec2)
    enddo
    Wbar_NLO = num/sigma
  end subroutine get_Wbar_NLO

  subroutine get_WS_NLO(i1,i2)
    !     NLO soft sector functions WS(i1,i2) = barS_i1 W(i1,i2)
    implicit none
    integer :: i,i1,i2,sec(2)
    double precision :: num,sigma
    include 'all_K1_sector_list.inc'
    num = sig2(i1,i2)
    sigma = 0d0
    do i=1,len
       sec=s_sector_list(i1,i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma + &
            sig2(sec(1),sec(2))
    enddo
    WS_NLO = num/sigma
  end subroutine get_WS_NLO

  subroutine get_WSbar_NLO(i1,i2)
    !     NLO sector functions WSbar(i1,i2)
    implicit none
    integer :: i,j,a,b,i1,i2,sec1,sec2
    double precision :: num,sigma
    include 'all_sector_list.inc'
    num = sig2(i1,i2)
    sigma = 0d0
    j=1
    do i=1,lensectors
       sec1=bar_indices(j)
       sec2=bar_indices(j+1)
       j=j+2
       if(sec1.ne.i1) cycle
       sigma = sigma + sig2(sec1,sec2)
    enddo
    WSbar_NLO = num/sigma
  end subroutine get_WSbar_NLO

  subroutine get_WC_NLO(i1,i2,i3,ir)
    !     NLO collinear sector functions WC(i1,i2) = barC_i1i2 W(i1,i2)
    implicit none
    integer :: i,i1,i2,i3,ir,sec(2)
    double precision :: num,sigma
    include 'all_K1_sector_list.inc'
    num = sig2(i1,ir)
    sigma = 0d0
    do i=1,len
       sec=c_sector_list(i1,i2,i3,i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma + sig2(sec(1),sec(2))
    enddo
    WC_NLO = num/sigma
  end subroutine get_WC_NLO

  subroutine get_WCbar_NLO(i1,i2,ir)
    !     NLO collinear sector functions WCbar(i1,i2) = barC_i1i2 Wbar(i1,i2)
    implicit none
    integer :: i,i1,i2,i3,ir,sec(2)
    double precision :: num,sigma
    num = sig2(i1,ir)
    sigma = sig2(i1,ir) + sig2(i2,ir)
    WCbar_NLO = num/sigma
  end subroutine get_WCbar_NLO

  subroutine get_WSC_NNLO(i,j,k,l,ir)
    implicit none
    integer :: i,j,k,l,ir
    integer :: legs(4),m,b,d
    double precision :: num,sigma
    double precision :: sum_b_T1,sum_d_T2,sum_d_T3
    num=sig2(i,j)**alpha_mod*sig2(k,l)/w(k,ir)
    ! by construction of legs none of the indices sum over i
    legs=[j,k,l,ir]
    ! sig2(i,b)**alpha
    sum_b_T1=0d0
    do m=1,4
       b=legs(m)
       if(any(legs(1:m-1)==b))cycle
       sum_b_T1=sum_b_T1+sig2(i,b)**alpha_mod
    enddo
    ! sig2(i,d) where d!=k
    sum_d_T2=0d0
    do m=1,4
       d=legs(m)
       if(d.eq.k)cycle
       ! let's say in 3445 the legs will be [4,4,5,6], so now we want to avoid double-summing the index
       ! should/might work for 3454, legs will be [4,5,4,6]
       if(any(legs(1:m-1)==d))cycle
       sum_d_T2=sum_d_T2+sig2(i,d)
    enddo
    ! sig2(i,d) where d!=l
    sum_d_T3=0d0
    do m=1,4
       d=legs(m)
       if(d.eq.l)cycle
       if(any(legs(1:m-1)==d))cycle
       sum_d_T3=sum_d_T3+sig2(i,d)
    enddo
    sigma=sum_b_T1*(sig2(k,l)/w(k,ir)+sig2(l,k)/w(l,ir))+sig2(k,l)**alpha_mod/w(k,ir)*sum_d_T2+sig2(l,k)**alpha_mod/w(l,ir)*sum_d_T3
    WSC_NNLO=num/sigma
  end subroutine get_WSC_NNLO

  subroutine get_WSS_SC_NNLO(a,b,c,d)
    !     NNLO double-soft soft-collinear sector functions WSS_SC(a,b,c,d)
    implicit none
    integer :: i,a,b,c,d,ic,sec(4)
    double precision :: num,sigma
    include 'all_K2_sector_list.inc'
    num = sig2(a,b)**alpha_mod*sig2(c,d)
    sigma=0d0
    if(b.eq.c) then
        ic=b
    elseif(b.eq.d)then
        ic=c
    endif
    do i=1,len
       sec=ss_sc_sector_list(a,ic,d,i,:)
       if(all(sec.eq.0))cycle
       sigma = sigma + sig2(sec(1),sec(2))**alpha_mod*sig2(sec(3),sec(4))
    enddo
    WSS_SC_NNLO = num/sigma
  end subroutine get_WSS_SC_NNLO

  subroutine get_WSC_CC_NNLO(a,b,c,d,ir)
    !     NNLO soft-collinear double-collinear sector functions WSC_CC(a,b,c,d)
    implicit none
    integer :: a,b,c,d,ir,ic
    double precision :: num,sigma
    num = (sig2(a,b)/w(a,ir))**alpha_mod*sig2(c,d)/w(c,ir)
    if(b.eq.c) then
        ic=b
    elseif(b.eq.d)then
        ic=c
    endif
    sigma = (sig2(a,ic)/w(a,ir))**alpha_mod*sig2(ic,d)/w(ic,ir) + (sig2(a,ic)/w(a,ir))**alpha_mod*sig2(d,ic)/w(d,ir) &
          + (sig2(a,d)/w(a,ir))**alpha_mod*sig2(ic,d)/w(ic,ir)  + (sig2(a,d)/w(a,ir))**alpha_mod*sig2(d,ic)/w(d,ir)  &
          + (sig2(ic,d)/w(ic,ir))**alpha_mod*sig2(a,d)/w(a,ir)  + (sig2(d,ic)/w(d,ir))**alpha_mod*sig2(a,ic)/w(a,ir)
    WSC_CC_NNLO = num/sigma
  end subroutine get_WSC_CC_NNLO

  subroutine get_WSS_SC_CC_NNLO(a,b,c,d,ir)
    !     NNLO double-soft soft-collinear double-collinear sector functions WSS_SC_CC(a,b,c,d)
    implicit none
    integer :: a,b,c,d,ir,ic
    double precision :: num,sigma
    num = (sig2(a,b)/w(a,ir))**alpha_mod*sig2(c,d)/w(c,ir)
    if(b.eq.c) then
        ic=b
    elseif(b.eq.d)then
        ic=c
    endif
    sigma =  (sig2(a,ic)/w(a,ir))**alpha_mod*sig2(ic,d)/w(ic,ir) + (sig2(a,d)/w(a,ir))**alpha_mod*sig2(ic,d)/w(ic,ir) &
          +  (sig2(ic,d)/w(ic,ir))**alpha_mod*sig2(a,d)/w(a,ir)
    WSS_SC_CC_NNLO = num/sigma
  end subroutine get_WSS_SC_CC_NNLO



  subroutine swap(x,y)
  integer :: x,y,t
  t=x; x=y; y=t
  end subroutine swap

end module sectors4_module
