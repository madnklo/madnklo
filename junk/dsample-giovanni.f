      subroutine setgrid(j,xo,a,itype)
c*************************************************************************
c     Presets the grid for a 1/(x-a)^itype distribution down to xo
c*************************************************************************
      implicit none
c
c     Constants
c
      include 'genps.inc'
c
c     Arguments
c
      integer j, itype                !grid number
      double precision  xo            !minimum value
      double precision  a             !offset for peak
c
c     Local
c
      integer i,k
      integer ngu, ngd
c
c     Global
c
      double precision   grid(2, ng, 0:maxinvar)
      common /data_grid/ grid

      logical            flat_grid
      common/to_readgrid/flat_grid                !Tells if grid read from file

c----- 
c  Begin Code
c-----
      if (flat_grid) then
         if (itype.gt.1) then
            write(*,'(a,i4,2e15.5,i4)') 'Setting grid',j,xo,a,itype
            if (a .ge. xo) then
               write(*,*) 'Can not integrate over singularity'
               write(*,*) 'Set grid',j,xo,a
               return
            endif
         else
            write(*,'(a,i4,1e15.5,i4)') 'Setting grid',j,xo,itype            
         endif
c     grid(2,1,j) = xo
         grid(2,ng,j)=xgmax
         if (itype .eq. 1) then
c
c     We'll use most for the peak, but save some for going down
c
            ngu = ng *0.9
            ngd = ng-ngu

            do i=1,ngu-1
c-------------------
c     tjs 6/30/2009; tjs & ja 2/25/2011
c     New form for setgrid
c-------------------
c               grid(2,i+ngd,j)=((1d0-a)/(xo-a))**(1d0-dble(i)/dble(ngu))
c               grid(2,i+ngd,j)=1d0/grid(2,i+ngd,j)+a
c               grid(2,i+ngd,j) = xo + ((dble(i)+xo-a)/(dble(ngu)+xo-a))**2
               grid(2,i+ngd,j) = xo**(1-dble(i)/dble(ngu))

            enddo
c
c     Now lets go down the other side
c
            grid(2,ngd,j) =  xo
            do i=1,ngd-1
c               grid(2,i,j) = ((1d0-a)/(xo-a))**(1d0-dble(i)/dble(ngd))
               grid(2,ngd-i,j) = xo-(grid(2,ngd+i,j)-xo)
               if (grid(2,ngd-i,j) .lt. -1d0) then
                  write(*,*) 'Error grid set too low',grid(2,ngd-i,j)
                  do k=1,ng
                     write(*,*) k,grid(2,k,j)
                  enddo
                  stop
               endif
            enddo
c
c     tjs, ja 2/25/11
c     Make sure sample all the way down to zero only if minimum positive
c     
            if (grid(2,1,j) .gt. 0) grid(2,1,j) = 0d0
c            write(*,*) "Adjusted bin 1 to zero"

         elseif (itype .eq. 2) then
            do i=2,ng-1
               grid(2,i,j)=(1d0/(xo-a))*(1d0-dble(i)/dble(ng))+
     $              (dble(i)/dble(ng))*(1d0/(1d0-a))
               grid(2,i,j)=1d0/grid(2,i,j)+a
            enddo         
         else
            write(*,*) 'No modification in setgrid',itype
         endif
         do i=1,ng
c             write(*,*) j,i,grid(2,i,j)
         enddo
         call sample_write_g(j,'_0')
      else
         write(*,*) 'No modification is setgrid, grid read from file'
      endif
      end

      subroutine sample_write_g(idim,cpost)
c**************************************************************************
c     Writes out grid in function form for dimension i with extension cpost
c     
c**************************************************************************
      implicit none
c
c     Constants
c
      include 'genps.inc'
c
c     Arguments
c
      integer idim
      character*(*) cpost
c
c     Local
c
      character*60 fname
      integer i
      double precision xo,yo
c
c     Global
c
      double precision   grid(2, ng, 0:maxinvar)
      common /data_grid/ grid

c-----
c  Begin Code
c-----
      return
      if (idim .lt. 1 .or. idim .gt.maxinvar) then
         write(*,*) 'Error invalid dimension in sample_write_f',idim
         return
      endif
      if (idim .lt. 10) then
         write(fname,'(a,i1,a,a)') 'g_',idim,cpost,'.dat'
      elseif (idim .lt. 100) then
         write(fname,'(a,i2,a,a)') 'g_',idim,cpost,'.dat'
      endif
      open(unit=21,file=fname,status='unknown',err=99)
      do i=1,ng-1
         xo = (grid(2,i,idim)+grid(2,i+1,idim))/2d0
         yo =1d0/(-grid(2,i,idim)+grid(2,i+1,idim))
         write(21,*) xo,yo
      enddo
      close(21)
      return
 99   write(*,*) 'Error opening file ',fname
      end
