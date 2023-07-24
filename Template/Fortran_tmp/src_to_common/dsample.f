      subroutine sample_get_x(wgt, x, j, ipole, xmin, xmax)
c************************************************************************
c     Returns maxdim random numbers between 0 and 1, and the wgt
c     associated with this set of points, and the iteration number
c     This routine chooses the point within the range specified by
c     xmin and xmax for dimension j in configuration ipole
c************************************************************************
      implicit none
c
c     Constants
c
      include 'genps.inc'
      include 'maxconfigs.inc'
c
c     Arguments
c
      double precision wgt, x, xmin, xmax
      integer j, ipole
c
c     Local
c
      integer  im, ip,ij,icount,it_warned
      double precision xbin_min,xbin_max,ddum(maxdim),xo,y
c
c     External
c
      double precision xbin
      external         xbin
c
c     Global
c
      double precision tmean, trmean, tsigma
      integer             dim, events, itm, kn, cur_it, invar, configs
      common /sample_common/
     .     tmean, trmean, tsigma, dim, events, itm, kn, cur_it, invar, configs

      double precision    grid(2, ng, 0:maxinvar)
      common /data_grid/ grid
      integer           Minvar(maxdim,lmaxconfigs)
      common /to_invar/ Minvar

      integer           ituple
      common /to_random/ituple

      double precision      spole(maxinvar),swidth(maxinvar),bwjac
      common/to_brietwigner/spole        ,swidth        ,bwjac

      integer nzoom
      double precision  tx(1:3,maxinvar)
      common/to_xpoints/tx, nzoom

      data ddum/maxdim*0d0/
      data icount/0/
      data it_warned/0/

      integer            lastbin(maxdim)
      common /to_lastbin/lastbin

c-----
c  Begin Code
c-----


      ituple=1                  !DEBUG


      
c$$$      if (it_warned .ne. cur_it) then
c$$$         icount=0
c$$$         it_warned = cur_it
c$$$      endif
c$$$      if (ituple .eq. 2) then   !Sobel generator
c$$$         print*,'Sorry Sobel generator disabled'
c$$$         stop
c$$$c         call sobel(ddum)
c$$$c         write(*,'(7f11.5)')(ddum(j)*real(ng),j=1,dim)
c$$$      endif
c$$$      if (ituple .eq. 1) then
c$$$c         write(*,*) 'Getting variable',ipole,j,minvar(j,ipole)
c$$$         xbin_min = xbin(xmin,minvar(j,ipole))
c$$$         xbin_max = xbin(xmax,minvar(j,ipole))
c$$$         if (xbin_min .gt. xbin_max-1) then
c$$$c            write(*,'(a,4e15.4)') 'Bad limits',xbin_min,xbin_max,
c$$$c     &           xmin,xmax
c$$$c            xbin_max=xbin_min+1d-10
c$$$            xbin_max = xbin(xmax,minvar(j,ipole))
c$$$            xbin_min = min(xbin(xmin,minvar(j,ipole)), xbin_max)
c$$$         endif

      
c$$$c
c$$$c     Line which allows us to keep choosing same x
c$$$c
c$$$
c$$$
c$$$
c$$$c         if (swidth(j) .ge. 0) then
c$$$         if (nzoom .le. 0) then
c$$$            call ntuple(ddum(j), xbin_min,xbin_max, j, ipole)
c$$$         else
c$$$c            write(*,*) 'Reusing num',j,nzoom,tx(2,j)
c$$$
c$$$            call ntuple(ddum(j),max(xbin_min,dble(int(tx(2,j)))),
c$$$     $           min(xbin_max,dble(int(tx(2,j))+1)),j,ipole)
c$$$
c$$$            if(max(xbin_min,dble(int(tx(2,j)))).gt.
c$$$     $           min(xbin_max,dble(int(tx(2,j))+1))) then
c$$$c               write(*,*) 'not good'
c$$$            endif
c$$$
c$$$c            write(*,'(2i6,4e15.5)') nzoom,j,ddum(j),tx(2,j),
c$$$c     $           max(xbin_min,dble(int(tx(2,j)))),
c$$$c     $           min(xbin_max,dble(int(tx(2,j))+1))
c$$$
c$$$c            ddum(j) = tx(2,j)                 !Use last value
c$$$
c$$$
c$$$         endif
c$$$         tx(1,j) = xbin_min
c$$$         tx(2,j) = ddum(j)
c$$$         tx(3,j) = xbin_max
c$$$      elseif (ituple .eq. 2) then
c$$$         if (ipole .gt. 1) then
c$$$            print*,'Sorry Sobel not configured for multi-pole.'
c$$$            stop
c$$$         endif
c$$$         ddum(j)=ddum(j)*dble(ng)
c$$$      else
c$$$         print*,'Error unknown random number generator.',ituple
c$$$         stop
c$$$      endif
c$$$
c$$$      im = ddum(j)
c$$$      if (im.ge.ng)then
c$$$         im = ng -1
c$$$         ddum(j) = ng
c$$$      endif
c$$$      if (im.lt.0) im = 0
c$$$      ip = im + 1
      ij = Minvar(j,ipole)
c------
c     tjs 3/5/2011  save bin used to avoid looking up when storing wgt
c------
c$$$      lastbin(j) = ip
c$$$c
c$$$c     New method of choosing x from bins
c$$$c
c$$$      if (ip .eq. 1) then         !This is in the first bin
c$$$         xo = grid(2, ip, ij)-xgmin
c$$$         x = grid(2, ip, ij) - xo * (dble(ip) - ddum(j))
c$$$      else           
c$$$         xo = grid(2, ip, ij)-grid(2,im,ij)
c$$$         x = grid(2, ip, ij) - xo * (dble(ip) - ddum(j))
c$$$      endif
c
c     Now we transform x if there is a B.W., S, or T  pole
c
      if (ij .gt. 0) then
c         write(*,*) "pole, width",ij,spole(ij),swidth(ij)
         if (swidth(ij) .gt. 0d0) then
c            write(*,*) 'Tranpole called',ij,swidth(ij)
            y = x                             !Takes uniform y and returns
            call transpole(spole(ij),swidth(ij),y,x,wgt) !x on BW pole or 1/x 
         endif
      endif
c
c     Simple checks to see if we got the right point note 1e-3 corresponds
c     to the fact that the grids are required to be separated by 1e-14. Since
c     double precision is about 18 digits, we expect things to agree to
c     3 digit accuracy.
c
c$$$      if (abs(ddum(j)-xbin(x,ij))/(ddum(j)+1d-22) .gt. 1e-3) then
c$$$         if (icount .lt. 5) then
c$$$            write(*,'(a,i4,2e14.6,1e12.4)')
c$$$     &           'Warning xbin not returning correct x', ij,
c$$$     &           ddum(j),xbin(x,ij),xo
c$$$         elseif (icount .eq. 5) then
c$$$            write(*,'(a,a)')'Warning xbin still not working well. ',
c$$$     &           'Last message this iteration.'
c$$$         endif
c$$$         icount=icount+1
c$$$      endif
c$$$      if (x .lt. xmin .or. x .gt. xmax) then
c$$$c         write(*,'(a,4i4,2f24.16,1e10.2)') 'Bad x',ij,int(xbin_min),ip,
c$$$c     &        int(xbin_max),xmin,x,xmax-xmin
c$$$      endif
c$$$
c$$$      wgt = wgt * xo * dble(xbin_max-xbin_min)
c      print*,'Returning x',ij,ipole,j,x
      end


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

      flat_grid=.true.
            
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
c         call sample_write_g(j,'_0')
      else
         write(*,*) 'No modification is setgrid, grid read from file'
      endif
      end


      double precision function xbin(y,j)
c**************************************************************************
c     Subroutine to determine which value y  will map to give you the
c     value of x when put through grid j.  That is what random number
c     do you need to be given to get the value x out of grid j and will be
c     between 0 < x < ng.
c**************************************************************************
      implicit none
c
c     Constants
c
      include 'genps.inc'
      double precision tol
      parameter       (tol=1d-12)
c
c     Arguments
c      
      double precision y
      integer j
c
c     Local
c
      integer i,jl,ju
      double precision x,xo
c
c     Global
c
      double precision    grid(2, ng, 0:maxinvar)
      common /data_grid/ grid
      double precision      spole(maxinvar),swidth(maxinvar),bwjac
      common/to_brietwigner/spole        ,swidth        ,bwjac
c
c     Data
c
      data spole,swidth/maxinvar*0d0,maxinvar*0d0/
c-----
c  Begin Code
c-----

      
      bwjac = 1d0
      if (j .gt. 0) then
         if (swidth(j) .gt. 0d0) then
            call  untranspole(spole(j),swidth(j),x,y,bwjac)
         else
            x=y
         endif
      else
         x=y
      endif
      if (x .eq. xgmax) then
         i=ng
         xbin = dble(ng)
      elseif (x .eq. xgmin) then
         xbin=0d0
      elseif(x .le. grid(2,1,j)) then
         i=1
         xo = grid(2,i,j)-xgmin
         xbin = dble(i)+(x-grid(2,i,j))/xo
      else
         jl = 1
         ju = ng
         do while (ju-jl .gt. 1)                    !Binary search
            i = (ju-jl)/2+jl
            if (grid(2,i,j) .le. x) then
               jl=i
            else
               ju=i
            endif
         enddo
         i=ju
         xo = grid(2,i,j)-grid(2,i-1,j)
         xbin = dble(i)+(x-grid(2,i,j))/xo
      endif
c      jbin=i
c      x = 
c      if (x+tol .gt. grid(2,i,j) .and. i .ne. ng) then
c         write(*,'(a,2e23.16,e9.2)') 'Warning in DSAMPLE:JBIN ',
c     &                x,grid(2,i,j),tol
c         x=2d0*grid(2,i,j)-x
c         jbin=i+1
c      endif
      end
