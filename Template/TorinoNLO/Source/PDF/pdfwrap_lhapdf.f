      subroutine pdfwrap
c     MadNkLO, hadron beams: initialise the LHAPDF set lhaid (run_card)
c     and take alpha_s(MZ) and the order of its running from the set.
c     The PDF sets are looked for in the directory PDFsets found by
c     FindPDFPath (bin/run_all.py links <output>/lib/PDFsets to the
c     data directory of the LHAPDF given in input/madnklo_config.txt),
c     then in the default paths of LHAPDF.
      implicit none
      include 'pdf.inc'
      include '../alfas.inc'
      real*8 zmass
      data zmass/91.188d0/
      character*150 LHAPath
      character*20 parm(20)
      double precision value(20)
      real*8 alphasPDF
      external alphasPDF
c
      call FindPDFPath(LHAPath)
      call SetPDFPath(LHAPath)
      value(1)=lhaid
      parm(1)='DEFAULT'
      call pdfset(parm,value)
      call GetOrderAs(nloop)
      nloop=nloop+1
      asmz=alphasPDF(zmass)
c
      return
      end


      subroutine FindPDFPath(LHAPath)
c     Directory PDFsets in the current directory or in lib/ of one of
c     the parent directories (up to 6 levels)
      implicit none
      character LHAPath*150,up*3
      data up/'../'/
      logical exists
      integer i
c
      LHAPath='./PDFsets'
      inquire(file=LHAPath, exist=exists)
      if(exists)return
      LHAPath='lib/PDFsets'
      inquire(file=LHAPath, exist=exists)
      if(exists)return
      do i=1,6
         LHAPath=up//LHAPath
         inquire(file=LHAPath, exist=exists)
         if(exists)return
      enddo
      write(*,*)'pdfwrap: no PDFsets directory found in lib/ of the ',
     &     'output (it is created by bin/run_all.py from ',
     &     'input/madnklo_config.txt)'
      stop 1
c
      return
      end
