      subroutine ct_log(label,val)
      implicit none
      character*(*) label
      double precision val
      logical consistency_check
      common/cconscheck/consistency_check
      character*8 ct_limit
      character*40 ct_tag
      double precision ct_lam
      integer iu_ct
      parameter(iu_ct=88)
      common/cctlog/ct_lam,ct_limit,ct_tag
      if(.not.consistency_check)return
      write(iu_ct,'(a,1x,es23.15e3,1x,a,1x,es23.15e3)')
     &     trim(ct_limit),ct_lam,trim(label),val
      end


      subroutine ct_log_begin(limstr,tag)
      implicit none
      character*(*) limstr,tag
      character*8 ct_limit
      character*40 ct_tag
      double precision ct_lam
      common/cctlog/ct_lam,ct_limit,ct_tag
      ct_limit=limstr
      ct_tag=tag
      ct_lam=0d0
      end


      subroutine ct_log_lam(lam)
      implicit none
      double precision lam
      logical consistency_check
      common/cconscheck/consistency_check
      character*8 ct_limit
      character*40 ct_tag
      double precision ct_lam
      integer iu_ct
      parameter(iu_ct=88)
      logical isopen,seen
      integer i,nseen,mxseen
      parameter(mxseen=64)
      character*60 fname,seen_files(mxseen)
      save nseen,seen_files
      data nseen/0/
      common/cctlog/ct_lam,ct_limit,ct_tag
      ct_lam=lam
      if(.not.consistency_check)return
      inquire(unit=iu_ct,opened=isopen)
      if(.not.isopen)then
         fname='consistency_'//trim(ct_limit)//'_'//trim(ct_tag)//'.log'
c        the first opening in this run replaces an old file; the later
c        ones (one per event) append to it
         seen=.false.
         do i=1,nseen
            if(seen_files(i).eq.fname)seen=.true.
         enddo
         if(seen)then
            open(unit=iu_ct,file=fname,status='old',position='append')
         else
            open(unit=iu_ct,file=fname,status='replace')
            if(nseen.lt.mxseen)then
               nseen=nseen+1
               seen_files(nseen)=fname
            endif
         endif
         write(iu_ct,'(a)')'# limit lambda label value'
      endif
      end


      subroutine ct_log_end()
      implicit none
      integer iu_ct
      parameter(iu_ct=88)
      logical isopen
      inquire(unit=iu_ct,opened=isopen)
      if(isopen)close(iu_ct)
      end
