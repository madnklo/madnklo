      subroutine ct_log(label,val)
      implicit none
      character*(*) label
      double precision val
      logical consistency_check
      common/cconscheck/consistency_check
      if(consistency_check) write(88,*) label, val
      end
