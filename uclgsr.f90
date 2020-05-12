      subroutine uclgsr(kdum,result,ierr)
      implicit none
      integer ierr
      real(4) result
      character*(*) kdum
!
      ierr=0
!
      read (5,*)result
!      write (6,*)'in uclgsr, result=',result
      return
      end
