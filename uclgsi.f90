      subroutine uclgsi(kdum,iresult,ierr)
      implicit none
      integer iresult, ierr
      character*(*) kdum
!
      ierr=0
!
      read (5,*)iresult
!      write (6,*)'in uclgsi, iresult=',iresult
      return
      end
