      subroutine uclgsr8(kdum,result,ierr) 
!
!     Name:  uclgsr8.f90
!     Description:
!        Read in a real*8 quantity from input string
!        mod of a xanlib canned routine
!        T. Kallman
!     Parameters:
!        Input:
!        kdum=input string
!        Output:
!        result=quantity 
!        ierr=error flag
!     Dependencies:  uclgsr
!     Called by: rread1
!
      implicit none 
      integer ierr 
      real(8) result 
      real(4) result4 
      character*(*) kdum 
!                                                                       
      ierr=0 
!                                                                       
      call uclgsr(kdum,result4,ierr) 
      result=result4 
!                                                                       
      return 
      end                                           
