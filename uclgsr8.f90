      subroutine uclgsr8(kdum,result,ierr) 
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
