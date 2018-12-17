      subroutine xwrite(string,ll) 
      implicit none 
      integer ll 
      character*120 string 
!                                                                       
!     Not used                                                          
      integer javi 
      javi=ll 
!      ll=javi                                                          
!                                                                       
      write (6,*)string 
      return 
      end                                           
