      subroutine dfact(n,x) 
!                                                                       
!     Name: dfact.f90  
!     Description:  
!       to calculate the factorial of an integer n.  the output x is          
!       the natural log of n factorial, in double precision.                  
!       author:  T. Kallman                                               
!     Parameters:
!        Input:
!        n=integer
!        Output:
!        x=output
!     Dependencies: none
!     Called by:  anl1
!                                                                       

      implicit none 
      integer n,i 
      real(8) x 
!                                                                       
      x=0. 
      if(n.eq.0) return 
      do i=1,n 
       x=x+log(float(i)) 
       enddo 
!                                                                       
      return 
      END                                           
