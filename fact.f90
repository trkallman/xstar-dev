      subroutine fact(n,x) 
!                                                                       
!     Name: dfact.f90  
!     Description:  
!       to calculate the factorial of an integer n.  the output x is          
!       the natural log of n factorial, in double precision.                  
!       nb redundamt.  same as dfact
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
      real(8) x 
      integer i,n 
!                                                                       
      x=0. 
      if(n.ne.0) then 
      do  i=1,n 
        x=x+log(float(i)) 
        enddo 
      endif 
!                                                                       
      return 
      END                                           
