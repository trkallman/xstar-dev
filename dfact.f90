      subroutine dfact(n,x) 
!                                                                       
! to calculate the factorial of an integer n.  the output x is          
! the natural log of n factorial, in double precision.                  
!     author:  T. Kallman                                               
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
