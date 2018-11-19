      subroutine fact(n,x) 
!                                                                       
! to calculate the factorial of an integer n.  the output x is          
! the natural log of n factorial.                                       
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
