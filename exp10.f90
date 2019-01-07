      function exp10(x) 
!                                                                       
!     Name: exp10.f90  
!     Description:  
!        calculates 10^x
!     Parameters:
!         Input: 
!         x=independent value
!         Output:
!         exp10=10^x
!     Dependencies:  none
!     Called by:  ucalc

      use globaldata
      implicit none 
      real(8) exp10,x 
!                                                                       
      exp10=exp(2.30259*x) 
!                                                                       
      return 
      END                                           
