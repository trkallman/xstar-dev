      real(8) function pow(x,y) 
!
!     Name: pow.f90  
!     Description:  
!       this function calculates exponentiation (x^y) with limits
!     Parmameters:
!          x=mantissa
!          y=exponent
!          Output:
!          pow=x^y
!     Dependencies:  expo
!     Called by: calc_maxwell_rates

      real(8) x,y,expo 
!      pow=x**y                                                         
      pow=expo(y*log(max(1.d-49,x))) 
      return 
      end                                           
