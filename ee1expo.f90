      real(8) function ee1expo(x) 
!                                                                       
!     Name: ee1expo.f90  
!     Description:  
!       this function computes the first exponential integral.             
!       same as ee1exp but with limit on exponential.
!     Parmameters:
!          x=independent variable
!          Output:
!          ee1exp=first exponential integral
!     Dependencies:  expo
!     Called by: ucalc     
!                                                                       
!     this routine computes the first exponential integral.             
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      real(8) x,expo 
!                                                                       
      if ( x.ge.1. ) then 
      ee1expo=(1./x)*(0.250621+x*(2.334733+x))/(1.68153+x*(3.330657+x)) 
      return 
      endif 
!                                                                       
      ee1expo = (-log(x)-0.57721566+                                    &
     &      x*(0.99999193+x*(-0.24991055+x*(0.05519968+                 &
     &      x*(-0.00976004+x*0.0010707857)))))*expo(x)                  
!                                                                       
      return 
      END                                           
