!      real*8 function voigte(vs,a)                                     
      real*8 function voigtelorentzian(vs,a) 
!                                                                       
!     fake version does lorentzian                                      
!                                                                       
!     the integral of this function over vs (-infinity -> infinity) is 1
!                                                                       
      real*8 a,vs 
!                                                                       
      voigtelorentzian=a/(vs*vs+a*a)/3.14 
!                                                                       
      Return 
      END                                           
