      real*8 function voigtedoppler(vs,a) 
!                                                                       
!     fake version does doppler                                         
!     the integral of this function over vs (-infinity -> infinity) is s
!                                                                       
      real*8 a,vs,expo 
!                                                                       
      voigtedoppler=expo(-vs*vs) 
!                                                                       
      return 
      end                                           
