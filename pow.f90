      real(8) function pow(x,y) 
      real(8) x,y,expo 
!      pow=x**y                                                         
      pow=expo(y*log(max(1.e-34,x))) 
      return 
      end                                           
