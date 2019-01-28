      real(8) function pescl(tau) 
!                                                                       
!     Name:  pescl.f90
!     Description:
!       this routine calculates escape probability for a line transition  
!       inputs: optical depths-- tau for a line transition                
!       Currently using expression from Kwan and Krolik 1982
!     Parameters:
!       Input: 
!       tau=line center optical depth
!       output:
!       pescl=escape probability
!     Dependencies:  none
!     Called by: ucalc
!                                                                       
      implicit none 
!                                                                       
      real(8) tau 
      real(8) pi,tauw,aa,bb 
!                                                                       
      data pi/3.1415927/ 
!                                                                       
      tauw=1.e5 
!     *** need to determine tauw from line profiles?***                 
      if(tau.lt.1.0) then 
        if(tau.lt.1.e-5) then 
          pescl=1.0 
          go to 10 
        end if 
        aa=2.0*tau 
        pescl=(1.0-exp(-aa))/aa 
        go to 10 
      end if 
      bb=0.5*sqrt(log(tau))/(1.0+tau/tauw) 
      pescl=1./(tau*sqrt(pi)*(1.2+bb)) 
   10 continue 
      pescl=pescl/2.0 
!                                                                       
      return 
      end                                           
