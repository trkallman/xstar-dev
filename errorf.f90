      subroutine errorf(dw,wgaus,er) 
!
!     Name: enxt.f90  
!     Description:  
!     This routine calculates the error function
!     Parameters:
!        Input:
!        dw=
!        wgaus=
!        Output:
!        er=error function
!     Called by: none
!     Dependencies: none
!
      implicit none 
      real(8) dw,wgaus,er 
      real(8) aerf1,aerf2,aerf3,te 
!                                                                       
!     calculates the error-function                                     
      aerf1=0.3480242 
      aerf2=-0.0958798 
      aerf3=0.7478556 
      te=1./(1.+0.47047*abs(dw)/wgaus) 
      er=1.-(aerf1*te+aerf2*te**2+aerf3*te**3)*exp(-(dw/wgaus)**2) 
      if (dw.lt.0.) er=0.5-0.5*er 
      if (dw.ge.0.) er=0.5+0.5*er 
      return 
      END                                           
