      subroutine remtms(ct) 
!
!     Name:  remtms.f90
!     Description:
!       calculate time using etime
!     Parameters:
!        ct=not used
!     
      implicit none 
!                                                                       
      real(8) ct 
      real(4) a(2), etime 
!                                                                       
       ct=etime(a) 
!       ct=0.                                                           
!      write (6,*)'in remtms:',ct,a                                     
!                                                                       
      return 
      end                                           
