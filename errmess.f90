      subroutine errmess(lun11,nlen,str1) 
!
!     Name: errmess.f90  
!     Description:  
!     This routine prints an error message
!     author: apec
!     Parameters:
!           input:
!           lun11=logical unit for printing
!           nlen=length of message string
!           str1(nlen)=string 
!     Dependencies:  none
!     Called by: calc_maxwell_rates, exintn
!
      character*(*) str1 
      integer lun11,nlen,ntmp
      ntmp=nlen
      write (lun11,*) 'in errmess:',str1 
      return 
      END                                           
