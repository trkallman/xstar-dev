      subroutine fitsclose(lun11,unit,status) 
!                                                                       
!     Name: fitsclose.f90  
!     Description:  
!        Close the file & release the unit number                          
!        author: T. Bridgman                                               
!
!     List of Parameters:
!           Input:
!           lun11=logical unit number for printing
!           unit=File unit number to clos                   
!           status= Returned status code                
!     Dependencies:  none
!     Called by:  xstar, writespectra, writespectra2, 
!                 writespectra3, writespectra4, pprint
!                                                                       

      implicit none 
      integer unit, status,lun11 
!                                                                       
!      write (6,*)'closing fits unit ',unit                             
      call ftclos(unit, status) 
      call ftfiou(unit, status) 
      call frelunx(unit) 
      close(unit) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      return 
      END                                           
