      subroutine fitsclose(lun11,unit,status) 
!                                                                       
!     Close the file & release the unit number                          
!     author: T. Bridgman                                               
!                                                                       
!     Parameters:                                                       
!        unit    integer            File unit number                    
!        status  integer            Returned status code                
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
