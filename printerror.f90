      subroutine printerror(lun11,status) 

!     Name:  printerror.f90
!     Description:
!       print out the fitsio error messages to the user                   
!       author:  T. Bridgman                                              
!     Parameters:
!       Input:
!       lun11=logical unit number
!       status=cfitsio status code
!     Dependencies: none
!     Called by: fheader, fparmlist, writespectra, writespectra2, 
!       writespectra3, writespectra4, rstepr, rstepr2, rstepr3, rstepr4, 
!       fstepr, fstepr2, fstepr3, fstepr4, pprint, fwrtascii, savd, unsavd
!                                                                       
      implicit none 
                                                                        
      integer status, lun11 
      character errtext*30,errmessage*80 
!                                                                       
!     check if status is ok (no error); if so, simply return            
      if (status .le. 0)return 
!                                                                       
!     get the text string which describes the error                     
      call ftgerr(status,errtext) 
      write (lun11,*)'fitsio error status =',status,': ',errtext 
                                                                        
!     read and print out all the error messages on the fitsio stack     
      call ftgmsg(errmessage) 
      do while (errmessage .ne. ' ') 
          write (lun11,*)errmessage 
          call ftgmsg(errmessage) 
      end do 
!                                                                       
      END                                           
