      subroutine printerror(lun11,status) 

      implicit none 
                                                                        
!     print out the fitsio error messages to the user                   
!     author:  T. Bridgman                                              
!                                                                       
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
