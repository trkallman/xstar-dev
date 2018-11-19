      subroutine deletefile(filename,status) 
!                                                                       
!     a simple little routine to delete a fits file                     
!     author:  T. Bridgman                                              
!                                                                       

      implicit none 
      integer status,unit,blocksize 
      character*(*) filename 
      character(50) kcom 
      character(1) ktmp 
      integer mm,ll,lenact 
!                                                                       
      data kcom/'rm -f                                             '/ 
!                                                                       
      mm=lenact(filename) 
                                                                        
      do ll=1,mm 
        read (filename(ll:ll),'(a1)')ktmp 
        write (kcom(6+ll:6+ll),'(a1)')ktmp 
        enddo 
      do ll=mm+7,50 
        write (kcom(ll:ll),'(a1)')' ' 
        enddo 
!       write (6,*)'executing:',kcom                                    
!      this is slow but it works                                        
!      call system(kcom)                                                
!      return                                                           
!                                                                       
!     simply return if status is greater than zero                      
!      write (6,*)'in deletefile',unit,filename,status                  
      if (status .gt. 0)return 
!                                                                       
!     get an unused logical unit number to use to open the fits file    
      call ftgiou(unit,status) 
      call getlunx(unit) 
!      write (6,*)'after ftgiou',unit,status                            
!                                                                       
!     try to open the file, to see if it exists                         
      call ftopen(unit,filename,1,blocksize,status) 
!      write (6,*)'after ftopen',unit,status                            
!                                                                       
      if (status .eq. 0)then 
!         file was opened;  so now delete it                            
          call ftdelt(unit,status) 
!      write (6,*)'after ftdelt 1',unit,status                          
      else if (status .eq. 103)then 
!         file doesn't exist, so just reset status to zero and clear err
          status=0 
          call ftcmsg 
      else 
!         there was some other error opening the file; delete the file a
          status=0 
          call ftcmsg 
          call ftdelt(unit,status) 
!      write (6,*)'after ftdelt 2',unit,status                          
      end if 
                                                                        
                                                                        
!     free the unit number for later reuse                              
!      call ftfiou(unit, status)                                        
      call frelunx(unit) 
      close(unit) 
!      write (6,*)'after ftfiou',unit,status                            
!                                                                       
      END                                           
