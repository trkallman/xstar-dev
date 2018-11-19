        subroutine getlunx(iounit) 
                                                                        
!       get an unallocated logical unit number                          
!                                                                       
!        O  (i) iounit - An unopened logical unit number                
!                                                                       
!       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996    
!       Copied with minor changes from the FITSIO routine ftgiou.       
                                                                        
        integer iounit 
                                                                        
        iounit=0 
        call lunlstx(iounit) 
      END                                           
