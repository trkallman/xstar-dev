        subroutine frelunx(iounit) 
                                                                        
!       free specified logical unit number; if iounit=-1, then free all 
!                                                                       
!        I  (i) iounit - The logical unit number to be freed            
!                                                                       
!       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996    
!       Copied with minor changes from the FITSIO routine ftfiou.       
                                                                        
        integer iounit 
                                                                        
        call lunlstx(iounit) 
      END                                           
