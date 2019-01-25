        subroutine frelunx(iounit) 
                                                                        
!     Name: frelunx.f90
!     Description:
!       this sub-routine 
!       free specified logical unit number; if iounit=-1, then free all 
!       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996    
!       Copied with minor changes from the FITSIO routine ftfiou.       
!     Parameters:                                        
!        Input:
!        I  (i) iounit - The logical unit number to be freed            
!     Dependencies:  lunlstx
!     Called by: deletefile, fitsclose
!
        integer iounit 
                                                                        
        call lunlstx(iounit) 
      END                                           
