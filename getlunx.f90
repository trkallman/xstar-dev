      subroutine getlunx(iounit) 
                                                                        
!     Name: getlunx.f90
!     Description:
!       get an unallocated logical unit number                          
!       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996    
!       Copied with minor changes from the FITSIO routine ftgiou.       
!     Parameters:
!        Output:
!        O  (i) iounit - An unopened logical unit number                
!     Dependencies: lunlstx                                    
!     Called by: deletefile, readtbl,fheader,rread1,xstar,xstarsetup
!                                                                        
        integer iounit 
                                                                        
        iounit=0 
        call lunlstx(iounit) 
      END                                           
