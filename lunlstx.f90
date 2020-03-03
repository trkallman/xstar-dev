        subroutine lunlstx(iounit) 
                                                                        
!     Name: lunlstx.f90
!     Description:
!       generic routine to manage logical unit numbers in the range 10-40
!       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996    
!       Copied with minor changes from the FITSIO routine ftxiou.       
!     Parameters:
!       Input:
!       I/O  (i) iounit - The logical unit number to be allocated/freed 
                                                                        
        integer iounit,i 
        integer array(400) 
        save array 
        data array/400*0/ 
                                                                        
!        write (6,*)'entering lunlstx',iounit                           
!                                                                       
        if (iounit .eq. 0)then 
!           get an unused logical unit number                           
            do 10 i=400,1,-1 
                                                                        
!        The following would be a more robust way of testing for        
!        an available unit number, however, this cannot work            
!        when building XANLIB using the IRAF/SPP version, because       
!        IRAF does not use Fortran I/O.                                 
!                                                                       
!                inquire(unit=iounit, exist=exists, opened=open)        
!                if(exists .and. .not. open)then                        
!                    array(iounit-9)=1                                  
!                    return                                             
!                end if                                                 
                                                                        
!               write (6,*)'i,array(i)',i,array(i)                      
               if (array(i) .eq. 0)then 
                     array(i)=1 
                     iounit=i+9 
!                     write (6,*)'allocating:',iounit,i                 
                     return 
                 end if 
   10       continue 
!           error: all units are allocated                              
            iounit=-1 
            call xaerror(                                               &
     &           'GETLUNx has no more available unit numbers.', 1)      
                                                                        
        else if (iounit .eq. -1)then 
!           deallocate all the unit numbers                             
            do 20 i=1,400 
                 array(i)=0 
   20       continue 
                                                                        
        else 
!            deallocat a specific unit number                           
             if (iounit .ge. 10 .and. iounit .le. 49)then 
!                write (6,*)'deallocating:',iounit,iounit-9             
                array(iounit-9)=0 
             end if 
        endif 
      END                                           
