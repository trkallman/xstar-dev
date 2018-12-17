                                                                        
      integer function lenact(cbuf) 
      implicit none 
      character cbuf*(*) 
!---                                                                    
! function to return the active length of a character string, not       
! counting any trailing blanks.  n.b. an all blank string will          
! return zero as the length.                                            
!---                                                                    
! cbuf    i    string whose length is to be measured.                   
!---                                                                    
! 1988-jun-13 - standard fortran version [aft]                          
!---                                                                    
      integer   i 
!---                                                                    
      do 190 i=len(cbuf),1,-1 
         if(cbuf(i:i).ne.' ') then 
            lenact=i 
            return 
         end if 
  190 continue 
      lenact=0 
      return 
      end                                           
