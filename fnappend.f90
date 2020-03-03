      subroutine fnappend(knam,nint) 
!
!     Name: fnappend.f90
!     Description:  
!       Add an integer to a string in elements 3 and 4
!     Parameters:
!        Input:
!        knam= input string (on input)
!        nint=integer to append
!        Output:
!        knam= output string (on output)
!     Dependencies: none
!     Called by:  xstar
!
      character(16) knam 
      if (nint.gt.9) then 
          write (knam(3:4),'(i2)')nint 
        else 
          write (knam(3:3),'(a1)')'0' 
          write (knam(4:4),'(i1)')nint 
        endif 
      return 
      end                                           
