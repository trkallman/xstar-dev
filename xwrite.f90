      subroutine xwrite(string,ll) 
!
!     Name:  xwrite.f90
!     Description:
!        Write a string to the standard output
!        Rewrite of a xanlib routine
!        Author:  T. Kallman
!     Parameters:
!        Input:
!        string=string to write
!        ll=not used
!     Dependencies: none
!     Called by: xstar
!
      implicit none 
      integer ll 
      character(120) string 
!                                                                       
!     Not used                                                          
      integer javi 
      javi=ll 
!      ll=javi                                                          
!                                                                       
      write (6,*)string 
      return 
      end                                           
