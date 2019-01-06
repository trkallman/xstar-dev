       subroutine calt67(temp,np1r,gamma) 
!                                                                       
!     Name: calt67.f90  
!     Description:  
!      Takes coefficients in data type 67 and returns effective collision  
!      strenghts for He-like ions according to Keenan, McCann, & Kingston 
!      (1987) eq. (2)                                                     
!      author: M. Bautista                                              
!     Parameters:
!        Input:
!        temp=temperature in K
!        np1r=pointer to real data
!        nrdt=number real data parameters
!        gamma=upsilon
!     Dependencies: none
!     called by:  ucalc
!                                                                       
      use globaldata
       implicit none 
!                                                                       
!                                                                       
       real(8) gamma,temp,tp 
       integer np1r 
!                                                                       
       tp=log10(temp) 
       gamma=masterdata%rdat1(np1r-1+1)+masterdata%rdat1(np1r-1+2)*tp   &
     &       +masterdata%rdat1(np1r-1+3)*tp*tp 
!                                                                       
       return 
      END                                           
