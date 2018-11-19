       subroutine calt68(temp,np1r,np1i,gamma) 
!                                                                       
!    Takes coefficients in data type 68 and returns effective collision 
!    strenghts for He-like ions according to Sanpson & Zhang.           
!      author: M. Bautista                                              
!                                                                       
      use globaldata
       implicit none 
!                                                                       
       real(8) z,tt,gamma,temp 
       integer np1r,np1i 
!                                                                       
       z=float(masterdata%idat1(np1i-1+3)) 
       tt=log10(temp/z/z/z) 
       gamma=masterdata%rdat1(np1r-1+1)+masterdata%rdat1(np1r-1+2)*tt   &
     &       +masterdata%rdat1(np1r-1+3)*tt*tt 
!                                                                       
       return 
      END                                           
