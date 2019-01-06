      subroutine calt72(temp,np1r,nrdt,rate,lun11,lpri) 
!                                                                       
!     Name: calt72.f90  
!     Description:  
!       Takes coefficients in data type 72 and returns capture rates        
!       (in s^-1) for DR through satellite levels considered explicitly.    
!       author: M. Bautista                                              
!     Parameters:
!        Input:
!        temp=temperature in K
!        np1r=pointer to beginning of real data
!        nrdt=number of reals
!        lpri=print switch
!        lun11=logical unit number
!        Output:
!        rate=rate
!     Dependencies: none
!     called by:  ucalc
!                                                                       
      use globaldata
       implicit none 
!                                                                       
!                                                                       
       real(8) temp,rate,dele,s,rtmp,temp4,ekt 
       integer np1r,nrdt 
       integer lpri,lun11 
!                                                                       
       dele=masterdata%rdat1(np1r-1+2) 
!     using safranova's expression                                      
       s=3.3e-11*(13.6/(0.861707*temp/1.e+4))**1.5 
!       s=4.141292e-22/(temp**1.5)                                      
       rtmp=masterdata%rdat1(np1r+2) 
       temp4=temp/1.e+4 
       ekt=0.861707*temp4 
       if (nrdt.lt.3) rtmp=1. 
       rate=s*exp(-dele/ekt)*(masterdata%rdat1(np1r-1+1)/1.e+13)*rtmp 
!     $          /2.                                                    
       if (lpri.gt.1) write (lun11,*)'in calt72:',dele,s,temp,rtmp,rate 
!                                                                       
       return 
      END                                           
