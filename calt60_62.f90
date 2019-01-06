       subroutine calt60_62(temp,m,idata,np1r,np1i,Upsilon) 
!                                                                       
!     Name: calt60_62.f90  
!     Description:  
!       This rutine takes the coefficients in data type 60 and 62 (reals     
!       as dtype and integers as itype) and return effective collision       
!       strengths according to fits by Callaway (1994).                      
!       "temp" is the temperature in Kelvin and "m" is the number of         
!       reals in dtype. "idata" is the data type either 60 or 62.            
!      author: M. Bautista                                              
!     Parameters:
!        Input:
!         temp=temperature in K                                              
!         m=number ofreals
!         idata=data type either 60 or 62
!         np1r=pointer to real data
!         np1i=pointer to integer data
!        Output:
!         Upsilon=upsilon
!     Dependencies: none
!     called by:  ucalc
!                                                                       
      use globaldata
        implicit none 
!                                                                       
        integer m 
!        real(8) dtype(m)                                                
        real(8) temp,Upsilon 
!        integer itype(m)                                               
        integer i,idata,np1i,np1r 
        real(8) t1,de,tmax,tt,rat 
!                                                                       
        t1=temp*6.33652e-6 
        if (temp.gt.1.e9) t1=6.33652e+3 
        de=1./float(masterdata%idat1(np1i))**2                          &
      &       -1./float(masterdata%idat1(np1i-1+2))**2 
        tmax=4.*de 
        tmax=1. 
        tt=t1 
        if (t1.gt.tmax) tt=tmax 
        if (idata.eq.60) then 
         rat=0. 
         do i=1,m-2 
          rat=rat+masterdata%rdat1(np1r-1+i+2)*(tt**(i-1)) 
         enddo 
         Upsilon=rat 
        else 
         rat=0. 
         do i=1,m-5 
          rat=rat+masterdata%rdat1(np1r-1+i+2)*(tt**(i-1)) 
         enddo 
         Upsilon=rat+masterdata%rdat1(np1r-1+m-2)                       &
     &         *log(masterdata%rdat1(np1r-1+m-1)*tt)                    &
     &         *exp(-masterdata%rdat1(np1r-1+m)*tt)                          
        endif 
!                                                                       
         if (t1.gt.tt) then 
          upsilon=Upsilon*(1.+log(t1/tmax)/(log(t1/tmax)+1.)) 
         endif 
!                                                                       
      return 
      END                                           
