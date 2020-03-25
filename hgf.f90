      subroutine hgf(ia,ib,ic,x,hyp) 
!                                                                       
!     Name: hgf.f90  
!     Description:  
!       subroutine hgf calculates the value, hyp, of the                  
!       hypergeometric fn at x for constants ia,ib,ic                     
!       real(8)  ser,hyp                                                   
!       author:  M. Bautista                                              
!     Parameters:
!         Input:
!         ia=first constant
!         ib=second constant
!         ic=third constant
!         x= independent variable
!         Output:
!         hyp=output
!     Dependencies: none
!     Called by: anl1
!                                                                       

      implicit none 
      integer ia,ib,ic,i,j,n 
      real(8) x,hyp,ser 
!                                                                       
      ser=1. 
      hyp=1. 
      i=-ia 
      j=-ib 
      i=min(i,j) 
      do 10 n=0,i 
      ser=ser*(ia+n)*(ib+n)*x/((n+1.)*(ic+n)) 
      hyp=hyp+ser 
   10  continue 
!         
      return 
      END                                           
