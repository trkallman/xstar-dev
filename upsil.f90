      function upsil(k,eij,c,p1,p2,p3,p4,p5,t) 
!
!     Name: upsil.f90
!     Description:
!       this routine calculates upsilons for Burgess and Tully            
!       author:  M. Bautista                                              
!     Parameters:
!       Input:                                                                  
!       t = electron temperature in Kelvin                                
!       p# = spline knot values                                           
!       c = abscissa scale parameter                                      
!       k = transition type                                               
!       eij = transition energy (Ryd)                                     
!       Output:
!       upsil=upsilon value
!     Dependencies:  splinem
!     Called by:  ucalc
!                                                                       
!
      use globaldata
      implicit none 
      real(8) e, eij, c, p1, p2, p3, p4, p5, t 
      real(8) y, splinem, upsil, x 
      integer k 
!                                                                       
                                                 !<<<<<< CORRECTED LINE 
       e=abs(t/(1.57888e5*eij)) 
       if ((k.eq.1).or.(k.eq.4.)) x=log((e+c)/c)/log(e+c) 
       if ((k.eq.2).or.(k.eq.3)) x=e/(e+c) 
       y=splinem(p1,p2,p3,p4,p5,x) 
       if (k.eq.1) y=y*log(e+2.71828) 
       if (k.eq.3) y=y/(e+1) 
       if (k.eq.4) y=y*log(e+c) 
       upsil=y 
!                                                                       
      return 
      END                                           
