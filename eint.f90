       subroutine eint(t,e1,e2,e3) 
!                                                                       
!     Name: ee1exp.f90  
!     Description:  
!       returns the values of the exponential integral function of order     
!       1, 2, and 3                                                          
!       author:  T. Kallman                                               
!     Parmameters:
!          t=independent variable
!          Output:
!          e1,e2,e3= exponential integrals
!     Dependencies:  expint
!     Called by: calt57, szirc, szirco, ucalc, calt73
!                                                                       
!                                                                       

       implicit none 
       real(8) t,e1,e2,e3,ss,expo 
!                                                                       
       e1=0. 
       e2=0. 
       e3=0. 
!       if (t.gt.50.) return                                            
       call expint(t,ss) 
       e1=ss/t/exp(t) 
       e2=exp(-t)-t*e1 
       e3=0.5*(expo(-t)-t*e2) 
       return 
      END                                           
