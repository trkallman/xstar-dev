       subroutine eint(t,e1,e2,e3) 
!                                                                       
!  returns the values of the exponential integral function of order     
!  1, 2, and 3                                                          
!     author:  T. Kallman                                               
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
