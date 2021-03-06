      subroutine szirco(nn,t,rz,cii) 
!                                                                       
!     Name:  szirco.f90
!     Description:
!          calculates ionization rate for ionization of hydrogen atoms
!          calculates electron impact ionizition rates from semiempirical    
!          formula (eq.35) from Smpson & Zhang (1988, ApJ 335, 516)          
!          author:  M. Bautista                                              
!
!     List of Parameters:
!          Input:
!          n: principal quantum number
!          t: temperature
!          rz: ion charge
!          rno:  principal quantum number
!          lpri: print switch
!          lun11: logical unit number for printing
!          Output:
!          cii:  collisional ionization rate
!
!      Dependencies:
!          none
!      Called by:
!          irc.f90
!
!     calculates electron impact ionizition rates from semiempirical    
!     formula (eq.35) from smpson & zhang (1988, apj 335, 516)          
!     author:  M. Bautista                                              
!                                                                       
!                                                                       

       implicit none 
!                                                                       
       real(8) abethe(11), hbethe(11), rbethe(11) 
       integer nn,i 
       real(8) t,rz,cii,boltz,eion,const,an,hn,rrn,tt,rn,yy,             &
     &      e1,e2,e3,term1,term2,term3                                  
!                                                                       
       data(abethe(i),i=1,11)/ 1.134, 0.603, 0.412, 0.313, 0.252,       &
     &       0.211, 0.181, 0.159, 0.142, 0.128, 1.307 /                 
       data(hbethe(i),i=1,11)/ 1.48, 3.64, 5.93, 8.32, 10.75, 12.90,    &
     &       15.05, 17.20, 19.35, 21.50, 2.15 /                         
       data(rbethe(i),i=1,11)/ 2.20, 1.90, 1.73, 1.65, 1.60, 1.56,      &
     &       1.54, 1.52, 1.52, 1.52, 1.52 /                             
!                                                                       
       boltz=1.38066e-16 
       eion=2.179874e-11 
       const=4.6513e-3 
!                                                                       
       if (nn.lt.11) then 
         an=abethe(nn) 
         hn=hbethe(nn) 
         rrn=rbethe(nn) 
       else 
         an=abethe(11)/float(nn) 
         hn=hbethe(11)*float(nn) 
         rrn=rbethe(11) 
       endif 
       tt= t*boltz 
       rn=float(nn) 
!       rz=float(nz)                                                    
       yy=rz*rz/(rn*rn)*eion/tt 
       call eint(yy,e1,e2,e3) 
       term1=e1/rn-(exp(-yy)-yy*e3)/(3.*rn) 
       term2=(yy*e2-2.*yy*e1+exp(-yy))*3.*hn/rn/(3.-rrn) 
       term3=(e1-e2)*3.36*yy 
       cii=const*sqrt(tt)*(rn**5)/(rz**4)*an*yy* (                      &
     &   term1+term2+term3)                                             
!       write (lun11,*)'in szirco:',nn,t,an,hn,rrn,rn,yy,e1,e2,e3,term1, 
!     $    term2,term3,cii                                              
!      give result                                                      
!                                                                       
      END                                           
