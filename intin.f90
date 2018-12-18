      subroutine intin(x1,x2,x0,t,ri2,ri3,lpri,lun11) 
!                                                                       
!     Name: intin.f90  
!     Description:  
!           Calculates the integrals needed by milne                   
!           Milne calculates the milne rate by dividing 
!           into intervals and assuming linear variation of
!           the cross section between boundaries
!           the integrals ri1,ri2 and ri3 are exponential integrals
!           over these intervals
!
!     List of Parameters:
!           Input: 
!           x1:  lower bound on x
!           x2:  upper bound on x
!           x0:  fiducial value of x
!           t:   temperature
!           lpri:  print switch
!           lun11: logical unit number for printing
!           Output:
!           ri2:  integral
!           ri3:  integral
!
!     Dependencies:
!           none
!
!     Called by:
!           milne.f90 
!
!     this routine does the integrals needed by milne                   
!     author:  M. Bautista                                              
!                                                                       
      implicit none 
!                                                                       
      real(8)  x1,x2,x0,t,ri2,ri3 
      real(8)  ryk,s1,s2,s0,del,rr 
      integer lpri,lun11 
!                                                                       
       ryk=7.2438d+15 
       s1=x1*ryk/t 
       s2=x2*ryk/t 
       s0=x0*ryk/t 
       del=ryk/t 
       if (lpri.gt.1)                                                   &
     &  write (lun11,*)'in intin:',s1,s2,s0,del                         
       if ((s1-s0).lt.90.) then 
!           ri2=dexpo(s0-s1)*(s1*s1+2.*s1+2.)-dexpo(s0-s2)*(s2*s2+2.*s2+
           ri2=exp(s0-s1)*((s1*s1+2.*s1+2.)                             &
     &        -exp(s1-s2)*(s2*s2+2.*s2+2.))/del/sqrt(del)               
           if (lpri.gt.1)                                               &
     &     write (lun11,*)'ri2=',ri2                                    
           if ((s0.lt.1.d-3).and.(s2.lt.1.d-3).and.(s1.lt.1.d-3))       &
     &       ri2=0.                                                     
         else 
           ri2=0. 
         endif 
!       ri2=ri2/(del**1.5)                                              
       if (lpri.gt.1)                                                   &
     &     write (lun11,*)'ri2=',ri2                                    
!       rr=dexpo(s0-s1)*(s1**3)-dexpo(s0-s2)*(s2**3)                    
       rr=exp(s0-s1)*((s1**3)-exp(s1-s2)*(s2**3)) 
       if (lpri.gt.1)                                                   &
     &     write (lun11,*)'rr=',rr                                      
       ri3=(rr/del/sqrt(del)+3.*ri2)/del 
       if (lpri.gt.1)                                                   &
     &  write (lun11,*)'in intin:',s1,s2,s0,del,ri2,rr,ri3              
       return 
      END                                           
