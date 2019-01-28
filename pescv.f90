      real(8) function pescv(tau) 
!                                                                       
!     Name:  pescv.f90
!     Description:
!       this routine calculates escape probability for an rrc
!       inputs: optical depths-- tau for a line transition                
!       Currently using trivial expression
!     Parameters:
!       Input: 
!       tau=line center optical depth
!       output:
!       pescl=escape probability
!     Dependencies:  none
!     Called by: ucalc
!                                                                       
      implicit none 
!                                                                       
      real(8) tau 
!      real(8) taubar
      real(8) eps 
!                                                                       
!     NB optically thin rrcs                                            
!     test for lte                                                      
      pescv=0.5 
!      return                                                           
!      if (e.lt.13.6) return                                            
!     fudge because of too much case b                                  
!      taubar=tau/200.                                                  
!      taubar=tau/5.                                                    
!      taubar=tau*100.                                                  
!      pescv=1./(taubar+1.)                                             
!      pescv=expo(-taubar)                                              
      pescv=exp(-tau) 
!     fudge because of numerical problem with large tau.                
      eps=1.e-12 
      pescv=max(pescv,eps) 
      pescv=pescv/2. 
!                                                                       
      return 
      end                                           
