      subroutine amcrs(n,l,temp,ic,z1,rm,ne,sum,ecm,psi,il,cn,          &
     &                 lun11,lpri)                                      
!                                                                       
!     Name: amcrs.f90  
!     Description:  
!      Calculates angular momentum changing collision rates using    
!      either the pengelly & seaton (1964) formula (amcol) or the impact     
!      parameter method of seaton (1962) (impact) if the energy levels are   
!      non-degenerate.  the ps routine is used if the ratio of amcol/impact  
!      is greater than 0.94 since this code is faster.  ** beware - there    
!      may be problems if ne is too large ( > 1.e+7).  pc1 will be used in   
!      amcol rather than pc3 and the change will not occur.                  
!      author: M. Bautista                                              
!
!     List of Parameters:
!           Input:
!           n = principal quantum number of initial state                     
!           l = orbital quantum number of initial state                       
!           temp = temperature in kelvin                                      
!           ic = ionic charge of target particle                              
!           z1 = charge of incident particle                                  
!           rm = mass of incident particle in units of electron mass me       
!           ne = electron number density                                      
!           sum = total spontaneous transition rate out of n,l                
!           ecm = energy difference between nl and nl-1                       
!           psi = see notes for defn                                          
!           il = flag to decide whether to use impact or amcol                
!           Output:
!           cn = transition rate for nl -> nl-1                               
!     Dependencies: amcol, velimp, impact
!     Called by:  ucalc
!                                                                       
      implicit none 
!                                                                       
      real(8) ne 
      real(8) temp,z1,rm,sum,ecm,psi,cr 
      real(8) en,dnl,rho1,rhom,cn,rat 
      integer n,l,ic,il,lun11,lpri 
!                                                                       
      en=real(n) 
      dnl=6.*z1/ic*z1/ic*n*n*(n*n-l*l-l-1) 
      rho1=0.72/sum 
      if(ecm.ne.0.) rho1=min(rho1,5.946e-12/ecm) 
      rhom=3.929e11*rho1*temp/sqrt(dnl)/rm 
      if(rhom.lt.10.) go to 30 
      call amcol(n,l,temp,ic,z1,rm,ne,sum,ecm,cn,lun11,lpri) 
!      if (lpri.ge.1)                                                    &
!     &  write (lun11,*)'after amcol:',n,l,temp,ic,z1,rm,ne,sum,ecm,cn
      cn=0 
! mab il=0                                                              
        il=0 
        if(ecm.ne.0.) then 
          if(il.eq.0) then 
          if (lpri.ge.1)                                                &
     &    write (lun11,*)'call impact 1',en,l,temp,ic,z1,rm,ecm,psi     
          call impact(en,l,temp,ic,z1,rm,ecm,psi,cr) 
          rat=cn/cr 
          cn=cr 
          if(rat.gt. 0.94) il=1 
          endif 
        endif 
!     go to 40                                                          
!                                                                       
   30  if(ecm.eq.0.) then 
      call velimp(n,l,temp,ic,z1,rm,ne,sum,cn) 
      else 
      if (lpri.ge.1)                                                    &
     &write (lun11,*)'call impact 2',en,l,temp,ic,z1,rm,ecm,psi         
      call impact(en,l,temp,ic,z1,rm,ecm,psi,cn) 
      endif 
                                                                        
!      if(ne.gt.1.e14) then                                             
!     call impact(en,l,temp,ic,z1,rm,ecm,psi,cn)                        
!      endif                                                            
!                                                                       
! 40    continue                                                        
!                                                                       
      return 
      END                                           
