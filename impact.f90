      subroutine impact(en,l,temp,ic,z1,rm,ecm,psi,cr) 
!                                                                       
!     Name: impact.f90  
!     Description:  
!     impact parameter collision cross-sections using the method of seaton:
!        Proceedings of the Physical Society, Volume 79, Issue 6, 
!          pp. 1105-1117 (1962).
!     author:  M. Bautista                                              
!     List of Parameters:
!           Input:
!           en:  principal quantum number
!           l: angular momentum quantum number
!           temp:  temperature in K
!           ic = ionic charge of target particle                              
!           z1 = charge of incident particle                                  
!           rm = mass of incident particle in units of electron mass me       
!           ecm:  energy splitting in cm^-1
!           psi = see notes for defn                                          
!           Output:
!           cr:  rate
!      Called by: amcrs
!      Dependencies:  impcfn
!
      implicit none 
!                                                                       
      real(8) en,temp,z1,rm,ecm,psi,cr 
      integer l,ic 
!                                                                       
      real(8)  b,xsi,phi,bo,xsw,phw,del 
      real(8) tk,fi,wo,ev,po,w,wi,ff,crinc,ric 
      integer inc,jm,j 
!                                                                       
      tk=8.617e-5*temp 
      inc=1 
      jm=90*inc 
      cr=0. 
      fi=0. 
      wo=0. 
      b=10. 
      ric=float(ic) 
      ev=abs(ecm)/8065.48 
!                                                                       
      po=(3.*en*en-l*(l+1))/2./ic 
!                                                                       
! strong coupling                                                       
!                                                                       
   21  del=b/100./inc 
      do 20 j=1,jm 
      b=b-del 
      call impcfn(b,xsi,phi) 
!     write (lun11,*)b,xsi,phi                                          
      w=ric*rm*ev/(b)*sqrt(2.*(xsi)*psi) 
      wi=w+ecm/8065.48/2. 
!     write (lun11,*)wi/tk                                              
      if(wi/tk.ge.100.) go to 13 
      if(wi.le.0.) go to 20 
!                                                                       
! weak coupling                                                         
!                                                                       
      bo=dble(po*ev/2./w*sqrt(wi*rm/13.60)) 
      call impcfn(bo,xsw,phw) 
!                                                                       
! the minimum of the weak and strong coupling x-sections is used        
      ff=min((xsi/2.+phi),(phw)) 
      ff=ff*exp(-wi/tk) 
!                                                                       
      crinc=(fi+ff)/2.*(wi-wo) 
      cr=crinc+cr 
      if(cr.lt.1.e-20) go to 20 
      fi=ff 
      wo=wi 
      if(crinc/cr.lt.1.e-5) go to 13 
!                                                                       
   20  continue 
      go to 21 
   13     cr=6.900e-5*z1*z1*sqrt(rm/temp)*psi*cr/tk 
!                                                                       
      return 
      END                                           
