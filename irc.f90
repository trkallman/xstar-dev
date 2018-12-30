      subroutine irc(n,t,rc,rno,se,lpri,lun11) 
!                                                                       
!     Name:  irc.f90
!     Description:
!          calculates ionization rate for ionization of hydrogen atoms
!
!     List of Parameters:
!          Input:
!          n: principal quantum number
!          t: temperature
!          rc: effective ion charge
!          rno: effective principal quantum number
!          lpri: print switch
!          lun11: logical unit number for printing
!          Output:
!          se:  collisional ionization rate
!
!      Dependencies:
!          szirc
!      Called by:
!           ucalc
!
!
! irc calculates the excitation rate, se [cm**3/s], for ionization      
! of hydrogen atoms from state n due to electron collisions, assuming   
! the continuum starts at level rno.  the                               
! energy loss rate, elost [ev*cm**3/s], is also determined.             
! cin is the 3-body recombination rate, determined from cni by          
! detailed balance.                                                     
! ref. johnson (1972)                                                   
! note: this does not make sense:  reference does not exist?
!     author:  m. bautista                                              
!                                                                       
      implicit none 
!                                                                       
      real(8) t,rc,rno,se 
      integer n 
      real(8) xo,yn,an,bn,rn,g0,g1,g2,zn,ey 
      real(8) ez 
      integer lpri,lun11 
!                                                                       
      if (lpri.ne.0) write (lun11,*)'in irc',                           &
     & n,t,rc,rno                                                       
                                                 ! mab                  
      if(rc.ne.1.) then 
       call szirc(n,t,rc,rno,se,lpri,lun11) 
       return 
      endif 
!                                                                       
      xo=1.-n*n/rno/rno 
      yn=xo*157803./(t*n*n) 
      if(n-2) 100,200,300 
  100  an=1.9603*n*(1.133/3./xo**3-0.4059/4./xo**4+0.07014/5./xo**5) 
      bn=2./3.*n*n/xo*(3.+2./xo-0.603/xo/xo) 
      rn=0.45 
      go to 400 
!                                                                       
  200  an=1.9603*n*(1.0785/3./xo**3-0.2319/4./xo**4+0.02947/5./xo**5) 
      bn=(4.-18.63/n+36.24/(n*n)-28.09/(n*n*n))/n 
      bn=2./3.*n*n/xo*(3.+2./xo+bn/xo/xo) 
      rn=0.653 
      go to 400 
!                                                                       
  300  g0=(0.9935+0.2328/n-0.1296/(n*n))/3./xo**3 
      g1=-(0.6282-0.5598/n+0.5299/(n*n))/(n*4.)/xo**4 
      g2=(0.3887-1.181/n+1.470/(n*n))/(n*n*5.)/xo**5 
      an=1.9603*n*(g0+g1+g2) 
      bn=(4.-18.63/n+36.24/(n*n)-28.09/(n*n*n))/n 
      bn=(3.+2./xo+bn/xo/xo)*2.*n*n/3./xo 
      rn=1.94*n**(-1.57) 
!                                                                       
  400  continue 
      rn=rn*xo 
      zn=rn+yn 
      call expint(yn,ey) 
      call expint(zn,ez) 
      se=an*(ey/yn/yn-exp(-rn)*ez/zn/zn) 
      ey=1.+1./yn-ey*(2./yn+1.) 
      ez=exp(-rn)*(1.+1./zn-ez*(2./zn+1.)) 
      se=se+(bn-an*log(2.*n*n/xo))*(ey-ez) 
      se=se*sqrt(t)*yn*yn*n*n*1.095e-10/xo 
      if (lpri.ne.0) write (lun11,*)'in irc',                           &
     & xo,yn,an,bn,rn,zn,ey,ez,se                                       
!      cii=se*n*n                                                       
! mab                                                                   
!     cii=cii/100.                                                      
!      cni=se*exp(-yn)                                                  
!      cin=4.144219e-16*n*n/t/sqrt(t)*se                                
!      elost=se*13.60/(n*n)                                             
!                                                                       
      return 
      END                                           
