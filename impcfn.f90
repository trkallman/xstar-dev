      subroutine impcfn(x,xsi,phi) 
!                                                                       
! data for functions used in the impact parameter method are generated  
! using polynomials fitted to seaton's (1962) values using least square 
!     author:  M. Bautista                                              
!                                                                       
      implicit none 
!                                                                       
      real(8)  a(6),b(6),x,xsi,phi,pi,y 
      integer n 
!                                                                       
      pi=2.*acos(0.) 
      a(1)=0.9947187 
      a(2)=0.6030883 
      a(3)=-2.372843 
      a(4)=1.864266 
      a(5)=-0.6305845 
      a(6)=8.1104480d-02 
      b(1)=0.2551543 
      b(2)=-0.5455462 
      b(3)=0.3096816 
      b(4)=4.2568920d-02 
      b(5)=-2.0123060d-02 
      b(6)=-4.9607030d-03 
!                                                                       
      if(x.gt.2.) go to 25 
      xsi=0. 
      phi=0. 
      do 20 n=1,6 
      xsi=xsi+a(n)*x**(n-1) 
      y=log(x) 
      phi=phi+b(n)*y**(n-1) 
   20  continue 
      if(x.eq.1.) phi=b(1) 
      if(x.lt.0.05) then 
      xsi=1.0+0.01917/0.05*x 
      y=log(1.1229/x) 
      phi=y+x*x/4.*(1.-2.*y*y) 
      endif 
      go to 30 
!                                                                       
   25  xsi=pi*x*exp(-2.*x)*(1.+0.25/x+1./32./x/x) 
      phi=pi/2.*exp(-2.*x)*(1.+0.25/x-3./32./x/x) 
!                                                                       
   30  continue 
!                                                                       
      return 
      END                                           
