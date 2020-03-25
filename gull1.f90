!-----------------------------------------------------------------------
      subroutine gull1(n,rs,gus,gls,lpri,lun11) 
!                                                                       
!     Name: gull1.f90  
!     Description:  
!        this subroutine calculates the value of |g(n,l;r,l')|**2              
!        given n the principal qn and r for all l=[o,n-1] and                  
!        l'=l+1 or l'=l-1.  ref burgess (1964), 
!        and brockelhurst (1971) MNRAS 153, 471.
!      author: m. bautista                                              
!      Parameters:
!           Input:
!           n=number of levels
!           rs=k/z^2, k=wavenumber of outgoing electron
!           lpri=print switch
!           lun11=logical unit number
!           Output:
!           gus=quantities needed for radial matrix element (eqn 3.17)
!           gls=quantities needed for radial matrix element (eqn 3.17)
!     Called by: hphotx
!     Dependencies:  fact
!                                                                       

      implicit none 
!                                                                       
      integer n,lpri,lun11 
      real(8)  cn,clu,cll,g0,gu(100),gl(100),fn,pi,s,r 
      real(8)  dn,dl 
      real(8) f1,rs 
      real(8) gus(100),gls(100) 
      integer n1,l 
!                                                                       
      dn=dfloat(n) 
      r=dble(rs) 
      pi=2.*acos(0.) 
      n1=2*n-1 
      call fact(n1,f1) 
      g0=0.5*log(pi/2.)+log(8.*dn)                                      &
     &   +dn*log(4.*dn)-dble(f1)                                        
        if(r.eq.0.) then 
        gu(n)=g0-2.*dn 
        else 
      s=sqrt(r) 
      gu(n)=g0-2.*datan(dfloat(n)*s)/s                                  &
     &  -0.5*log(1.-exp(-2.*pi/s))                                      
        end if 
      gu(n)=exp(gu(n)) 
!                                                                       
      fn=1.d-300/gu(n) 
      gu(n)=gu(n)*fn 
!                                                                       
      if(n.eq.1) go to 40 
      gu(n-1)=(2.*dn-1.)*(1.+dn*dn*r)*dn*gu(n) 
      gl(n)=(1.+dn*dn*r)*gu(n)/(2.*dn) 
      gl(n-1)=(2.*dn-1.)*(4.+(dn-1.)*(1.+dn*dn*r))*gl(n) 
!                                                                       
      do 10 l=n-1,3,-1 
      dl=dfloat(l) 
      gu(l-1)=(4.*dn*dn                                                 &
     &  -4.*dl*dl+dl*(2.*dl-1.)*(1.+dn*dn*r))*gu(l)                     
      gu(l-1)=gu(l-1)-4*dn*dn*(dn-dl)                                   &
     &        *(dn+dl)*(1+(dl+1)*(dl+1)*r)*gu(l+1)                      
      gl(l-1)=(4*dn*dn-4*(dl-1)*(dl-1)                                  &
     &        +(dl-1)*(2*dl-1)*(1+dn*dn*r))*gl(l)                       
      gl(l-1)=gl(l-1)-4*dn*dn*(dn-dl)                                   &
     &        *(dn+dl)*(1+(dl-1)*(dl-1)*r)*gl(l+1)                      
   10  continue 
      gl(1)=0. 
      gu(1)=(4.*dn*dn-16.+6.*(1.+dn*dn*r))*gu(2) 
      gu(1)=gu(1)-4.*dn*dn*(dn-2.)*(dn+2.)*(1.+9.*r)*gu(3) 
!                                                                       
      cn=log(dn)-dn*log(4.*dn*dn)                                       &
     &    -(2.*dn+4.)*log(1.+dn*dn*r)                                   
      gu(1)=cn+log(1.+r)+2.*log(gu(1))-2.*log(fn) 
      clu=cn+log(1.+r) 
      cll=cn 
!                                                                       
      do 30 l=1,n-1 
      dl=dfloat(l) 
      clu=clu+log(4.*dn*dn*(dn-dl)                                      &
     &    *(dn+dl)*(1.+(dl+1.)*(dl+1.)*r))                              
      cll=cll+log(4.*dn*dn*(dn-dl)                                      &
     &    *(dn+dl)*(1.+(dl-1.)*(dl-1.)*r))                              
      gu(l+1)=clu+2.*log(gu(l+1))-2.*log(fn) 
      gl(l+1)=cll+2.*log(gl(l+1))-2.*log(fn) 
   30  continue 
      go to 60 
!                                                                       
   40  gl(1)=0. 
      gu(1)=2.*log(gu(1))-log(4.)                                       &
     &      -5.*log(1.+r)-2.*log(fn)                                    
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'in gull1:',r,fn,gu(1),gu(2),gl(1),gl(2)          
! converts results to single precision to give in retudn                
       do l=1,100 
        gus(l)=(gu(l)) 
        gls(l)=(gl(l)) 
       enddo 
                                                                        
   60  return 
      END                                           
