      subroutine hphotx(ener,ic,nq,xsec,lun11,lpri) 
!                                                                       
!     Name: hphotx.f90  
!     Description:  
!        this subroutine calculates the photoionization cross section 
!       in the hydrogenic approximation, using expressions from 
!       Brocklehurst 1971 MNRAS 153 471
!        real(8)  en,cons,r,rk,theta1,theta2,gu,gl                       
!       author:  M. Bautista                                              
!     Parameters:
!       Input:
!       ic=ion charge                                                 
!       np=principal quantum number                                   
!       ll=angular momentum number                                    
!       Output:
!       ener=photon energy in ryds with respect to the ionization
!         threshold                                                     
!       xsec=array containing the cross section in mb (18^{-18} cm
!         for all l=[0,nq-1]                                            
!     Dependencies:  gull1
!     Called by:  ucalc
!                                                                       
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      real(8) gu(100),gl(100),xsec(100) 
      real(8) cons,ener,en,r,rk,theta1,theta2 
      integer ic,nq,lun11,lpri,lm 
!                                                                       
      cons=.54492*acos(0.) 
!        write (lun11,*)'in hphotx:',ener,ic,nq,cons                    
        en=ener 
!        r=sqrt(en)                                                     
        r=sqrt(en) 
        rk=r/float(ic*ic) 
        if (lpri.ne.0)                                                  &
     &   write (lun11,*)'before call gull1:',nq,rk,r,en                 
        call gull1(nq,rk*rk,gu,gl,lpri,lun11) 
!        call gull1(nq,rk,gu,gl,lun11)                                  
       do lm=0,nq-1 
        theta1=(1.+nq*nq*rk*rk)*exp(gu(lm+1)) 
        theta2=(1.+nq*nq*rk*rk)*exp(gl(lm+1)) 
        if (lpri.ne.0)                                                  &
     &   write (lun11,*)'after call gull1:',lm,gu(lm+1),gl(lm+1),       &
     &           theta1,theta2                                          
        xsec(lm+1)=cons*((lm+1)*theta1+lm*theta2)/(2.*lm+1.) 
        xsec(lm+1)=float(nq*nq)/float(ic*ic)*xsec(lm+1) 
       enddo 
!                                                                       
      return 
      END                                           
