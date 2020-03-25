      subroutine pexs(nmin,kdim,zc,eion,far,gam,scal,                &
     &                e,axs,ierr,lpri,lun11)                            
!                                                                       
!     Name:  pexs.f90
!     Description:
!       Compute photoexcitation cross-section assuming one                
!       Rydberg serie converging to a threshold.                          
!       Author: P. Palmeri 2005
!     Parameters:
!       Input:
!       nmin = starting princ. quant. num. of the serie                   
!       zc = effective charge, Z-Ne+1                                     
!       eion = threshold energy in Ry                                     
!       far = oscillator strength of the serie member with                
!            n=nmin                                                     
!       gam = resonance width in Ry                                       
!       e = external energy grid in Ry                                    
!       kdim = dimension of the external energy grid                      
!       scal = scaling factor                                             
!       Output:
!       axs = cross section                                               
!       ierr = error indicator (=0: OK; .ne.0:error)                      
!     Dependencies:  none
!     Called by:  ucalc
!                                                                       
      implicit none 
!                                                                       
      integer nmax,nmin 
      real(8)  pi 
!                                                                       
      parameter(pi=3.14159,nmax=30) 
!                                                                       
      real(8)  x(nmax),a(nmax),e(*),axs(*) 
      real(8)  zc,eion,far,gam,scal 
      integer lpri,ierr,lun11,kdim,im 
      integer nres,jmin,jmax,jres,jj,i,ii,ij,kk,n 
      real(8)  del,xmin,xres,axtp,res 
!                                                                       
      data x,a/nmax*0.,nmax*0./ 
!                                                                       
<<<<<<< HEAD
      if (lpri.gt.0) write (lun11,*)'in pexs',nmin,kdim,zc,          &
=======
      if (lpri.ne.0) write (lun11,*)'in pexs',nmin,kdim,zc,          &
>>>>>>> 2d75308c63b9789458ce092c697c7853fcdde44a
     &        eion,far,gam,scal                                         
!                                                                       
      ierr=0 
      if(nmin.ge.nmax) then 
       ierr=1 
       return 
      endif 
      nres=nmax 
      do 10 n=nmin,nmax 
!                                                                       
!   energy of the resonance ...                                         
!                                                                       
       x(n)=-(zc/dble(n))**2. 
!                                                                       
!   area of the resonance ...                                           
!                                                                       
       a(n)=8.06725*far*dble(nmin**3)/dble(n**3) 
!                                                                       
!   search for unresolved limit in term of member ...                   
!                                                                       
       if(n.gt.nmin) then 
        del=x(n)-x(n-1) 
         res=gam/2. 
        if(del.gt.res) nres=n 
       endif 
                                                                        
        if (lpri.gt.0) write (lun11,*)n,x(n),a(n),del,res 
                                                                        
   10 continue 
!                                                                       
!   define shifted energy range ...                                     
!                                                                       
      xmin=x(nmin)-30.*gam 
      xres=x(nres) 
      jmin=1 
      jmax=kdim 
      jres=jmax 
      do 20 i=1,kdim 
        axs(i)=0. 
       e(i)=e(i)-eion 
       im=max(1,i-1) 
       if(i.gt.1.and.e(im).le.xmin.and.e(i).gt.xmin)                    &
     &    jmin=i-1                                                      
       if(i.gt.1.and.e(im).le.xres.and.e(i).gt.xres)                    &
     &    jres=i-1                                                      
       if(i.gt.1.and.e(im).lt.0..and.e(i).ge.0.)                        &
     &    jmax=i-1                                                      
   20 continue 
      if(jmin.eq.jmax) jmax=jmin+1 
      if (lpri.gt.0) write (lun11,*)'jmin,jres:',jmin,jres,             &
     &              jmax,nmin,nmax                                      
      do 30 ii=nmin,nmax 
       do 3000 jj=jmin,jres 
!                                                                       
!   constant-width Lorentzian resonances ...                            
!                                                                       
        axtp=a(ii)/pi*gam/2.                                            &
     &  /((e(jj)-x(ii))**2+(gam/2.)**2)                                 &
     &  +a(ii)/pi*gam/2./((abs(e(jj))-x(ii))**2                         &
     &   +(gam/2.)**2)                                                  
!                                                                       
!   near-threshold pill-up (oscill. strength conservation)              
!                                                                       
        axs(jj)=axs(jj)+axtp 
      if ((lpri.gt.0).and.(axs(jj).gt.1.d-24))                          &
     &  write (lun11,*)'30 loop',jj,e(jj),axtp,axs(jj)                  
3000  continue 
   30 continue 
!                                                                       
!   near-threshold extrapolation ...                                    
!                                                                       
      do ij=jres+1,jmax 
        axs(ij)=axs(jres) 
        if (lpri.gt.0) write (lun11,*)ij,axs(ij) 
      enddo 
!                                                                       
!   scaling of the xs ...                                               
!                                                                       
      do 40 kk=1,kdim 
       axs(kk)=scal*axs(kk) 
!                                                                       
!   return to the "usual" energy grid ...                               
!                                                                       
       e(kk)=e(kk)+eion 
       if ((lpri.gt.0).and.(axs(kk).gt.1.d-24))                         &
     &  write (lun11,*)kk,e(kk),axs(kk)                                 
   40 continue 
!                                                                       
!                                                                       
      return 
      END                                           
