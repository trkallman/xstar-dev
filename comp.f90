      subroutine comp(lpri,lun11,epi,ncn2,bremsa,cmp1,cmp2) 
!                                                                       
!                                                                       
!     this subroutine computes the heating - cooling due to compton     
!     scattering.                                                       
!     non-relativistic version                                          
!     the rate is returned in the common block coheat.                  
!                                                                       
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      real(8) bremsa(ncn),epi(ncn),cmp1,cmp2 
      integer lpri,lun11,ncn2,lprisv,numcon,i
      real(8) sigth,ery,alpha,beta,fac1,fac3,delt,fac4,ebar,c1,c2,fac2,  &
     &  tmp1,tmp2
      
!                                                                       
!                                                                       
      data c2/8.219e-06/ 
!                                                                       
      lprisv=lpri 
!      lpri=3                                                           
      if (lpri.ge.1) write (lun11,*)'in comp' 
!                                                                       
      sigth = 6.65e-25 
      c1=1.95639e-6 
      tmp1 = 0. 
      tmp2 = 0. 
      c2 = 0. 
!                                                                       
!     due to continuum.                                                 
      fac1 = sigth*bremsa(1)*epi(1)*(1.-c2*epi(1)) 
      fac3 = sigth*bremsa(1)*4. 
      numcon=ncn2 
      do 100 i = 2,numcon 
         delt = epi(i) - epi(i-1) 
         fac2 = sigth*bremsa(i)*epi(i)*(1.-c2*epi(i)) 
         tmp1 = tmp1 + (fac1+fac2)*delt/2. 
         fac1 = fac2 
         fac4 = sigth*bremsa(i)*4. 
         tmp2 = tmp2 + (fac3+fac4)*delt/2. 
         fac3 = fac4 
         if ( lpri.gt.2 ) write (lun11,99001) i,epi(i),bremsa(i),       &
     &                           fac1,fac3,tmp1,tmp2                    
  100 continue 
!                                                                       
      ebar = tmp1*4./(1.e-30+tmp2) 
      if ( lpri.gt.2 ) write (lun11,*) 'ebar=',ebar 
!                                                                       
!                                                                       
      if (lpri.gt.2)  write (lun11,*)c1,tmp1,tmp2 
      cmp1 = c1*tmp1 
      cmp2 = c1*tmp2 
      if (lpri.gt.2) write (lun11,*)cmp1,cmp2 
      lpri=lprisv 
!                                                                       
      return 
99001 format (' ',i4,6e12.4) 
      end                                           
