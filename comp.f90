      subroutine comp(lpri,lun11,epi,ncn2,bremsa,cmp1,cmp2) 
!                                                                       
!                                                                       
!     Name: comp.f90  
!     Description:  
!       this  computes the heating - cooling due to compton    
!       scattering. 
!       non-relativistic version                                          
!       Output is cmp1, cmp2, to be used as follows:
!         htcomp = cmp1*xnx*ergsev 
!         clcomp = ekt*cmp2*xnx*ergsev 
!       where htcomp, clcomp are heating, cooling rates in erg cm^-3 s^-1
!       xnx is electron number density, ergsev=1.602197e-12, ekT=kT in eV
!       author:  T. Kallman (from xstar1)                                 
!     Parameters:
!           Input:
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           t: temperature in 10^4K
!           lun11=logical unit number for printing
!           lpri=print switch
!           Output:
!           cmp1=compton heating coefficient
!           cmp2=compton cooling coefficient
!     Dependencies: cmpfnc
!     Called by:  not currently called
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
