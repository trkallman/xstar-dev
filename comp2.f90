      subroutine comp2(lpri,lun11,epi,ncn2,bremsa,t,cmp1,cmp2)      
!                                                                       
!     Name: comp2.f90  
!     Description:  
!       this  computes the heating - cooling due to compton    
!       scattering.
!       relativistic version, using rates from I. Khabibulin (private comm.)
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
!     Called by:  func
!                                                                       
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      real(8) bremsa(ncn),epi(ncn) 
      real(8) t,cmp1,cmp2,emc2,sigth0,tmp1,                              &
     &     ekt,xx,sxx,zrmstp,eee,ee,sum1,sum2,sum3,tmp1o,eeeo,          &
     &     eeo,ans,cfake,hfake,cohc,cmpfnc                              
      integer lpri,lun11,ncn2,lprisv,numcon,kl 
!                                                                       
      data emc2/5.11e+5/,sigth0/6.65e-25/ 
!                                                                       
      lprisv=lpri 
!      lpri=2                                                           
      if (lpri.ge.1) write (lun11,*)'in comp2' 
!                                                                       
      sigth0 = 6.65e-25 
      tmp1 = 0. 
!                                                                       
      ekt = t*0.861707 
      xx = emc2/(ekt+1.e-10) 
      sxx = 1./xx 
      zrmstp = bremsa(1) 
      eee = epi(1) 
      ee = eee/emc2 
      tmp1 = zrmstp*cmpfnc(ee,sxx,lun11,lpri) 
      sum1 = 0. 
      sum2 = 0. 
      sum3 = 0. 
      numcon=ncn2 
      do kl = 2,numcon 
         tmp1o = tmp1 
         eeeo = eee 
         eeo = ee 
         eee = epi(kl) 
         ee = eee/emc2 
         zrmstp = bremsa(kl) 
         tmp1 = zrmstp*cmpfnc(ee,sxx,lun11,lpri) 
         sum1 = sum1 + (tmp1+tmp1o)*(eee-eeeo)/2. 
         sum2 = sum2 + (bremsa(kl)+bremsa(kl-1))*(eee-eeeo)/2. 
         sum3 = sum3 + (bremsa(kl)*ee+bremsa(kl-1)*eeo)*(eee-eeeo)/2. 
         if (lpri.ne.0) write (lun11,*)kl,eee,ee,zrmstp,tmp1,sum1 
         enddo 
      ans = sum1 
      cfake=sum2*sigth0 
      hfake=sum3*sigth0 
      cohc = -ans*sigth0 
      cmp1=hfake 
      cmp2=(-cohc+hfake)/ekt 
                                                                        
                                                                        
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'cmp1,cmp2:',cmp1,cmp2,cfake,hfake                
!                                                                       
      lpri=lprisv 
!                                                                       
      return 
      end                                           
