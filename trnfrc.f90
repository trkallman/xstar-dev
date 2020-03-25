      subroutine trnfrc(lpri,lun11,ldir,                             &
     &      r,xpxcol,xpx,                                               &
     &      epi,ncn2,zremsz,dpthc,opakc,                                &
     &      zrems,bremsa,bremsint)                               
!                                                                       
!     Name: trnfrc.f90  
!     Description:  
!           Calculates local mean intensity
!           Currently assumes single stream.
!
!     List of Parameters:
!     Input:
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           ldir:  direction: 1=outward, -1=inward
!           r: radius in nebula (cm)
!           xpxcol:  column density (cm^-2)
!           xpx: H number density (cm^-3)
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           dpthc(2,ncn): optical depth in continuum bins 
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           zrems(4,ncn):  master spectrum array.  (erg/s/erg/10^38)
!     Output:
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!
!     Dependencies: none
!     Called by:  xstar
!
!     this routine calculates continuum transfer                        
!     author:  T. Kallman (from xstar1)                                 
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8) epi(ncn),zremsz(ncn),dpthc(2,ncn),bremsa(ncn) 
      real(8) bremsint(ncn),opakc(ncn),                                  &
     &          zrems(5,ncn)
      integer lpri,lun11,ldir,ncn2,jkp,jk,ncnm 
      real(8) r,xpxcol,xpx,r19,fpr2,rmax,sumtmp 
!                                                                       
      r19=r/(1.e+19) 
      fpr2=(12.56)*r19*r19 
      ncnm=ncn2-1 
      bremsa(ncn2)=0. 
      bremsa(ncnm)=0. 
      bremsint(ncn2)=0. 
      bremsint(ncnm)=0. 
      rmax=xpxcol/xpx 
      if (lpri.gt.0) write (lun11,*)'in trnfrc:',rmax,xpxcol,xpx 
      do 1 jkp=1,ncnm 
         jk=ncnm+1-jkp 
!                                                                       
!        for outward only                                               
!                                                                       
         if (ldir.lt.0) then 
             bremsa(jk)=zrems(1,jk)/fpr2 
!     $                      +flinel(jk)                                
           else 
             bremsa(jk)=zremsz(jk)*exp(-dpthc(1,jk))/fpr2 
           endif 
         sumtmp=(bremsa(jk)+bremsa(jk+1))*(epi(jk+1)-epi(jk))/2. 
         bremsint(jk)=bremsint(jk+1)+sumtmp*(1.602197e-12) 
         if (lpri.gt.0) write (lun11,*)jk,epi(jk),dpthc(1,jk),          &
     &       zremsz(jk),opakc(jk),                                      &
     &       zremsz(jk)*exp(-dpthc(1,jk))/fpr2,bremsa(jk)               &
     &       ,bremsint(jk)                                              
    1    continue 
!                                                                       
      return 
      END                                           
