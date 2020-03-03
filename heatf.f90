      subroutine heatf(lpri,lun11,                                      &
     &       t,r,delr,xee,xpx,                                          &
     &       epi,ncn2,                                                  &
     &       ncsvn,                                                     &
     &       brcems,cmp1,cmp2,httot,cltot,httot2,cltot2,hmctot,         &
     &             htcomp,clcomp,clbrems)                
!                                                                       
!     Name: heatf.f90  
!     Description:  
!           Adds Compton and brems heating, cooling, 
!           to total rates from calc_rates_level
!
!     List of Parameters:
!     Input:
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           t: temperature in 10^4K
!           r:  radius in nebula (cm)
!           cfrac:  covering fraction (affects line and continuum 
!                forward-backward ratio
!           delr: thickness of current spatial zone (cm)
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           abel(nl):  element abundances relative to H=1
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           nlev:  number of levels in ion
!           cmp1:  coefficient for use in compton heating calculation
!           cmp2:  coefficient for use in compton cooling calculation
!           brcems(ncn):  bremsstrahlung emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!       Output:
!           httot: total heating rate (erg s^-1 cm^-3) 
!           cltot: total cooling rate (erg s^-1 cm^-3) 
!           hmctot: 2*(heating-cooling)/(heating+cooling)
!           cllines:  total cooling rate due to lines (erg s^-1 cm^-3) 
!           clcont:   total cooling rate due to recombination (erg s^-1 cm^-3) 
!           htcomp: total heating rate due to compton (erg s^-1 cm^-3) 
!           clcomp: total cooling rate due to compton (erg s^-1 cm^-3) 
!           clbrems: total cooling rate due to bremsstrahlung (erg s^-1 cm^-3) 
!           
!        Dependencies:  Calls drd
!        Called by: calc_hmc_all
!
!     this routine calculates heating and cooling.                      
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!     energy bins                                                       
      real(8) epi(ncn) 
!     state variables                                                   
      real(8) r,t,xpx,delr,delrl 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) xee 
      integer ncn2,lpri,lun11 
      integer ncsvn 
      character(1) kblnk 
      real(8) brcems(ncn)
      real(8) cmp1,cmp2,cltot,                                          &
     &     hmctot,htcomp,clcomp,clbrems,ekt,                            &
     &     fac,                                                         &
     &     xnx,httot,httot2,cltot2
      real(8) ergsev
      real(8) tmp2,tmp2o 
      integer lskp 
      integer kl,                                                       &
     &    numcon
      integer                                                           &
     &     lprisv
!                                                                       
!                                                                       
      data kblnk/' '/ 
      data ergsev/1.602197e-12/ 
!                                                                       
      lprisv=lpri 
!                                                                       
      xnx=xpx*xee 
      if (lpri.ge.1) lpri=2 
      if (lpri.gt.1) write (lun11,*)'in heatf',httot,cltot,delr,r 
      if (lpri.gt.1) write (lun11,*)ncsvn 
      numcon=ncn2 
!                                                                       
!     comment these out to implement scattering                         
      clbrems=0. 
      lskp=1 
      tmp2=0. 
      do kl=1,numcon 
        tmp2o=tmp2 
        tmp2 =brcems(kl) 
        if ( kl.ge.2 ) clbrems=clbrems+(tmp2+tmp2o)                     &
     &                *(epi(kl)-epi(kl-lskp))*ergsev/2.                 
        enddo 
!                                                                       
      delrl=delr 
!                                                                       
      fac=delrl 
      ekt = t*(0.861707) 
      htcomp = cmp1*xnx*ergsev 
      clcomp = ekt*cmp2*xnx*ergsev 
      httot=httot+htcomp 
      cltot=cltot+clcomp+clbrems 
      httot2=httot2+htcomp 
      cltot2=cltot2+clcomp+clbrems 
      if (lpri.gt.1) write (lun11,9953)htcomp,clcomp,cmp1,cmp2,         &
     &   clbrems,httot,cltot                                            
      hmctot=2.*(httot-cltot)/(1.e-37+httot+cltot) 
      if (lpri.gt.1) write (lun11,*)hmctot 
 9953 format (1h , ' compton heating, cooling=',8e12.4) 
      lpri=lprisv 
!                                                                       
                                                                        
      return 
      end                                           
