      subroutine xstarcalc(lpri2,lnerrd,nlimdt,                         &
     &       lpri,lprid,lun11,tinf,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       xilev,bilev,rnist,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                           
                                                                        
!     Name: xstarcalc.f90  
!     Description:  
!           Calculates all relevant quantities for one spatial zone.
!           First calls dsec for thermal equilibrium calculation (if needed)
!           Then calls func and funcsyn
!
!     List of Parameters:
!           Input: 
!           lpri2: local print switch.  Turns on printing for 
!                 the func and funcsyn calls only
!           lnerrd: thermal equilibrium error switch
!           nlimdt: thermal equilibrium solver iteration limit
!           lpri:  print switch
!           lprid: thermal equilibrium solver (dsec) print switch
!           lun11: logical unit number for printing
!           tinf:  temperature lower limit
!           vturbi:  turbulent speed (km/s)
!           critf: threshold value for ion fraction to be included in 
!                   level population calculation
!           t: temperature in 10^4K
!           trad: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!           r:  radius in nebula (cm)
!           delr: thickness of current spatial zone (cm)
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           ababs(nl):  element abundances relative to H=1
!           cfrac:  covering fraction (affects line and continuum 
!                forward-backward ratio
!           p:  pressure in dynes/cm^2
!           lcdd: constant pressure switch, 1=constant pressure 
!                      0=constant density
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!       Output:
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!           tau0(2,nnnl):  line optical depths
!           tauc(2,nnml):  rrc optical depths
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           ntotit:  number of iterations for thermal equilibrium
!           xii(nni):  ion fractions, xiin(1)=H, xiin(2)=He0, xiin(3)=He+ etc
!           rrrt(nni): total recombination rates for each ion (s^-1)
!           pirt(nni): total photoionization rates for each ion(s^-1)
!           htt(nni): total heating rate for each ion (approximate) 
!                       (erg s^-1 cm^-3)
!           cll(nni): total cooling rate for each ion (approximate) 
!           httot: total heating rate (erg s^-1 cm^-3) 
!           cltot: total cooling rate (erg s^-1 cm^-3) 
!           hmctot: 2*(heating-cooling)/(heating+cooling)
!           elcter:  charge conservation error
!           xilev(nnml):  level populations (relative to parent element)
!           bilev(nnml):  level departure coefficient
!           rnist(nnml):  lte level populations (relative to parent element)
!           elum(nnnl):  line luminosities (erg/s/10^38)
!           rcem(2,nnnl):  line emissivities  (erg cm^-3 s^-1) /10^38
!                  inward and outward
!           oplin(nnnl):  line opacities  (cm^-1)
!           rccemis(2,ncn): continuum emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!                  inward and outward
!           brcems(ncn):  bremsstrahlung emissivity (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           opakcont(ncn):  continuum opacities lines excluded (cm^-1)
!           cemab(nnml):  rrc emissivities (erg cm^-3 s^-1) 
!           cabab(nnml):  total energy absorbed by rrc (erg cm^-3 s^-1) 
!           opakab(nnml):  rrc opacities (cm^-1)
!           fline(2,nnnl):  line emissivity (net radiative)
!              (erg cm^-3 s^-1) 
!           flinel(ncn):  line emissivity binned into continuum bins 
!              (erg cm^-3 s^-1 erg^-1)
!     Dependencies: dsec,func,funcsyn
!     Called by:  xstar
!
      use globaldata
      implicit none 
!                                                                       
!     global xstar data
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
      real(8) fline(2,nnnl),flinel(ncn) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     level populations                                                 
      real(8) xilev(nnml),bilev(nnml),rnist(nnml)
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating and cooling                                               
      real(8) htt(nni),cll(nni) 
      real(8) rrrt(nni),pirt(nni) 
!     element abundances                                                
      real(8) ababs(nl) 
!                                                                       
!     state variables                                                   
      real(8) p,r,t,xpx,delr 
!     heating-cooling variables                                         
      real(8) httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,          &
     &     clcont,hmctot                                                
      real(8) trad 
      real(8) cfrac,critf,vturbi,xee,tinf 
      integer lcdd,ncn2 
!     variables associated with thermal equilibrium solution            
      integer ntotit,lnerrd 
!     switches                                                          
      integer lprid,lpri,nlimdt 
!     strings for atomic data read                                      
      integer nlsvn,ncsvn,lun11,np2,lprisv,lpri2 
!                                                                       
!                                                                       
      lprisv=lpri 
      lpri=0 
      if (nlimdt.ne.0) then 
        call dsec(lnerrd,nlimdt,                                        &
     &       lpri,lprid,lun11,tinf,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &         cllines,clcont,htcomp,clcomp,clbrems,                    &
     &       xilev,bilev,rnist,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                        
        endif 
!       do ll=1,nnml                                                    
!         xilev(ll)=0.                                                  
!         enddo                                                         
!                                                                       
      if (lpri2.eq.1) lpri=1 
      call func(lpri,lun11,vturbi,critf,                                &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       xilev,bilev,rnist)
      call funcsyn(lpri,lun11,vturbi,critf,                             &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,                                                       &
     &       xilev,bilev,rnist,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                           
!                                                                       
       lpri=lprisv 
!                                                                       
!                                                                       
      return 
      end                                           
