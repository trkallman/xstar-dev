      subroutine xstarcalc(lpri2,lnerrd,nlimdt,                         &
     &       lpri,lprid,lun11,tinf,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,zeta,              &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       epim,ncn2m,bremsam,                                        &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,        &
     &       elcter,                                                    &
     &       cllines,clcont,htcomp,clcomp,clbrems,htfreef,              &
     &       httot2,cltot2,                                             &
     &       xilevg,bilevg,rnisg,                                       &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel,elin,errc)                           
                                                                        
!     Name: xstarcalc.f90  
!     Description:  
!           Calculates all relevant quantities for one spatial zone.
!           First calls dsec for thermal equilibrium calculation (if needed)
!           Then calls calc_hmc_all and calc_emis_all
!
!     List of Parameters:
!           Input: 
!           lpri2: local print switch.  Turns on printing for 
!                 the func and calc_emis calls only
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
!           abel(nl):  element abundances relative to H=1
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
!           xilevg(nnml):  level populations (relative to parent element)
!           bilevg(nnml):  level departure coefficient
!           rnisg(nnml):  lte level populations (relative to parent element)
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
!     Dependencies: dsec,calc_hmc_all,calc_emis_all
!     Called by:  xstar
!
      use globaldata
      implicit none 
!                                                                       
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,ndl) 
        integer:: ilev(10,ndl),nlpt(ndl),iltp(ndl) 
        character(1) :: klev(100,ndl) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
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
      real(8) epim(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
      real(8) bremsam(ncn)
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     level populations                                                 
      real(8) xilevg(nnml),bilevg(nnml),rnisg(nnml)
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating and cooling                                               
      real(8) htt(nni),cll(nni) 
      real(8) htt2(nni),cll2(nni) 
      real(8) rrrt(nni),pirt(nni) 
!     element abundances                                                
      real(8) abel(nl) 
      real(8) elin(nnnl)
      real(8) errc(nnml)
!     limits on ion indeces vs element
      integer mml(nl),mmu(nl)
!                                                                       
!     state variables                                                   
      real(8) p,r,t,xpx,delr 
!     heating-cooling variables                                         
      real(8) httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,         &
     &     clcont,hmctot,httot2,cltot2,zeta,htfreef
      real(8) trad 
      real(8) cfrac,critf,vturbi,xee,tinf 
      integer lcdd,ncn2,ncn2m
!     variables associated with thermal equilibrium solution            
      integer ntotit,lnerrd 
!     switches                                                          
      integer lprid,lpri,nlimdt 
!     strings for atomic data read                                      
      integer nlsvn,ncsvn,lun11,np2,lprisv,lpri2 
!                                                                       
!      write (lun11,*)'in xstarcalc',lun11,lpri,lprid
!
      call bremsmap(bremsa,bremsam,bremsint,epi,epim,ncn2,ncn2m,        &
     &        lpri2,lun11)       
!                                                                
      lprisv=lpri 
      lpri=0 
      if (nlimdt.ne.0) then 
        call dsec(lnerrd,nlimdt,                                        &
     &       lpri,lprid,lun11,tinf,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,zeta,              &
     &       mml,mmu,                                                   &
     &       epim,ncn2m,bremsam,bremsint,                               &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,elcter, &
     &         cllines,clcont,htcomp,clcomp,clbrems,htfreef,            &
     &       httot2,cltot2,                                             &
     &       xilevg,bilevg,rnisg)
        endif 
!                                                                       
!      lpri2=-1
      call calc_hmc_all(lpri2,lun11,vturbi,critf,                       &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,zeta,              &
     &       mml,mmu,                                                   &
     &       epim,ncn2m,bremsam,bremsint,                               &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,elcter, &
     &       cllines,clcont,htcomp,clcomp,clbrems,htfreef,              &
     &       httot2,cltot2,                                             &
     &       xilevg,bilevg,rnisg)  
!                                                                       
      call calc_emisab_all(lpri2,lun11,vturbi,critf,                    &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,                   &
     &       mml,mmu,                                                   &
     &       epim,ncn2m,bremsam,bremsint,                               &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,xilevg,bilevg,rnisg,                                   &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab)                         
!
      call calc_emis_all(lpri2,lun11,vturbi,critf,                      &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,                   &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,xilevg,bilevg,rnisg,                                   &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab,elin,errc)                         
!                                                                       
       lpri=lprisv 
!                                                                       
!                                                                       
      return 
      end                                           
