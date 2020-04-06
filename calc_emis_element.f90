      subroutine calc_emis_element(ml_element,lpri,lun11,xeltp,        &
     &       vturbi,critf,temperature,trad,radius,delr,xee,xpx,         &
     &       abel,cfrac,pressure,lcdd,                                  &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,xileve,bileve,rnise,ipmat,                             &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab)                         

!                                                                       
!     Name: calc_emis_element.f90  
!     Description:  
!           Calculates emissivities and opacitiesfor 
!           one element.  
!
!     List of Parameters:
!           Input:
!           ml_element: poiter to element record
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           vturbi: turbulent speed in km/s
!           critf:  critical fractional abundance for ion inclusion
!           t: temperature in 10^4K
!           trad: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!           r:  radius in nebula (cm)
!           delr: thickness of current spatial zone (cm)
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           xh1:  H+ number density (cm^-3)
!           xh0:  neutral H number density (cm^-3)
!           cfrac:  covering fraction (affects line and continuum 
!                forward-backward ratio
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!           tau0(nnnl):  line optical depths
!           tauc(nnml):  rrc optical depths
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           xileve(nd):  level populations (relative to parent element)
!           bileve(nd):  departure coefficients for levels
!           rnise(nd): lte level populations
!           nmat: index into rates array corresponding to first level
!           also uses variables from globaldata
!           Output:
!           rcem(2,nnnl):  line emissivities  (erg cm^-3 s^-1) /10^38
!                  inward and outward
!           rccemis(2,ncn): continuum emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!                  inward and outward
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           opakcont(ncn):  continuum opacities lines excluded (cm^-1)
!           cemab(nnml):  rrc emissivities (erg cm^-3 s^-1) 
!           
!        Dependencies:  Calls ucalc,drd
!        Called by: 
!                  
!     this routine calculates rates affecting level populations         
!     author: T. Kallman                                                
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
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     continuum emissivities                                            
      real(8) brcems(ncn) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) xileve(nd),rnise(nd),bileve(nd)
      real(8) xilevi(nd),rnisi(nd),bilevi(nd)
      real(8) xii(nni)
!     line opacities                                                    
      real(8) oplin(nnnl) 
      real(8) rcem(2,nnnl)
!     element abundances                                                
!     state variables                                                   
      real(8) radius,temperature,xpx,delr,pressure
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee 
      integer ncn2,lpri,lun11,lfpi,np2,lcdd
      integer nlsvn,ncsvn 
      real(8) tsq,xh1,xh0,cfrac,critf
      integer ipmat
!     continuum flux                                                    
      real(8) tauc(2,nnml) 
!     element abundances                                                
      real(8) abel(nl) 
!     limits on ion indeces vs element
      integer mml(nl),mmu(nl)
!                                                                       
      character(1) kblnk 
      real(8) xeltp
      integer np1i,np1r,np1k 
      integer nlev,jk,klion,                                            &
     &        ltyp,mm,nsp,jkk_ion,                                      &
     &        lrtyp,lcon,nrdt,nidt,nkdt,                                &
     &        ml_ion,ml_element,ml_ion_data_type,                       &
     &        ml_element_test
!                                                                       
      data kblnk/' '/ 
!
      save kblnk
!                                                                       
!                                                                       
      if (lpri.ge.1)                                                    &
     &  write (lun11,901)temperature,xee,xpx,lcdd,pressure,delr                              
901   format (1x,'    in calc_emis_element, inputs:',3(1pe10.3),        &
     &i6,2(1pe10.3))
!                                                                       
!     print element information
      call drd(ltyp,lrtyp,lcon,                                         &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_element,             &
     &            0,lun11)                                        
      jk=masterdata%idat1(np1i)
      if (lpri.ge.1) then 
        write (lun11,9339)ml_element,                                   &
     &    (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)
        endif 
 9339 format (1x, ' element:',i12,1x,12(1a1)) 
!
!     lfpi value: photoionization and recombination, no opacities       
      lfpi=2 
!                                                                       
      tsq=sqrt(temperature) 
      ipmat=0
      nsp=1
!                                                                       
      ml_ion_data_type=12
!     step thru ions
      ml_ion=derivedpointers%npfirst(ml_ion_data_type)
      do while (ml_ion.ne.0)
!
!       test if element belongs to parent of ion
        ml_element_test=derivedpointers%npar(ml_ion)
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)'ml_element_test=',ml_element_test,ml_element
        if (ml_element_test.eq.ml_element) then
!
!         get ion index
          call drd(ltyp,lrtyp,lcon,                                     &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_ion,                 &
     &            0,lun11)                                        
          klion=masterdata%idat1(np1i)
          jkk_ion=masterdata%idat1(np1i+nidt-1)
          nlev=derivedpointers%nlevs(jkk_ion)
!
!         test for ion abundancs
!          if (xii(jkk_ion).gt.critf) then
          if ((klion.ge.mml(jk)).and.(klion.le.mmu(jk))) then
!
!           get level data                                          
            if (lpri.gt.1) write (lun11,*)'nlev=',nlev
            if (lpri.gt.1) write (lun11,*)'ipmat=',ipmat
            do mm=1,nlev
!             get level pointer                                     
              bilevi(mm)=bileve(mm+ipmat)
              rnisi(mm)=rnise(mm+ipmat)
              xilevi(mm)=xileve(mm+ipmat)
              if (lpri.gt.1)                                            &
     &         write (lun11,*)'before calc_emis_ion',                   &
     &         mm,rnise(mm),xilevi(mm),bilevi(mm)
              enddo
            call calc_emis_ion(ml_ion,lpri,lun11,xeltp,                 &
     &       vturbi,critf,temperature,trad,radius,delr,                 &
     &       xee,xpx,xh0,xh1,cfrac,pressure,lcdd,                       &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &                   leveltemp,                                     &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xilevi,bilevi,rnisi,                                       &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab)                         
!
!           end of test for ion abundancs
            endif
!
          ipmat=ipmat+nlev-1
!
!         end of test if element belongs to parent of ion
          endif
!
!       end of loop over ions
        ml_ion=derivedpointers%npnxt(ml_ion)
        enddo 
!
      if (lpri.gt.0) write (lun11,*)'returning from calc_emis_element'
!
      return 
      end                                           
