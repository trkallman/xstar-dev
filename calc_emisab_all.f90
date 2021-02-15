      subroutine calc_emisab_all(lpri,lun11,vturbi,critf,               &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,                   &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,xilevg,bilevg,rnisg,                                   &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab)                         
                                                                        
!     List of Parameters:
!           Input:
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           vturbi: turbulent speed in km/s
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
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!           tau0(nnnl):  line optical depths
!           tauc(nnml):  rrc optical depths
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           xilevg(nnml):  level populations (global)
!           bilevg(nnml):  level departure coefficient
!           rnisg(nnml):  lte level populations (global)
!           Output:
!           rcem(2,nnnl):  line emissivities  (erg cm^-3 s^-1) /10^38
!                  inward and outward
!           rccemis(2,ncn): continuum emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!                  inward and outward
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           opakcont(ncn):  continuum opacities lines excluded (cm^-1)
!           cemab(nnml):  rrc emissivities (erg cm^-3 s^-1) 
!           also uses variables from globaldata
!           
!        Dependencies:  Calls 
!        Called by: xstar, dsec
!
!     this routine steps through data and calculates                    
!     new version attempts to avoid rates for unabundant ions           
!     author: T. Kallman                                                
!                                                                       
!     with data structures designed for Lucy's iterative method         
!       nsup is a pointer from level n to superlevel N                  
!                                                                       
!     no longer calls full calc_emis_ion in main loop.                          
!     calc_emis_ion calls moved to calc_emis_all
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
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
      real(8) rccemis(2,ncn)
!     level populations                                                 
      real(8) xilevg(nnml),rnisg(nnml),bilevg(nnml)
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) tauc(2,nnml) 
      real(8) xii(nni)
      real(8) abel(nl)
!     limits on ion indeces vs element
      integer mml(nl),mmu(nl)
!                                                                       
      real(8) xileve(nd),rnise(nd),bileve(nd)
      integer lpri,lun11,lcdd,ncn2,np2,ncsvn,nlsvn 
      real(8) vturbi,critf,t,trad,r,delr,xee,xpx,cfrac,p
      real(8) xh1,xh0,htfreef
      real(8) xeltp
      integer nlev,klion,jkk,ipmat,ltyp,jkk_ion,                        &
     &     lrtyp,lcon,nrdt,nidt,nkdt,ll,jk,mm,mmtmp,                    &
     &     nnz,mlm,np1i,np1r,np1k,lprisv,ml                 
      integer ml_element_data_type,ml_element                    
      integer np1ki,nkdti,ml_ion,ml_ion_data_type,ml_element_test
      integer ipmatmax
      integer nnzz,nnnn
      real(8) sigth,xnx
!                      
      data sigth/6.65e-25/                                                 
      save sigth
!            
      lprisv=lpri 
      if (lpri.ge.1)                                                    &
     &  write (lun11,901)t,xee,xpx,lcdd,p,delr,cfrac                              
901   format (1x,'in calc_emisab_all, inputs:',3(1pe10.3),i6,3(1pe10.3))
       if (lcdd.ne.1)                                                   &
     &   xpx = p/1.38e-12/max(t,1.d-24)                                 
       xnx=xpx*xee

!                                                                       
      xh0=xpx*xilevg(1)*abel(1) 
      xh1=xpx*(1.-xilevg(1))*abel(1) 
!     NB this is a fudge for testing
!      xin1=6.358e-04
!      xh0=xpx*xin1*abel(1) 
!      xh1=xpx*(1.-xin1)*abel(1) 
!     NB testing
      xh1=0.
      xh0=0.

!
       do ll=1,nnml 
         cemab(1,ll)=0.
         cemab(2,ll)=0.
         cabab(ll)=0.
         opakab(ll)=0.
         enddo
       do ll=1,nnnl 
         rcem(1,ll)=0. 
         rcem(2,ll)=0. 
         oplin(ll)=0. 
         enddo
!
!     step thru elements                         
      ml_element_data_type=11 
      ml_element=derivedpointers%npfirst(ml_element_data_type) 
      do while (ml_element.ne.0) 
!
!       get data for this element
        call drd(ltyp,lrtyp,lcon,                                       &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,ml_element,                    &
     &     0,lun11)                                               
!        call dprinto(ltyp,lrtyp,lcon,                                  &
!     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)   
!
!       test if abundance
        xeltp=0. 
        jk=masterdata%idat1(np1i)
        nnz=masterdata%idat1(np1i+nidt-2)
        if (jk.gt.0) xeltp=abel(jk) 
        if (xeltp.gt.1.d-24) then 
!
          if (lpri.gt.1)                                                &
     &        write (lun11,902)jk,ml_element,nnz,                       &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,min(8,nkdt))
902           format (1x,'  element:',3(i12,1x),8(1a1))
!
          ipmat=0
!
!         step thru ions and fill arrays of old populations
!         test if first ion pointer
          ml_ion_data_type=12
          ml_ion=derivedpointers%npfirst(ml_ion_data_type)
          do while (ml_ion.ne.0) 
!
            if (lpri.gt.1)                                              &
     &         write (lun11,*)'ml_ion=',ml_ion
!
!           test if element belongs to parent of ion
            ml_element_test=derivedpointers%npar(ml_ion)
            if (lpri.gt.1)                                              &
     &      write (lun11,*)'ml_element_test=',ml_element_test,ml_element
!
!           test if belongs to element
            if (ml_element_test.eq.ml_element) then
!
!               retrieve ion name from kdati                              
                mlm=ml_ion
                call drd(ltyp,lrtyp,lcon,                               &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                    &
     &            0,lun11)                                        
                klion=masterdata%idat1(np1i)
!
                jkk=masterdata%idat1(nidt+np1i-1)
                if (lpri.gt.1)                                          &
     &            write (lun11,903)jkk,ml_ion,                          &
     &               (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)
903             format (1x,'      ion:',2(i12,1x),8(1a1))
!
!               get level data                                          
                if (lpri.gt.1) write (lun11,*)'jkk=',jkk
                nlev=derivedpointers%nlevs(jkk)
                if (lpri.gt.1) write (lun11,*)'nlev=',nlev
                do mm=1,nlev
!                 get level pointer                                     
                  mmtmp=derivedpointers%npilev(mm,jkk) 
                  rnise(mm+ipmat)=rnisg(mmtmp)
                  bileve(mm+ipmat)=bilevg(mmtmp)
                  xileve(mm+ipmat)=xilevg(mmtmp)
!                  write (lun11,*)'mapping populations:',mm,mm+ipmat,    &
!     &                 mmtmp,xileve(mm+ipmat)
                  enddo
                if (lpri.gt.1) then 
                  write (lun11,*)'  index level             energy',    &
     &              '  stat.wt. population   LTE ' 
                  do mm=1,nlev 
                    write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20)&
     &                       ,leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),&
     &                       xileve(mm+ipmat),rnise(mm+ipmat)                                   
 9022               format (2x,i4,1x,20a1,4(1pe10.3)) 
                    enddo 
                  endif
                ipmat=ipmat+nlev-1
!
!               end of test if belongs to element
                endif
!
!           end of step thru ions
            ml_ion=derivedpointers%npnxt(ml_ion)
            ml_element_test=0
            if (ml_ion.gt.0)ml_element_test=derivedpointers%npar(ml_ion)
            enddo
!
!         now fill rate matrix for element
          call calc_emisab_element(ml_element,lpri,lun11,xeltp,         &
     &       vturbi,critf,t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,      &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,xileve,bileve,rnise,ipmat,                             &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab)                         
!                                                                       
!
!         end of test if abundance
          endif 
!
!       end of step thru elements
        if  (ml_element.ne.0)                                           &
     &          ml_element=derivedpointers%npnxt(ml_element) 
        enddo 
!
!                                                                       
      if (lpri.gt.1) write (lun11,*)'leaving calc_emisab_all' 
!                                                                       
      lprisv=lpri 
!         
      return 
      end                                           
