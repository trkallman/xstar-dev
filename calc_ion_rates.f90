      subroutine calc_ion_rates(ml_ion,lpri,lun11,                      &
     &                   vturbi,t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,    &
     &                   epi,ncn2,bremsa,bremsint,                      &
     &                   leveltemp,                                     &
     &                   tau0,tauc,                                     &
     &                   np2,ncsvn,nlsvn,                               &
     &                   pirti,rrrti,                                   &
     &                   rnisi,nlev)

!                                                                       
!     Name: calc_ion_rates.f90  
!     Description:  
!           Calculates  rates affecting level ion fractions for 
!           one ion.  Puts results in pirti,rrrti
!           formerly:  func1
!
!     List of Parameters:
!           Input:
!           jkk: index of ion in xstar scheme 1=H0, 432=Zn29+
!           kl:  index of ion relative element: 1=neutral, n=hydrogenic
!           lpriz: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           vturbi: turbulent speed in km/s
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
!           rnisi:  lte level populations  relative to ion 
!           nlev:  number of levels for the ion
!           also uses variables from globaldata
!           Output:
!           pirti:  total ionization rate
!           rrrti:  total recombination rate
!           
!        Dependencies:  Calls ucalc,drd
!        Called by: calc_hmc_all
!                  
!     this routine calculates rates affecting level populations         
!     author: T. Kallman                                                
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,nd) 
        integer:: ilev(10,nd),nlpt(nd),iltp(nd) 
        character(1) :: klev(100,nd) 
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
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     element abundances                                                
!     state variables                                                   
      real(8) r,t,xpx,delr 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee 
      integer ncn2,lpri,lun11,lfpi,np2 
      integer nlsvn,ncsvn 
      character(49) kdesc2 
      real(8) tsq,ans1,ans2,xh1,xh0,cfrac 
      real(8) abund1,abund2,ptmp1,ptmp2,ans3,ans4,ans5,ans6,opakb1 
      integer idest1,idest2,idest3,idest4
!     continuum flux                                                    
      real(8) tauc(2,nnml) 
!                                                                       
      character(1) kblnk 
      real(8) rnisi(nd)
      real(8) rnisse(nd) 
      real(8) tau1,tau2,ptmp
      real(8) pirti,rrrti
      integer np1i,np1r,np1k 
      integer nlev,ml,lpriu,mm,lprisv,                                  &
     &        llo,lup,ltyp,jkk_ion,                                     &
     &        lrtyp,lcon,nrdt,nidt,nkdt,                                &
     &        ml_data,ml_ion,ml_data_type,ml_data_par
!                                                                       
      data kblnk/' '/ 
!                                                                       
      lprisv=lpri
!      if (lpri.ge.1) lpri=2
!
      if (lpri.ge.1)                                                    &
     &  write (lun11,901)t,xee,xpx,delr                              
901   format (1x,'      in calc_ion_rates, inputs:',                    &
     &           4(1pe10.3))
!                                                                       
!     lfpi value: photoionization only
      lfpi=1
!                                                                       
      tsq=sqrt(t) 
      rrrti=0.
      pirti=0.
!                                                                       
!     get ion index
      call drd(ltyp,lrtyp,lcon,                                         &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_ion,                 &
     &            0,lun11)                                        
      jkk_ion=masterdata%idat1(np1i+nidt-1)
      if (lpri.ge.1)                                                    &
     &            write (lun11,903)jkk_ion,ml_ion,                      &
     &               (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)
903             format (1x,'      ion:',2(i12,1x),8(1a1))
!
      call calc_rates_level_lte(jkk_ion,lpri,lun11,t,xee,xpx,           &
     &          leveltemp,nlev)
      nlev=derivedpointers%nlevs(jkk_ion)
      if (lpri.gt.1) then 
        write (lun11,*)'      nlev=',nlev
        write (lun11,*)'  index level             energy',              &
     &              '  stat.wt.                LTE ' 
        do mm=1,nlev 
          write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20)          &
     &       ,leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),                &
     &        rnisi(mm)                                   
 9022     format (2x,i4,1x,20a1,2(1pe10.3),10x,1pe10.3) 
          enddo 
        endif
!
!     find the rates affecting this element
!     step thru data types
      ml_data_type=0 
      do while (ml_data_type.lt.ntyp) 
!
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)'ml_data_type=',ml_data_type
        ml_data_type=ml_data_type+1 
!
!       loop over data
        ml_data=derivedpointers%npfi(ml_data_type,jkk_ion) 
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)'ml_data=',ml_data_type,ml_ion,jkk_ion,ml_data  
        ml_data_par=0
        if (ml_data.ne.0) ml_data_par=derivedpointers%npar(ml_data)
        if (lpri.gt.1)                                                  &
     &      write (lun11,*)'ml_data_par=',ml_data,ml_data_par,ml_ion
        do while ((ml_data.ne.0).and.(ml_data_par.eq.ml_ion)) 
!
!         step thru records of this type                   
          call drd(ltyp,lrtyp,lcon,                                     &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,ml_data,                  &
     &          0,lun11)     
          idest1=0 
          if (lrtyp.eq.7) idest1=masterdata%idat1(np1i+nidt-2) 
          if ((lrtyp.eq.1).or.(lrtyp.eq.15).or.(lrtyp.eq.8)             &
     &       .or.(lrtyp.eq.42)                                          &
     &       .or.(lrtyp.eq.6).or.((lrtyp.eq.7).and.(idest1.eq.1))) then
!           calculate rates                                         
            lpriu=min(1,lpri)
            abund1=0. 
            abund2=0.
            ptmp2=0.5
            ptmp1=0.5
            tau1=0.
            tau2=0.
            ptmp=(ptmp1+ptmp2) 
            call ucalc(ltyp,lrtyp,ml_data,lcon,jkk_ion,vturbi,cfrac,    &
     &                nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,          &
     &                ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,  &
     &                abund1,abund2,ptmp1,ptmp2,xpx,opakb1,             &
     &                opakc,opakcont,rccemis,lpriu,kdesc2,              &
     &                r,delr,t,trad,tsq,xee,xh1,xh0,                    &
     &                epi,ncn2,bremsa,bremsint,                         &
     &                leveltemp,                                        &
     &                rnisi,rnisse,nlev,lfpi,lun11,                     &
     &                np2,ncsvn,nlsvn)               
            if ((lrtyp.eq.1).or.(lrtyp.eq.15).or.                       &
     &          (lrtyp.eq.42).or.                                       &
     &          ((lrtyp.eq.7).and.(idest1.eq.1))) then
              pirti=pirti+ans1
              endif
            if ((lrtyp.eq.8).or.(lrtyp.eq.6)) then
              rrrti=rrrti+ans1
              endif
            if ((lpri.ge.1))                                            &
     &        write (lun11,9004)jkk_ion,lrtyp,ltyp,idest1,              &
     &        idest2,llo,lup,ml_data,ans1,ans2,rrrti,pirti
 9004       format(7x,7i6,i12,' ion  ',4(1pe10.3))                                  
            endif
!
!         end of loop over data
          ml_data=derivedpointers%npnxt(ml_data) 
          ml_data_par=0
          if (ml_data.ne.0) ml_data_par=derivedpointers%npar(ml_data)
          enddo 
!
!       end of loop over data types
        enddo 
!
      lpri=lprisv
!
      return 
      end                                           
