      subroutine calc_hmc_ion(ml_ion,lpri,lun11,                        &
     &                   vturbi,t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,    &
     &                   epi,ncn2,bremsa,bremsint,                      &
     &                   tau0,tauc,                                     &
     &                   np2,ncsvn,nlsvn,                               &
     &                   pirti,rrrti,                                   &
     &                   rnisi,                                         &
     &                   ajisi,cjisi,cjisi2,indbi,nindbi,nlev)

!                                                                       
!     Name: calc_hmc_ion.f90  
!     Description:  
!           Calculates all rates affecting level populations for 
!           one element.  Puts results into array used in abundance calculation
!           formerly:  func2
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
!           ajisi(2,ndb) Matrix of rates (s^-1) 
!             ajisb(1,n)=forward rate,ajisb(2,n)=reverse rate
!           cjisi(ndb): Matrix of cooling rates (erg s^-1)
!           indbi(2,ndb): Index array for ajisb, cjisb, relative to ion
!           nindb current length of ajisb
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
      integer idest1,idest2,idest3,idest4,nindbi 
!     continuum flux                                                    
      real(8) tauc(2,nnml) 
!                                                                       
      character(1) kblnk 
      real(8) rnisi(nd)
      real(8) xilevi(nd)
      real(8) rnisse(nd) 
      real(8) ajisi(2,ndb),cjisi(ndb),cjisi2(ndb) 
      integer indbi(2,ndb) 
      real(8) tau1,tau2,airtmp,e1,e2,pescl,pescv,ptmp
      real(8) pirti,rrrti
      integer np1i,np1r,np1k 
      integer nlev,ml,lpriu,mm,lprisv,                                  &
     &        llo,lup,ltyp,jkk_ion,                                     &
     &        lrtyp,lcon,nrdt,nidt,nkdt,kkkl,                           &
     &        ml_data,ml_ion,ml_data_type,ml_data_par
!                                                                       
      data kblnk/' '/ 
!                                                                       
      lprisv=lpri
!      if (lpri.ge.1) lpri=2
!
      if (lpri.ge.1)                                                    &
     &  write (lun11,901)t,xee,xpx,delr                              
901   format (1x,'      in calc_hmc_ion, inputs:',4(1pe10.3))
!                                                                       
!     lfpi value: photoionization and recombination, no opacities       
      lfpi=2 
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
      call calc_rates_level_lte(jkk_ion,lpri,lun11,t,xee,xpx,nlev)
      nlev=derivedpointers%nlevs(jkk_ion)
      if (lpri.ge.1) then 
        write (lun11,*)'      nlev=',nlev
        write (lun11,*)'  index level             energy',              &
     &              '  stat.wt. population   LTE ' 
        do mm=1,nlev 
          write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20)          &
     &       ,leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),                &
     &        xilevi(mm),rnisi(mm),rnisse(mm)                                   
 9022     format (2x,i4,1x,20a1,5(1pe10.3)) 
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
          if (lpri.gt.1)                                                &
     &      write (lun11,*)'ml_data_par=',ml_data,ml_data_par,ml_ion
        do while ((ml_data.ne.0).and.(ml_data_par.eq.ml_ion)) 
!
!           step thru records of this type                   
            call drd(ltyp,lrtyp,lcon,                                   &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_data,                &
     &            0,lun11)     
            if (.not.((lrtyp.eq.1).and.(ltyp.eq.53)).and.(lrtyp.ne.15)  &
     &          .and.(lrtyp.ne.8))   then
!            if ((lrtyp.eq.7).or.(lrtyp.eq.5).or.(lrtyp.eq.40).or.       &
!     &        (lrtyp.eq.3).or.(lrtyp.eq.23).or.(lrtyp.eq.4).or.         &
!     &        (lrtyp.eq.14).or.(lrtyp.eq.9).or.(lrtyp.eq.41).or.        &
!     &        (lrtyp.eq.42)) then
!           calculate rates                                         
!            lpriu=min(1,lpri)
            lpriu=lpri
            abund1=0. 
            abund2=0.
            ptmp2=0.5
            ptmp1=0.5
            tau1=0.
            tau2=0.
            if (lrtyp.eq.4) then
              kkkl=derivedpointers%nplini(ml) 
              if ((kkkl.gt.0).and.(kkkl.le.nnnl)) then
                tau1=tau0(1,kkkl) 
                tau2=tau0(2,kkkl) 
                ptmp1=pescl(tau1)*(1.-cfrac)                            
                ptmp2=pescl(tau2)*(1.-cfrac)                         &
     &              +2.*pescl(tau1+tau2)*cfrac 
!                ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac     
                endif
              endif
            if (lrtyp.eq.7) then
              kkkl=derivedpointers%npconi2(ml_data) 
              if ((kkkl.gt.0).and.(kkkl.le.ndat2)) then 
                 tau1=tauc(1,kkkl) 
                 tau2=tauc(2,kkkl) 
                 ptmp1=pescv(tau1)*(1.-cfrac) 
                 ptmp2=pescv(tau2)*(1.-cfrac)                        &
     &                  +2.*pescv(tau1+tau2)*cfrac                      
                 endif
               endif
            ptmp=(ptmp1+ptmp2) 
!           note that radiative rates and emissivities will     
!           have the escape probabilities in them             
!           when calculated in ucalc                          
            call ucalc(ltyp,lrtyp,ml_data,lcon,jkk_ion,vturbi,cfrac,    &
     &                nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,          &
     &                ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,  &
     &                abund1,abund2,ptmp1,ptmp2,xpx,opakb1,             &
     &                opakc,opakcont,rccemis,lpriu,kdesc2,              &
     &                r,delr,t,trad,tsq,xee,xh1,xh0,                    &
     &                epi,ncn2,bremsa,bremsint,                         &
     &                rnisi,rnisse,nlev,lfpi,lun11,                     &
     &                np2,ncsvn,nlsvn)               
!           this Statement prevents double counting of PI from  
!           excited levels.  All rate type 1 data should not h
!           lower level associated with it, but some does.    
            if ((lrtyp.eq.1).and.(idest1.ne.1)) ans1=0. 
!            if (lrtyp.eq.1) then 
!               ans2=0. 
!               ans4=0. 
!               endif 
            if ((lrtyp.eq.7).or.(lrtyp.eq.1)) then
               if (idest1.eq.1) pirti=pirti+ans1
               rrrti=rrrti+ans2
               endif
            llo=idest1 
            lup=idest2 
            if ((idest1.gt.0).and.(idest2.gt.0).and.                    &
     &         (llo.le.nd).and.(lup.le.nd)) then                 
               e1=leveltemp%rlev(1,idest1) 
               e2=leveltemp%rlev(1,idest2) 
               if ((lrtyp.ne.7).and.(lrtyp.ne.41)) then
                 if ((e1/(1.d-24+e2)-1.).lt.1.e-8)  then 
                     lup=idest2 
                     llo=idest1 
                     if ((lpri.gt.1).and.(lrtyp.eq.4))                  &
      &                write (lun11,*)'not switching energies',         &
      &                 idest1,idest2,llo,lup,e1,e2
                   else 
                     lup=idest1 
                     llo=idest2 
                     if ((lpri.gt.1).and.(lrtyp.eq.4))                  &
      &                write (lun11,*)'switching energies',             &
      &                 idest1,idest2,llo,lup,e1,e2
                   endif 
                 endif
              airtmp=ans2 
              nindbi=nindbi+1 
              if (nindbi.gt.ndb) stop 'array indexing error' 
              ajisi(1,nindbi)=ans1 
              ajisi(2,nindbi)=airtmp 
              cjisi(nindbi)=0. 
              cjisi2(nindbi)=0. 
              indbi(1,nindbi)=lup
              indbi(2,nindbi)=llo
              if (lpri.gt.1)                                           &
     &           write (lun11,*)nindbi,indbi(1,nindbi),                &
     &           indbi(2,nindbi),ajisi(1,nindbi),ajisi(2,nindbi),      &
     &           cjisi(nindbi)            
              nindbi=nindbi+1 
              if (nindbi.gt.ndb) stop 'array indexing error' 
              ajisi(1,nindbi)=airtmp 
              ajisi(2,nindbi)=ans1 
              cjisi(nindbi)=0. 
              cjisi2(nindbi)=0. 
              indbi(1,nindbi)=llo
              indbi(2,nindbi)=lup
              if (lpri.gt.1)                                            &
                write (lun11,*)nindbi,indbi(1,nindbi),indbi(2,nindbi),  &
     &             ajisi(1,nindbi),ajisi(2,nindbi),cjisi(nindbi)              
              nindbi=nindbi+1 
              if (nindbi.gt.ndb) stop 'array indexing error' 
              ajisi(1,nindbi)=-ans1 
              ajisi(2,nindbi)=-ans1 
!              cjisi(nindbi)=-ans3
              cjisi(nindbi)=ans4*xpx
              cjisi2(nindbi)=ans6*xpx
              indbi(1,nindbi)=llo
              indbi(2,nindbi)=llo
              if (lpri.gt.1)                                            &
                write (lun11,*)nindbi,indbi(1,nindbi),indbi(2,nindbi),  &
                       ajisi(1,nindbi),ajisi(2,nindbi),cjisi(nindbi)            
              nindbi=nindbi+1 
              if (nindbi.gt.ndb) stop 'array indexing error' 
              ajisi(1,nindbi)=-airtmp 
              ajisi(2,nindbi)=-airtmp 
              cjisi(nindbi)=-ans3*xpx
              cjisi2(nindbi)=-ans5*xpx
              indbi(1,nindbi)=lup
              indbi(2,nindbi)=lup
              if (lpri.gt.1)                                            &
                write (lun11,*)nindbi,indbi(1,nindbi),indbi(2,nindbi),  &
                       ajisi(1,nindbi),ajisi(2,nindbi),cjisi(nindbi)  
              if ((lpri.ge.1))                                          &
     &           write (lun11,9044)jkk_ion,lrtyp,ltyp,idest1,           &
     &                  idest2,llo,lup,ml_data,ans1,ans2,ans3,ans4
!     &                  idest2,llo,lup,ml_data,ans1,ans2,ans3,ans4,     &
!     &                    ans5,ans6,eth,kkkl,                           &
!     &                    ptmp1,ptmp2                              
! 9004            format(7x,7i6,i8,' level',7(1pe10.3),i6,7(1pe10.3))                                  
 9044            format(7x,7i6,i12,' level',4(1pe10.3))
              endif 
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
