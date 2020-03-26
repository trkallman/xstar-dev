      subroutine calc_emis_ion(ml_ion,lpri,lun11,xeltp,                 &
     &       vturbi,critf,t,trad,r,delr,xee,xpx,xh0,xh1,cfrac,p,lcdd,   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xilevi,bilevi,rnisi,                                       &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab)                         
!                                                                       
!                                                                       
!     Name: calc_emis_ion.f90  
!     Description:  
!           Calculates emissivities and opacitiesfor 
!           one ion.  
!
!     List of Parameters:
!           Input:
!           ml_ion: poiter to ion record
!           lpri: print switch, 1=on, 0=off
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
!           rniss:  lte level populations
!           rnisse:  lte level populations with exponential removed
!           nmat: index into rates array corresponding to first level
!           nlev:  number of levels for the ion
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
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) xilevi(nd),rnisi(nd),bilevi(nd)
      real(8) rnisse(nd)
!     element abundances                                                
!     state variables                                                   
      real(8) r,t,xpx,delr 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee 
      integer ncn2,lpri,lun11,lfpi,np2,lcdd
      integer nlsvn,ncsvn 
      character(49) kdesc2 
      real(8) tsq,ans1,ans2,xh1,xh0,cfrac,critf,p,xeltp
      real(8) abund1,abund2,ptmp1,ptmp2,ans3,ans4,ans5,ans6,opakb1 
      integer idest1,idest2,idest3,idest4
!     continuum flux                                                    
      real(8) tauc(2,nnml) 
!                                                                       
      character(1) kblnk 
      real(8) tau1,tau2,e1,e2,pescl,pescv,eth 
      real(8) htt,cll,ergsev,rcemm
      integer np1i,np1r,np1k 
      integer nlev,ml,lpriu,lprisv,                                     &
     &        llo,lup,ltyp,jkk_ion,                                     &
     &        lrtyp,lcon,nrdt,nidt,nkdt,kkkl,jkkl,                      &
     &        ml_data,ml_ion,ml_data_type,ml_data_par,mm
!                                                                       
      data kblnk/' '/ 
      data ergsev/1.602197e-12/
!     
      lprisv=lpri
!      if (lpri.ge.1) lpri=2
!
      if (lpri.ge.1)                                                    &
     &  write (lun11,901)t,xee,xpx,lcdd,p,delr                              
901   format (1x,'      in calc_emis_ion, inputs:',3(1pe10.3),          &
     &i6,2(1pe10.3))
!                                                                       
!     lfpi value: photoionization and recombination, no opacities       
      lfpi=2 
!                                                                       
      tsq=sqrt(t) 
!                                                                       
!     get ion index
      call drd(ltyp,lrtyp,lcon,                                         &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_ion,                 &
     &            0,lun11)                                        
      jkk_ion=masterdata%idat1(np1i+nidt-1)
      if (lpri.ge.1)                                                    &
     &            write (lun11,903)jkk_ion,ml_ion,                      &
     &               (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)
903             format (1x,'      ion:',i12,1x,i12,1x,8(1a1))
!
      call calc_rates_level_lte(jkk_ion,lpri,lun11,t,xee,xpx,           &
     &              leveltemp,nlev)
!
      if (lpri.gt.1) then 
                  write (lun11,*)'      index level        energy',     &
     &              ' stat.wt. population LTE population' 
        do mm=1,nlev 
          write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20)          &
     &       ,leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),                &
     &       xilevi(mm),rnisi(mm)                                    
 9022     format (4x,i4,1x,20a1,4(1pe10.3)) 
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
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_data,                &
     &            0,lun11)     
!
          if (ml_data_type.eq.7) then
            idest1=masterdata%idat1(np1i+nidt-2) 
            idest2=nlev+masterdata%idat1(np1i-1+nidt-3)-1 
            kkkl=derivedpointers%npconi2(ml_data) 
            if ((kkkl.ne.0).and.(kkkl.le.ndat2)                         &
     &         .and.(idest1.gt.0)) then                                     
              llo=idest1
              lup=idest2
              eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
              abund1=xilevi(llo)*xeltp 
              abund2=xilevi(lup)*xeltp 
              tau1=tauc(1,kkkl) 
              tau2=tauc(2,kkkl) 
              ptmp1=pescv(tau1)*(1.-cfrac)                                    
              ptmp2=pescv(tau2)*(1.-cfrac)                           &
     &            +2.*pescv(tau1+tau2)*cfrac 
              lpriu=lpri
              call ucalc(                                               &
     &            ltyp,lrtyp,ml_data,lcon,jkk_ion,vturbi,cfrac,         &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,              &
     &            ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,      &
     &            abund1,abund2,ptmp1,ptmp2,xpx,opakab(kkkl),           &
     &            opakc,opakcont,rccemis,lpriu,kdesc2,                  &
     &            r,delr,t,trad,tsq,xee,xh1,xh0,                        &
     &            epi,ncn2,bremsa,bremsint,                             &
     &            leveltemp,                                            &
     &            rnisi,rnisse,nlev,lfpi,lun11,                         &
     &            np2,ncsvn,nlsvn)               
              cabab(kkkl)=ans3 
              cabab(kkkl)=cabab(kkkl)*abund1*xpx 
              cemab(1,kkkl)=ptmp1*ans4/(ptmp1+ptmp2) 
              cemab(2,kkkl)=ptmp2*ans4/(ptmp1+ptmp2) 
              cemab(1,kkkl)=cemab(1,kkkl)*abund2*xpx 
              cemab(2,kkkl)=cemab(2,kkkl)*abund2*xpx 
              cll=cll+cemab(1,kkkl)+cemab(2,jkkl) 
              htt=htt+abund1*ans3 
              if (lpri.ge.1)                                            &
     &            write (lun11,9002)jkk_ion,lrtyp,ltyp,                 &
     &              idest1,idest2,llo,lup,ml_data,ans1,ans2,            &
     &              ans3,ans4,cemab(1,kkkl)+cemab(2,kkkl),              &
     &              opakab(kkkl),eth,                                   &
     &              kkkl,htt,cll,ptmp1,ptmp2                           
 9002             format (7x,7i6,i12,' h-c ',                           &
     &              7(1pe10.3),i12,2(1pe10.3),4(1pe10.3))                    
              endif
            endif
!
          if ((ml_data_type.eq.4).or.                                   &
     &         (ml_data_type.eq.14).or.                                 &
     &         (ml_data_type.eq.9)) then
            if ((ml_data_type.eq.4).or.(ml_data_type.eq.9)) then
              idest1=masterdata%idat1(np1i) 
              idest2=masterdata%idat1(np1i+1) 
              jkkl=derivedpointers%nplini(ml_data) 
              tau1=tau0(1,jkkl) 
              tau2=tau0(2,jkkl) 
              endif
            if (ml_data_type.eq.14) then
              idest1=masterdata%idat1(np1i-1+nidt-3) 
              idest2=masterdata%idat1(np1i+nidt-3) 
              jkkl=0
              tau1=0.
              tau2=0.
              endif
            if ((idest1.gt.0).and.(idest2.gt.0).and.                    &
     &        (idest1.lt.nlev).and.(idest2.lt.nlev)) then             
              e1=leveltemp%rlev(1,idest1) 
              e2=leveltemp%rlev(1,idest2) 
              if (e1.lt.e2) then 
                  lup=idest2
                  llo=idest1
                else 
                  lup=idest1
                  llo=idest2
                endif 
              abund1=xilevi(llo)*xpx*xeltp 
              abund2=xilevi(lup)*xpx*xeltp 
              ptmp1=pescl(tau1)*(1.-cfrac)                            
              ptmp2=pescl(tau2)*(1.-cfrac)                           &
     &            +2.*pescl(tau1+tau2)*cfrac 
!             ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac     
!              we need to call ucalc again because rcem               
!              already has the abundance in from func3p               
              lpriu=lpri
              call ucalc(                                               &
     &             ltyp,lrtyp,ml_data,lcon,jkk_ion,vturbi,cfrac,        &
     &             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,             &
     &             ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,     &
     &             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,                &
     &             opakc,opakcont,rccemis,lpriu,kdesc2,                 &
     &             r,delr,t,trad,tsq,xee,xh1,xh0,                       &
     &             epi,ncn2,bremsa,bremsint,                            &
     &             leveltemp,                                           &
     &             rnisi,rnisse,nlev,lfpi,lun11,                        &
     &             np2,ncsvn,nlsvn)               
             if (ml_data_type.eq.14) then
               rcemm=abund2*ans4 
               rccemis(2,3)=rccemis(2,3)+                               &
     &                  rcemm/(epi(4)-epi(3)+1.e-24)/ergsev/12.56       
!               if (lpri.gt.0) write (lun11,*)jkkl,tau0(1,jkkl),       
                cll=cll+rcemm 
!                clcont=clcont+rcemm 
                endif
              if (ml_data_type.eq.4) then
  !             note need to change from v2.55 owing to change in order ans3 ans4
                rcem(1,jkkl)=-abund2*ans3*ptmp1/(ptmp1+ptmp2) 
                rcem(2,jkkl)=-abund2*ans3*ptmp2/(ptmp1+ptmp2) 
                oplin(jkkl)=opakb1*abund1 
                cll=cll+rcem(1,jkkl)+rcem(2,jkkl) 
                htt=htt+abund1*ans3 
                endif
              if ((lpri.ge.1))                                          &
     &               write (lun11,9002)jkk_ion,lrtyp,                   &
     &                 ltyp,idest1,idest2,                              &
     &                 llo,lup,ml_data,ans1,ans2,ans3,ans4,             &
     &                 rcem(1,jkkl)+rcem(2,jkkl),oplin(jkkl),           &
     &                 masterdata%rdat1(np1r),jkkl,cll,htt              &
     &                  ,ptmp1,ptmp2                                    
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
      if (lpri.gt.0) write (lun11,*)'returning from calc_emis_ion'
      lpri=lprisv
!
      return 
      end                                           
