      subroutine calc_emis_ion(ml_ion,lpri,lun11,xeltp,                 &
     &       vturbi,critf,t,trad,r,delr,xee,xpx,xh0,xh1,cfrac,p,lcdd,   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xileve,bileve,rnise,ipmat,                                 &
     &       rcem,oplin,brcems,rccemis,opakc,opakcont,cemab,            &
     &       cabab,opakab,elin,errc,nlbin,ncbin)                         
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
      real(8) xileve(nd),rnise(nd),bileve(nd)
!     element abundances                                                
!     state variables                                                   
      real(8) r,t,xpx,delr 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee 
      integer ncn2,lpri,lun11,lfpi,np2,lcdd
      integer nlsvn,ncsvn 
      integer ipmat
      character(49) kdesc2 
      real(8) tsq,ans1,ans2,xh1,xh0,cfrac,critf,p,xeltp
      real(8) abund1,abund2,ptmp1,ptmp2,ans3,ans4,ans5,ans6,opakb1 
      integer idest1,idest2,idest3,idest4
!     continuum flux                                                    
      real(8) tauc(2,nnml) 
      real(8) elin(nnnl)
      real(8) errc(nnml)
      integer nlbin(nrank,ncn),ncbin(nrank,ncn)
!                                                                       
      character(1) kblnk 
      real(8) tau1,tau2,e1,e2,pescl,pescv,eth,ener
      real(8) htt,cll,ergsev,rcemm,cemtmp1,cemtmp2
      integer np1i,np1r,np1k 
      integer nlev,ml,lpriu,lprisv,                                     &
     &        llo,lup,ltyp,jkk_ion,                                     &
     &        lrtyp,lcon,nrdt,nidt,nkdt,kkkl,jkkl,                      &
     &        ml_data,ml_ion,ml_data_type,ml_data_par,mm
      integer nnzz,nnnn,nbinc,nb1
!                                                                       
      data kblnk/' '/ 
      save kblnk
      data ergsev/1.602197e-12/
      save ergsev
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
      nnzz=masterdata%idat1(np1i+1)
      nnnn=nnzz-masterdata%idat1(np1i)+1
      jkk_ion=masterdata%idat1(np1i+nidt-1)
      if (lpri.ge.1)                                                    &
     &            write (lun11,903)jkk_ion,ml_ion,                      &
     &               (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)
903             format (1x,'      ion:',i12,1x,i12,1x,8(1a1))
!
      call calc_rates_level_lte(jkk_ion,lpri,lun11,t,xee,xpx,           &
     &              nnzz,nnnn,leveltemp,nlev)
!
      if (lpri.gt.1) then 
                  write (lun11,*)'      index level        energy',     &
     &              ' stat.wt. population LTE population' 
        do mm=1,nlev 
          write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20)          &
     &       ,leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),                &
     &       xileve(mm+ipmat),rnise(mm+ipmat)                                    
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
!         test if rate type 7
          if (ml_data_type.eq.7) then
!
!           test if continuum pointer ok
            idest1=masterdata%idat1(np1i+nidt-2) 
            idest2=nlev+masterdata%idat1(np1i-1+nidt-3)-1 
            kkkl=derivedpointers%npconi2(ml_data) 
!            if (lpri.ne.0) write (lun11,*)'testing 7:',kkkl,idest1,     &
!     &        errc(kkkl)
            if ((kkkl.ne.0).and.(kkkl.le.ndat2)                         &
     &         .and.(idest1.gt.0)) then 
!
              if ((errc(kkkl).gt.epi(1)).and.(errc(kkkl).lt.epi(ncn2))) &
     &           then                                    
!
                llo=idest1+ipmat
                lup=idest2+ipmat
!                eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
                eth=12398.54/errc(kkkl)
                nb1=nbinc(eth,epi,ncn2) 
                if (lpri.gt.1) write (lun11,*)'testing rrc in func3',   &
     &            kkkl,errc(kkkl),nb1,ncbin(1,nb1),ncbin(2,nb1)
!
!               test to see if in strong rrc list
                mm=1
                if (lpri.ge.2) write (lun11,*)'test to see if in list', &
     &              mm,kkkl,ncbin(mm,nb1)
                do while ((ncbin(mm,nb1).ne.0).and.                     &
     &            (ncbin(mm,nb1).ne.kkkl).and.(mm.lt.nrank))
                  mm=mm+1
                  if (lpri.gt.1) write (lun11,*)mm,kkkl,ncbin(mm,nb1)
                  enddo
!               nb a fudge for testing
!               ncbin(mm,nb1)=kkkl        
                if ((ncbin(mm,nb1).eq.kkkl)                             &
     &           .or.(ncbin(1,nb1).eq.9999999)) then
!
                  abund1=xileve(llo)*xeltp 
                  abund2=xileve(lup)*xeltp 
                  tau1=tauc(1,kkkl) 
                  tau2=tauc(2,kkkl) 
                  ptmp1=pescv(tau1)*(1.-cfrac)                                    
                  ptmp2=pescv(tau2)*(1.-cfrac)                          &
     &            +2.*pescv(tau1+tau2)*cfrac 
                  lpriu=lpri
                  ans1=0.
                  ans2=0.
                  ans3=0.
                  ans4=0.
                  ans5=0.
                  ans6=0.
                  call ucalc(                                           &
     &              ltyp,lrtyp,ml_data,lcon,jkk_ion,vturbi,cfrac,       &
     &              nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,            &
     &              ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,    &
     &              abund1,abund2,ptmp1,ptmp2,xpx,opakab(kkkl),         &
     &              opakc,opakcont,rccemis,lpriu,kdesc2,                &
     &              r,delr,t,trad,tsq,xee,xh1,xh0,                      &
     &              epi,ncn2,bremsa,bremsint,                           &
     &              leveltemp,                                          &
     &              nlev,lfpi,lun11,                                    &
     &              np2,ncsvn,nlsvn)             
                  if (lpri.ge.1)                                        &
     &            write (lun11,9002)jkk_ion,lrtyp,ltyp,                 &
     &              idest1,idest2,llo,lup,ml_data,ans1,ans2,            &
     &              ans3,ans4,cemab(1,kkkl)+cemab(2,kkkl),              &
     &              opakab(kkkl),eth,                                   &
     &              kkkl,htt,cll,ptmp1,ptmp2,abund1,abund2                                   
 9002             format (7x,7i6,i12,' h-c ',                           &
     &              7(1pe10.3),i12,2(1pe10.3),6(1pe10.3))               
!
!                 end of test if rrc in strong list
                  endif     
!
!              end of test if rrc energy in range
               endif
!
!             end of test if continuum pointer ok
              endif
!
!           end of test if rate type 7
            endif
!
!         test if rate type 9
          if (ml_data_type.eq.9) then
!
!           test if pointer ok
            idest1=masterdata%idat1(np1i) 
            idest2=masterdata%idat1(np1i+1) 
!            write (lun11,*)ml,ltyp,lrtyp,idest1,idest2,jkkl,nlev
            if ((masterdata%rdat1(np1r).gt.0.01)                        &
     &        .and.(idest1.gt.0).and.(idest2.gt.0).and.                 &
     &        (idest1.lt.nlev).and.(idest2.lt.nlev)) then                 
!
                  llo=idest1+ipmat
                  lup=idest2+ipmat
                  eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
                  nb1=nbinc(eth,epi,ncn2) 
!
                  abund1=xileve(llo)*xeltp*xpx
                  abund2=xileve(lup)*xeltp*xpx
                  tau1=0.
                  tau2=0.
                  ptmp1=                                                &
     &             pescl(tau1)*(1.-cfrac)                                    
                  ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac 
                  lpriu=lpri
                  ans1=0.
                  ans2=0.
                  ans3=0.
                  ans4=0.
                  ans5=0.
                  ans6=0.
                  call ucalc(                                           &
     &              ltyp,lrtyp,ml_data,lcon,jkk_ion,vturbi,cfrac,       &
     &              nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,            &
     &              ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,    &
     &              abund1,abund2,ptmp1,ptmp2,xpx,opakab(kkkl),         &
     &              opakc,opakcont,rccemis,lpriu,kdesc2,                &
     &              r,delr,t,trad,tsq,xee,xh1,xh0,                      &
     &              epi,ncn2,bremsa,bremsint,                           &
     &              leveltemp,                                          &
     &              nlev,lfpi,lun11,                                    &
     &              np2,ncsvn,nlsvn)             
                  cemtmp1=abund2*ptmp1*ans4 
                  cemtmp2=abund2*ptmp2*ans4 
                  if (lpri.ge.1)                                        &
     &            write (lun11,9002)jkk_ion,lrtyp,ltyp,                 &
     &              idest1,idest2,llo,lup,ml_data,ans1,ans2,            &
     &              ans3,ans4,cemtmp1+cemtmp2,opakb1,eth,               &
     &              jkkl,htt,cll,ptmp1,ptmp2,abund1,abund2                                   
!
!                 end of test if continuum pointer ok
                  endif
!
!           end of test if rate type 9
            endif
!
!         test if rate type 4
          if ((ml_data_type.eq.4).or.                                   &
     &         (ml_data_type.eq.9)) then
!
!           check to see if pointer is sensible
            idest1=masterdata%idat1(np1i) 
            idest2=masterdata%idat1(np1i+1) 
            jkkl=derivedpointers%nplini(ml_data) 
            if ((jkkl.ne.0).and.(jkkl.le.ndat2)                         &
     &         .and.(idest1.gt.0)) then 
!
!             check to see if pointer is in strong line list
              ener=12398.41/(elin(jkkl)+1.d-36)
              nb1=nbinc(ener,epi,ncn2)
              mm=1
              if (lpri.gt.1) write (lun11,*)mm,jkkl,nlbin(mm,nb1)
              do while ((nlbin(mm,nb1).ne.0)                            &
     &          .and.(nlbin(mm,nb1).ne.jkkl).and.(mm.lt.nrank))
                if (lpri.gt.1) write (lun11,*)mm,jkkl,nlbin(mm,nb1)
                mm=mm+1
                enddo
!             nb a fudge for testing
!             nlbin(mm,nb1)=jkkl        
              if (lpri.gt.1) write (lun11,*)mm,jkkl,nlbin(mm,nb1)
              if ((nlbin(mm,nb1).eq.jkkl)                               &
     &          .or.(nlbin(1,nb1).eq.9999999)) then
!
                tau1=tau0(1,jkkl) 
                tau2=tau0(2,jkkl) 
                e1=leveltemp%rlev(1,idest1) 
                e2=leveltemp%rlev(1,idest2) 
                if (e1.lt.e2) then 
                  lup=idest2+ipmat
                  llo=idest1+ipmat
                  else 
                  lup=idest1+ipmat
                  llo=idest2+ipmat
                  endif 
                abund1=xileve(llo)*xpx*xeltp 
                abund2=xileve(lup)*xpx*xeltp 
!                if (lpri.ne.0) write (lun11,*)'abund1,abund2:',         &
!     &            abund1,abund2,llo,lup,xileve(llo),xileve(lup),        &
!     &           idest1,idest2,ipmat,xpx,xeltp
                ptmp1=pescl(tau1)*(1.-cfrac)                            
                ptmp2=pescl(tau2)*(1.-cfrac)                            &
     &            +2.*pescl(tau1+tau2)*cfrac 
!               ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac     
!               we need to call ucalc again because rcem               
!               already has the abundance in from func3p               
                lpriu=lpri
                ans1=0.
                ans2=0.
                ans3=0.
                ans4=0.
                ans5=0.
                ans6=0.
                call ucalc(                                             &
     &             ltyp,lrtyp,ml_data,lcon,jkk_ion,vturbi,cfrac,        &
     &             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,             &
     &             ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,     &
     &             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,                &
     &             opakc,opakcont,rccemis,lpriu,kdesc2,                 &
     &             r,delr,t,trad,tsq,xee,xh1,xh0,                       &
     &             epi,ncn2,bremsa,bremsint,                            &
     &             leveltemp,                                           &
     &             nlev,lfpi,lun11,                                     &
     &             np2,ncsvn,nlsvn)               
                if ((lpri.ge.1))                                        &
     &               write (lun11,9002)jkk_ion,lrtyp,                   &
     &                 ltyp,idest1,idest2,                              &
     &                 llo,lup,ml_data,ans1,ans2,ans3,ans4,             &
     &                 rcem(1,jkkl)+rcem(2,jkkl),oplin(jkkl),           &
     &                 masterdata%rdat1(np1r),jkkl,cll,htt              &
     &                  ,ptmp1,ptmp2,abund1,abund2                                    
!
!               end of test for in strong line list
                endif
!
!             end of test of pointer sensible
              endif
!
!           end of test if rate type
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
      if (lpri.gt.0)                                                    &
     &  write (lun11,*)'     returning from calc_emis_ion'
      lpri=lprisv
!
      return 
      end                                           
