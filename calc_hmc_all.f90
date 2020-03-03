      subroutine calc_hmc_all(lpri,lun11,vturbi,critf,                  &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,zeta,              &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xiin,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,elcter,&
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       httot2,cltot2,                                             &
     &       xilevg,bilevg,rnisg)
                                                                        
!                                                                       
!     Name: calc_hmc_all.f90  
!     Description:  
!           Master routine which steps throuch elements and ions, 
!           calls routines which calculate rates,
!           calculate ion fractions, level populations, heating cooling
!           this is called from within the heating=cooling loop
!           formerly func
!
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
!           Output:
!           xiin(nni):  ion fractions, xiin(1)=H, xiin(2)=He0, xiin(3)=He+ etc
!           rrrts(nni): total recombination rates for each ion (s^-1)
!           pirts(nni): total photoionization rates for each ion(s^-1)
!           htt(nni): total heating rate for each ion (approximate) 
!                       (erg s^-1 cm^-3)
!           cll(nni): total cooling rate for each ion (approximate) 
!           httot: total heating rate (erg s^-1 cm^-3) 
!           cltot: total cooling rate (erg s^-1 cm^-3) 
!           hmctot:  (httot-cltot)*2./(httot+cltot)
!           elcter:  charge conservation error (relative to H)
!           cllines:  total cooling rate due to lines (erg s^-1 cm^-3) 
!           clcont:  total cooling rate due to continuum (erg s^-1 cm^-3) 
!           cllines:  total cooling rate due to lines (erg s^-1 cm^-3) 
!           htcomp:  compton heating rate (erg s^-1 cm^-3) 
!           clcomp:  compton cooling rate (erg s^-1 cm^-3) 
!           clbrems:  bremsstrahlung cooling rate (erg s^-1 cm^-3) 
!           xilevg(nnml):  level populations relative to all elements
!                (note that xilevt which is relative to parent element
!                      is used internally)
!           rnisg(nnml): lte level populations relative to all elements
!                (note that rnist which is relative to parent element
!                      is used internally)
!           bilevg(nnml): departure coefficient relative to all elements
!                (note that bilevt which is relative to parent element
!                      is used internally)
!           also uses variables from globaldata
!           
!        Dependencies:  Calls calc_ion_rates,calc_rates_level,
!                   calc_num_level,calc_level_rates_lte,istruc,
!                   msolvelucy,chisq,comp2,bremem,heatf
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
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn)
!     level populations                                                 
      real(8) xilevg(nnml),rnisg(nnml),bilevg(nnml)
!     ion abundances                                                    
      real(8) xiin(nni) 
      real(8) rrrt(nni),pirt(nni)
      real(8) rrrti(nl),pirti(nl)
      real(8) tauc(2,nnml) 
      real(8) htt(nni),cll(nni) 
      real(8) htt2(nni),cll2(nni) 
!     element abundances                                                
      real(8) abel(nl) 
!     limits on ion indeces vs element
      integer mml(nl),mmu(nl)
!                                                                       
      real(8) xii(nl)
      real(8) xileve(nd),rnise(nd),bileve(nd) 
      integer lpri,lun11,lcdd,ncn2,np2,ncsvn,nlsvn 
      real(8) vturbi,critf,t,trad,r,delr,xee,xpx,cfrac,p,zeta,          &
     &     hmctot,elcter,cllines,clcont,htcomp,clcomp,clbrems           
      real(8) xh1,xh0,httot,cltot,httot2,cltot2 
      real(8) cmp1,cmp2,                                                &
     &     enelec,                                                      &
     &     xeltp,                                                       &
     &     xisum,cl,ht,cl2,ht2                                             
      integer nlev,klion,                                               &
     &     jkk,ipmat,ltyp,                                              &
     &     lrtyp,lcon,nrdt,nidt,nkdt,ipmatsv,                           &
     &     jk,mm,                                                       &
     &     mlion,mmtmp,                                                 &
     &     nnz,mlm,np1i,np1r,np1k,lprisv                 
      integer ml_element_data_type,ml_element           
      integer ml_ion_data_type,ml_ion,ml_element_test,ml
!                                                                       
!            
      lprisv=lpri 
      if (lpri.ge.1)                                                    &
     &  write (lun11,901)t,xee,xpx,lcdd,p,delr                              
901   format (1x,'in calc_hmc_all, inputs:',3(1pe10.3),i6,2(1pe10.3))
       if (lcdd.ne.1)                                                   &
     &   xpx = p/1.38e-12/max(t,1.d-24)                                 
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
      httot=0.
      cltot=0.
      httot2=0.
      cltot2=0.
      enelec=0.
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
!        write (lun11,*)'jk=',jk,nnz,abel(jk)
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
            ml_element_test=derivedpointers%npar(ml_ion)
            enddo
!
!         now fill rate matrix for element
          call calc_hmc_element(ml_element,lpri,lun11,                  &
     &                 critf,vturbi,t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,&
     &                 zeta,mml,mmu,                                    &
     &                 epi,ncn2,bremsa,bremsint,                        &
     &                 tau0,tauc,                                       &
     &                 np2,ncsvn,nlsvn,                                 &
     &                 rnise,bileve,xileve,cl,ht,xii,rrrti,pirti)
!
!                                                                       
!         step thru levels and store populations
!         test if first ion pointer
          ipmatsv=0
          xisum=0.
          ml_ion_data_type=12
          ml_ion=derivedpointers%npfirst(ml_ion_data_type)
          do while  (ml_ion.ne.0)
!
            if (lpri.gt.1)                                              &
     &         write (lun11,*)'ml_ion=',ml_ion
!
!           test if element belongs to parent of ion
            ml_element_test=derivedpointers%npar(ml_ion)
            if (lpri.gt.1)                                              &
     &      write (lun11,*)'ml_element_test=',ml_element_test,ml_element
            if (ml_element_test.eq.ml_element) then
!
!               retrieve ion name from kdati                              
                mlm=ml_ion
                call drd(ltyp,lrtyp,lcon,                               &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                    &
     &            0,lun11)                                        
                klion=masterdata%idat1(np1i)
                jkk=masterdata%idat1(nidt+np1i-1)
                xiin(jkk)=xii(klion)
                rrrt(jkk)=rrrti(klion)
                pirt(jkk)=pirti(klion)
                if (lpri.gt.1)                                          &
     &           write (lun11,*)'saving xii:',klion,jkk,xii(klion),     &
     &               rrrti(klion),pirti(klion)
                xisum=xisum+xii(klion)
                enelec=enelec+xii(klion)*float(klion-1)*xeltp
                if (lpri.ge.1)                                          &
     &            write (lun11,903)jkk,ml_ion,                          &
     &               (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)
!
!               get level data                                          
                if (lpri.gt.1) write (lun11,*)'jkk=',jkk
                nlev=derivedpointers%nlevs(jkk)
                do mm=1,nlev-1
!                 get level pointer                                     
                  mmtmp=derivedpointers%npilev(mm,jkk) 
                  rnisg(mmtmp)=rnise(mm+ipmatsv)
                  bilevg(mmtmp)=bileve(mm+ipmatsv)
                  xilevg(mmtmp)=xileve(mm+ipmatsv)
                  enddo
                mmtmp=derivedpointers%npilev(nlev,jkk) 
                if (mmtmp.gt.nnml) stop 'mmtmp error' 
                xilevg(mmtmp)=xileve(ipmatsv+nlev) 
                rnisg(mmtmp)=rnise(ipmatsv+nlev) 
                bilevg(mmtmp)=xilevg(mmtmp)/(rnisg(mmtmp)+1.d-48)
!
                if (lpri.ge.1) then 
                  write (lun11,*)'level populations:' 
                  call calc_rates_level_lte(jkk,lpri,lun11,t,xee,xpx,   &
     &                       nlev)
                  do mm=1,nlev 
                    write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20)&
     &                       ,leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),&
     &                       xileve(mm+ipmatsv),rnise(mm+ipmatsv)        
                    enddo 
                  endif 
!
                ipmatsv=ipmatsv+nlev-1
!
!               end of test if element belongs to parent of ion
                endif
!
!           end of step thru ions
            ml_ion=derivedpointers%npnxt(ml_ion)
            ml_element_test=derivedpointers%npar(ml_ion)
            enddo
!
!
          cltot=cltot+cl*xeltp 
          httot=httot+ht*xeltp 
          cltot2=cltot2+cl2*xeltp 
          httot2=httot2+ht2*xeltp 
          htt(jk)=ht*xeltp 
          cll(jk)=cl*xeltp 
          htt2(jk)=ht2*xeltp 
          cll2(jk)=cl2*xeltp 
!
          enelec=enelec+max(0.d0,1.-xisum)*float(nnz)*xeltp
!
!         end of test if abundance
          endif 
!
!       end of step thru elements
        if  (ml_element.ne.0)                                           &
     &          ml_element=derivedpointers%npnxt(ml_element) 
        enddo 
!
      elcter=-enelec+xee
      if (lpri.ne.0) write (lun11,*)'elcter=',elcter,enelec,xee

      call comp2(lpri,lun11,epi,ncn2,bremsa,t,cmp1,cmp2)             
!
!     nonrelativistic compton                                           
!     call comp(lpri,lun11,epi,ncn2,bremsa,cmp1,cmp2)                 
!     ferland compton                                                   
!      call comp3(lpri,lun11,epi,ncn2,bremsa,cmp1,cmp2)               
!      call freef(lpri,lun11,epi,ncn2,t,xpx,xee,opakc)                  
      call bremem(lpri,lun11,xee,xpx,t,epi,ncn2,brcems,opakc) 
      call heatf(lpri,lun11,                                            &
     &       t,r,delr,xee,xpx,                                          &
     &       epi,ncn2,                                                  &
     &       ncsvn,                                                     &
     &       brcems,cmp1,cmp2,httot,cltot,httot2,cltot2,hmctot,         &
     &             htcomp,clcomp,clbrems)                
!                                                                       
      if (lpri.gt.1) write (lun11,*)'leaving calc_hmc_all' 
!                                                                       
      lprisv=lpri 
!         
!                                                                       
      return 
      end                                           
