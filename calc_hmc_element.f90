      subroutine calc_hmc_element(ml_element,lpri,lun11,                &
     &                 critf,vturbi,t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,&
     &                 zeta,mml,mmu,                                    &
     &                 epi,ncn2,bremsa,bremsint,                        &
     &                 leveltemp,                                       &
     &                 tau0,tauc,                                       &
     &                 np2,ncsvn,nlsvn,                                 &
     &                 rnise,bileve,xileve,cl,ht,xii,rrrti,pirti)

!                                                                       
!     Name: calc_hmc_element.f90  
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
!           rnist:  lte level populations relative to element
!           nmat: index into rates array corresponding to first level
!           also uses variables from globaldata
!           Output:
!           ajise(2,ndb) Matrix of rates (s^-1) 
!             ajisb(1,n)=forward rate,ajisb(2,n)=reverse rate
!           cjise(2,ndb): Matrix of cooling rates (erg s^-1)
!           indbe(2,ndb): Index array for ajisb, cjisb
!           ipmat: total number of levels in element
!           nindbe current length of ajisb
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
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     element abundances                                                
!     state variables                                                   
      real(8) r,t,xpx,delr 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee,critf,zeta
      integer ncn2,lpri,lun11,lfpi,np2 
      integer nsup(nd) 
      integer nlsvn,ncsvn 
      integer nindbe
!     continuum flux                                                    
      real(8) tauc(2,nnml) 
!     ion fractions
      real(8) xii(nl)
!     limits on ion indeces vs element
      integer mml(nl),mmu(nl)
      real(8) xh0,xh1,cfrac,tsq
!                                                                       
      real(8) cl,ht,cl2,ht2
      real(8) xileve(nd),rnise(nd),bileve(nd) 
      real(8) x(nd) 
      character(1) kblnk 
      real(8) rnisi(nd)
!     allocating these variables adds significant overhead
!      real(8), allocatable, dimension(:,:) :: ajisi
!      real(8), allocatable, dimension(:) :: cjisi,cjisi2 
!      real(8), allocatable, dimension(:,:) :: ajise
!      real(8), allocatable, dimension(:) :: cjise,cjise2
!      integer, allocatable, dimension(:,:) :: indbi
!      integer, allocatable, dimension(:,:) :: indbe
       real(8) ajisi(2,ndb)
       real(8) ajise(2,ndb)
       real(8) cjisi(ndb)
       real(8) cjisi2(ndb)
       real(8) cjise(ndb)
       real(8) cjise2(ndb)
      integer indbi(2,ndb)
      integer indbe(2,ndb)
      real(8) pirti(nl),rrrti(nl),xitp(nl+1),pirttmp,rrrttmp
      real(8) xitpul,xitpll
      real(8) tt1,tt2
      integer np1i,np1r,np1k 
      integer lfl,lfu
      integer nlev,jk,nnz,nindbi,nindbio,                               &
     &        ltyp,mm,nsp,kk,                                           &
     &        lrtyp,lcon,nrdt,nidt,nkdt,                                &
     &        ml_ion,ml_element,ml_ion_data_type,                       &
     &        ml_element_test,klion,jkk_ion,ipmatsv,                    &
     &        nit,nit2,nit3,nitmx,nitmx2,lprim,ipmat2,ipmat
!                                                                       
      data kblnk/' '/ 
!                                                                       
!      allocate(ajisi(2,ndb))
!      allocate(ajise(2,ndb))
!      allocate(cjisi(ndb))
!      allocate(cjisi2(ndb))
!      allocate(cjise(ndb))
!      allocate(cjise2(ndb))
!      allocate(indbi(2,ndb))
!      allocate(indbe(2,ndb)) 
!                                                                       
      if (lpri.ge.1)                                                    &
     &  write (lun11,901)t,xee,xpx,delr                              
901   format (1x,'    in calc_hmc_element, inputs:',4(1pe10.3))
!
      if (lpri.gt.1) then
        write (lun11,*)'populations'
        do mm=1,nd
          write (lun11,*)mm,xileve(mm)
          enddo 
        endif
!
!     set up super levels
      if (lpri.gt.1) write (lun11,*)'zeroing:',ipmat
      nsp=1
      do mm=1,nd
        nsup(mm)=0 
        enddo 
!
!
!     lfpi value: photoionization and recombination, no opacities       
      lfpi=2 
!                                                                       
      tsq=sqrt(t) 
!                                                                       
!     print element information
      call drd(ltyp,lrtyp,lcon,                                         &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_element,             &
     &            0,lun11)                                        
      jk=masterdata%idat1(np1i)
      nnz=masterdata%idat1(np1i+nidt-2)
      if (lpri.ge.1) then
          write (lun11,902)jk,ml_element,nnz,                           &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)
902           format (1x,'  element:',3(i12,1x),8(1a1))
          endif 
!
!     first pass for total rates
!     step thru ions
      ml_ion_data_type=12
      ml_ion=derivedpointers%npfirst(ml_ion_data_type)
      ipmat=0
      nindbi=0
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
!
!         get level data                                          
          nlev=derivedpointers%nlevs(jkk_ion)
          if (lpri.gt.1) write (lun11,*)'nlev=',nlev
          if (lpri.gt.1) write (lun11,*)'ipmat=',ipmat
          do mm=1,nlev
!           get level pointer                                     
!            rnisi(mm)=rnise(mm+ipmat)
            rnisi(mm)=0.
            if (lpri.gt.1)                                              &
     &       write (lun11,*)'before calc_ion_rates',                    &
     &         mm,rnise(mm+ipmat)
            enddo
          ipmat=ipmat+nlev-1
!
          nindbio=nindbi
          call calc_ion_rates(ml_ion,lpri,lun11,                        &
     &                   vturbi,t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,    &
     &                   epi,ncn2,bremsa,bremsint,                      &
     &                   leveltemp,                                     &
     &                   tau0,tauc,                                     &
     &                   np2,ncsvn,nlsvn,                               &
     &                   pirttmp,rrrttmp,                               &
     &                   rnisi,nlev)
          pirti(klion)=pirttmp
          rrrti(klion)=rrrttmp
          if (lpri.gt.1) write (lun11,*)klion,pirttmp,rrrttmp
!
!         end of test if element belongs to parent of ion
          endif
!
!       end of loop over ions
        ml_ion=derivedpointers%npnxt(ml_ion)
        enddo 
!
!     calculate ion fractions
      call istruc(pirti,rrrti,xitp,nnz,lpri,lun11) 
!
!     find limits on ions
      lfu=0
      lfl=0
      mmu(jk)=nnz+2
      mml(jk)=0
      mm=0
      xitpul=0.
      xitpll=0.
      do while ((mm.le.nnz).and.((lfu.eq.0).or.(lfl.eq.0)))
        mm=mm+1
        mmu(jk)=mmu(jk)-1
        mml(jk)=mml(jk)+1
        if ((xitpul.lt.critf).and.(xitp(mmu(jk)).ge.critf)              &
     &     .and.(lfu.eq.0)) then
          lfu=mmu(jk)
          endif
        if ((xitpll.lt.critf).and.(xitp(mml(jk)).ge.critf)              &
     &     .and.(lfl.eq.0)) then
          lfl=mml(jk)
          endif
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)'searching for ion limits:',xitp(mml(jk)),     &
     &    xitp(mmu(jk)),xitpll,xitpul,mml(jk),mmu(jk),mm,lfl,lfu
        xitpll=xitp(mml(jk))
        xitpul=xitp(mmu(jk))
        enddo
      mmu(jk)=lfu
      mml(jk)=lfl
      if (lfu.eq.0) mmu(jk)=nnz+1
      if (lfl.eq.0) mml(jk)=1
      mml(jk)=max(1,mml(jk)-1)
      mmu(jk)=min(nnz,mmu(jk)+1)
      if (lpri.ge.1) write (lun11,*)'ion limits:',mml(jk),mmu(jk),mm
!
      call levwkelement(ml_element,lpri,ipmatsv,t,xee,xpx,leveltemp,    &
     &    lun11,rnise,mml(jk),mmu(jk))    
!
!     second pass for ion by ion rates
!     step thru ions
      ipmat=0
      ipmat2=0
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
!
          nlev=derivedpointers%nlevs(jkk_ion)
!
!         at the ion level we work only with ions between the limits
!         ipmat points to the ions in the element arrays, and so includes
!         ions outside the limits
!         ipmat2 is the pointer to the ion level population arrays...
          if ((klion.ge.mml(jk)).and.(klion.le.mmu(jk))) then
!
!           get level data                                          
            if (lpri.gt.1) write (lun11,*)'nlev=',nlev
            if (lpri.gt.1) write (lun11,*)'ipmat=',ipmat
            if (lpri.gt.1) write (lun11,*)'ipmat2=',ipmat2
            do mm=1,nlev
!             get level pointer                                     
              rnisi(mm)=rnise(mm+ipmat)
              if (lpri.gt.1)                                            &
     &         write (lun11,*)'before calc_hmc_ion',                    &
     &         mm,rnise(mm+ipmat),xileve(mm+ipmat)
              enddo
!
           nindbio=nindbi
           call calc_hmc_ion(ml_ion,lpri,lun11,                         &
     &                   vturbi,t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,    &
     &                   epi,ncn2,bremsa,bremsint,                      &
     &                   leveltemp,                                     &
     &                   tau0,tauc,                                     &
     &                   np2,ncsvn,nlsvn,                               &
     &                   pirttmp,rrrttmp,                               &
     &                   rnisi,                                         &
     &                   ajisi,cjisi,cjisi2,indbi,nindbi,nlev)
!
            if (lpri.gt.0) write (lun11,*)'   pirttmp,rrrtmp:',         &
     &           pirttmp,rrrttmp
!
!           map to element matrix
            do mm=nindbio+1,nindbi
              do kk=1,2
                ajise(kk,mm)=ajisi(kk,mm)
                indbe(kk,mm)=indbi(kk,mm)+max(0,ipmat2)
                if (lpri.gt.1) write (lun11,*)mm,kk,indbe(kk,mm),       &
     &                       indbi(kk,mm),ajisi(kk,mm),ipmat2
                enddo
              cjise(mm)=cjisi(mm)
              cjise2(mm)=cjisi2(mm)
              enddo
!
!           set up superlevel pointers                   
            if (lpri.ge.1) write (lun11,*)'superlevel pointers',        &
     &           nlev,nsp,ipmat2
            nsup(1+ipmat2)=nsp 
            if (nlev.gt.2) then
              nsp=nsp+1 
              do mm=2,nlev-1
                nsup(mm+ipmat2)=nsp 
                if (lpri.ge.1) write (lun11,*)mm,mm+ipmat2,nsp
                enddo 
              endif 
            nsp=nsp+1 
!
            do mm=1,nlev
              x(mm+ipmat2)=xileve(mm+ipmat)
              if (lpri.gt.1)                                            &
     &         write (lun11,*)'mapping 2:',mm,mm+ipmat2,mm+ipmat,       &
     &                 x(mm+ipmat2)
              enddo
!
            ipmat2=ipmat2+nlev-1
!
!           end of test for ion in range
            endif
!                                                                       
          ipmat=ipmat+nlev-1
          if (lpri.ge.1) write (lun11,*)'ipmat=',ipmat,nsp,nlev
!
!         end of test if element belongs to parent of ion
          endif
!
!       end of loop over ions
        ml_ion=derivedpointers%npnxt(ml_ion)
        enddo 
!
      nindbe=nindbi
!
!     more superlevels
      nsup(ipmat2+1)=nsp 
      x(ipmat2+1)=0. 
!     nb this statement may be needed
      ipmat2=ipmat2+1 

!     now calculate populations
      call remtms(tt1) 
      nitmx=40 
      nitmx2=40 
      lprim=0
!       if (lpri.gt.0) lprim=3                                    
!      if (lpri.gt.0) lprim=1                                    
      if (lprim.ne.0) then
        write (lun11,*)'before msolvelucy'
        do mm=1,ipmat
          write (lun11,*)mm,x(mm)
          enddo 
        endif
      call msolvelucy(ajise,cjise,cjise2,                               &
     &          indbe,nindbe,nsup,nsp,ipmat2,                           &
     &          x,ht,cl,ht2,cl2,nit,nit2,nit3,nitmx,nitmx2,lun11,lprim)    
!      if (lpri.ge.1) then
!        write (lun11,*)'after msolvelucy'
!        do mm=1,ipmat
!          write (lun11,*)mm,x(mm)
!          enddo 
!        endif
      if (lpri.ge.1)                                                    &
     &       call chisq(ajise,indbe,nindbe,                             &
     &                    ipmat2,x,lun11,lpri)                           
      call remtms(tt2) 
      if (lpri.gt.0)                                                    &
     &  write (lun11,981)zeta,t,ipmat2,jk,abs(tt2-tt1),nit,nit2,nit3,   &
     &       ht,cl         
  981   format (1x,'after msolvelucy',2(1pe11.3),2i6,(1pe11.3),         &
     &           3i8,2(1pe11.3)) 
      if (lpri.eq.-1)                                                   &
     &  write (13,981)zeta,t,ipmat2,jk,abs(tt2-tt1),nit,nit2,nit3,      &
     &        ht,cl         
!
!     now map back to element abundances
!     step thru ions
      ml_ion=derivedpointers%npfirst(ml_ion_data_type)
      ipmat=0
      ipmat2=0
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
          xii(klion)=0.
!
!         get level data                                          
          nlev=derivedpointers%nlevs(jkk_ion)
          if (lpri.gt.1) write (lun11,*)'nlev=',nlev
          if (lpri.gt.1) write (lun11,*)'ipmat=',ipmat
          if (lpri.gt.1) write (lun11,*)'klion=',klion,mml(jk),mmu(jk)
!
!         use calculated abundances if in range
          if ((klion.ge.mml(jk)).and.(klion.le.mmu(jk))) then
!
!           map to element matrix
            do mm=1,nlev
              xileve(mm+ipmat)=x(mm+ipmat2) 
              bileve(mm+ipmat)=x(mm+ipmat2)/(1.e-37+rnise(mm+ipmat))
              if (mm.lt.nlev)                                           &
     &           xii(klion)=xii(klion)+x(mm+ipmat2)
              if (lpri.gt.1)                                            &
     &         write (lun11,*)'saving populations',mm,mm+ipmat,         &
     &          mm+ipmat2,xileve(mm+ipmat),rnise(mm+ipmat),             &
     &          bileve(mm+ipmat),xii(klion)
              enddo
            if (lpri.gt.1)                                              &
     &       write (lun11,*)'ipmat,ipmat2,nlev:',ipmat,ipmat2,nlev
!
            ipmat2=ipmat2+nlev-1
!
            else
!
            do mm=1,nlev
              xileve(mm+ipmat)=0.
              bileve(mm+ipmat)=0.
              if (lpri.gt.1)                                            &
     &         write (lun11,*)'not saving populations',mm,mm+ipmat,     &
     &           xileve(mm),rnise(mm),bileve(mm)

              enddo
!
!           end of test for ion in range
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
      if (lpri.gt.0) write (lun11,*)'leaving calc_hmc_element'
!
!      deallocate(ajisi)
!      deallocate(ajise)
!      deallocate(cjisi)
!      deallocate(cjisi2)
!      deallocate(cjise)
!      deallocate(cjise2)
!      deallocate(indbi)
!      deallocate(indbe) 
!                                                                       
      return 
      end                                           
