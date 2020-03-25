      subroutine calc_rates_level_lte(jkk,lpri,lun11,t,xee,xpx,         &
     &              leveltemp,nlev)
!                                                                       
!     Name: calc_rates_level_lte.f90  
!     Description:  
!           Calculates quantities associated with level data
!           for one ion.  Also unpacks level data into temporary 
!           arrays in the module globaldata
!           formerly:  func2l
!
!     List of Parameters:
!           Input:
!           jkk: index of ion in xstar scheme 1=H0, 432=Zn29+
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           vturbi: turbulent speed in km/s
!           t: temperature in 10^4K
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           Output:
!           rniss:  lte level populations
!           rnisse:  lte level populations with exponential removed
!           nlev:  number of levels
!     Dependencies: levwk, drd
!     Called by: calc_hmc_all, fstepr2,fstepr3,calc_emis_all,pprint
!
!     this routine calculates rates affecting level populations         
!     author: T. Kallman                                                
!     note that in this routine rniss indeces are relative to ground
!                                                                       
                                                                        
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,nd) 
        integer:: ilev(10,nd),nlpt(nd),iltp(nd) 
        character(1) :: klev(100,nd) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
      integer np1i,np1r,np1k 
      character(1) kblnk 
!                                                                       
      integer nlevmx,mltype,ml,mllz,nlev,nidt,nrdt,nkdt,lun11,          &
     &        lpri,ltyp,lrtyp,lcon,jkk,mm,lk,mlpar,mlm,lpril                  
      real(8) xee,xpx,t
!                                                                       
      data kblnk/' '/ 
!                                                                       
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in calc_rates_level_lte, inputs:',t,            &
     &          xee,xpx                                                 
!                                                                       
!     now find level data                                               
!     step thru types                                                   
      nlevmx=0 
      mltype=13 
      ml=derivedpointers%npfi(mltype,jkk) 
      if (lpri.gt.1) write (lun11,*)'jkk=',                             &
     &      jkk,ml,derivedpointers%npar(ml) 
      mllz=derivedpointers%npar(ml) 
!     step thru records of this type                                    
      mlpar=derivedpointers%npar(ml) 
      do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
         mlm=ml 
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         nlev=masterdata%idat1(np1i+nidt-2) 
         nlevmx=max(nlevmx,nlev) 
         if ((nlev.gt.0).and.(nlev.le.nd)) then 
           leveltemp%nlpt(nlev)=ml 
           leveltemp%iltp(nlev)=ltyp 
 9101      format (1x,'level quantities:',4i6,4(1pe12.5),3i6,8a1) 
           if (lpri.gt.1) write (lun11,9101)                            &
     &       ml,nlev,ltyp,lrtyp,(masterdata%rdat1(np1r+mm-1),mm=1,4),   &
     &       masterdata%idat1(np1i),masterdata%idat1(np1i+1),           &
     &       masterdata%idat1(np1i+2),                                  &
     &       (masterdata%kdat1(np1k+mm-1),mm=1,8)                    
           do  lk=1,nrdt 
             leveltemp%rlev(lk,nlev)=masterdata%rdat1(np1r+lk-1) 
             enddo 
           do lk=1,nidt 
             leveltemp%ilev(lk,nlev)=masterdata%idat1(np1i+lk-1) 
             enddo 
           do lk=1,nkdt 
             leveltemp%klev(lk,nlev)=masterdata%kdat1(np1k+lk-1) 
             enddo 
           do lk=nkdt+1,100 
             leveltemp%klev(lk,nlev)=kblnk 
             enddo 
           endif 
         ml=derivedpointers%npnxt(ml) 
         if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
         enddo 
      nlev=nlevmx 
      lpril=lpri
!                                                                       
      return 
      END                                           
