      subroutine numlevel(jkk,                                       &
     &       nlev) 
!                                                                       
!     Name: numlevel.f90  
!     Description:  
!            counts the levels for each ion                       
!     List of Parameters:
!           Input:
!           jkk: index of ion in xstar scheme 1=H0, 432=Zn29+
!           also uses variables from globaldata
!           Output:
!           nlev:  number of levels
!     Dependencies:  Calls drd
!     Called by: calc_hmc_all, calc_emis_all
!
!     author: T. Kallman                                                
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      integer nlevmx,mltype,ml,mllz,nlev,                               &
     &        jkk,ltyp,lrtyp,lcon,mlpar,lun11,mlm                       
      integer  nrdt,np1r,nidt,np1i,nkdt,np1k
!                                                                       
!          now find level data                                          
!          step thru types                                              
           nlevmx=0 
           mltype=13 
           ml=derivedpointers%npfi(mltype,jkk) 
           mllz=derivedpointers%npar(ml) 
!          step thru records of this type                               
           mlpar=derivedpointers%npar(ml) 
           do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
              mlm=ml 
              call drd(ltyp,lrtyp,lcon,                                 &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                      &
     &          0,lun11)                                          
              nlevmx=nlevmx+1 
              nlev=masterdata%idat1(np1i+nidt-2) 
              ml=derivedpointers%npnxt(ml) 
              if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
              enddo 
           nlev=nlevmx 
!                                                                       
      return 
      END                                           
