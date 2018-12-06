      subroutine func1(jkk,kl,nnz,lpri,lun11,vturbi,                    &
     &       t,trad,r,delr,xee,xpx,xh1,xh0,                             &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       np2,ncsvn,nlsvn,                                           &
     &       rniss,rnisse,nlev,pirt,rrrt2)                                                
!     Name: func1.f90  
!     Description:  
!           Calculates rates affecting ionization balance for 
!           one ion.  These rates are photoionization rate and 
!           recombination rate
!
!     List of Parameters:
!           Input:
!           jkk: index of ion in xstar scheme 1=H0, 432=Zn29+
!           kl:  index of ion relative element: 1=neutral, n=hydrogenic
!           nnz: atomic number of element
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
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           also uses variables from globaldata
!           Output:
!           pirt(nni):  photoionization rate coefficients (s^-1).  
!                  Only the element  corresponding to the current ion is filled.
!           rrrt(nni):  recombination rate coefficients (s^-1).  
!                  Only the element  corresponding to the current ion is filled.
!           
!        Dependencies:  Calls ucalc,drd
!                  
!                                                                       
!     this routine calculates rates affecting ion balance               
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
      real(8) rcdum(2,ncn) 
      real(8) fline(2,nnnl) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     element abundances                                                
!     state variables                                                   
      real(8) r,t,xpx,delr 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee 
      integer ncn2,lpri,lun11,np2 
      integer nlsvn,ncsvn 
      real(8) rniss(nd),rnisse(nd)
      real(8) pirt(31) 
      character(49) kdesc2 
      real(8) tsq,ans1,ans2,xh1,xh0,rrrt2 
      real(8) abund1,abund2,ptmp1,ptmp2,ans3,ans4,opakb1 
      integer idest1,idest2,idest3,idest4 
      integer np1i,np1r,np1k 
      integer nnzp,mltype,ml,mllz,nlev,mlpar,                           &
     &  ltyp,lrtyp,lcon,nrdt,nidt,nkdt,mlrdesc,llo,lup,                 &
     &  nnz,mm,jkk,lpriu,kl,mlm,lfpi                                 
!                                                                       
      if (lpri.ne.0)                                                    &
     &  write (lun11,*)'in func1, inputs:',t,                           &
     &         xee,xpx,bremsa(1),bremsa(10)                             
!                                                                       
!     lfpi value:  calculate photoionization rates only                 
      lfpi=1 
!                                                                       
!     zero temporaries                                                  
      tsq=sqrt(t) 
      nnzp=nnz+1 
      do mm=1,nnzp 
        pirt(mm)=0. 
        enddo 
      rrrt2=0. 
!     now find all the rates affecting this ion                         
!     step thru types                                                   
      mltype=1 
      do while (mltype.lt.ntyp) 
        mlrdesc=mltype 
        ml=derivedpointers%npfi(mltype,jkk) 
        if (((mlrdesc.eq.6).or.(mlrdesc.eq.1).or.(mlrdesc.eq.8)         &
     &    .or.(mlrdesc.eq.15).or.(mlrdesc.eq.7).or.(mlrdesc.eq.42))     &
     &    .and.(ml.ne.0)) then                                          
          mllz=derivedpointers%npar(ml) 
          mlpar=derivedpointers%npar(ml) 
          do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
!           step thru records of this type                              
            if (masterdata%nptrs(3,ml).eq.mlrdesc) then 
              mlm=ml-1 
              call drd(ltyp,lrtyp,lcon,                                 &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                      &
     &          0,lun11)                                          
!              if (lpri.ne.0)  call dprinto(ltyp,lrtyp,lcon,            
!     $          nrdt,np1r,nidt,np1i,nkdt,np1k,                         
!     $          masterdata%rdat1,idat1,kdat1,lun11)                           
!             calculate rates                                           
              lpriu=lpri 
              abund1=0. 
              abund2=0. 
              ptmp1=0. 
              ptmp2=0. 
              idest1=0 
              if (lrtyp.eq.7) idest1=masterdata%idat1(np1i+nidt-2) 
              if (((lrtyp.eq.1).and.(ltyp.ne.53))                       &
     &             .or.((lrtyp.eq.7).and.(idest1.eq.1))                 &
     &             .or.(lrtyp.eq.15).or.(lrtyp.eq.42)) then             
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,              &
     &             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,             &
     &             ans3,ans4,idest1,idest2,idest3,idest4,               &
     &             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,                &
     &             opakc,opakcont,rcdum,fline,lpriu,kdesc2,             &
     &             r,delr,t,trad,tsq,xee,xh1,xh0,                       &
     &             epi,ncn2,bremsa,bremsint,                            &
     &             rniss,rnisse,nlev,lfpi,lun11,                        &
     &             np2,ncsvn,nlsvn)                
                 if (idest1.eq.1) then 
                   llo=idest3 
                   lup=idest4-idest3 
                   pirt(kl+lup)=pirt(kl+lup)+ans1 
                   if (lpri.ge.1)                                       &
     &              write (lun11,9001)jkk,lrtyp,ltyp,llo,lup+idest3,    &
     &                      ml,ans1,pirt(kl+lup)                        
 9001               format (1x,6i6,' ion ',1pe10.3,14x,1pe10.3) 
                   endif 
                 endif 
              if ((lrtyp.eq.6).or.(lrtyp.eq.8)) then 
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,              &
     &             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,             &
     &             ans3,ans4,idest1,idest2,idest3,idest4,               &
     &             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,                &
     &             opakc,opakcont,rcdum,fline,lpriu,kdesc2,             &
     &             r,delr,t,trad,tsq,xee,xh1,xh0,                       &
     &             epi,ncn2,bremsa,bremsint,                            &
     &             rniss,rnisse,nlev,lfpi,lun11,                        &
     &             np2,ncsvn,nlsvn)                
                 rrrt2=rrrt2+ans1 
!                Commented... llo has a funny value...    jg            
                 llo=1 
                 if (lpri.ge.1)                                         &
     &              write (lun11,9002)jkk,lrtyp,ltyp,llo,ml,            &
     &                    ans1,rrrt2                                    
 9002            format (1x,4i6,6x,i6,' ion ',14x,                      &
     &                              1pe10.3,14x,1pe10.4)                
                 endif 
              endif 
            ml=derivedpointers%npnxt(ml) 
            if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
            enddo 
          endif 
        mltype=mltype+1 
        enddo 
!                                                                       
!                                                                       
      return 
      end                                           
