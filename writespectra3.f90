      subroutine writespectra3(lun11,lpri,nparms,parname,partype,    &
     &       parval,parcomm,atcredate,epi,ncn2,dpthc,dpthcont,          &
     &       np2,                                                       &
     &       elum,zrems,zremsz,kmodelname,nloopctl)             
                                                                        
!                                                                       
!     Name: writespectra3.f90  
!     Description:  
!           Writes out continuum into the file xout_cont1.fits
!           Not including binned lines
!
!     List of Parameters:
!     Input:
!           lun11: logical unit number for printing
!           lpri: print switch, 1=on, 0=off
!           nparms: number of input parameters
!           parname(nparms): names of input parameters
!           partype(nparms): types of input parameters
!           parval(nparms): values of input parameters
!           parcomm(nparms): comments of input parameters
!           atcredate:  atomic data file creation date (string length 63)
!           t:  temperature (10^4 K)
!           vturbi:  ion turbulent speed (km/s)
!           epi(ncn):  continuum energy bins (eV)
!           ncn2:  number of continuum energy bins
!           dpthc(2,ncn): optical depth in continuum bins 
!           np2: atomic data parameter, number of records in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           zrems(4,ncn):  master spectrum array.  (erg/s/erg/10^38)
!           zremsz(ncn):  input spectrum  (erg/s/erg/10^38)
!           kmodelname:  model name 
!           nloopcntl:  loop control variable
!
!     Dependencies: none
!     Called by:  xstar
!
!     Write extension containing the spectra for                        
!     this particular model.                                            
!                                                                       
!     Modifications:                                                    
!       04/01/1999,WTB: Disabled appending loop control value to        
!               extension name due to changes in xstar2xspec design     
!                                                                       
!     author:  T. Bridgman                                              
!                                                                       
      use globaldata
      implicit none 
      integer ncn2 
!                                                                       
!     passed parameters                                                 
      character(30) kmodelname 
      integer nparms, nloopctl, lun11 
      character(20) parname(55) 
      character(10) partype(55) 
      real(8) parval(55) 
      character(30) parcomm(55) 
!     line luminosities                                                 
      real(8) elum(3,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     the atomic data creation date                                     
      character(63) atcredate 
!     continuum lum                                                     
      real(8) zrems(5,ncn),zremsz(ncn) 
!     continuum optical depths                                          
      real(8) dpthc(2,ncn),dpthcont(2,ncn)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: zrtmp
      REAL(4), DIMENSION(:), ALLOCATABLE :: rtmp
      character(16) knam,klabs(5),kunits(5),kform(5),kblnk16 
      character(30) extname 
      integer unit,istatus, kl,lpri 
      integer ll, numcon
      integer np2, tbcol(5), nrows, rowlen, kk
      integer frow, felem, colnum, tfields, status, verbose,mm 
!                                                                       
!     Not used                                                          
      real(8) javir 
      integer javi 
!      character(80) javik                                              
!                                                                       
      data kblnk16/'                '/ 
                                                                        
      ALLOCATE(zrtmp(5,ncn))
      ALLOCATE(rtmp(ncn))
!
      javir=epi(1) 
!      epi(1)=javir                                                     
      javir=dpthc(1,1) 
      javi=ncn2 
      javi=np2 
!      np2=javi                                                         
      javi=derivedpointers%npfirst(1) 
      javi=derivedpointers%nplini(1) 
      javi=derivedpointers%npcon(1) 
      javi=derivedpointers%npconi(1) 
      javi=derivedpointers%npilev(1,1) 
      javi=derivedpointers%npilevi(1) 
      javi=derivedpointers%npconi2(1) 
                                                                        
      javi=masterdata%idat1(1) 
      javir=masterdata%rdat1(1) 
!      javik=kdat1(1)                                                   
      javi=masterdata%nptrs(1,1) 
      javi=derivedpointers%npar(1) 
      javi=derivedpointers%npnxt(1) 
      javi=derivedpointers%npfi(1,1) 
      javi=derivedpointers%nplin(1) 
      javir=elum(1,1) 
                                                                        
!                                                                       
                                                                        
      verbose=lpri 
!                                                                       
!     open and prepare the fits file for spectral data                  
      if(verbose.gt.0) write (lun11,*)'writespectra3: opening header'&
     &  ,kmodelname                                                      
      knam='xout_cont1.fits' 
      call fheader(unit,knam,atcredate,kmodelname,istatus) 
      if(istatus.gt.0) call printerror(lun11,istatus) 
                                                                        
!     write extension of parameter values                               
      if(verbose.gt.0)                                                  &
     &  write (lun11,*)'writespectra: write parameter list'             
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,&
     &               parcomm,nloopctl,istatus,lun11)                    
      if(istatus.gt.0) call printerror(lun11,istatus) 
      if(verbose.gt.0)                                                  &
     &  write (lun11,*)'writespectra: building data tables'             
                                                                        
!     build spectra data tables                                         
      numcon=ncn2 
      do ll=1,ncn2 
        zrtmp(4,ll)=0. 
        zrtmp(5,ll)=0. 
        enddo 
      do kl=1,numcon 
         zrtmp(4,kl)=zrtmp(4,kl)+zrems(4,kl) 
         zrtmp(5,kl)=zrtmp(5,kl)+zrems(5,kl) 
         zrtmp(3,kl)=zremsz(kl)*exp(-dpthcont(1,kl)) 
!         write (lun11,968)kl,epi(kl),zremsz(kl),                       
!     $          zrtmp1(kl),zrtmp2(kl)                                  
         zrtmp(2,kl)=zremsz(kl) 
         zrtmp(1,kl)=epi(kl) 
         enddo 
                                                                        
!     write the spectral data to the extension                          
      do mm=1,5 
        kunits(mm)=kblnk16 
        klabs(mm)=kblnk16 
        kform(mm)=kblnk16 
        enddo 
      klabs(1)='energy          ' 
      kform(1)='E11.3' 
      kunits(1)='eV' 
      klabs(2)='incident        ' 
      kform(2)='E11.3' 
      kunits(2)='erg/s/erg' 
      klabs(3)='transmitted     ' 
      kform(3)='E11.3' 
      kunits(3)='erg/s/erg' 
      klabs(4)='emit_inward     ' 
      kform(4)='E11.3' 
      kunits(4)='erg/s/erg' 
      klabs(5)='emit_outward    ' 
      kform(5)='E11.3' 
      kunits(5)='erg/s/erg' 
!     build extension name                                              
      extname='XSTAR_SPECTRA' 
!      if(nloopctl.gt.0) then                                           
!          write(ktmp2,'(i4.4)')nloopctl                                
!          extname='xstar_spectra_' // ktmp2                            
!          endif                                                        
      if(verbose.gt.0)                                                  &
     &  write (lun11,*)'writespectra: writing spectral data'            
!      call writespectra(unit,ktmp1,zrtmp,5,999,ncn,                    
!     $                 klabs,kform,kunits)                             
                                                                        
!     append a new empty extension onto the end of the primary array    
      status=0 
      call ftcrhd(unit,status) 
      if(verbose.gt.0)                                                  &
     &   write (lun11,*)'writespectra: writing header table'            
                                                                        
      tfields=5 
      nrows=ncn2 
      rowlen=0 
      do mm=1,5 
      tbcol(mm)=0 
      enddo 
                                                                        
!     write the required header parameters for the ascii table          
      status=0 
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,   &
     &            extname,status)                                       
      if (status .gt. 0)call printerror(lun11,status) 
      status=0 
!                                                                       
!     map each column to a 1-d array before writing to the file         
      do kk=1,tfields 
        if(verbose.gt.0)                                                &
     &    write (lun11,*)'writespectra: building column ',kk            
        frow=1 
        felem=1 
        colnum=kk 
        do ll=1,nrows 
          rtmp(ll)=sngl(zrtmp(kk,ll))
          enddo 
        if(verbose.gt.0)                                                &
     &    write (lun11,*)'writespectra: writing column ',kk             
        status=0 
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status) 
        if (status .gt. 0)call printerror(lun11,status) 
        enddo 
                                                                        
!     compute checksums                                                 
      if(verbose.gt.0) write (lun11,*)'writespectra:writingchecksum' 
      call ftpcks(unit,status) 
!     check for any error, and if so print out error messages           
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if(verbose.gt.0) write (lun11,*)'writespectra: closing file' 
      call fitsclose(lun11,unit,istatus) 
!
      DEALLOCATE(zrtmp)
      DEALLOCATE(rtmp)
!                                                                       
      return 
      end                                           
