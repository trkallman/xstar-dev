      subroutine writespectra(lun11,lpri,lwri,nparms,                   &
     &       parname,partype,parval,parcomm,atcredate,                  &
     &       t,vturbi,epi,ncn2,dpthc,                                   &
     &       nlsvn,                                                     &
     &       elinb,elum,zrems,zremsz,kmodelname,nloopctl)             
!                                                                       
!     Name: writespectra.f90  
!     Description:  
!           Writes out spectrum binned in continuum bins 
!              including emission lines into the file xout_spect1.fits
!
!     List of Parameters:
!     Input:
!           lun11: logical unit number for printing
!           lpri: print switch, 1=on, 0=off
!           atcredate:  atomic data file creation date (string length 63)
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
!           nlsvn: atomic data parameter, number of lines in atomic database
!           elum(2,nnnl):  line luminosities (erg/s/10^38)
!           zrems(4,ncn):  master spectrum array.  (erg/s/erg/10^38)
!           zremsz(ncn):  input spectrum  (erg/s/erg/10^38)
!           kmodelname:  model name 
!           nloopcntl:  loop control variable
!
!     Dependencies: binemis
!     Called by:  xstar
!     
!     Write extension containing the spectra for                        
!     this particular model.                                            
!                                                                       
!     Modifications:                                                    
!       04/01/1999,WTB: Disabled appending loop control value to        
!               extension name due to changes in xstar2xspec design     
!       051/17/2003 TK added auger damping                              
!     author:  T. Bridgman                                              
!       3/2/2017:  changed to reduce memory.  zrtmp is now a temporary used 
!               to store zrems during binemis call.
!                                                                       
      use globaldata
!                                                                       
      implicit none 
!                                                                       
!     passed parameters                                                 
      character(30) kmodelname 
      integer nparms, nloopctl, lun11 
      character(20) parname(55) 
      character(10) partype(55) 
      real(8) parval(55) 
      character(30) parcomm(55) 
!     line luminosities                                                 
      real(8) elum(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum lum                                                     
      real(8) zrems(5,ncn),zremsz(ncn) 
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
      real(8) elinb(nnnl)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: zrtmp
      REAL(4), DIMENSION(:), ALLOCATABLE :: rtmp
      character(16) knam,klabs(6),kunits(6),kform(6),kblnk16 
      character(30) extname 
      integer unit,istatus 
      integer nlsvn, ll 
      integer tbcol(6), nrows, rowlen, kk 
      integer frow, felem, colnum, tfields, status, verbose,mm
      real(8) eliml, elimh 
      real(8) vturbi 
!     the atomic data creation date                                     
      character(63) atcredate 
!                                                                       
! jg                                                                    
      real(8) t,xlum 
                                                                        
      integer lpri, ncn2,lwri 
!                                                                       
!     Not used                                                          
      integer javi 
!
      data kblnk16/'                '/ 
!                                                                       
       save kblnk16

      ALLOCATE(rtmp(ncn))
      ALLOCATE(zrtmp(5,ncn))
!                                                                       
      javi=derivedpointers%npfirst(1) 
      javi=derivedpointers%nplini(1) 
      javi=derivedpointers%npcon(1) 
      javi=derivedpointers%npconi(1) 
      javi=derivedpointers%npilev(1,1) 
      javi=derivedpointers%npilevi(1) 
      javi=derivedpointers%npconi2(1) 
                                                                        
!                                                                       
!                                                                       
      verbose=lpri 
      eliml=0.1 
      elimh=1.d+5 
      elimh=min(elimh,8.9d+4) 
!                                                                       
!     open and prepare the fits file for spectral data                  
      if(verbose.gt.0) write (lun11,*)'writespectra: opening header',&
     &  kmodelname                                                      
      knam='xout_spect1.fits' 
      call fheader(unit,knam,atcredate,kmodelname,istatus) 
      if(istatus.gt.0) call printerror(lun11,istatus) 
!                                                                       
!                                                                       
!     write extension of parameter values                               
      if(verbose.gt.0)                                                  &
     &     write (lun11,*)'writespectra: write parameter list'          
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,&
     &               parcomm,nloopctl,istatus,lun11)                    
      if(istatus.gt.0) call printerror(lun11,istatus) 
      if(verbose.gt.0)                                                  &
     &  write (lun11,*)'writespectra: building data tables'             
                                                                        
      xlum=parval(10) 
      do ll=1,ncn2
        do mm=1,5
          zrtmp(mm,ll)=zrems(mm,ll)
          enddo
        enddo
      if (lwri.ge.0)                                                    &
     & call binemis(lun11,lpri,xlum,                                    &
     &       t,vturbi,epi,ncn2,dpthc,                                   &
     &       nlsvn,                                                     &      
     &       elinb,eliml,elimh,elum,zrems,zremsz)
!                                                                       
!     write the spectral data to the extension                          
      do mm=1,6 
        kunits(mm)=kblnk16 
        klabs(mm)=kblnk16 
        kform(mm)=kblnk16 
        enddo 
      klabs(1)='energy          ' 
      kform(1)='E13.5' 
      kunits(1)='eV' 
      klabs(2)='incident        ' 
      kform(2)='E13.5' 
      kunits(2)='erg/s/erg' 
      klabs(3)='transmitted     ' 
      kform(3)='E13.5' 
      kunits(3)='erg/s/erg' 
      klabs(4)='emit_inward     ' 
      kform(4)='E13.5' 
      kunits(4)='erg/s/erg' 
      klabs(5)='emit_outward    ' 
      kform(5)='E13.5' 
      kunits(5)='erg/s/erg' 
      klabs(6)='scattered       ' 
      kform(6)='E13.5' 
      kunits(6)='erg/s/erg' 
!     build extension name                                              
      extname='XSTAR_SPECTRA' 
!      if(nloopctl.gt.0) then                                           
!          write(ktmp2,'(i4.4)')nloopctl                                
!          extname='xstar_spectra_' // ktmp2                            
!          endif                                                        
      if(verbose.gt.0)                                                  &
     &   write (lun11,*)'writespectra: writing spectral data'           
                                                                        
!     append a new empty extension onto the end of the primary array    
      status=0 
      call ftcrhd(unit,status) 
      if(verbose.gt.0)                                                  &
     &    write (lun11,*)'writespectra: writing header table'           
                                                                        
      tfields=5
      nrows=ncn2 
      rowlen=0 
      do mm=1,6 
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
     &      write (lun11,*)'writespectra: building column ',kk          
        frow=1 
        felem=1 
        colnum=kk 
        do ll=1,nrows 
          if (kk.gt.1) then
          rtmp(ll)=sngl(zrems(kk-1,ll))
          else
          rtmp(ll)=sngl(epi(ll))
          endif
          enddo 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra: writing column ',kk            
        status=0 
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status) 
        if (status .gt. 0)call printerror(lun11,status) 
        enddo 
                                                                        
!     compute checksums                                                 
      if(verbose.gt.0) write (lun11,*)'writespectra: writingchecksum' 
      status=0 
      call ftpcks(unit,status) 
!     check for any error, and if so print out error messages           
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if(verbose.gt.0) write (lun11,*)'writespectra: closing file' 
      call fitsclose(lun11,unit,istatus) 
!                                                                       
      do ll=1,ncn2
        do mm=1,5
          zrems(mm,ll)=zrtmp(mm,ll)
          enddo
        enddo
!
      DEALLOCATE(rtmp)
      DEALLOCATE(zrtmp)
!                                                                       
      return 
      end                                           
