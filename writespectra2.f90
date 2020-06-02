      subroutine writespectra2(lun11,lpri,nparms,parname,partype,    &
     &       parval,parcomm,atcredate,epi,ncn2,dpthc,                   &
     &       np2,nlsvn,                                                 &
     &       elum,tau0,kmodelname,nloopctl)                     
                                                                        
!                                                                       
!     Name: writespectra2.f90  
!     Description:  
!           Writes out lines into the file xout_lines1.fits
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
!           elum(2,nnnl):  line luminosities (erg/s/10^38)
!           tau0(2,nnnl):  line optical depths
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
!                                                                       
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
      real(8) elum(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     the atomic data creation date                                     
      character(63) atcredate 
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
      real(4) rtmp
      character(16) knam,klabs(9),kunits(9),kform(9),kblnk16 
      character(30) extname 
      integer unit,istatus, nilin, nkdt,nidt,lcon,lrtyp,ltyp,ml,status 
      integer nlsvn, ln, lnn, nrdt,mllz
      integer tbcol(9), nrows, rowlen 
      integer np2, kk 
      integer frow, felem, colnum, tfields, verbose,mm 
      real(8) eliml, elimh, elmmtpp,elin 
      integer ntptr
      character(10) kion
      character(20) klevl,klevu
      integer lpri,lpril 
      integer jkk, nlev 
      integer nlplmx,nilin2,nlpl,lmm,kltmpn,kltmpo,                     &
     &         llo,lup,llofnd,lupfnd,kk4,                               &
     &         k,kl2,lm,kk2,mlpar,mlm,np1i,np1k,np1r                    
      real(8) elcomp 
      character(20) ktmp2,kblnk20 
!     Database manipulation quantities                                  
      character(1) kblnk,kdtmp(200) 
      integer kltmp(1000) 
      real(4) elsv
      logical done 
!                                                                       
!     Not used                                                          
      real(8) javir 
      integer javi 
!                                                                       
      data kblnk/' '/ 
      data kblnk16/'                '/ 
      data kblnk20/'                    '/ 
!                                                                       
       save kblnk,kblnk16,kblnk20

      if (nlsvn.le.5) return 
                                                                        
      javir=epi(1) 
      epi(1)=javir 
      javir=dpthc(1,1) 
      javi=ncn2 
      np2=javi 
      javi=derivedpointers%npfirst(1) 
      javi=derivedpointers%nplini(1) 
      javi=derivedpointers%npcon(1) 
      javi=derivedpointers%npconi(1) 
      javi=derivedpointers%npilev(1,1) 
      javi=derivedpointers%npilevi(1) 
      javi=derivedpointers%npconi2(1) 
!                                                                       
      verbose=lpri 
      eliml=0.1 
      elimh=1.0e10 
!                                                                       
!     open and prepare the fits file for spectral data                  
      if(verbose.gt.0) write (lun11,*)'writespectra2: opening header'&
     &  ,kmodelname                                                      
      knam='xout_lines1.fits' 
      call fheader(unit,knam,atcredate,kmodelname,istatus) 
      if(istatus.gt.0) call printerror(lun11,istatus) 
                                                                        
!     write extension of parameter values                               
      if(verbose.gt.0)                                                  &
     & write (lun11,*)'writespectra2: write parameter list'             
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,&
     &               parcomm,nloopctl,istatus,lun11)                    
      if(istatus.gt.0) call printerror(lun11,istatus) 
      if(verbose.gt.0)                                                  &
     &  write (lun11,*)'writespectra2: building data tables'            
!                                                                       
!     build spectra data tables                                         
      if (verbose.gt.0) write (lun11,*)' ' 
      kltmpo=0 
      lpril=verbose
      if (verbose.gt.0)                                                 &
     &  write (lun11,*)'emission line luminosities (erg/sec/10**38))'   
      nlplmx=600 
      eliml=0.1 
      elimh=1.0e10 
!     find the strongest lines.                                         
      do  lm=1,nlplmx 
       kltmp(lm)=0 
       enddo 
      nlpl=1 
      do lnn=1,nlsvn 
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml 
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        elin=abs(masterdata%rdat1(np1r)) 
        if ((lrtyp.ne.14).and.(lrtyp.ne.9)) then 
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          nilin2=masterdata%idat1(np1i-1+nidt) 
          elmmtpp=(elum(2,ln)+elum(1,ln))/2. 
          if (verbose.gt.0)                                             &
     &         write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml              
          if ((ln.gt.0).and.(ln.lt.nnnl)                                &
     &         .and.(elin.ge.eliml).and.(elin.le.elimh)                 &
     &         .and.(elin.le.8.9e+6)                                    &
     &         .and.(elmmtpp.gt.1.e-36)                                 &
     &         .and.(nilin2.gt.0).and.(nilin2.le.nni))                  &
     &           then                                                   
!                                                                       
            lmm=0 
            elcomp=1.e+10 
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp)) 
              lmm=lmm+1 
              kl2=kltmp(lmm) 
              elcomp=0. 
              if (kl2.gt.0)                                             &
     &          elcomp=(elum(2,kl2)+elum(1,kl2))/2.                     
              enddo 
!                                                                       
            if (verbose.gt.0)                                           &
     &       write (lun11,8516)ln,elin,elmmtpp                          
 8516       format (1h ,i4,2e12.4) 
            kltmpo=ln 
            do  k=lmm,min(nlplmx,nlpl) 
              if ((lpril.ne.0).and.(kltmp(k).ne.0))                     &
     &         write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo           
              kltmpn=kltmp(k) 
              kltmp(k)=kltmpo 
              kltmpo=kltmpn 
              enddo 
             nlpl=min(nlplmx,nlpl+1) 
            if (verbose.gt.0)                                           &
     &       write (lun11,*)'done with 557 loop',lmm                     
            endif 
          endif 
        enddo 
      if (nlpl.gt.0) kltmp(nlpl)=kltmpo 
!      nlpl=nlpl-1                                                      
!
!     write the spectral data to the extension                          
      do mm=1,9 
        kunits(mm)=kblnk16 
        klabs(mm)=kblnk16 
        kform(mm)=kblnk16 
        enddo 
      klabs(1)='index           ' 
      kform(1)='I6' 
      kunits(1)='  ' 
      klabs(2)='ion             ' 
      kform(2)='A9' 
      kunits(2)=' ' 
      klabs(3)='lower_level     ' 
      kform(3)='A20' 
      kunits(3)='  ' 
      klabs(4)='upper_level     ' 
      kform(4)='A20' 
      kunits(4)='  ' 
      klabs(5)='wavelength      ' 
      kform(5)='E11.3' 
      kunits(5)='A' 
      klabs(6)='emit_inward     ' 
      kform(6)='E11.3' 
      kunits(6)='erg/s/10**38' 
      klabs(7)='emit_outward    ' 
      kform(7)='E11.3' 
      kunits(7)='erg/s/10**38' 
      klabs(8)='depth_inward    ' 
      kform(8)='E11.3' 
      kunits(8)='  ' 
      klabs(9)='depth_outward   ' 
      kform(9)='E11.3' 
      kunits(9)='  ' 
!     build extension name                                              
      extname='XSTAR_LINES' 
!      if(nloopctl.gt.0) then                                           
!          write(ktmp2,'(i4.4)')nloopctl                                
!          extname='xstar_spectra_' // ktmp2                            
!          endif                                                        
      if(verbose.gt.0)                                                  &
     &   write (lun11,*)'writespectra2: writing spectral data'          
                                                                        
!     append a new empty extension onto the end of the primary array    
      status=0 
      call ftcrhd(unit,status) 
      if(verbose.gt.0)                                                  &
     &    write (lun11,*)'writespectra2: writing header table'          
                                                                        
      tfields=9 
      nrows=nlpl 
      rowlen=0 
      do mm=1,9 
        tbcol(mm)=0 
        enddo 
!                                                                        
!     write the required header parameters for the ascii table          
      status=0 
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,   &
     &            extname,status)                                       
      if (status .gt. 0)call printerror(lun11,status) 
      status=0 
!                                                                       
!     map each column to a 1-d array before writing to the file         
      if (verbose.gt.0)                                                 &
     &    write (lun11,959)                                             
  959 format (1x,'index, ion, wavelength, transmitted, reflected') 
      kk2=0 
      do  kk=1,nlpl 
        ln=kltmp(kk) 
        if (ln.ne.0) then 
          ml=derivedpointers%nplin(ln) 
          klevl=kblnk20 
          klevu=kblnk20 
          if (ml.ne.0) then 
            if (verbose.gt.0)                                           &
     &        write (lun11,*)'   ',ln,ml                                
            mlm=ml 
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            llo=masterdata%idat1(np1i) 
            lup=masterdata%idat1(np1i+1) 
            elsv=sngl(abs(masterdata%rdat1(np1r)))
            nilin=derivedpointers%npar(ml) 
            mlm=nilin
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            do mm=1,nkdt 
              kdtmp(mm)=masterdata%kdat1(np1k-1+mm) 
              enddo 
            do mm=nkdt+1,10 
              kdtmp(mm)=kblnk 
              enddo 
            write(kion,'(10a1)')(kdtmp(mm),mm=1,10) 
            done=.false. 
            jkk=1 
            do while (.not.done) 
              ml=derivedpointers%npfi(13,jkk) 
              if (ml.ne.0) then 
                if ((derivedpointers%npar(ml).eq.nilin).or.(jkk.gt.nni))&
     &            done=.true.                                           
                endif 
              jkk=jkk+1 
              enddo 
            if (jkk.gt.nni) ml=0 
            if (ml.ne.0) then 
              mllz=derivedpointers%npar(ml) 
              mlpar=derivedpointers%npar(ml) 
              lupfnd=0 
              llofnd=0 
              do while ((ml.ne.0).and.(mlpar.eq.mllz)                   &
     &           .and.((llofnd.ne.1).or.(lupfnd.ne.1)))                 
                mlm=ml 
                call drd(ltyp,lrtyp,lcon,                               &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                    &
     &            0,lun11)                                        
                nlev=masterdata%idat1(np1i+nidt-2) 
                if (lpri.gt.0)                                          &
     &            call dprinto(ltyp,lrtyp,lcon,                      &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)
                if (lpri.gt.0)                                          &
     &            write (lun11,*)nlev,llo,lup,llofnd,lupfnd             
                if (nlev.eq.llo) then 
                  do mm=1,20 
                    if (mm.le.nkdt) then 
                       write(ktmp2(mm:mm),'(a1)')                      &
     &                       masterdata%kdat1(np1k-1+mm) 
                      else 
                        write(ktmp2(mm:mm),'(a1)')kblnk 
                      endif 
                    enddo 
                  klevl=ktmp2 
                  llofnd=1 
                  endif 
                if (nlev.eq.lup) then 
                  do mm=1,20 
                    if (mm.le.nkdt) then 
                       write(ktmp2(mm:mm),'(a1)')                       &
     &                       masterdata%kdat1(np1k+mm-1) 
                      else 
                        write(ktmp2(mm:mm),'(a1)')kblnk 
                      endif 
                    enddo 
                  klevu=ktmp2 
                  lupfnd=1 
                  endif 
                ml=derivedpointers%npnxt(ml) 
                mlpar=0 
                if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                enddo 
              kk2=kk2+1 
              kk2=min(kk2,600)
              ntptr=ln 
              frow=kk2
              felem=1 
              nrows=1
              kk4=1 
              colnum=kk4
              call ftpclj(unit,colnum,frow,felem,nrows,ntptr,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,ntptr,status
              kk4=2 
              colnum=kk4 
              call ftpcls(unit,colnum,frow,felem,nrows,kion,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,kion,status
              kk4=3 
              colnum=kk4 
              call ftpcls(unit,colnum,frow,felem,nrows,klevl,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,klevl,status
              kk4=4 
              colnum=kk4 
              call ftpcls(unit,colnum,frow,felem,nrows,klevu,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,klevu,status
              kk4=5 
              colnum=kk4
              call ftpcle(unit,colnum,frow,felem,nrows,elsv,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,elsv,status
              kk4=6 
              colnum=kk4 
              rtmp=sngl(elum(1,ntptr))
              call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,rtmp,status
              kk4=7 
              colnum=kk4 
              rtmp=sngl(elum(2,ntptr)) 
              call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,rtmp,status
              kk4=8 
              colnum=kk4
              rtmp=sngl(tau0(1,ntptr))
              call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,rtmp,status
              kk4=9 
              colnum=kk4
              rtmp=sngl(tau0(2,ntptr))
              call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status) 
              if (status .gt. 0)call printerror(lun11,status) 
              if (lpri.gt.0)                                            &
     &         write (lun11,*)unit,colnum,frow,felem,nrows,rtmp,status
              if (verbose.gt.0) then 
                write (lun11,*)ml,nilin,derivedpointers%npar(ml) 
                write (lun11,9955)kk,ln,(kdtmp(mm),mm=1,9),elsv,        &
     &               elum(1,ln),elum(2,ln)                              
                write (lun11,*)klevu
                write (lun11,*)klevl
                endif 
 9955           format (1x,2i8,1x,9a1,3(1pe11.3)) 
              endif 
            endif 
          endif 
        enddo 
!      if (nlpl.le.0) return                                            
!                                                                       
      nlpl=kk2 
      nlpl=max(nlpl,1) 
!                                                                       
                                                                        
                                                                        
!     compute checksums                                                 
      if(verbose.gt.0) write (lun11,*)'writespectra2:writingchecksum' 
      status=0 
      call ftpcks(unit,status) 
!     check for any error, and if so print out error messages           
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if(verbose.gt.0) write (lun11,*)'writespectra2: closing file' 
      call fitsclose(lun11,unit,istatus) 
!                                                                       
!                                                                       
      return 
      end                                           
