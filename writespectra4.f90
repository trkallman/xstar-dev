      subroutine writespectra4(lun11,lpri,nparms,parname,partype,parval,&
     &       parcomm,atcredate,epi,ncn2,dpthc,abel,                     &
     &       np2,ncsvn,                                                 &
     &       elumab,tauc,kmodelname,nloopctl)                   
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
      real*8 parval(55) 
      character(30) parcomm(55) 
      real*8 elumab(2,nnml),tauc(2,nnml) 
!     energy bins                                                       
      real*8 epi(ncn) 
!     the atomic data creation date                                     
      character(63) atcredate 
!     continuum optical depths                                          
      real*8 dpthc(2,ncn) 
      real*8 abel(nl) 
      character(8) kdtmpi(nnml),kdtmp8 
      character(20) kdtmpl(nnml)
      character(1) kblnk 
      character(16) knam,klabs(8),kunits(8),kform(8),kblnk16 
      character(30) extname 
      character(20) ktmp20
      integer unit,istatus,kl,nkdt,nidt,lcon,lrtyp,ltyp,ml,nkdti 
      integer nrdt 
      integer np2, tbcol(8), nrows, rowlen, kk 
      integer frow, felem, colnum, tfields, status, verbose,mm 
      real*8 eth,xeltp 
      integer lpril,lpri,klel,mlel,jk,mt2,mllel,nnz,jkk,klion,mlion,    &
     &        mlleltp,nlevmx,mltype,mllz,nlev,lk,kkkl,idest1,           &
     &        kksv,mlpar,mlm,np1k,np1ki,np1i,np1r, ncsvn
      REAL*4, DIMENSION(:), ALLOCATABLE :: rsv1
      REAL*4, DIMENSION(:), ALLOCATABLE :: rsv2
      REAL*4, DIMENSION(:), ALLOCATABLE :: rsv3
      REAL*4, DIMENSION(:), ALLOCATABLE :: rsv4
      REAL*4, DIMENSION(:), ALLOCATABLE :: rsv5
      INTEGER, DIMENSION(:), ALLOCATABLE :: ntptr
!      real*4 rsv1(nnml),rsv2(nnml),rsv3(nnml),rsv4(nnml),rsv5(nnml) 
!      integer ntptr(nnml) 
!                                                                       
!     Not used                                                          
      real*8 javir 
      integer javi 
!                                                                       
      data kblnk/' '/ 
      data kblnk16/'                '/ 
                                                                        
      ALLOCATE(rsv1(nnml))
      ALLOCATE(rsv2(nnml))
      ALLOCATE(rsv3(nnml))
      ALLOCATE(rsv4(nnml))
      ALLOCATE(rsv5(nnml))
      ALLOCATE(ntptr(nnml))
!
      javi=lpri 
      javir=epi(1) 
!      epi(1)=javir                                                     
      javir=dpthc(1,1) 
      javi=ncn2 
      javi=np2 
!      np2=javi                                                         
      javi=derivedpointers%nplini(1) 
      javi=derivedpointers%nplin(1) 
      javi=derivedpointers%npcon(1) 
      javi=derivedpointers%npconi(1) 
      javi=derivedpointers%npilev(1,1) 
      javi=derivedpointers%npilevi(1) 
                                                                        
!                                                                       
                                                                        
      lpril=lpri 
      verbose=lpri 
!     open and prepare the fits file for spectral data                  
      if(verbose.gt.0) write (lun11,*)'writespectra4: opening header',  &
     &  kmodelname                                                      
      knam='xout_rrc1.fits' 
      call fheader(unit,knam,atcredate,kmodelname,istatus) 
      if(istatus.gt.0) call printerror(lun11,istatus) 
                                                                        
!     write extension of parameter values                               
      if(verbose.gt.0)                                                  &
     &     write (lun11,*)'writespectra4: write parameter list'         
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,   &
     &               parcomm,nloopctl,istatus,lun11)                    
      if(istatus.gt.0) call printerror(lun11,istatus) 
      if(verbose.gt.0)                                                  &
     &  write (lun11,*)'writespectra4: building data tables'            
                                                                        
!     build spectra data tables                                         
      lpril=verbose 
!     print 500 strongest recombination continua                        
!      write (lun11,*)'recombination continuum luminosities',           
!     $  '(erg/sec/10**38))'                                            
!      write (lun11,*)'ion, level, energy (eV), RRC luminosity '        
!     lpril is flag for printing debug information                      
!      initialize line counter                                          
      kksv=0 
      jkk=0 
!      First look for element data (jk is element index)                
        klel=11 
        mlel=derivedpointers%npfirst(klel) 
        jk=0 
        do while (mlel.ne.0) 
          jk=jk+1 
          mt2=mlel-1 
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                          &
     &      0,lun11)                                              
          if (nidt.gt.0) then 
            mllel=masterdata%idat1(np1i-1+nidt) 
            xeltp=masterdata%rdat1(np1r) 
            xeltp=abel(mllel) 
            nnz=masterdata%idat1(np1i) 
            if (lpril.ne.0)                                             &
     &        write (lun11,*)'element:',jk,mlel,mllel,nnz,              &
     &           (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)                  
!           ignore if the abundance is small                            
            if (xeltp.lt.1.e-10) then 
                jkk=jkk+nnz 
              else 
!               now step thru ions (jkk is ion index)                   
                klion=12 
                mlion=derivedpointers%npfirst(klion) 
                jkk=0 
                kl=0 
                do while ((mlion.ne.0).and.(kl.lt.nnz)) 
                  jkk=jkk+1 
!                 retrieve ion name from kdati                          
                  mlm=mlion-1 
                  call drd(ltyp,lrtyp,lcon,                             &
     &              nrdt,np1r,nidt,np1i,nkdti,np1ki,mlm,                &
     &              0,lun11)                                      
!                 if not accessing the same element, skip to the next el
                  mlleltp=masterdata%idat1(np1i+nidt-2) 
                  if (mlleltp.eq.mllel) then 
                    kl=kl+1 
                    if (lpril.ne.0)                                     &
     &                write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,     &
     &                   (masterdata%kdat1(np1k-1+mm),mm=1,nkdti)         
!                   now find level data                                 
!                   step thru types                                     
                    nlevmx=0 
                    mltype=13 
                    ml=derivedpointers%npfi(mltype,jkk) 
                    mllz=derivedpointers%npar(ml) 
                    mlpar=derivedpointers%npar(ml) 
!                   step thru records of this type                      
                    do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
                      mlm=ml-1 
                      call drd(ltyp,lrtyp,lcon,                         &
     &                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,              &
     &                  0,lun11)                                  
                      nlev=masterdata%idat1(np1i+nidt-2) 
                      nlevmx=max(nlevmx,nlev) 
                      if ((nlev.gt.0).and.(nlev.le.nd)) then 
                        if (lpril.ne.0)                                 &
     &                    write (lun11,*)'level quantities:',           &
     &                    ml,nlev,ltyp,lrtyp,masterdata%rdat1(np1r),    &
     &                    masterdata%rdat1(np1r+1)  
                        do  lk=1,nrdt 
                          leveltemp%rlev(lk,nlev)                       &
     &                            =masterdata%rdat1(np1r-1+lk) 
                          enddo 
                        do lk=1,nkdt 
                          leveltemp%klev(lk,nlev)                       &
     &                            =masterdata%kdat1(np1k-1+lk) 
                          enddo 
                        do lk=nkdt+1,20 
                          leveltemp%klev(lk,nlev)=kblnk 
                          enddo 
                        endif 
                      ml=derivedpointers%npnxt(ml) 
                      if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                      enddo 
                    nlev=nlevmx 
                    mltype=7 
                    ml=derivedpointers%npfi(mltype,jkk) 
                    mllz=derivedpointers%npar(ml) 
                    mlpar=derivedpointers%npar(ml) 
                    do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
!                     step thru records of this type                    
                      mlm=ml-1 
                      call drd(ltyp,lrtyp,lcon,                         &
     &                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,              &
     &                  0,lun11)                                  
                      kkkl=derivedpointers%npconi2(ml) 
                      idest1=masterdata%idat1(np1i+nidt-2) 
                      if ((kkkl.gt.0).and.(kkkl.le.ndat2)               &
     &                  .and.((elumab(1,kkkl).gt.1.e-36)                &
     &                  .or.(elumab(2,kkkl).gt.1.e-36))) then           
                        kksv=kksv+1 
                        eth=leveltemp%rlev(4,idest1)                    &
     &                     -leveltemp%rlev(1,idest1) 
                        ntptr(kksv)=kkkl 
                        rsv1(kksv)=eth 
                        rsv2(kksv)=elumab(1,kkkl) 
                        rsv3(kksv)=elumab(2,kkkl) 
                        rsv4(kksv)=tauc(1,kkkl) 
                        rsv5(kksv)=tauc(2,kkkl) 
                        do mm=1,nkdti 
                           write (kdtmp8(mm:mm),'(a1)')                 &
     &                           masterdata%kdat1(np1ki-1+mm) 
                          enddo 
                        do mm=nkdti+1,8 
                          write (kdtmp8(mm:mm),'(a1)')kblnk 
                          enddo 
                        kdtmpi(kksv)=kdtmp8 
                        do lk=1,20 
                          write (ktmp20(lk:lk),'(a1)')                  &
     &                        leveltemp%klev(lk,idest1) 
                          enddo 
                        kdtmpl(kksv)=ktmp20 
                        if (lpril.ne.0)                                 &
     &                   write (lun11,*)jkk,idest1,                     &
     &                   eth,elumab(1,kkkl),elumab(2,kkkl)              
                        if (lpril.ne.0)                                 &
     &                   write (lun11,9293)kdtmpi(kksv),                &
     &                      (leveltemp%klev(lk,idest1),lk=1,20),eth,    &
     &                          elumab(1,kkkl)                          
 9293                   format(1x,20a1,20a1,2(1pe11.3)) 
                        endif 
                      ml=derivedpointers%npnxt(ml) 
                      if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                      enddo 
                    endif 
!                 Go to next ion                                        
                  mlion=derivedpointers%npnxt(mlion) 
                  enddo 
              endif 
            endif 
          if (mlel.ne.0) mlel=derivedpointers%npnxt(mlel) 
!         Go to next element                                            
          enddo 
                                                                        
!     write the spectral data to the extension                          
      do mm=1,8 
        kunits(mm)=kblnk16 
        klabs(mm)=kblnk16 
        kform(mm)=kblnk16 
        enddo 
      klabs(1)='index          ' 
      kform(1)='I6' 
      kunits(1)='  ' 
      klabs(2)='ion            ' 
      kform(2)='A9' 
      kunits(2)='  ' 
      klabs(3)='level          ' 
      kform(3)='A20' 
      kunits(3)='  ' 
      klabs(4)='energy         ' 
      kform(4)='E11.3' 
      kunits(4)='eV' 
      klabs(5)='emit_outward    ' 
      kform(5)='E11.3' 
      kunits(5)='erg/s' 
      klabs(6)='emit_inward     ' 
      kform(6)='E11.3' 
      kunits(6)='erg/s' 
      klabs(7)='depth_outward   ' 
      kform(7)='E11.3' 
      kunits(7)='  ' 
      klabs(8)='depth_inward    ' 
      kform(8)='E11.3' 
      kunits(8)='  ' 
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
                                                                        
      tfields=8 
      nrows=kksv 
      if (lpril.ne.0) write (6,*)'nrows=',nrows 
      rowlen=0 
      do mm=1,8 
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
      kk=1 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpclj(unit,colnum,frow,felem,nrows,ntptr,status) 
        if (status .gt. 0)call printerror(lun11,status) 
      kk=2 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpcls(unit,colnum,frow,felem,nrows,kdtmpi,status) 
        if (status .gt. 0)call printerror(lun11,status) 
      kk=3 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpcls(unit,colnum,frow,felem,nrows,kdtmpl,status) 
        if (status .gt. 0)call printerror(lun11,status) 
      kk=4 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpcle(unit,colnum,frow,felem,nrows,rsv1,status) 
        if (status .gt. 0)call printerror(lun11,status) 
      kk=5 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpcle(unit,colnum,frow,felem,nrows,rsv2,status) 
        if (status .gt. 0)call printerror(lun11,status) 
      kk=6 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpcle(unit,colnum,frow,felem,nrows,rsv3,status) 
        if (status .gt. 0)call printerror(lun11,status) 
      kk=7 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpcle(unit,colnum,frow,felem,nrows,rsv4,status) 
        if (status .gt. 0)call printerror(lun11,status) 
      kk=8 
        if(verbose.gt.0)                                                &
     &      write (lun11,*)'writespectra2: building column ',kk         
        frow=1 
        felem=1 
        colnum=kk 
        if(verbose.gt.0)                                                &
     &     write (lun11,*)'writespectra2: writing column ',kk           
        status=0 
        call ftpcle(unit,colnum,frow,felem,nrows,rsv5,status) 
        if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     compute checksums                                                 
      if(verbose.gt.0) write (lun11,*)'writespectra: writing checksum' 
      status=0 
      call ftpcks(unit,status) 
!     check for any error, and if so print out error messages           
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if(verbose.gt.0) write (lun11,*)'writespectra: closing file' 
      call fitsclose(lun11,unit,istatus) 
!                                                                       
                                                                        
      DEALLOCATE(rsv1)
      DEALLOCATE(rsv2)
      DEALLOCATE(rsv3)
      DEALLOCATE(rsv4)
      DEALLOCATE(rsv5)
      DEALLOCATE(ntptr)
!                                                                       
      return 
      end                                           
