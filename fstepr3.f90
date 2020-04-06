      subroutine fstepr3(unit,hdunum,radin,radout,rdel,t,prs,abel,   &
     &             xcol,xee,xpx,xi,                                     &
     &             np2,                                                 &
     &             cemab,cabab,opakab,tauc,                             &
     &             lun11,lpri,status)                                   
!                                                                       
!     Name: fstepr3.f90
!     Description:
!       Write rrc quantities for each radial zone to an individual 
!       extension of the file xoxx_detal3.fits
!       Append a FITS extension binary table containing                   
!       nrhs columns and at most nrhdimj rows                             
!       author: T. Bridgman                                               
!     Parameters:                               
!        Input:                        
!        unit    integer            File unit number                    
!        hdunum  integer            Number of last HDU written          
!        radin   real(8)               inner radius of shell             
!        radout  real(8)               outer radius of shell             
!        delr    real(8)               thickness of shell
!        temp    real(8)               temperature of shell in 10^4K          
!        pres    real(8)               pressure in shell                 
!        cemab(nnml):               rrc emissivities (erg cm^-3 s^-1) 
!        cabab(nnml):               total energy absorbed by 
!                                      rrc (erg cm^-3 s^-1)
!        opakab(nnml):              rrc opacities (cm^-1)
!        tauc(2,nnml):              rrc optical depths
!        lun11                      logical unit number for printing
!        lpri                       print switch
!        Output:
!        status  integer            Returned status code                
!     Dependencies:  none
!     called by:  xstar       
!                                                                       
!                                                                       
      use globaldata
      implicit none 
      integer mllz,mllz2 
                                                                        
      integer nptmpdim 
      parameter (nptmpdim=500000) 
!                                                                       
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,ndl) 
        integer:: ilev(10,ndl),nlpt(ndl),iltp(ndl) 
        character(1) :: klev(100,ndl) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
      real(4) rtmp 
      real(8) radin, radout,rdel, t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, np2
!     line opacities                                                    
      real(8) tauc(2,nnml) 
      real(8) cemab(2,nnml),opakab(nnml),cabab(nnml) 
      real(8) abel(nl) 
                                                                        
                                                                        
!     Internal work areas                                               
      integer nlev 
      real(4), dimension(:), allocatable :: rwrk1,rwrk3,rwrk4,rwrk5,  &
     &  rwrk6,rwrk7,rwrk8                 
      integer, dimension(:), allocatable :: ntptr,ntptr2
      character(20), dimension(:), allocatable :: klevl,klevu
      character(8), dimension(:), allocatable :: kion
      integer, dimension(:), allocatable :: ilevlo,ilevup 
      real(4), dimension(:), allocatable :: elsv
      integer tfields,varidat 
      character(16) ttype(12),tform(12),tunit(12) 
      integer colnum,frow,felem,hdutype,ll, ltyp 
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpri 
      integer jkk,nidti 
      integer ml,nlpl,                                                  &
     &         nlines,                                                  &
     &         kk,mlm,kkkl,idest1,idest2,lk,                            &
     &         mlel,mlion,mllel,mlleltp,mlpar,mt2,mltype,nkdti,         &
     &         jk,kl,klel,klion,nnz,np1ki,mmlv             
      integer np1i,np1r,np1k 
      real(8) eth,xeltp 
      character(33) extname 
      character(20) kblnk20 
      character(20) ktmp20 
      character(8) ktmp8
!     needed for upper level search                                     
      integer jkk3,nlevp,ndtmp,iltmp,lcon2,lrtyp2,ltyp2,                &
     &         np1r2,nrdt2,np1i2,nidt2,np1k2,nkdt2                      
      real(8) ett 
      integer nnzz,nnnn
                                                                        
!     Database manipulation quantities                                  
      integer nkdt 
      character(1) kblnk
                                                                        
      data kblnk/' '/ 
      data kblnk20/'                    '/ 
!                                                                       
      data tform/'1J','1J','1E','8A','20A','20A','1E','1E','1E',        &
     & '1E','1E','1E'/                                                  
                                                                        
      data ttype/'rrc index','level index','energy','ion',              &
     & 'lower_level','upper_level','emis_inward',                       &
     & 'emis_outward','integrated absn','opacity',                      &
     &  'tau_in','tau_out'/                                             
                                                                        
      data tunit/' ',' ','ev',' ',' ',' ','erg/cm^3/s',                 &
     & 'erg/cm^3/s','erg/cm^3/s','/cm',' ',' '/                         

       save kblnk,kblnk20,tform,ttype,tunit
                                                                        
     allocate(rwrk1(nptmpdim))
     allocate(rwrk3(nptmpdim))
     allocate(rwrk4(nptmpdim))
     allocate(rwrk5(nptmpdim))
     allocate(rwrk6(nptmpdim))
     allocate(rwrk7(nptmpdim))
     allocate(rwrk8(nptmpdim))
     allocate(ntptr(nptmpdim))
     allocate(ntptr2(nptmpdim))
     allocate(klevl(nptmpdim))
     allocate(klevu(nptmpdim))
     allocate(kion(nptmpdim))
     allocate(ilevlo(nptmpdim))
     allocate(ilevup(nptmpdim))
     allocate(elsv(nptmpdim))
!
      varidat=0 
!                                                                       
!                                                                       
                                                                        
!     Move to the last HDU (hdunum) in the file                         
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Moving to end-of-FITS file'              
      call ftmahd(unit,hdunum,hdutype,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     append a new empty extension after the last HDU                   
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr3: Create the new extension'               
      call ftcrhd(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
!                                                                       
!     Extracting data from the Atomic Database here                     
!                                                                       
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'in fstepr3 '                                     
!                                                                       
!     First look for element data (jk is element index)                 
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      jk=0 
      kk=0 
      jkk=0 
!                                                                       
!     step through elements                                             
      do while (mlel.ne.0) 
!                                                                       
!       get element data                                                
        jk=jk+1 
        mt2=mlel 
        call drd(ltyp,lrtyp,lcon,                                       &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                        &
     &        0,lun11)                                            
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=masterdata%rdat1(np1r) 
        xeltp=abel(mllel) 
        nnz=masterdata%idat1(np1i) 
        if (lpri.ge.1)                                                  &
     &        write (lun11,*)'element:',jk,mlel,mllel,nnz,              &
     &           (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)                    
!                                                                       
!       ignore if the abundance is small                                
        if (xeltp.lt.1.e-10) then 
            jkk=jkk+nnz 
          else 
!                                                                       
!           now step thru ions (jkk is ion index)                       
            klion=12 
            mlion=derivedpointers%npfirst(klion) 
            jkk=0 
            kl=0 
            do while ((mlion.ne.0).and.(kl.lt.nnz)) 
              jkk=jkk+1 
!                                                                       
!             retrieve ion name from kdati                              
              mlm=mlion 
              call drd(ltyp,lrtyp,lcon,                                 &
     &            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,                 &
     &            0,lun11)                                        
              nnzz=masterdata%idat1(np1i+1)
              nnnn=nnzz-masterdata%idat1(np1i)+1
!                                                                       
!             if not accessing the same element, skip to the next elemen
              mlleltp=masterdata%idat1(np1i+nidti-2) 
              if (mlleltp.eq.mllel) then 
!                                                                       
                kl=kl+1 
                if (lpri.ge.1)                                          &
     &            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,         &
     &               (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)            
!                                                                       
!               now find level data                                     
                call calc_rates_level_lte(jkk,lpri,lun11,t,xee,xpx,     &
     &              nnzz,nnnn,leveltemp,nlev)
!                                                                       
!               now step through rate type 7 data                       
                mltype=7 
                ml=derivedpointers%npfi(mltype,jkk) 
                mllz=0 
                if (ml.ne.0) mllz=derivedpointers%npar(ml) 
                mlpar=0 
                if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
!                                                                       
!                 get rrc data                                          
                  kkkl=derivedpointers%npconi2(ml) 
                  if (lpri.gt.0) write (lun11,*)kkkl,ml,idest1,         &
     &                    cemab(1,kkkl),cemab(2,kkkl)                   
!                                                                       
!                 test for non-zero rrc data                            
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)                   &
     &                .and.((cemab(1,kkkl).gt.1.e-36)                   &
     &                .or.(cemab(2,kkkl).gt.1.e-36)                     &
     &                .or.(cabab(kkkl).gt.1.e-36)                       &
     &                .or.(opakab(kkkl).gt.1.e-36))) then               
!                                                                       
!                                                                       
!                   increment buffer counter                            
                    kk=kk+1 
                    if (kk.ge.nptmpdim) then 
                      write (lun11,*)'buffer overflow in fstepr33' 
                      kk=nptmpdim 
                      endif 
!                                                                       
!                   get rrc  data                                       
                    mlm=ml 
                    call drd(ltyp,lrtyp,lcon,                           &
     &                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                &
     &                0,lun11)                                    
                    idest1=masterdata%idat1(np1i+nidt-2) 
                    nlevp=nlev 
                    idest2=nlevp+masterdata%idat1(np1i-1+nidt-3)-1 
!                                                                       
!                   label for lower level                               
                    do lk=1,20 
                      write (ktmp20(lk:lk),'(a1)')                      &
     &                        leveltemp%klev(lk,idest1) 
                      enddo 
                    klevl(kk)=ktmp20 
!                                                                       
!                   label for upper level                               
                    write (ktmp20(1:20),'(a20)')'continuum           ' 
                    klevu(kk)=ktmp20 
!                                                                       
!                   ion label                                           
                    do lk=1,nkdti 
                       write (ktmp8(lk:lk),'(a1)')                      &
     &                       masterdata%kdat1(np1ki+lk-1) 
                      enddo 
                    do lk=nkdti+1,8 
                      write (ktmp8(lk:lk),'(a1)')kblnk 
                      enddo 
!                                                                       
                    eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
                    ett=eth 
!                                                                       
                    go to 9009
!                   get upper level data                                
                    if (idest2.gt.nlevp) then 
                      jkk3=jkk+1 
                      if (lpri.ge.1)                                    &
     &                  write (lun11,*)'upper level',                   &
     &                     jkk3,ndtmp,nlevp,idest2          
                      ndtmp=derivedpointers%npfi(13,jkk3) 
                      if (lpri.ge.1)                                    &
     &                  write (lun11,*)jkk3,ndtmp,nlevp,idest2          
                      if ((ndtmp.le.0).or.(ndtmp.gt.np2))               &
     &                    stop 'ndtmp error' 
                      mllz2=derivedpointers%npar(ndtmp) 
                      iltmp=0 
                      do while ((ndtmp.ne.0).and.(ndtmp.le.np2).and.    &
     &                    (iltmp.ne.(idest2-nlevp+1)).and.              &
     &                    (derivedpointers%npar(ndtmp).eq.mllz2))              
                        mlm=ndtmp 
                        call drd(ltyp2,lrtyp2,lcon2,                    &
     &                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,      &
     &                    0,lun11)                                
                        iltmp=masterdata%idat1(np1i2+nidt2-2) 
                        if (lpri.ge.1) write (lun11,*)nidt2,iltmp,ndtmp 
                        ndtmp=derivedpointers%npnxt(ndtmp) 
                        if ((ndtmp.le.0).or.(ndtmp.gt.np2))             &
     &                          stop 'ndtmp error' 
                        enddo 
!                     NB fix to excited level PI and rec                
                      ett=ett+masterdata%rdat1(np1r2) 
                      eth=ett 
                      if (lpri.gt.1)                                    &
     &                  write (lun11,*) ndtmp,iltmp,idest2,ett          
!                     label for lower level                             
                      ktmp20=kblnk20 
                      do lk=1,nkdt2 
                         write (ktmp20(lk:lk),'(a1)')                   &
      &                        masterdata%kdat1(np1k2+lk-1) 
                        enddo 
                      klevu(kk)=ktmp20 
                      endif 
9009                  continue
!                                                                       
!                   other data                                          
                    kion(kk)=ktmp8 
                    elsv(kk)=sngl(eth)
                    ilevlo(kk)=idest1 
                    ilevup(kk)=idest2 
                    ntptr(kk)=kkkl 
                    mmlv=derivedpointers%npilev(idest1,jkk) 
                    ntptr2(kk)=mmlv 
                    if (lpri.ge.1)                                      &
     &                  write (lun11,981)kkkl,eth,idest1,               &
     &                    cemab(1,kkkl),cemab(2,kkkl)                   
  981                 format (1x,i6,1pe11.3,i6,6(1pe11.3)) 
!                                                                       
!                   done with this rrc                                  
                    endif 
!                                                                       
!                 end of loop over rrcs                                 
                  ml=derivedpointers%npnxt(ml) 
                  if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                  enddo 
!                                                                       
!               end of test for element                                 
                endif 
!                                                                       
!             Go to next ion                                            
              mlion=derivedpointers%npnxt(mlion) 
              enddo 
!                                                                       
!         end of test for non-zero element abund                        
          endif 
!                                                                       
        mlel=derivedpointers%npnxt(mlel) 
!       Go to next element                                              
        enddo 
!                                                                       
      nlpl=max(nlpl,1) 
!                                                                       
                                                                        
!     End of atomic database extraction                                 
!----------------------------------------------------------------       
!     define parameters for the binary table (see the above data stateme
      nrows=kk 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'before header write'                             
      tfields=12 
!     Build extension name                                              
      extname='XSTAR_RADIAL' 
                                                                        
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr3: Write table headers'                    
!     write the required header parameters for the binary table         
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,         &
     &              varidat,status)                                     
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr3: Add some more keywords'                 
                                                                        
!     Write some model parameters in the extension header               
      call ftpcom(unit,'***********************************',status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftpcom(unit,'Model Keywords',status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     Write values to 3 decimal places                                  
      rtmp=sngl(radin)
      call ftpkye(unit,'RINNER',rtmp,3,'[cm] Inner shell radius',       &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=sngl(radout)
      call ftpkye(unit,'ROUTER',rtmp,3,'[cm] Outer shell radius',       &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=sngl(rdel)
      call ftpkye(unit,'RDEL',rtmp,3,'[cm] distance from face',         &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=sngl(t)
      call ftpkye(unit,'TEMPERAT',rtmp,3,'[10**4K] Shell Temperature',  &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=sngl(prs) 
      call ftpkye(unit,'PRESSURE',rtmp,3,'[dynes/cm**2] Shell Pressure',&
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      rtmp=sngl(xcol)
      call ftpkye(unit,'COLUMN',rtmp,3,'[/cm**2] Column ',              &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=sngl(xee)
      call ftpkye(unit,'XEE',rtmp,3,'electron fraction',                &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      rtmp=sngl(xpx)
      call ftpkye(unit,'DENSITY',rtmp,3,'[/cm**3] Density',             &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      rtmp=sngl(xi)
      call ftpkye(unit,'LOGXI',rtmp,3,                                  &
     & '[erg cm/s] log(ionization parameter)',status)                   
      if (status .gt. 0)call printerror(lun11,status) 

      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after header write'                              
!-------------------------------------------------------------------    
!     Step through the columns and write them to the file               
!                                                                       
!     set 'global' parameters for writing FITS columns                  
      frow=1 
      felem=1 
                                                                        
                                                                        
!     column  1  (continuum index)                                      
      colnum=1 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum,nlpl             
      nlines=kk 
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  2  (level index)                                          
      colnum=2 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum,nlpl             
      nlines=kk 
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr2,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  3  (wavelength)                                           
      colnum=3 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,elsv,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
                                                                        
!     column  4  (Ion)                                                  
      colnum=4 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcls(unit,colnum,frow,felem,nlines,kion,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  5 (lower Level Designation)                               
      colnum=5 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr3: Writing Column ',colnum                 
      call ftpcls(unit,colnum,frow,felem,nlines,klevl,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  6 (Level Designation)                                     
      colnum=6 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcls(unit,colnum,frow,felem,nlines,klevu,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
                                                                        
      do ll=1,nlines 
         rwrk1(ll)=0. 
         if (ntptr(ll).ne.0) then 
           rwrk3(ll)=sngl(cemab(1,ntptr(ll)))
           rwrk4(ll)=sngl(cemab(2,ntptr(ll)))
           rwrk5(ll)=sngl(cabab(ntptr(ll)))
           rwrk6(ll)=sngl(opakab(ntptr(ll)))
           rwrk7(ll)=sngl(tauc(1,ntptr(ll)))
           rwrk8(ll)=sngl(tauc(2,ntptr(ll)))
           endif 
         enddo 
!                                                                       
!                                                                       
!     column  7                                                         
      colnum=7 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk3,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!                                                                       
!     column  8                                                         
      colnum=8 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk4,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!                                                                       
!     column  9                                                         
      colnum=9 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk5,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!                                                                       
!     column  10                                                        
      colnum=10 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk6,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!                                                                       
!     column  11                                                        
      colnum=11 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr3: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk7,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  12                                                        
      colnum=12 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr3: Writing Column ',colnum                 
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk8,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
                                                                        
!----------------------------------------------------------------       
!     Compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
     deallocate(rwrk1)
     deallocate(rwrk3)
     deallocate(rwrk4)
     deallocate(rwrk5)
     deallocate(rwrk6)
     deallocate(rwrk7)
     deallocate(rwrk8)
     deallocate(ntptr)
     deallocate(ntptr2)
     deallocate(klevl)
     deallocate(klevu)
     deallocate(kion)
     deallocate(ilevlo)
     deallocate(ilevup)
     deallocate(elsv)
 !
      return 
      end                                           
