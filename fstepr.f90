      subroutine fstepr(unit,hdunum,radin,radout,rdel,t,pres,        &
     &                xcol,xee,xpx,xi,                                  &
     &                xilev,rnist,                                      &
     &                lun11,lpri,status)                                
!                                                                       
!                                                                       
!     Name: fstepr.f90
!     Description:
!       Write level populations for each radial zone to an individual 
!       extension of the file xoxx_detail.fits
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
!        xilev   real(nrhdimj)       Fractional level population array  
!        rnist   real(nrhdimj)       LTE level populations
!        lun11                      logical unit number for printing
!        lpri                       print switch
!        Output:
!        status  integer            Returned status code                
!     Dependencies:  none
!     called by:  savd
!                                                        
      use globaldata
      implicit none 
!                                                                       
!     Allocation for passed parameters                                  
      real(8) xilev(nnml),rnist(nnml)
      real(8) radin, radout,rdel, t, pres, xcol,xee,xpx,xi
      real(4) rtmp 
      integer unit,hdunum, nrows, status
                                                                        
      real(4), dimension(:), allocatable :: rwrk1,rwrk2, elev
      integer, dimension(:), allocatable :: ntptr,natomic,mllev,nupper
      character(10), dimension(:), allocatable :: kion
      character(20), dimension(:), allocatable :: klevt
      integer tfields,varidat 
      character(16) ttype(9),tform(9),tunit(9) 
      integer colnum,frow,felem,hdutype, klel, mlel, jk, ltyp 
      integer lrtyp, lcon, nrdt, nidt, mmlv, mm, lun11, lpril,lpri 
      integer mllel, klion, mlion, jkk, kl
      integer mt2, mlleltp, nnz, nions 
      character(43) extname 
!     Database manipulation quantities                                  
      real(8)  xeltp 
      integer  nkdt 
      integer j,nkdti,np1ki 
      integer nlev 
      integer mm2,mmtmp,kkkl,lk,mlm 
      integer np1i,np1r,np1k
      real(8) eth 
      character(10) kdtmp 
                                                                        
      data tform/'1J','1I','1E','8A','1I','20A','1E','1E','1I'/ 
      data ttype/'index','ion_index','e_excitation','ion',              &
     & 'atomic_number','ion_level','population','lte',                  &
     &  'upper index'/                                                  
      data tunit/' ',' ','eV',' ',' ',' ',' ',' ',' '/ 
!                                                                       
      allocate(rwrk1(nnml))
      allocate(rwrk2(nnml))
      allocate(elev(nnml))
      allocate(ntptr(nnml))
      allocate(natomic(nnml))
      allocate(mllev(nnml))
      allocate(nupper(nnml))
      allocate(kion(nnml))
      allocate(klevt(nnml))
!
      lpril=lpri 
      varidat=0 
!                                                                       
      status=0 
!                                                                       
!     Move to the last HDU (hdunum) in the file                         
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Moving to end-of-FITS file'               
      call ftmahd(unit,hdunum,hdutype,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!                                                                       
!     append a new empty extension after the last HDU                   
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Create the new extension'                 
      call ftcrhd(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     Extracting data from the Atomic Database here                     
!                                                                       
!                                                                       
!     lpril is flag for printing debug information                      
       nions=0 
      if (lpril.ne.0) then 
        write (lun11,*)'raw data' 
        do j=1,nnml 
          if (xilev(j).gt.1.e-37)                                       &
     &     write (lun11,*)j,xilev(j)                                    
          enddo 
        endif 
!                                                                       
!     initialize line counter                                           
      mmlv=0 
!     First look for element data (jk is element index)                 
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      jkk=0 
      jk=0 
      do while (mlel.ne.0) 
!                                                                       
!       get element data                                                
       jk=jk+1 
        mt2=mlel 
        call drd(ltyp,lrtyp,lcon,                                       &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                           &
     &     0,lun11)                                               
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=masterdata%rdat1(np1r) 
        nnz=masterdata%idat1(np1i) 
        if (lpril.ne.0)                                                 &
     &    write (lun11,*)'element:',jk,mlel,mllel,nnz,                  &
     &    (masterdata%kdat1(np1k-1+mm),mm=1,nkdt),xeltp              
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
!                                                                       
              jkk=jkk+1 
!             retrieve ion name from kdati                              
              mlm=mlion 
              call drd(ltyp,lrtyp,lcon,                                 &
     &            nrdt,np1r,nidt,np1i,nkdti,np1ki,mlm,                  &
     &            0,lun11)                                        
                                                                        
!             if not accessing the same element, skip to the next elemen
              mlleltp=masterdata%idat1(np1i+nidt-2) 
              if (mlleltp.eq.mllel) then 
!                                                                       
                kl=kl+1 
                if (lpril.ne.0)                                         &
     &              write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,       &
     &                 (masterdata%kdat1(np1ki-1+mm),mm=1,nkdti)          
!                                                                       
!               get level data                                          
                call calc_rates_level_lte(jkk,lpri,lun11,t,xee,xpx,     &
     &              nlev)
!                                                                       
!               step thru levels                                        
                do mm2=1,nlev 
!                                                                       
!                 get level pointer                                     
                  mmtmp=derivedpointers%npilev(mm2,jkk) 
                  if (mmtmp.ne.0) then 
                    kkkl=mmtmp 
                    mmlv=mmtmp 
!                                                                       
!                   test for level pop                                  
                    if (xilev(kkkl).gt.1.d-34) then 
!                                                                       
!                     get data                                          
                      eth=leveltemp%rlev(1,mm2) 
                      nions=nions+1 
                      mllev(nions)=masterdata%idat1(np1i+nidt-2) 
!                     Note that rwrk1 must be written to the file before
!                     it is overwritten in subsequent columns           
                      rwrk1(nions)=sngl(xilev(mmlv))
                      rwrk2(nions)=sngl(rnist(mmlv))
                      elev(nions)=sngl(eth)
                      ntptr(nions)=kkkl 
                      natomic(nions)=nnz 
                      nupper(nions)=mm2 
                      do mm=1,nkdti 
                         write (kdtmp(mm:mm),'(a1)')                     &
     &                         masterdata%kdat1(np1ki-1+mm) 
                        enddo 
                      do mm=nkdti+1,9 
                        write (kdtmp(mm:mm),'(a1)')' ' 
                        enddo 
                      kion(nions)=kdtmp 
                      write(klevt(nions),'(20a1)')                      &
     &                     (leveltemp%klev(mm,mm2),mm=1,20)                    
                      if (lpri.ne.0) then 
                        write (lun11,*)nions,xilev(mmlv),               &
     &                       masterdata%rdat1(np1r),nnz,mmlv,kkkl            
                        write (lun11,9296)kkkl,                         &
     &                      (masterdata%kdat1(np1i-1+mm),mm=1,20),      &
     &                      (leveltemp%klev(lk,mm2),lk=1,20),           &
     &                      eth,xilev(kkkl),rnist(kkkl)
 9296                   format (1x,i6,1x,(40a1),7(1pe13.5)) 
                        endif 
!                                                                       
!                     end of test for level pop                         
                      endif 
!                                                                       
!                   end of test for level pointer                       
                    endif 
!                                                                       
!                 end of step thru levels                               
                  enddo 
!                                                                       
!               end of test for element                                 
                endif 
!                                                                       
!             Go to next ion                                            
              mlion=derivedpointers%npnxt(mlion) 
              enddo 
!                                                                       
!           end of test for abundance                                   
            endif 
!                                                                       
        mlel=derivedpointers%npnxt(mlel) 
!       Go to next element                                              
        enddo 
                                                                        
!                                                                       
                                                                        
!     End of atomic database extraction                                 
!----------------------------------------------------------------       
!     define parameters for the binary table (see the above data stateme
      nrows=nions 
      tfields=9 
!     Build extension name                                              
      extname='XSTAR_RADIAL' 
                                                                        
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Write table headers'                      
!     write the required header parameters for the binary table         
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,         &
     &              varidat,status)                                     
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Add some more keywords'                   
                                                                        
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
                                                                        
      rtmp=sngl(pres)
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
                                                                        
!-------------------------------------------------------------------    
!     Step through the columns and write them to the file               
!                                                                       
!     set 'global' parameters for writing FITS columns                  
      frow=1 
      felem=1 
                                                                        
!     column  1  (Line number)                                          
      colnum=1 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpclj(unit,colnum,frow,felem,nions,ntptr,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  2 (Level number of this ion)                              
      colnum=2 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpclj(unit,colnum,frow,felem,nions,mllev,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  3  (Energy)                                               
      colnum=3 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpcle(unit,colnum,frow,felem,nions,elev,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  4  (Ion)                                                  
      colnum=4 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpcls(unit,colnum,frow,felem,nions,kion,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  5  (Atomic Number)                                        
      colnum=5 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpclj(unit,colnum,frow,felem,nions,natomic,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  6 (Level Designation)                                     
      colnum=6 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpcls(unit,colnum,frow,felem,nions,klevt,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
!     column 7 (Level population)                                       
!     rwrk1 can be safely overwritten after this step                   
                                                                        
      colnum=7 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpcle(unit,colnum,frow,felem,nions,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
      colnum=8 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpcle(unit,colnum,frow,felem,nions,rwrk2,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  9 (upper level index)                                     
      colnum=9 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr: Writing Column ',colnum                   
      call ftpclj(unit,colnum,frow,felem,nions,nupper,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
!     Compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      deallocate(rwrk1)
      deallocate(rwrk2)
      deallocate(elev)
      deallocate(ntptr)
      deallocate(natomic)
      deallocate(mllev)
      deallocate(nupper)
      deallocate(kion)
      deallocate(klevt)

      return 
      end                                           
