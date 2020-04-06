      subroutine fstepr4(unit,hdunum,radin,radout,rdel,t,prs,        &
     &             xcol,xee,xpx,xi,                                     &
     &             epi,ncn2,zrems,dpthc,opakc,rccemis,                  &
     &             lun11,lpri,status)                                   
!                                                                       
!                                                                       
!     Name: fstepr4.f90
!     Description:
!       Write continuum quantities for each radial zone to an individual 
!       extension of the file xoxx_detal4.fits
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
!        epi(ncn)                   energy grid (eV)
!        ncn2                       number of energy points
!        zrems(5,ncn)               radiation field
!        dpthc(2,ncn)               continuum optical depths
!        opakc(ncn)                 continuum opacity
!        rccemis                    continuum emissivity
!        lun11                      logical unit number for printing
!        lpri                       print switch
!        Output:
!        status  integer            Returned status code                
!     Dependencies:  none
!     called by:  xstar       
!                                                                       
      use globaldata
      implicit none 
                                                                        
      integer nptmpdim 
      parameter (nptmpdim=ncn) 
!                                                                       
!     Allocation for passed parameters                                  
      real(4) rtmp 
      real(8) radin, radout,rdel, t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn) 
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn) 
      real(8) zrems(5,ncn) 
      integer ncn2 
                                                                        
      integer, dimension(:), allocatable :: ntptr
      real(4), dimension(:), allocatable :: rwrk1
      integer tfields,varidat 
      character(16) ttype(12),tform(12),tunit(12) 
      integer colnum,frow,felem,hdutype,izcol
      integer mm, lun11, lpri 
      integer nlines 
      character(33) extname 
                                                                        
!     Database manipulation quantities                                  
      character(1) kblnk
                                                                        
      data kblnk/' '/ 
!                                                                       
      data tform/'1J','1E','1E','1E','1E','1E','1E','1E','1E','1E',     &
     &'1E','1E'/
                                                                        
      data ttype/'index','energy','zrems(1)','zrems(2)','zrems(3)',     &
     &  'zrems(4)','zrems(5)','opacity','emis out','emis in',           &
     &  'fwd dpth','bck dpth'/                                                                                                                              
      data tunit/' ','ev','erg/s','erg/s','erg/s','erg/s',              &
     &   'erg/s','/cm','erg/cm**3/s','erg/cm**3/s',' ',' '/
!                                                                        
       save kblnk,tform,ttype,tunit
!
      allocate(ntptr(nptmpdim))
      allocate(rwrk1(nptmpdim))
!
      varidat=0 
!                                                                       
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'in fstepr4 input hdu',hdunum                     
!                                                                       
!     Move to the last HDU (hdunum) in the file                         
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Moving to end-of-FITS file'              
      call ftmahd(unit,hdunum,hdutype,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     append a new empty extension after the last HDU                   
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr4: Create the new extension'               
      call ftcrhd(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
!                                                                       
!     Extracting data from the Atomic Database here                     
!                                                                       
!                                                                       
                                                                        
!     End of atomic database extraction                                 
!----------------------------------------------------------------       
!     define parameters for the binary table (see the above data stateme
      nrows=ncn2 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'before header write'                             
      tfields=12
!     Build extension name                                              
      extname='XSTAR_RADIAL' 
                                                                        
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr4: Write table headers'                    
!     write the required header parameters for the binary table         
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,         &
     &              varidat,status)                                     
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'fstepr4: Add some more keywords'                 
                                                                        
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
                                                                        
      do mm=1,ncn2 
        ntptr(mm)=mm 
        enddo 
                                                                        
!     column  1  (continuum index)                                      
      colnum=1 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Writing Column ',colnum                  
      nlines=ncn2 
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      do mm=1,ncn2 
        rwrk1(mm)=sngl(epi(mm))
        enddo 
!                                                                       
!     column  2 energy                                                  
      colnum=2 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
!                                                                       
!     zrems
      do izcol=1,5
        do mm=1,ncn2 
          rwrk1(mm)=sngl(zrems(izcol,mm))
          enddo 
        colnum=2+izcol 
        call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
        enddo
!
!     column  8 opacity                                                 
      do mm=1,ncn2 
        rwrk1(mm)=sngl(opakc(mm))
        enddo 
!                                                                       
      colnum=8 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
!     column  9 emiss outward
      do mm=1,ncn2 
        rwrk1(mm)=sngl(rccemis(1,mm)) 
        enddo 
!                                                                       
      colnum=9 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
!     column  10 emiss inward
      do mm=1,ncn2 
        rwrk1(mm)=sngl(rccemis(2,mm))
        enddo 
!                                                                       
      colnum=10
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
!                                                                       
      do mm=1,ncn2 
        rwrk1(mm)=sngl(dpthc(1,mm))
        enddo 
!                                                                       
!     column  11 depth forward                                           
      colnum=11 
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
!                                                                       
      do mm=1,ncn2 
        rwrk1(mm)=sngl(dpthc(2,mm)) 
        enddo 
!                                                                       
!     column  12 depth backward                                          
      colnum=12
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'fstepr4: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
                                                                        
                                                                        
                                                                        
!----------------------------------------------------------------       
!     Compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!
      deallocate(ntptr)
      deallocate(rwrk1)
!                                                                        
      return 
      end                                           
