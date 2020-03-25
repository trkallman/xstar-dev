      subroutine rstepr4(unit,hdunum,radin,radout,rdel,t,prs,        &
     &             xcol,xee,xpx,xi,                                     &
     &             zrems,dpthc,opakc,rccemis,                           &
     &             lun11,lpri,status)                                   
!                                                                       
!     Name: rstepr4.f90
!     Description:
!       read continuum quantities for each radial zone to an individual 
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
!                                                                       
      implicit none 
                                                                        
!                                                                       
      real(4) rtmp 
      real(8) radin, radout,rdel, t, prs, xcol,xee,xpx,xi 
      integer unit,hdunum, nrows, status
!     continuum opacities                                               
      real(8) opakc(ncn) 
      real(8) rccemis(2,ncn)
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
!     continuum emissivities                                            
      real(8) zrems(5,ncn) 
                                                                        
      character(16) ttype(12),tform(12),tunit(12) 
      integer colnum,felem,hdutype
      integer mm, lun11,lpri 
      character(20) kcom 
                                                                        
!     Database manipulation quantities                                  
      integer nelems,nullj,nkeys,irow2,nspace,nhdu,izcol
      real(4) anynull
      character(1) kblnk,nullstr 
                                                                        
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
      status=0
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'in rstepr4 ',hdunum                              
!                                                                       
      call FTGHDN(unit, nhdu) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu                               
      call ftmahd(unit,1,hdutype,status) 
      call FTGHDN(unit, nhdu) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu                               
!                                                                       
!     Move to the appropriate HDU (hdunum) in the file                  
      mm=hdunum 
!                                                                       
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'rstepr4: Moving to extension',mm                  
      call ftmahd(unit,mm,hdutype,status)
      if (status.ne.0) return 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)unit,mm,hdutype,status                            
      if (status .gt. 0)call printerror(lun11,status) 
      call FTGHDN(unit, nhdu) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu                               
                                                                        
!     Determine the number of keywords in the header                    
      nkeys=0 
      call ftghsp(unit,nkeys,nspace,status) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftghsp:',unit,nkeys,nspace,status          
!                                                                       
!     Read each 80-character keyword record, and print it out           
      call ftgkyj(unit,'NAXIS2',nrows,kcom,status) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkyj:',nrows,kcom,status                 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'RINNER',rtmp,kcom,status) 
      radin=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',radin,kcom,status                  
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'ROUTER',rtmp,kcom,status) 
      radout=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',radout,kcom,status                 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'RDEL',rtmp,kcom,status) 
      rdel=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',rdel,kcom,status                   
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'TEMPERAT',rtmp,kcom,status) 
      t=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',t,kcom,status                      
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'PRESSURE',rtmp,kcom,status) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',prs,kcom,status                    
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'COLUMN',rtmp,kcom,status) 
      xcol=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye, xcol=',xcol,kcom,status            
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'XEE',rtmp,kcom,status) 
      xee=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',xee,kcom,status                    
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'DENSITY',rtmp,kcom,status) 
      xpx=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',xpx,kcom,status                    
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'LOGXI',rtmp,kcom,status) 
      xi=rtmp 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'after ftgkye',xi,kcom,status                     
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
      felem=1 
      nelems=1 
      nullstr=' ' 
      nullj=0 
      do irow2=1,nrows 
      if ((lpri.gt.0).and.(irow2.eq.499))                               &
     &   write (lun11,*)'row=',irow2                                    
        do izcol=1,5
          colnum=2+izcol 
          call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,           &
     &                  rtmp,anynull,status)                            
          zrems(izcol,irow2)=rtmp 
          enddo
        colnum=8 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        opakc(irow2)=rtmp 
        if ((lpri.gt.0).and.(irow2.eq.499)) write (lun11,*)rtmp
        colnum=9 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        rccemis(1,irow2)=rtmp 
        colnum=10
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        rccemis(2,irow2)=rtmp 
        colnum=11
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        dpthc(1,irow2)=rtmp 
        colnum=12 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        dpthc(2,irow2)=rtmp 
        enddo 
!                                                                       
                                                                        
!----------------------------------------------------------------       
!     Compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      return 
      end                                           
