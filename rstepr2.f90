      subroutine rstepr2(unit,hdunum,radin,radout,rdel,t,prs,        &
     &             xcol,xee,xpx,xi,                                     &
     &             rcem,oplin,tau0,                                     &
     &             lun11,lpri,status)                                   
!                                                                       
!                                                                       
!     Name: rstepr2.f90
!     Description:
!       Read line quantities for each radial zone to an individual 
!       extension of the file xoxx_detal2.fits
!       Append a FITS extension binary table containing                   
!       nrhs columns and at most nrhdimj rows                             
!       author: T. Bridgman                                               
!     Parameters:                               
!        Input:                        
!        unit    integer            File unit number                    
!        hdunum  integer            Number of last HDU written          
!        lun11                      logical unit number for printing
!        lpri                       print switch
!        Output:
!        radin   real(8)               inner radius of shell             
!        radout  real(8)               outer radius of shell             
!        delr    real(8)               thickness of shell
!        temp    real(8)               temperature of shell in 10^4K          
!        pres    real(8)               pressure in shell                 
!        rcem(2,nrhdimj)              line emissivities
!        oplin(nrhdimj)             line opacities
!        tau0(2,nrhdimj)            line optical depths
!        status  integer            Returned status code                
!     Dependencies:  none
!     called by:  xstar       
!                                                                       
      use globaldata
!                                                                       
      implicit none 
!                                                                       
!     Allocation for passed parameters                                  
      real(8) tau0(2,nnnl), rcem(2,nnnl) 
      real(4) rtmp 
      real(8) radin, radout,rdel, t, prs, xcol,xee,xpx,xi 
      integer unit,hdunum, nrows, status
!     line opacities                                                    
      real(8) oplin(nnnl) 
                                                                        
!     Internal work areas                                               
      integer ntptr 
      character(16) ttype(5),tform(5),tunit(5) 
      integer colnum,felem,hdutype
      integer mm, lun11, lpri 
      character(20) kcom 
      integer nhdu 
                                                                        
!     Database manipulation quantities                                  
      integer nelems,nullj,nkeys,irow2,nspace 
      real(4) anynull
      character(1) kblnk,nullstr 
                                                                        
      data kblnk/' '/ 
!                                                                       
      data tform/'1J','1E','1E','1E','1E'/ 
                                                                        
      data ttype/'index','energy','opacity','fwd dpth',                 &
     & 'bck dpth'/                                                      
                                                                        
      data tunit/' ','ev','/cm',' ',' '/ 
                                                                        
       save kblnk,tform,ttype,tunit
!                                                                       
!                                                                       
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'in rstepr2 ',hdunum                              
                                                                        
      call FTGHDN(unit, nhdu) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu                               
      call ftmahd(unit,1,hdutype,status) 
      call FTGHDN(unit, nhdu) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu,unit                               
!                                                                       
!     Move to the appropriate HDU (hdunum) in the file                  
      mm=hdunum 
!                                                                       
      if (lpri.gt.0)                                                    &
     & write(lun11,*)'rstepr2: Moving to extension',mm                  
      call ftmahd(unit,mm,hdutype,status) 
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
!                                                                       
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
      if (lpri.gt.0)                                                    &
     &   write (lun11,*)'row=',irow2                                    
        colnum=1 
        call ftgcvj(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  ntptr,anynull,status)                           
        colnum=6 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        rcem(1,ntptr)=rtmp 
        colnum=7 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        rcem(2,ntptr)=rtmp 
        colnum=8 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        oplin(ntptr)=rtmp 
        colnum=9 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        tau0(1,ntptr)=rtmp 
        colnum=10 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        tau0(2,ntptr)=rtmp 
        enddo 
!                                                                       
                                                                        
!----------------------------------------------------------------       
!     Compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      return 
      end                                           
