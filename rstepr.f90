      subroutine rstepr(unit,hdunum,radin,radout,rdel,t,prs,            &
     &             xcol,xee,xpx,xi,                                     &
     &             xilev,rnist,                                         &
     &             lun11,lpri,status)                                   
!                                                                       
!     Reads a FITS extension binary table containing                    
!     nrhs columns and at most nrhdimj rows                             
!     author: T. Kallman                                                
!                                                                       
!     Parameters:                                                       
!        unit    integer            File unit number                    
!        hdunum  integer            Number of last HDU written          
!        radin   real(4)               inner radius of shell             
!        radout  real(4)               outer radius of shell             
!                                   nb but now it is delr in the call   
!        t    real(4)               temperature of shell                 
!        prs    real(4)               pressure in shell                  
!        nrhdimj  integer            Maximum number of rows             
!        idat1   integer(nidat1)    Needed by the atomic database       
!        rdat1   real(nidat1)       Needed by the atomic database       
!        kdat1   char*nidat1        Needed by the atomic database       
!        nptrs                      Needed by the atomic database       
!        npnxt                      Needed by the atomic database       
!        npfi                       Needed by the atomic database       
!        npfirst                    Needed by the atomic database       
!        npcon                      Needed by the atomic database       
!        npconi                     Needed by the atomic database       
!        npcon2                     Needed by the atomic database       
!        xilev   real(nrhdimj)       Fractional level population array  
!        cemab   real(2,nrhdimj)     Recombination emission             
!        opakab  real(nrhdimj)       Opacity                            
!        tauc    real(2,nrhdimj)     Optical depth                      
!        poptol  real(4)               Tolerance for population level    
!        nloopctl integer           Loop control variable               
!        nzone   integer            Pass number through iteration proces
!        status  integer            Returned status code                
      use globaldata
!                                                                       
      implicit none 
                                                                        
!     Allocation for passed parameters                                  
      real(8) xilev(nnml),rnist(nnml)     
      real(4) rtmp 
      real(8) radin, radout, rdel,t, prs, xcol,xee,xpx,xi 
      integer unit,hdunum, nrows, status
                                                                        
!     Internal work areas                                               
      integer ntptr 
      character(16) ttype(5),tform(5),tunit(5) 
      integer colnum,felem,hdutype
      integer mm, lun11, lpri 
      character(20) kcom 
                                                                        
!     Database manipulation quantities                                  
      integer nelems,nullj,nkeys,irow2,nspace 
      real(4) anynull
      character(1) kblnk,nullstr 
      integer nhdu 
                                                                        
      data kblnk/' '/ 
!                                                                       
      data tform/'1J','1E','1E','1E','1E'/ 
                                                                        
      data ttype/'index','energy','opacity','fwd dpth',                 &
     & 'bck dpth'/                                                      
                                                                        
      data tunit/' ','ev','/cm',' ',' '/ 
                                                                        
!                                                                       
!                                                                       
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'in rstepr ',hdunum                               
!                                                                       
                                                                        
      call FTGHDN(unit, nhdu) 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu                               
      call ftmahd(unit,1,hdutype,status) 
      call FTGHDN(unit, nhdu) 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu,unit                               
!                                                                       
!     Move to the appropriate HDU (hdunum) in the file                  
      mm=hdunum 
!                                                                       
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'rstepr2: Moving to extension',mm                  
      call ftmahd(unit,mm,hdutype,status) 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)unit,mm,hdutype,status                            
      if (status .gt. 0)call printerror(lun11,status) 
      call FTGHDN(unit, nhdu) 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'current hdu ',nhdu                               
                                                                        
!     Determine the number of keywords in the header                    
      nkeys=0 
      call ftghsp(unit,nkeys,nspace,status) 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftghsp:',unit,nkeys,nspace,status          
!                                                                       
!                                                                       
!     Read each 80-character keyword record, and print it out           
      call ftgkyj(unit,'NAXIS2',nrows,kcom,status) 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkyj:',nrows,kcom,status                 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'RINNER',rtmp,kcom,status) 
      radin=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',radin,kcom,status                  
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'ROUTER',rtmp,kcom,status) 
      radout=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',radout,kcom,status                 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'RDEL',rtmp,kcom,status) 
      rdel=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',rdel,kcom,status                   
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'TEMPERAT',rtmp,kcom,status) 
      t=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',t,kcom,status                      
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'PRESSURE',rtmp,kcom,status) 
      prs=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',prs,kcom,status                    
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'COLUMN',rtmp,kcom,status) 
      xcol=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye, xcol=',xcol,kcom,status            
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftgkye(unit,'XEE',rtmp,kcom,status) 
      xee=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',xee,kcom,status                    
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'DENSITY',rtmp,kcom,status) 
      xpx=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',xpx,kcom,status                    
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      call ftgkye(unit,'LOGXI',rtmp,kcom,status) 
      xi=rtmp 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after ftgkye',xi,kcom,status                     
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
      felem=1 
      nelems=1 
      nullstr=' ' 
      nullj=0 
      do irow2=1,nrows 
        if (lpri.ne.0)                                                  &
     &   write (lun11,*)'row=',irow2                                    
        colnum=1 
        call ftgcvj(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  ntptr,anynull,status)                           
        colnum=7 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
!       note that here we switch the inward and outward                 
        xilev(ntptr)=rtmp 
        colnum=8 
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,             &
     &                  rtmp,anynull,status)                            
        rnist(ntptr)=rtmp 
        enddo 
!                                                                       
                                                                        
!----------------------------------------------------------------       
!     Compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      return 
      end                                           
