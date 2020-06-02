      subroutine fheader(unit,knam,atcredate,mdlname,status) 
                                                                        
!                                                                       
!   File Name:    fheader.f                                             
!   Author:       W.T. Bridgman                                         
!   Date:         January 1999                                          
!   Abstract:     Routines for writing a stardard primary FITS          
!                 header for XSTAR output files.                        
!                                                                       
!     Create a FITS file with an empty primary header                   
!     nrhs columns and nrows row                                        
!                                                                       
!     Parameters:                                                       
!        unit    integer            File unit number                    
!        knam    char*16            File name to create                 
!        mdlname char*30            Model name for this run             
!        status  integer            Returned status code                
!                                                                       

      implicit none 
!                                                                       
      character(16) knam, filename 
      character(30) mdlname 
!     the atomic data creation date                                     
      character(63) atcredate 
      integer unit, status 
      integer lun11 
                                                                        
      integer bitpix,naxis,naxes(2),group 
      logical simple,extend 
      integer blocksize 
                                                                        
      status=0 
!                                                                       
      filename=knam 
                                                                        
!     Delete the file if it already exists, so we can recreate it       
      call deletefile(filename,status) 
!      write (6,*)'in fheader after deletefile',filename,status         
                                                                        
!     Get an unused Logical Unit Number to use to open the FITS file    
!      call ftgiou(unit,status)                                         
      call getlunx(unit) 
!      write (6,*)'in fheader after getlun',unit,status                 
                                                                        
!     open the FITS file, with write access                             
      lun11=6 
!      write (lun11,*)'opening fits unit ',unit, filename               
      blocksize=1 
      call ftinit(unit,filename,blocksize,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     try to open the file, to see if it exists                         
!      call ftopen(unit,filename,1,blocksize,status) 
!      write (lun11,*)'after ftopen',unit,status                            
!      if (status .gt. 0)stop                                           
                                                                        
                                                                        
!     initialize parameters for primary array                           
      simple=.true. 
      bitpix=16 
      naxis=0 
      naxes(1)=0 
      naxes(2)=0 
      extend=.true. 
      group=1 
!     write the required primary header keywords                        
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,group,extend,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     now add additional keywords                                       
      call ftpcom(unit,'***********************************',status) 
      call ftpkys(unit,'CREATOR','XSTAR version 2.56e',                 &
     & 'Program which generated this file',status)                      
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     Extract the system date                                           
      call ftpdat(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     Save run-specific information                                     
      call ftpkys(unit,'MODEL',mdlname,'Model name for this run',status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftpkys(unit,'ATDATA',atcredate,                              &
     &  'Atomic data creation date',status)                             
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
      return 
      END                                           
