      subroutine fparmlist(unit,hdunum,mdlname,npar,parname,partype, &
     &                    parval,parcomm,nloopctl,status,lun11)         
!                                                                       
!     Name:  fparmlist.f90
!     Description:
!        Write the input parameters to a fits file header
!        author: T. Bridgman                                               
!     Parameters:                                                       
!        Input:
!        unit    integer            file unit number                    
!        hdunum  integer            number of last hdu written          
!        mdlname char*30            model name for this run             
!        npar    integer            number of parameters passed         
!        parname char*20(999)       parameter name                      
!        partype char*10(999)       parameter type                      
!        parval  real(999)          parameter values converted to reals 
!        parcomm char*30(999)       parameter comments & string values  
!        nloopctl integer           loop control parameter              
!     Output:
!        status  integer            returned status code                
!     Dependencies:  None (fitsio)
!     Called by: xstar, writespectra,writespectra2,writespectra3,writespectra4
!
      implicit none 
!     passed parameters                                                 
      character(30) mdlname 
      integer unit, status, hdunum, npar, nloopctl 
      character(20) parname(55) 
      character(10) partype(55) 
      real(8) parval(55) 
      real(4) parval4(55) 
      character(30) parcomm(55) 
!     parameter info                                                    
                               !jg                                      
      integer idat1(6000000) 
      integer lun11
      integer tfields,nrows,varidat 
      character(16) ttype(5),tform(5),tunit(5) 
      integer colnum,frow,felem,hdutype 
      integer ll,mm 
      character(30) extname 
      character(4) ktmp2 
!                                                                       
      data tform/'1I','20A','1E','10A','30A'/ 
      data ttype/'index','parameter','value','type','comment'/ 
      data tunit/' ',' ',' ',' ',' '/ 
!                                                                       
      nrows=npar 
      varidat=0 
!                                                                       
      do mm=1,55 
        parval4(mm)=sngl(parval(mm))
        enddo 
!                                                                       
!     move to the last hdu (hdunum) in the file                         
      call ftmahd(unit,hdunum,hdutype,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     append a new empty extension after the last hdu                   
      call ftcrhd(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     define parameters for the binary table (see the above data stateme
      tfields=5 
!                                                                       
!     build extension name                                              
      extname='PARAMETERS' 
      if(nloopctl.gt.0) then 
          write(ktmp2,'(i4.4)')nloopctl 
!          extname='parameters_' // ktmp2                               
          endif 
!                                                                       
!     write the required header parameters for the binary table         
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,         &
     &              varidat,status)                                     
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     save run-specific information                                     
      call ftpkys(unit,'MODEL',mdlname,'model name for this run',status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     set 'global' parameters for writing fits columns                  
      frow=1 
      felem=1 
!                                                                       
!     column  1  (index)                                                
      colnum=1 
      do ll=1,nrows 
         idat1(ll)=ll 
         enddo 
      call ftpclj(unit,colnum,frow,felem,nrows,idat1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     column  2  (parameter name)                                       
      colnum=2 
      call ftpcls(unit,colnum,frow,felem,nrows,parname,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  3  (parameter value)                                      
      colnum=3 
      call ftpcle(unit,colnum,frow,felem,nrows,parval4,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  4 (parameter type)                                        
      colnum=4 
      call ftpcls(unit,colnum,frow,felem,nrows,partype,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  5 (parameter comment)                                     
      colnum=5 
      call ftpcls(unit,colnum,frow,felem,nrows,parcomm,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
!     compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      return 
      END                                           
