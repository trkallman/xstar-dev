      subroutine fwrtascii(unit,extname,rdati,ncol,                  &
     &                      nidat1j,klabs, kform, kunits,lun11)         
!                                                                        
!     Name: fwrtascii.f90  
!     Description:  
!       write an ascii table extension containing                         
!       ncol columns and nidat1 rows                                      
!       modifications:                                                    
!        1998/12/17, wtb: fix fits keyword format problem.  enhanced    
!                    parameter list for more flexibility.               
!        1999/01/04, wtb: added file creation date, model name, creator 
!                    code & checksum                                    
!        1999/01/25, wtb: convert this routine so it just writes an     
!                    ascii table extension.                             
!       author: T. Bridgman                                               
!     parameters:                                                       
!        unit    integer            file unit number                    
!        extname char*30            name of the ascii extension         
!        rdati   real(ncol*nidat1j)  data array                         
!        ncol    integer            number of columns                   
!        nrhdim  integer            maximum number of rows & columns    
!        nidat1j  integer            actual number of rows              
!        klabs   char*16(ncol)      column labels                       
!        kform   char*16(ncol)      column numeric format               
!        kunits  char*15(ncol)      column units                        
!     Dependencies:  none
!     Called by: pprint.f90
!                                                                       
!                                                                       

      implicit none 
      integer nrhmx,nrhmx1 
                                                                        
      parameter (nrhmx1=999) 
      parameter (nrhmx=3999) 
                                                                        
!     passed parameters                                                 
      real(8) rdati(nrhmx1,*) 
      real(4) rdat(nrhmx) 
      character(16) klabs(*), kform(*), kunits(*) 
      integer ncol, nidat1j 
!      character(30) extname  !jg                                       
      character(10) extname 
                                                                        
      integer, dimension(:), allocatable :: tbcol
      integer unit, status, tfields, nrows, rowlen, verbose,lun11 
      integer felem,frow,colnum,kk,ll 
!                                                                       
      allocate(tbcol(nrhmx))
!
      status=0 
      verbose=0 
      tfields=ncol 
      nrows=nidat1j 
      rowlen=0 
      tbcol(1)=0 
!     append a new empty extension onto the end of the primary array    
      call ftcrhd(unit,status) 
                                                                        
      if(verbose.gt.0) write(6,*)'fwrtascii: writing header table' 
!     write the required header parameters for the ascii table          
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,   &
     &            extname,status)                                       
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
!     map each column to a 1-d array before writing to the file         
      do kk=1,tfields 
        if(verbose.gt.0) write(6,*)'fwrtascii: building column ',kk 
        frow=1 
        felem=1 
        colnum=kk 
        do ll=1,nidat1j 
          rdat(ll)=sngl(rdati(kk,ll)) 
          enddo 
        if(verbose.gt.0) write(6,*)'fwrtascii: writing column ',kk 
        call ftpcle(unit,colnum,frow,felem,nrows,rdat,status) 
        enddo 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     compute checksums                                                 
      if(verbose.gt.0) write(6,*)'fwrtascii: writing checksum' 
      call ftpcks(unit,status) 
!     check for any error, and if so print out error messages           
      if (status .gt. 0)call printerror(lun11,status) 
!
      deallocate(tbcol)
!                                                                       
      END                                           
