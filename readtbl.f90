      subroutine readtbl(np1r,np1i,np1k,                             &
     &       np2,                                                       &
     &       filename,credate,lpri,lun11)                           
!                                                                       
!     Name:  readtbl.f90
!     Description:
!       reads in atomic data                                              
!       written by Ke Zhang, Nov. 9, 2001                               
!     Parameters:
!       Input:
!          filename=file name
!          lpri=print switch
!          lun11=logical unit number
!       Output:
!          credate=creation date
!          np1r=number of reals
!          np1i=number of integers
!          np1k=number of characters
!          np2=number of records
!    Dependencies:  none
!    Called by: xstarsetup
!                                                                       
       use globaldata
       implicit none 
!                                                                       
      integer nidatt 
      parameter (nidatt=nrdat1) 
!                                                                       
      integer lpri,lun11,np1r,np1i,np1k,np2 
      real(4), dimension(:), allocatable ::  rdat14
      integer, dimension(:), allocatable :: ntptr
      integer mm 
                                                                        
      real(4) nulle 
      integer nullj 
      character filename*256,nullstr*1 
      character credate*63, comment*50 
      logical anynull 
                                                                        
      integer status,unit,readwrite,blocksize,hdutype 
      integer row,col,j,i,lenact 
                                                                        
      allocate(rdat14(nidatt))
      allocate(ntptr(nidatt))

      status=0 
      nullstr=' ' 
      nullj=0 
      nulle=0. 
              
      if (lpri.ne.0) write (lun11,*)'in readtbl' 
!                                                                       
! get an unused logical unit number and open the fits file              
!      call ftgiou(unit,status)                                         
      call getlunx(unit) 
      readwrite=0 
      call ftopen(unit,filename,readwrite,blocksize,status) 
                                                                        
! Read the primary header & get the file creation date                  
      call ftmahd(unit,1,hdutype,status) 
      call ftgkys(unit,'DATE',credate,comment,status) 
      if(status .gt.0) call printerror(lun11,status) 
      write(lun11,*)'Atomic Data Version: ',credate(1:lenact(credate)) 
                                                                        
      col=1 
      row=1 
                                                                        
! move to the next extension : POINTERS                                 
      call ftmrhd(unit,1,hdutype,status) 
                                                                        
! read LENGTH keywords: total records #                                 
      call ftgkyj(unit,'LENGTH',np2,comment,status) 
                                                                        
! read POINTERS data                                                    
      call ftgcvj(unit,col,row,1,10*np2,nullj,                          &
     &               ntptr,anynull,status)                              
      do i=1,np2 
        do j=1,10 
          masterdata%nptrs(j,i)=ntptr(10*(i-1)+j) 
        enddo 
      enddo 
      write (lun11,*)'in readtbl:' 
      write (lun11,*)'number of pointers=',np2 
                                                                        
                                                                        
! move to the next extension : REALS                                    
      call ftmrhd(unit,1,hdutype,status) 
                                                                        
! read LENGTH keywords: total reals #                                   
      call ftgkyj(unit,'LENGTH',np1r,comment,status) 
      write (lun11,*)'number of reals=',np1r 
                                                                        
! read REAL(4) data                                                      
      call ftgcve(unit,col,row,1,np1r,nulle,                            &
     &               rdat14,anynull,status)                             
      do mm=1,nrdat1 
        masterdata%rdat1(mm)=rdat14(mm) 
        enddo 
                                                                        
                                                                        
! move to the next extension : INTEGERS                                 
      call ftmrhd(unit,1,hdutype,status) 
                                                                        
! read LENGTH keywords: total integers #                                
      call ftgkyj(unit,'LENGTH',np1i,comment,status) 
      write (lun11,*)'number of integers=',np1i 
                                                                        
! read INTEGER data                                                     
      call ftgcvj(unit,col,row,1,np1i,nullj,                            &
     &             masterdata%idat1,anynull,status)                       
                                                                        
                                                                        
! move to the next extension : CHARS                                    
      call ftmrhd(unit,1,hdutype,status) 
                                                                        
! read LENGTH keywords: total chars #                                   
      call ftgkyj(unit,'LENGTH',np1k,comment,status) 
      write (lun11,*)'number of characters=',np1k 
                                                                        
! read CHAR data                                                        
      call ftgcvb(unit,col,row,1,np1k,nullstr,                          &
     &             masterdata%kdat1,anynull,status)                    
                                                                        
                                                                        
! close the file and free the unit number                               
      call ftclos(unit, status) 
!      call ftfiou(unit, status)                                        
      close(unit) 
                                                                        
! check for any error, and if so print out error messages               
      if (status .gt. 0) call printerror(lun11,status) 
                                                                        
      deallocate(rdat14)
      deallocate(ntptr)

      return 
      END                                           
