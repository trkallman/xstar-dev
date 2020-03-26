      subroutine drd(ltyp,lrtyp,lcon,lrdat,np1r,lidat,np1i,lkdat,np1k,  &
     &                 np2,lpri,lun11)                            
!                                                                       
!     Name: drd.f90  
!     Description:  
!     this  routine gets pointers for one record of the database.
!     author:  T. Kallman                                               
!     List of Parameters:
!           Input:
!           np2: record of database
!           lpri= print switch
!           lun11=logical unit number for printout
!           Output:
!           ltyp=data type
!           lrtyp=rate type
!           lcon=continuation switch
!           lrdat=number of reals
!           np1r=pointer to first element of master real array to be printed
!           lidat=number of integers
!           np1i=pointer to first element of master integer array to be printed
!           lkdat=number of chars
!           np1k=pointer to first element of master char array to be printed
!     Dependencies:  none
!     Called by:  ucalc, pprint, setptrs, writespectra2, deleafnd
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
!                                                                       
      integer nrd,lcon,ltyp,lrtyp,lrdat,lidat,                          &
     &        lkdat,np2,lpri,lun11,np1,np1r,np1i,np1k
!                                                                       
      if ( lpri.gt.0 ) write (lun11,*) 'in drd, np2=' , np2,            &
     &                                  ntyp                            
!      if ((ltyp.le.0).or.(ltyp.gt.ntyp))                               
!     $    stop 'data typing error'                                     
      nrd = 0 
      lcon=1 
!        this is an old convention, now deprecated
!        np2 = np2 + 1 
        nrd = nrd + 1 
        np1 = masterdata%nptrs(1,np2) 
        ltyp = masterdata%nptrs(2,np2) 
        lrdat = masterdata%nptrs(5,np2) 
        lidat = masterdata%nptrs(6,np2) 
        lkdat = masterdata%nptrs(7,np2) 
        lrtyp = masterdata%nptrs(3,np2) 
        lcon = masterdata%nptrs(4,np2) 
        np1r = masterdata%nptrs(8,np2) 
        np1i = masterdata%nptrs(9,np2) 
        np1k = masterdata%nptrs(10,np2) 
        if ( lpri.gt.0 ) write (lun11,*) 'in  drd:' , np2 , np1, ltyp,  &
     &                                 lrtyp , lrdat , lidat            
        if ( lpri.gt.0 ) write (lun11,99001) lkdat , lcon , np1r ,np1i, &
     &                        np1k                                      
      lcon = 0 
      if ( lpri.gt.0 ) write (lun11,*) 'leaving drd' , np2 
!                                                                       
      return 
99001 format (8x,5i8) 
      END                                           
