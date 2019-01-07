      subroutine dprint(ltyp,lrtyp,lcon,                                &
     &  lrdat,rdat,lidat,idat,lkdat,kdat,                               &
     &  np1r,np1i,np1k,np2,                                             &
     &  lpri,lun11)                             
!                                                                       
!     Name: dprint.f90  
!     Description:  
!     this  routine puts one record of the database into the database arrays
!     author:  T. Kallman                                               
!     List of Parameters:
!           Input:
!           ltyp=data type
!           lrtyp=rate type
!           lcon=continuation switch
!           lrdat=number of reals
!           rdat=array of real data
!           lidat=number of integers
!           idat=array of integer data
!           lkdat=number of characters
!           kdat=array of character data
!           lpri=print switch
!           lun11=logical unit number for printout
!           Output (updated by routine):
!           np1r=pointer to first element of master databas real 
!                array to be filled
!           np1i=pointer to first element of master databas integer
!                array to be filled
!           np1k=pointer to first element of master databas character
!                array to be filled
!           np2=number of records in database
!     Dependencies: none
!     Called by: not called by xstar routines themselves; included for reference
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      integer nptmpdim 
      parameter (nptmpdim=200000) 
!                                                                       
      integer ltyp,lrtyp,lcon,lrdat,lidat,lkdat,                        &
     &  np1r,np1i,np1k,np2,ml,lpri,lun11,lprisv                         
      real(8) rdat(nptmpdim) 
      integer idat(nptmpdim) 
      character(1) kdat(nptmpdim) 
!!                                                                      
      lprisv=lpri 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'in dprint, np2=',np2                             
      if (np2.ge.ndat2) then 
          write (6,*) 'data index error' 
          return 
          endif 
      masterdata%nptrs(1,np2)=1 
      masterdata%nptrs(2,np2)=ltyp 
      masterdata%nptrs(3,np2)=lrtyp 
      masterdata%nptrs(4,np2)=lcon 
      masterdata%nptrs(5,np2)=lrdat 
      masterdata%nptrs(6,np2)=lidat 
      masterdata%nptrs(7,np2)=lkdat 
      masterdata%nptrs(8,np2)=np1r 
      masterdata%nptrs(9,np2)=np1i 
      masterdata%nptrs(10,np2)=np1k 
      if (lpri.ne.0) then 
        write (lun11,*)'in dprint:',np2,ltyp,lrtyp,lrdat,lidat,lkdat 
        write (lun11,*)'          ',lcon,np1r,np1i,np1k 
        endif 
      np2=np2+1 
      if (lrdat.gt.0) then 
        do ml=1,lrdat 
           masterdata%rdat1(np1r)=rdat(ml) 
           np1r=np1r+1 
           enddo 
        endif 
      if (lidat.gt.0) then 
        do  ml=1,lidat 
           masterdata%idat1(np1i)=idat(ml) 
           np1i=np1i+1 
           enddo 
        endif 
      if (lkdat.eq.0) then 
        do ml=1,lkdat 
            masterdata%kdat1(np1k)=kdat(ml) 
            np1k=np1k+1 
            enddo 
        endif 
!                                                                       
      if ((np1k.gt.nkdat1).or.(np1i.gt.nidat1).or.(np1r.gt.nrdat1)      &
     &   .or.(np2.gt.ndat2)) then                                       
        write (lun11,*)'dprint index error,',np1k,np1i,np1r,np2 
        return 
        endif 
!                                                                       
      lpri=lprisv 
!                                                                       
      return 
      end                                           
