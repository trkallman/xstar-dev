      subroutine dprint(ltyp,lrtyp,lcon,                                &
     &  lrdat,rdat,lidat,idat,lkdat,kdat,                               &
     &  np1r,np1i,np1k,np2,                                             &
     &  lpri,lun11)                             
!                                                                       
!     this  routine prints one element of the database                  
!     author:  T. Kallman                                               
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
