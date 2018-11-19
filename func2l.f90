      subroutine func2l(jkk,lpri,lun11,t,xee,xpx,                       &
     &              rniss,rnisse,nlev)
!                                                                       
!     this routine calculates rates affecting level populations         
!     author: T. Kallman                                                
!     note that in this routine rniss indeces are relative to ground
!                                                                       
                                                                        
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      integer np1i,np1r,np1k 
      character(1) kblnk 
!                                                                       
      integer nlevmx,mltype,ml,mllz,nlev,nidt,nrdt,nkdt,lun11,          &
     &        lpri,ltyp,lrtyp,lcon,jkk,mm,lk,mlpar,mlm                  
      real(8) xee,xpx,t,bb 
      real(8) rniss(nd),rnisse(nd) 
!                                                                       
      data kblnk/' '/ 
!                                                                       
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in func2l, inputs:',t,                          &
     &          xee,xpx                                                 
!                                                                       
!     now find level data                                               
!     step thru types                                                   
      nlevmx=0 
      mltype=13 
      ml=derivedpointers%npfi(mltype,jkk) 
      if (lpri.gt.1) write (lun11,*)'jkk=',                             &
     &      jkk,ml,derivedpointers%npar(ml) 
      mllz=derivedpointers%npar(ml) 
!     step thru records of this type                                    
      mlpar=derivedpointers%npar(ml) 
      do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
         mlm=ml-1 
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         nlev=masterdata%idat1(np1i+nidt-2) 
         nlevmx=max(nlevmx,nlev) 
         if ((nlev.gt.0).and.(nlev.le.nd)) then 
           leveltemp%nlpt(nlev)=ml 
           leveltemp%iltp(nlev)=ltyp 
 9101      format (1x,'level quantities:',4i6,4(1pe12.5),3i6,8a1) 
           if (lpri.gt.1) write (lun11,9101)                            &
     &       ml,nlev,ltyp,lrtyp,(masterdata%rdat1(np1r+mm-1),mm=1,4),   &
     &       masterdata%idat1(np1i),masterdata%idat1(np1i+1),           &
     &       masterdata%idat1(np1i+2),                                  &
     &       (masterdata%kdat1(np1k+mm-1),mm=1,8)                    
           do  lk=1,nrdt 
             leveltemp%rlev(lk,nlev)=masterdata%rdat1(np1r+lk-1) 
             enddo 
           do lk=1,nidt 
             leveltemp%ilev(lk,nlev)=masterdata%idat1(np1i+lk-1) 
             enddo 
           do lk=1,nkdt 
             leveltemp%klev(lk,nlev)=masterdata%kdat1(np1k+lk-1) 
             enddo 
           do lk=nkdt+1,100 
             leveltemp%klev(lk,nlev)=kblnk 
             enddo 
           endif 
         ml=derivedpointers%npnxt(ml) 
         if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
         enddo 
      nlev=nlevmx 
      call levwk(rniss,rnisse,bb,lpri,nlev,t,xee,xpx,lun11)      
!                                                                       
      return 
      END                                           
