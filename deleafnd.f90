      subroutine deleafnd(jkk,lup,                                   &
     &   delea,lfnd,lpri,lun11)                           
!                                                                       
!     Name: deleafnd.f90  
!     Description:  
!       finds the damping parameter for a transition
!     Parameters:
!           Input:
!           jkk=ion number
!           lup=level index for upper level
!           lpri=print switch
!           lun11=logical unit number
!           Output:
!           delea=damping parameter in s^-1
!           lfnd=found flag
!      Dependencies: none
!      Called by:  ucalc
!
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      integer jkk,mllz,lup,lfnd,iltmp,lcon,ltyp,lrtyp,                  &
     &   nrdt,np1r,nidt,np1i,nkdt,np1k,nilin,                           &
     &   lpri,lun11,nitmp,iion,mlm,ndtmp                    
      real(8) delea 
!                                                                       
!
!     find associated type 86 data                                      
!                                                                       
!     this is not needed                                                
      go to 9092 
      iion=0 
!      nilin=npar(ml)                                                   
      nilin=derivedpointers%npfirst(13) 
      mlm=nilin 
      call drd(ltyp,lrtyp,lcon,                                         &
     &           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                     &
     &           0,lun11)                                         
      if (lpri.ge.1)                                                    &
     &   write (lun11,*)'searching for ion',jkk,mlm,                    &
     &   masterdata%idat1(np1i+nidt-1)                                       
       do while ((masterdata%idat1(np1i+nidt-1).ne.jkk)                 &
     &     .and.(iion.lt.nni))                                          
        iion=iion+1 
        nitmp=derivedpointers%npfi(13,iion) 
        mlm=nitmp 
        call drd(ltyp,lrtyp,lcon,                                       &
     &           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                     &
     &           0,lun11)                                         
        if (lpri.ge.1)                                                  &
     &    write (lun11,*)iion,masterdata%idat1(np1i+nidt-1),            &
     &           nitmp                                                  
        enddo 
 9092   continue 
      iion=jkk 
!     get damping parameter for iron auger damped lines                 
      ndtmp=derivedpointers%npfi(41,iion) 
      if (lpri.gt.1) write (lun11,*)'ndtmp=',iion,ndtmp 
      if (ndtmp.eq.0) then 
          lfnd=0 
          return 
        else 
          if (lpri.gt.1)                                                &
     &    write (lun11,*)'  found ion',lup,ndtmp                        
          mllz=derivedpointers%npar(ndtmp) 
          mlm=ndtmp 
          call drd(ltyp,lrtyp,lcon,                                     &
     &           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                     &
     &           0,lun11)                                         
          iltmp=masterdata%idat1(np1i+1) 
          do while ((ndtmp.ne.0).and.(lup.ne.iltmp)                     &
     &            .and.(derivedpointers%npar(ndtmp).eq.mllz))                  
            mlm=ndtmp 
            call drd(ltyp,lrtyp,lcon,                                   &
     &             nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                   &
     &             0,lun11)                                       
            if (lpri.ge.2)                                              &
     &        call dprinto(ltyp,lrtyp,lcon,                          &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)      
            iltmp=masterdata%idat1(np1i+1) 
            if (lpri.gt.1)                                              &
     &       write (lun11,*)'   ',nidt,                                 &
     &       iltmp,ndtmp,lup                                            
            ndtmp=derivedpointers%npnxt(ndtmp) 
            enddo 
        endif 
      if (lpri.gt.1) write (lun11,*)'lup,iltmp',                        &
     &                     lup,iltmp                                    
      if (lup.eq.iltmp) then 
          lfnd=1 
          delea=masterdata%rdat1(np1r+2)*(4.136e-15) 
          if (lpri.gt.1) write (lun11,*)masterdata%rdat1(np1r+2),delea 
        else 
          lfnd=0 
        endif 
                                                                        
!                                                                       
      return 
      END                                           
