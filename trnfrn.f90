      subroutine trnfrn(lpri,lun11,                                     &
     &       nlsvn,ncsvn,ncn2,                                          &
     &       zrems,zremso,elumab,elumabo,elum,elumo)                    
!                                                                       
!     this routine updates escaping continua and lines                  
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
                                                                        
      real(8) elumab(2,nnml),elumabo(2,nnml) 
      real(8) zrems(5,ncn),zremso(5,ncn) 
      real(8) elum(3,nnnl),elumo(3,nnnl) 
      integer numcon,ncn2,kl,ll,jkk,jk,lpri,lun11,                      &
     &     ncsvn,nlsvn                                                  
!                                                                       
!     transfer continuum                                                
      if (lpri.ge.1) write (lun11,*)'in trnfrn' 
      numcon=ncn2 
      do kl=1,numcon 
        do ll=1,4 
          zremso(ll,kl) = zrems(ll,kl) 
          enddo 
        enddo 
!                                                                       
!     transfer lines                                                    
      do jkk=1,nlsvn 
        jk=jkk 
        do ll=1,2 
          elumo(ll,jk)=elum(ll,jk) 
          if (lpri.ge.1) write (lun11,*)jk,elum(1,jk),elumo(1,jk) 
          enddo 
        enddo 
!                                                                       
!     transfer RRCs                                                     
      do jkk=1,ncsvn 
        jk=jkk 
        do ll=1,2 
          elumabo(ll,jk)=elumab(ll,jk) 
          if (lpri.ge.1) write (lun11,*)jk,elumab(1,jk),elumabo(1,jk) 
          enddo 
        enddo 
!                                                                       
!                                                                       
      return 
      end                                           
