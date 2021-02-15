      subroutine trnfrn(lpri,lun11,                                  &
     &       nlsvn,ncsvn,ncn2,                                          &
     &       zrems,zremso,elumab,elumabo,elum,elumo)                    
!                                                                       
!     Name: trnfrn.f90  
!     Description:  
!           updates global continua and lines quantities
!
!     List of Parameters:
!     Input:
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!     Output:
!           zrems(4,ncn):  master spectrum array.  (erg/s/erg/10^38)
!           zrems(4,ncn):  old master spectrum array.  (erg/s/erg/10^38)
!           elumab(2,nnml):  rrc luminosities (erg s^-1)/10^38 
!           elumabo(2,nnml):  old rrc luminosities (erg s^-1)/10^38 
!           elum(2,nnnl):  line luminosities (erg/s/10^38)
!           elum(2,nnnl):  old line luminosities (erg/s/10^38)
!
!     Dependencies: none
!     Called by:  xstar
!
!     this routine updates escaping continua and lines                  
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
                                                                        
      real(8) elumab(2,nnml),elumabo(2,nnml) 
      real(8) zrems(5,ncn),zremso(5,ncn) 
      real(8) elum(2,nnnl),elumo(2,nnnl) 
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
