      subroutine init(lunlog,bremsa,bremsint,tau0,dpthc,dpthcont,tauc,  &
     &   xii,rrrt,pirt,htt,cll,httot,cltot,                             &
     &   cllines,clcont,htcomp,clcomp,clbrems,                          &
     &   xilev,rcem,oplin,rccemis,brcems,opakc,opakscatt,               &
     &   cemab,cabab,opakab,elumab,elumabo,elum,elumo,                  &
     &   zrems,zremso,fline,flinel)                                     
!                                                                       
!     this routine initializes everything                               
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
!     line luminosities                                                 
      real(8) elum(3,nnnl),elumo(3,nnnl) 
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
      real(8) fline(2,nnnl),flinel(ncn) 
!     continuum lum                                                     
      real(8) zrems(5,ncn)
      real(8) zremso(5,ncn)
!     continuum optical depths                                          
      real(8) dpthc(2,ncn),dpthcont(2,ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakscatt(ncn) 
!     level populations                                                 
      real(8) xilev(nnml) 
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) elumab(2,nnml),elumabo(2,nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating/cooling                                                   
      real(8) htt(nni),cll(nni) 
      real(8) rrrt(nni),pirt(nni) 
      real(8) httot,cltot,cllines,clcont,htcomp,clcomp,clbrems 
      integer i,lunlog 
!!                                                                      
!                                                                       
       do i = 1,nnml 
         elumab(1,i)=0. 
         elumab(2,i)=0. 
         elumabo(1,i)=0. 
         elumabo(2,i)=0. 
         cabab(i)=0. 
         cemab(1,i)=0. 
         cemab(2,i)=0. 
         opakab(i)=0. 
!         xilev(i)=1.                                                   
         xilev(i)=0. 
         tauc(1,i) = 0. 
         tauc(2,i) =0. 
         enddo 
!                                                                       
      httot=0. 
      cltot=0. 
      cllines=0. 
      clcont=0. 
      htcomp=0. 
      clcomp=0. 
      clbrems=0. 
!                                                                       
      do i = 1,ncn 
         rccemis(1,i)=0. 
         rccemis(2,i)=0. 
         brcems(i)=0. 
         flinel(i)=0. 
         zrems(1,i)=0. 
         zrems(2,i)=0. 
         zrems(3,i)=0. 
         zrems(4,i)=0. 
         zremso(1,i)=0. 
         zremso(2,i)=0. 
         zremso(3,i)=0. 
         zremso(4,i)=0. 
         bremsint(i)=0. 
         bremsint(i)=0. 
         bremsa(i)=0. 
         dpthc(1,i) = 0. 
         dpthc(2,i)=0. 
         dpthcont(1,i) = 0.
         dpthcont(2,i)=0.
!         dpthc(2,i)=1.e+10                                             
         opakc(i)=0. 
         opakscatt(i)=0. 
         enddo 
       do  i = 1,nnnl 
         fline(1,i)=0. 
         fline(2,i)=0. 
         rcem(1,i)=0. 
         rcem(2,i)=0. 
         elum(1,i)=0. 
         elum(2,i)=0. 
         elum(3,i)=0. 
         elumo(1,i)=0. 
         elumo(2,i)=0. 
         elumo(3,i)=0. 
!         tau0(2,i)=1.e+20                                              
!         tau0(1,i) = 1.e+20                                            
         tau0(2,i)=0. 
         tau0(1,i) = 0. 
         oplin(i)=0. 
         enddo 
                                                                        
!      write (lunlog,*)'NB no backward escape'                          
!      write (tmpst,*)'NB no backward escape'                           
!      call xwrite(tmpst,10)                                            
      do i=1,nni 
         xii(i)=0. 
         htt(i)=0. 
         cll(i)=0. 
         rrrt(i)=0. 
         pirt(i)=0. 
         enddo 
!                                                                       
!                                                                       
      return 
      END                                           
