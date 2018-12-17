      subroutine xstarcalc(lpri2,lnerrd,nlimdt,                         &
     &       lpri,lprid,lun11,tinf,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       xilev,bilev,rnist,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                           
                                                                        
      use globaldata
      implicit none 
!                                                                       
!     global xstar data
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
      real(8) fline(2,nnnl),flinel(ncn) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     level populations                                                 
      real(8) xilev(nnml),bilev(nnml),rnist(nnml)
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating and cooling                                               
      real(8) htt(nni),cll(nni) 
      real(8) rrrt(nni),pirt(nni) 
!     element abundances                                                
      real(8) ababs(nl) 
!                                                                       
!     state variables                                                   
      real(8) p,r,t,xpx,delr 
!     heating-cooling variables                                         
      real(8) httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,          &
     &     clcont,hmctot                                                
      real(8) trad 
      real(8) cfrac,critf,vturbi,xee,tinf 
      integer lcdd,ncn2 
!     variables associated with thermal equilibrium solution            
      integer ntotit,lnerrd 
!     switches                                                          
      integer lprid,lpri,nlimdt 
!     strings for atomic data read                                      
      integer nlsvn,ncsvn,lun11,np2,lprisv,lpri2 
!                                                                       
!                                                                       
      lprisv=lpri 
      lpri=0 
      if (nlimdt.ne.0) then 
        call dsec(lnerrd,nlimdt,                                        &
     &       lpri,lprid,lun11,tinf,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &         cllines,clcont,htcomp,clcomp,clbrems,                    &
     &       xilev,bilev,rnist,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                        
        endif 
!       do ll=1,nnml                                                    
!         xilev(ll)=0.                                                  
!         enddo                                                         
!                                                                       
      if (lpri2.eq.1) lpri=1 
      call func(lpri,lun11,vturbi,critf,                                &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       xilev,bilev,rnist)
      call funcsyn(lpri,lun11,vturbi,critf,                             &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,                                                       &
     &       xilev,bilev,rnist,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                           
!                                                                       
       lpri=lprisv 
!                                                                       
!                                                                       
      return 
      end                                           
