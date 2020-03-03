      subroutine stpcut(ldirt,lpri,lun11,                            &
     &      ncsvn,nlsvn,                                                &
     &      epi,ncn2,opakc,opakcont,oplin,opakab,delr,                  &
     &      dpthc,dpthcont,tau0,tauc)                               
!                                        
!     Name: stpcut.f90
!     Description:                               
!       this routine updates optical depths
!       author:  T. Kallman                                               
!     Paramters:
!       Input:
!       ldirt=direction switch
!       lpri=print switch
!       lun11=logical unit number for printing
!       ncsvn=number of rrcsin atomic database
!       nlsvn=number of rrcsin atomic database
!       epi(ncn)=continum bins (eV)
!       ncn2=length of epi
!       opakc(ncn)=continuum opacity (cm^-1)
!       opakcont(ncn)=continuum opacity, continuum only(cm^-1)
!       oplin(nnnl)=line opacities (cm^-1)
!       opakab(nnml)=rrc opacities at threshold (cm^-1)
!       delr=step size (cm)
!       t=temperature/10^4K
!       Output:
!       dpthc(2,ncn)=continuum optical depths
!       dpthcont(2,ncn)=continuum optical depths, continuum only
!       tau0(2,nnnl)=line optical depths
!       tauc(2,nnml)=rrc optical depths
!     Dependencies: none
!     called by:  xstar
!                                                                       
      implicit none 
!                                                                       
      include './PARAM' 
      integer nbtpp 
      parameter (nbtpp=10000) 
!                                                                       
      real(8) epi(ncn),opakc(ncn),dpthc(2,ncn),dpthcont(2,ncn)
      real(8) tauc(2,nnml),opakcont(ncn)
      real(8) opakab(nnml) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
      real(8) optmp(ncn) 
      integer ldirt,lpri,lun11,ncsvn,nlsvn,ncn2
      integer llk,i,lind
      real(8) delr,dpthmx,optpp,dpthtmp,dpthtmpcont
!                                                                       
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'in stpcut:',delr,nlsvn                           
!                                                                       
      do llk=1,ncn2 
        optmp(llk)=0. 
        enddo 
!                                                                       
!     calculate continuum depths                                        
      if (lpri.ne.0) write (lun11,*)'calculating depths in stpcut' 
      lind=1 
      if (ldirt.gt.0) lind=2 
      dpthmx=0. 
      optpp=0. 
      do  i = 1,ncn2 
!         opakc(i)=opakc(i)+optmp(i)                                    
         optpp=min(optmp(i),1.e+3/(1.e-24+delr)) 
         dpthtmp=(opakc(i)+optpp)*delr 
         dpthtmpcont=opakcont(i)*delr
         if (lpri.ne.0) write (lun11,*)i,epi(i),opakc(i),optmp(i),      &
     &        dpthc(lind,i),dpthtmp                                     
         dpthc(lind,i) = dpthc(lind,i) + dpthtmp 
         dpthcont(lind,i) = dpthcont(lind,i) + dpthtmpcont
         if (dpthtmp.gt.dpthmx) then 
           dpthmx=dpthtmp 
           endif 
         enddo 
!                                                                       
!     calculate line depths                                             
      do  i = 1,nlsvn 
         tau0(lind,i) = tau0(lind,i) + oplin(i)*delr 
         if (lpri.ne.0)                                                 &
     &    write (lun11,*)i,lind,oplin(i),delr,tau0(lind,i)              
         enddo 
!                                                                       
!     calculate level depths                                            
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'in stpcut:',delr,lind                            
      do i = 1,ncsvn 
         tauc(lind,i) = tauc(lind,i) + opakab(i)*delr 
         if (lpri.ne.0)                                                 &
     &    write (lun11,*)i,opakab(i),tauc(lind,i)                       
         enddo 
!                                                                       
!                                                                       
      return 
      END                                           
