      subroutine stpcut(ldirt,lpri,lun11,vturbi,                        &
     &      np2,ncsvn,nlsvn,                                            &
     &      epi,ncn2,opakc,opakcont,oplin,opakab,delr,t,                &
     &      dpthc,dpthcont,tau0,tauc,eliml,elimh)                               
!                                                                       
!     this routine updates.  calculates depths, etc.                    
!     author:  T. Kallman                                               
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
      integer ldirt,lpri,lun11,np2,ncsvn,nlsvn,ncn2
      integer llk,i,lind
      real(8) vturbi,eliml,elimh,delr,t,dpthmx,optpp,dpthtmp,dpthtmpcont
!                                                                       
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'in stpcut:',delr,nlsvn                           
!                                                                       
      do llk=1,ncn2 
        optmp(llk)=0. 
        enddo 
!                                                                       
!     calculate continuum depths                                        
      if (lpri.ne.0) write (lun11,*)'calculating depths in stpcuta' 
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