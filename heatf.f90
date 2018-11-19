      subroutine heatf(jkk,lpri,lun11,                                  &
     &       t,r,cfrac,delr,xee,xpx,abel,                               &
     &       epi,ncn2,bremsa,                                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       nlev,                                                      &
     &       brcems,cmp1,cmp2,httot,cltot,hmctot,                       &
     &             cllines,clcont,htcomp,clcomp,clbrems)                
!                                                                       
!     this routine calculates heating and cooling.                      
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn) 
!     state variables                                                   
      real(8) r,t,xpx,delr,delrl 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) xee 
      integer ncn2,lpri,lun11 
      integer nlsvn,ncsvn 
!     level populations                                                 
      real(8) abel(nl) 
      character(1) kblnk 
      real(8) brcems(ncn)
      real(8) xeltp,cllines,clcont,cmp1,cmp2,cltot,                      &
     &     hmctot,htcomp,clcomp,clbrems,etst,ekt,                       &
     &     fac,optpp,optp2,tautmp,cfrac,                                &
     &     xnx,httot
      real(8) ener,ergsev,tmpscat
      real(8) tmp2,tmp2o 
      integer lskp 
      integer jk,klel,kl,                                               &
     &     klion,mlleltp,mllel,mlel,mlion,mt2,nilin,nnz,numcon,         &
     &     nblin                                                        
      integer np1i,np1r,np1k,np1ki 
      integer nlev,nlevmx,mltype,ml,mllz,jkk,ltyp,                      &
     &     lrtyp,lcon,nrdt,nkdt,lk,kkkl,                                &
     &     lprisv,mm,nbinc,mlpar,mlm,np2
!                                                                       
!                                                                       
      data kblnk/' '/ 
      data ergsev/1.602197e-12/ 
!                                                                       
      lprisv=lpri 
!                                                                       
      xnx=xpx*xee 
      if (lpri.ge.1) lpri=2 
      if (lpri.gt.1) write (lun11,*)'in heatt',httot,cltot,delr,r 
      if (lpri.gt.1) write (lun11,*)ncsvn 
      numcon=ncn2 
!                                                                       
!     comment these out to implement scattering                         
      clbrems=0. 
      lskp=1 
      tmp2=0. 
      do kl=1,numcon 
        tmp2o=tmp2 
        tmp2 =brcems(kl) 
        if ( kl.ge.2 ) clbrems=clbrems+(tmp2+tmp2o)                     &
     &                *(epi(kl)-epi(kl-lskp))*ergsev/2.                 
        enddo 
!                                                                       
      delrl=delr 
!                                                                       
      fac=delrl 
      ekt = t*(0.861707) 
      htcomp = cmp1*xnx*ergsev 
      clcomp = ekt*cmp2*xnx*ergsev 
      httot=httot+htcomp 
      cltot=cltot+clcomp+clbrems 
      if (lpri.ge.1) write (lun11,9953)htcomp,clcomp,cmp1,cmp2,         &
     &   clbrems,httot,cltot                                            
      hmctot=2.*(httot-cltot)/(1.e-37+httot+cltot) 
      if (lpri.ge.1) write (lun11,*)hmctot 
 9953 format (1h , ' compton heating, cooling=',8e12.4) 
      lpri=lprisv 
!                                                                       
                                                                        
      return 
      end                                           
