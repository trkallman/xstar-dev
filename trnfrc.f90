      subroutine trnfrc(lpri,lun11,ldir,                                &
     &      r,xpxcol,xpx,                                               &
     &      epi,ncn2,zremsz,dpthc,opakc,                                &
     &      zrems,bremsa,bremsint)                               
!                                                                       
!     this routine calculates continuum transfer                        
!     author:  T. Kallman (from xstar1)                                 
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8) epi(ncn),zremsz(ncn),dpthc(2,ncn),bremsa(ncn) 
      real(8) bremsint(ncn),opakc(ncn),                                  &
     &          zrems(5,ncn)
      integer lpri,lun11,ldir,ncn2,jkp,jk,ncnm 
      real(8) r,xpxcol,xpx,r19,fpr2,rmax,sumtmp 
!                                                                       
      r19=r/(1.e+19) 
      fpr2=(12.56)*r19*r19 
      ncnm=ncn2-1 
      bremsa(ncn2)=0. 
      bremsa(ncnm)=0. 
      bremsint(ncn2)=0. 
      bremsint(ncnm)=0. 
      rmax=xpxcol/xpx 
      if (lpri.ne.0) write (lun11,*)'in trnfrc:',rmax,xpxcol,xpx 
      do 1 jkp=1,ncnm 
         jk=ncnm+1-jkp 
!                                                                       
!        for outward only                                               
!                                                                       
         if (ldir.lt.0) then 
             bremsa(jk)=zrems(1,jk)/fpr2 
!     $                      +flinel(jk)                                
           else 
             bremsa(jk)=zremsz(jk)*exp(-dpthc(1,jk))/fpr2 
           endif 
         sumtmp=(bremsa(jk)+bremsa(jk+1))*(epi(jk+1)-epi(jk))/2. 
         bremsint(jk)=bremsint(jk+1)+sumtmp*(1.602197e-12) 
         if (lpri.ne.0) write (lun11,*)jk,epi(jk),dpthc(1,jk),          &
     &       zremsz(jk),opakc(jk),                                      &
     &       zremsz(jk)*exp(-dpthc(1,jk))/fpr2,bremsa(jk)               &
     &       ,bremsint(jk)                                              
    1    continue 
!                                                                       
      return 
      END                                           