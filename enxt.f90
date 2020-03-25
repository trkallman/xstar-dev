      subroutine enxt(eth,nb1,lpri,epi,ncn2,t,lfast,lun11,           &
     &                  jk,nskp,nphint,lrcalc)                          
!                                                                       
!     Name: enxt.f90  
!     Description:  
!     This routine finds next energy bin for photoionizaion rate integrations 
!     author: T. Kallman                                                
!     Parameters:
!           input:
!           eth=photoionization threshold energy (eV)
!           nb1=index of threshold energy
!           lpri=print switch
!           epi(ncn)=energy grid (eV)
!           ncn2=length of epi
!           t=temperature (10^4 K)
!           lfast=mode switch:  1,2 --> nskp=1, nphint corresponds to 40 keV; 
!                               3 --> 16 steps to nphint, 
!                                 nphint corresponds to max(eth+3kT,3*eth);
!                               other --> nskp=1, nphint correspnds to 10 keV
!           lun11=logical unit number for printing
!           jk=current energy bin index
!           Output:
!           nskp=number of bins to skip for next index
!           nphint=maximum index for integration
!           lrcalc=rate calculation switch
!     Dependencies: nbinc
!     Called by:  ucalc
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8) epi(ncn) 
      real(8) ergsev,bk,tm,t,eth,bktm,exptst,epii 
      integer nb1,ncn2,lfast,lun11,jk,nphint,lrcalc,lpri 
      integer nskp,numcon2,nbinc,numcon3,nskp1,numcon,nskp2 
!                                                                       
      data ergsev/1.602197e-12/ 
      data bk/1.38062e-16/ 
!                                                                       
       if (lpri.gt.2)                                                   &
     &  write (lun11,*)'in enxt:',eth,nb1,t,lfast,jk,lpri,           &
     &                    epi(1),epi(ncn2),ncn2                         
      tm=t*1.e4 
      bktm=bk*tm/ergsev 
      if (lfast.le.2) then 
         numcon2=max(2,ncn2/50) 
         nphint=ncn2-numcon2 
         nskp=1 
         nskp2=1 
      elseif (lfast.eq.3) then 
         nphint=nbinc(max(3.*eth,eth+3.*bktm),epi,ncn2) 
         nphint=max(nphint,nb1+1) 
         nskp=max(1,int((nphint-nb1)/16)) 
         nskp2=nskp 
      else 
        nphint=nbinc(1.d+4,epi,ncn2) 
        nskp=1 
        nskp2=1 
        endif 
      nskp1=nskp 
      epii=epi(jk) 
      exptst=(epii-eth)/bktm 
      if (exptst.lt.3.) then 
         lrcalc=1 
         nskp=nskp1 
       else 
         lrcalc=0 
         nskp=nskp2 
       endif 
       nphint=max(nphint,nb1+nskp) 
       numcon=ncn2 
       numcon2=max(2,ncn2/50) 
       numcon3=numcon-numcon2 
       nphint=min(nphint,numcon3) 
       if (lpri.gt.2)                                                   &
     &  write (lun11,*)'in enxt:',eth,nb1,t,lfast,jk,nskp,           &
     &   nphint,lrcalc                                                  
!                                                                       
      return 
      END                                           
