      subroutine calt77(lpri,lun11,temp,den,                            &
     &                       np1r,np1i,cul,clu)             
!                                                                       
!  This rutine takes the coefficients in data type 77 (dtype77 reals    
!  and itype77 integers) and returns the collisional transition rates   
!  (in s-1) from the superlevel (cul) and to the superlevel (clu) from  
!  the spectroscopic level given by itype77(3).                         
!  The wavelength for the transition is also given in wav               
!  temp, den, and ic are the temperature, electron density              
!  and effective charge of the ion respectively.                        
!      author: M. Bautista                                              
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8) temp,den,cul,clu,rne,rte,div,rm,rec1,rec2,                 &
     &         gg,wav,rec,xt                                            
      integer nden,ntem,in,it,kt1,nll,k,nl1,nl2,il,                      &
     &        lpri,lun11,np1r,np1i                                      
!                                                                       
      rne=log10(den) 
      rte=log10(temp) 
      nden=masterdata%idat1(np1i-1+1) 
      ntem=masterdata%idat1(np1i-1+2) 
      if (rne.gt.masterdata%rdat1(np1r-1+nden)) then 
!       print*,'DENSITY TOO HIGH AT CALT77'                             
!       print*,'z=',ic,' temp=',temp,' Ne=',den,nden,rdat1(np1r-1+nden) 
       rne=masterdata%rdat1(np1r-1+nden) 
      endif 
      if (rte.gt.(masterdata%rdat1(np1r-1+nden+ntem)+1.)) then 
!       print*,'TEMPERATURE TOO HIGH AT CALT77'                         
!       print*,'z=',ic,' temp=',temp,' Ne=',den                         
         rte=(masterdata%rdat1(np1r-1+nden+ntem)+1.) 
      endif 
      if (rte.lt.(masterdata%rdat1(np1r-1+nden+1)-1.)) then 
       rte=masterdata%rdat1(np1r-1+nden+1)-1. 
      endif 
!                                                                       
      wav=masterdata%rdat1(np1r-1+nden*ntem+nden+ntem+1) 
      if (rne.le.masterdata%rdat1(np1r-1+1)) then 
       in=1 
      else 
       in=0 
    5  in=in+1 
       if (rne.ge.masterdata%rdat1(np1r-1+in+1).and.in.lt.nden) goto 5 
      endif 
      if (rte.lt.masterdata%rdat1(np1r-1+nden+1)) then 
       it=1 
      else 
       it=0 
    6  it=it+1 
       if (it.ge.ntem) then 
        it=ntem-1 
       else 
        if (rte.ge.masterdata%rdat1(np1r-1+nden+it+1)) goto 6 
       endif 
      endif 
!                                                                       
      kt1=nden+ntem+(in-1)*ntem+it 
      div=masterdata%rdat1(np1r-1+nden+it+1)                            &
     &      -masterdata%rdat1(np1r-1+nden+it) 
      rm=(masterdata%rdat1(np1r-1+kt1+1)                                &
     &      -masterdata%rdat1(np1r-1+kt1))/(div+1.d-36) 
      rec1=masterdata%rdat1(np1r-1+kt1)                                 &
     &      +rm*(rte-masterdata%rdat1(np1r-1+nden+it)) 
      kt1=kt1+ntem 
      rm=(masterdata%rdat1(np1r-1+kt1+1)-masterdata%rdat1(np1r-1+kt1))  &
     &    /(masterdata%rdat1(np1r-1+nden+it+1)-                         &
     &    masterdata%rdat1(np1r-1+nden+it)+1.d-36)                     
      rec2=masterdata%rdat1(np1r-1+kt1)                                 &
     &      +rm*(rte-masterdata%rdat1(np1r-1+nden+it)) 
!                                                                       
      rm=(rec2-rec1)/(masterdata%rdat1(np1r-1+in+1)                     &
     &      -masterdata%rdat1(np1r-1+in)+1.d-36) 
      rec=rec1+rm*(rne-masterdata%rdat1(np1r-1+in)) 
      if (lpri.gt.1) write (lun11,*)'in calt77:',                       &
     & temp,den,nden,ntem,rte,rne,wav,in,it,div,rm,                     &
     & rec1,rec2,rec                                                    
      cul=10.**rec 
!                                                                       
      nll=masterdata%idat1(np1i-1+3) 
      k=0 
    7 k=k+1 
      nl1=k*(k-1)/2+1 
      nl2=(k+1)*k/2+1 
      if (nll.ge.nl2) goto 7 
      il=nll-nl1 
      gg=float(2*il+1)*2. 
      xt=1.43817e+8/wav/temp 
      if (xt.lt.100) then 
       clu=cul*exp(-xt)/gg 
      else 
       clu=0.e0 
      endif 
!                                                                       
!                                                                       
      return 
      END                                           
