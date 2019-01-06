      subroutine calt71(temp,den,ic,m,np1r,np1i,wav,aij,                 &
     &                  lun11,lpri)                                     
!                                                                       
!     Name: calt71.f90  
!     Description:  
!      This routine takes the coefficients in data type 71 (dtype71 reals    
!      in itype71 integers) and returns the radiative transition prbability 
!      (in s-1) from the superlevels to the spectroscopic level given by    
!      itype71(3).                                                          
!      The wavelength for the transition is also given in wav               
!      temp, den, and ic are the temperature, electron density              
!      and effective charge of the ion respectively.                        
!      author: M. Bautista                                              
!     Parameters:
!        Input:
!        temp=temperature in K
!        den=density in cm^-3
!        ic=ion charge
!        m=?
!        np1r=pointer to beginning of real data
!        np1i=pointer to beginning of integer data
!        lpri=print switch
!        lun11=logical unit number
!        Output:
!        wav=wavelength of transition (A)
!        aij=transition probability (s^-1)
!     Dependencies: none
!     called by:  ucalc
!                                                                       
      use globaldata
       implicit none 
!                                                                       
      integer m 
      real(8) wav,aij,temp,den,rne,rte,dtmp,rm,rec1,                     &
     &     rec,rec2                                                     
      integer lun11,lpri,nden,ntem,in,it,kt1,ic 
      integer javi,np1r,np1i 
!                                                                       
      javi=m 
      m=javi 
!                                                                       
      rne=log10(den) 
      rte=log10(temp) 
      nden=masterdata%idat1(np1i-1+1) 
      ntem=masterdata%idat1(np1i-1+2) 
      if (lpri.gt.1) write (lun11,*)'in calt71:',nden,ntem 
                                                                        
      if (nden.eq.1 .and. ntem.eq.1) then 
        if (masterdata%rdat1(np1r-1+3).gt.30.) then 
          dtmp=log10(masterdata%rdat1(np1r-1+3)) 
        else 
          dtmp=masterdata%rdat1(np1r-1+3) 
        endif 
       wav=masterdata%rdat1(np1r-1+4) 
       aij=10.**dtmp 
!       aij=min(aij,1.e+12)                                             
       if (lpri.gt.1) write (lun11,*)'early return',aij,wav 
       return 
      endif 
      if (rne.gt.masterdata%rdat1(np1r-1+nden)) then 
        if (lpri.gt.1) then 
          write (lun11,*)'DENSITY TOO HIGH AT CALT71' 
          write (lun11,*)'z=',ic,' temp=',temp,' Ne=',den,nden,         &
     &           masterdata%rdat1(np1r-1+nden)                             
          endif 
        rne=min(rne,masterdata%rdat1(np1r-1+nden)) 
      endif 
      if (rte.gt.(masterdata%rdat1(np1r-1+nden+ntem)+1.)) then 
       rte=masterdata%rdat1(np1r-1+nden+ntem)+1. 
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
      kt1=nden+ntem+(in-1)*ntem+it 
      rm=(masterdata%rdat1(np1r-1+kt1+1)-masterdata%rdat1(np1r-1+kt1))  &
     &   /(masterdata%rdat1(np1r-1+nden+it+1)-                          &
     &    masterdata%rdat1(np1r-1+nden+it))                            
      rec1=masterdata%rdat1(np1r-1+kt1)                                 &
     &   +rm*(rte-masterdata%rdat1(np1r-1+nden+it)) 
      kt1=kt1+ntem 
      rm=(masterdata%rdat1(np1r-1+kt1+1)-masterdata%rdat1(np1r-1+kt1))  &
     &    /(masterdata%rdat1(np1r-1+nden+it+1)-                         &
     &    masterdata%rdat1(np1r-1+nden+it))                                  
      rec2=masterdata%rdat1(np1r-1+kt1)                                 &
     &    +rm*(rte-masterdata%rdat1(np1r-1+nden+it)) 
!                                                                       
      rm=(rec2-rec1)                                                    &
     &   /(masterdata%rdat1(np1r-1+in+1)-masterdata%rdat1(np1r-1+in)) 
      rec=rec1+rm*(rne-masterdata%rdat1(np1r-1+in)) 
      aij=10.**rec 
!      aij=min(aij,1.e+12)                                              
      if (lpri.gt.1) write (lun11,*)'late return',rm,rec2,              &
     &       rec1,rec,aij,wav                                           
!                                                                       
      return 
      END                                           
