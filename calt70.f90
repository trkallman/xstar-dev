      subroutine calt70(temp,den,eth,ic,m,np1r,np1i,                    &
     &                  nx,xe,xs,rec,al,lun11,lpri)                     
!                                                                       
!     Name: calt70.f90  
!     Description:  
!      This routine takes the coefficients in data type 70 (dtype70 reals   
!      in itype70 integers) and returns the recombination rate (in s-1cm-3) 
!      and the correstpondent phot. x-section for the superlevel. m is the  
!      dimension of dtype70. nx is the number of points in the x-section    
!      xe() contains the photon energy in Ry and xx() is the x-section      
!      in Mb.                                                               
!      temp, den, and ic are the temperature, electron density              
!      and effective charge of the ion respectively.                        
!      eth is the threshold energy for the superlevel in Ry.                
!      new revision (2017) quadratic density dependence tk
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
!        nx=number of energy points in cross section
!        xe=energy grid for cross section
!        xs=cross section (cm^2)
!        rec=recombination rate coefficient (cm^3 s^-1)
!        al=milne recombination rate coefficient (cm^3 s^-1)
!     Dependencies: milne
!     called by:  ucalc
!                                                                       
                                                                        
      use globaldata
       implicit none 
      integer nptmpdim 
      parameter (nptmpdim=200000) 
!                                                                       
      integer m 
      real(8) xe(*),xs(*),rne,rte,rme,rmt,                              &
     &      rec1,rec2,rec3,rec,scale,al,crit,temp,den,eth,dt            
      real(8) x1,x2,x3,y1,y2,y3,aa,bb,cc,denom,rmt1,rmt3 
      integer nden,ntem,nxs,in,it,kt1,kt3,i1,imax,nx,                   &
     &      lun11,lpri,lprim,ic,i,np1r,np1i,in1,in2,in3                 
!                                                                       
!     alpf: fitting coef. for hydrogenic recombination n=12,l=0         
!      dimension alpf(3)                                                
!      data alpf/-7.1094841E-02,-9.0274535E-02,-14.26129/               
!                                                                       
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in calt70:',temp,den,eth,ic,m,                  &
     &          masterdata%rdat1(np1r),masterdata%idat1(np1i-1+1)      
      rne=log10(den) 
      rte=log10(temp) 
      nden=masterdata%idat1(np1i-1+1) 
      ntem=masterdata%idat1(np1i-1+2) 
      nxs=masterdata%idat1(np1i-1+3) 
      if (nden.gt.1) then 
!     nb changing the data                                              
      masterdata%rdat1(np1r+1)=min(masterdata%rdat1(np1r+1),8.d0) 
      if (rne.gt.masterdata%rdat1(np1r-1+nden)) then 
!       write (lun11,*)'DENSITY TOO HIGH AT SUPREC'                     
!       write (lun11,*)'z=',ic,' temp=',temp,' Ne=',den                 
!       return                                                          
       rne=min(rne,masterdata%rdat1(np1r-1+nden)) 
      endif 
      if (rne.le.masterdata%rdat1(np1r-1+1)) then 
       in=1 
      else 
       in=int(rne/masterdata%rdat1(np1r-1+nden)*nden)-1 
       if (in.ge.nden) in=in-1 
    5  in=in+1 
       if (in.lt.nden .and. rne.ge.masterdata%rdat1(np1r-1+in+1)) goto 5 
       if (rne.lt.masterdata%rdat1(np1r-1+in)) then 
        in=in-2 
        goto 5 
       endif 
      endif 
      else 
       in=1 
      endif 
      if (rte.lt.masterdata%rdat1(np1r-1+nden+1)) then 
       it=1 
      else 
         dt=(masterdata%rdat1(np1r-1+nden+ntem)                         &
     &         -masterdata%rdat1(np1r-1+nden+1))/float(ntem) 
       it=int((rte-masterdata%rdat1(np1r-1+nden+1))/dt) 
    6  it=it+1 
       if (it.ge.ntem) then 
        it=ntem-1 
       else 
        if (rte.ge.masterdata%rdat1(np1r-1+nden+it+1)) goto 6 
        if (rte.lt.masterdata%rdat1(np1r-1+nden+it)) then 
         it=it-2 
         goto 6 
        endif 
       endif 
      endif 
      kt1=nden+ntem+(in-1)*ntem+it 
      rmt=(masterdata%rdat1(np1r-1+kt1+1)-masterdata%rdat1(np1r-1+kt1)) &
     &    /(masterdata%rdat1(np1r-1+nden+it+1)-                         &
     &    masterdata%rdat1(np1r-1+nden+it))                                  
      rec1=masterdata%rdat1(np1r-1+kt1)                                 &
     &      +rmt*(rte-masterdata%rdat1(np1r-1+nden+it)) 
      if (nden.gt.1) then 
         if ((nden.gt.2).and.(in.gt.1).and.(in.lt.nden)) then 
!             now we implement quadratic                                
              in1=in-1 
              in2=in 
              in3=in+1 
              y2=rec1 
              x2=masterdata%rdat1(np1r-1+in2) 
              kt3=nden+ntem+(in3-1)*ntem+it 
              rmt3=(masterdata%rdat1(np1r-1+kt3+1)                      &
     &              -masterdata%rdat1(np1r-1+kt3))                      &
     &         /(masterdata%rdat1(np1r-1+nden+it+1)                     &
     &             -masterdata%rdat1(np1r-1+nden+it))         
              rec3=masterdata%rdat1(np1r-1+kt3)                         &
     &              +rmt3*(rte-masterdata%rdat1(np1r-1+nden+it)) 
              y3=rec3 
              x3=masterdata%rdat1(np1r-1+in3) 
              kt1=nden+ntem+(in1-1)*ntem+it 
              rmt1=(masterdata%rdat1(np1r-1+kt1+1)                      &
     &           -masterdata%rdat1(np1r-1+kt1))                         &
     &         /(masterdata%rdat1(np1r-1+nden+it+1)                     &
     &          -masterdata%rdat1(np1r-1+nden+it))         
              rec1=masterdata%rdat1(np1r-1+kt1)                         &
     &           +rmt1*(rte-masterdata%rdat1(np1r-1+nden+it)) 
              y1=rec1 
              x1=masterdata%rdat1(np1r-1+in1) 
              denom=((x1*x1-x2*x2)*(x1-x3)-(x1*x1-x3*x3)*(x1-x2)) 
              aa=((y1-y2)*(x1-x3)-(y1-y3)*(x1-x2))/denom 
              bb=-((y1-y2)*(x1*x1-x3*x3)-(y1-y3)*(x1*x1-x2*x2))/denom 
              cc=y1-aa*x1*x1-bb*x1 
              rec=aa*rne*rne+bb*rne+cc 
              if (lpri.gt.1) write (lun11,*)'quadratic:',x1,x2,x3,      &
     &             y1,y2,y3,aa,bb,cc,denom,rec                          
           else 
              kt1=kt1+ntem 
              rmt=(masterdata%rdat1(np1r-1+kt1+1)                       &
     &              -masterdata%rdat1(np1r-1+kt1))                      &
     &         /(masterdata%rdat1(np1r-1+nden+it+1)-                    &
     &         masterdata%rdat1(np1r-1+nden+it))                             
              rec2=masterdata%rdat1(np1r-1+kt1)                         &
     &              +rmt*(rte-masterdata%rdat1(np1r-1+nden+it)) 
              rme=(rec2-rec1)                                           &
     &         /(masterdata%rdat1(np1r-1+in+1)                          &
     &             -masterdata%rdat1(np1r-1+in)) 
              rec=rec1+rme*(rne-masterdata%rdat1(np1r-1+in)) 
            endif 
        else 
          rec=rec1 
        endif 
      if (lpri.gt.1) write (lun11,*)nden,ntem,it,in,kt1,                &
     &     masterdata%rdat1(np1r-1+kt1),rme,rmt,rec2,rec1,rec              
      rec=10.**rec 
!                                                                       
      i1=ntem*nden+ntem+nden 
      do i=1,nxs 
       xe(i)=masterdata%rdat1(np1r-1+i1+(i-1)*2+1) 
       xs(i)=masterdata%rdat1(np1r-1+i1+(i-1)*2+2) 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)i,xe(i),xs(i)                                    
      enddo 
      lprim=0 
      call milne(temp,nxs,xe,xs,eth,al,lun11,lprim) 
      scale=rec/(1.e-24+al) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in calt70:',rec,al,scale,xs(1),nxs,eth           
      crit=1.e-6 
      imax=nxs 
      do i=1,nxs 
       xs(i)=xs(i)*scale 
       xs(i)=min(xs(i),1000000.d0) 
       if (xs(i).gt.xs(1)*crit) imax=i 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)i,xe(i),xs(i)                                    
      enddo 
      nxs=imax 
      nx=nxs 
!                                                                       
      return 
      END                                           
