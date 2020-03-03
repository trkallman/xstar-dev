      subroutine calt70(temp8,den8,eth8,ic,m,np1r,np1i,              &
     &                  nx,xe8,xs8,rec8,al8,lun11,lpri,ierror)
!                                                                       
!  This rutine takes the coefficients in data type 70 (dtype70 reals    
!  in itype70 integers) and returns the recombination rate (in s-1cm-3) 
!  and the correstpondent phot. x-section for the superlevel. m is the  
!  dimension of dtype70. nx is the number of points in the x-section    
!  xe() contains the photon energy in Ry and xx() is the x-section      
!  in Mb.                                                               
!  temp, den, and ic are the temperature, electron density              
!  and effective charge of the ion respectively.                        
!  eth is the threshold energy for the superlevel in Ry.                
!                                                                       
! revised on 26 March 2019 by M. Bautista
!  ierror will returns 1 when recombination to the superlevels approaches
!  0. This indicates that the density may be too high for propper account of
!  recombination reduction.
!  Recombination rate coefficients are input in linear scale but interpolated in
! log scale.
!  does a cleaner split of rates by density and temperature mesh before
! interpolation
!
!  INTERPOLATE INTEMPERATURE BEFORE DENSITY INTERPOLATION TO AVOID DISCONTINUITIES
!  Use second order polynomial interpolation in T and linear interpolation in Ne
!                                                                       
! ********************************************************************* 
                                                                        
      use globaldata
      real(8) temp8,den8,rec8,al8,eth8,xe8(100),xs8(100)
      dimension dtype70(m),itype70(11)
! alpf: fitting coef. for hydrogenic recombination n=12,l=0             
      dimension alpf(3) 
      real temps(200),dens(200),rcoefs(200,200) 
      data alpf/-7.1094841E-02,-9.0274535E-02,-14.26129/ 
      ierror=0 
      in=0
      it=0
      if (lpri.gt.1) write (lun11,*)'in calt70, m=',m
      do ll=1,m
        dtype70(ll)=sngl(masterdata%rdat1(np1r-1+ll))
        if (lpri.gt.1) write (lun11,*)ll,dtype70(ll)
        enddo
      do ll=1,11
        itype70(ll)=masterdata%idat1(np1i-1+ll)
        enddo
      den=sngl(den8)
      temp=sngl(temp8)
      eth=sngl(eth8)
! split data                                                            
      nden=itype70(1) 
      ntem=itype70(2) 
      dmin=dtype70(1) 
      dmax=dtype70(nden) 
      tmin=dtype70(nden+1) 
      tmax=dtype70(nden+ntem) 
      nxs=itype70(3) 
      if (lpri.gt.1) write (lun11,*)'densities',nden,ntem
      do i=1,nden 
       dens(i)=dtype70(i) 
        if (lpri.gt.1) write (lun11,*)i,dens(i)
      enddo 
      if (lpri.gt.1) write (lun11,*)'temperatures'
      do i=1,ntem 
       temps(i)=dtype70(i+nden) 
        if (lpri.gt.1) write (lun11,*)i,temps(i)
      enddo 
      m=nden+ntem 
      do it=1,ntem 
       do id=1,nden 
        m=m+1 
!       this is a fudge due to some confusion 
!          about log vs no log
        if (dtype70(m).gt.-1.e-31) then
          rcoefs(it,id)=log10(dtype70(m)+1.e-30)  
          else
          rcoefs(it,id)=dtype70(m)
          endif
        if (lpri.gt.1) write (lun11,*)it,id,m,rcoefs(it,id)
       enddo 
      enddo 
                                                                        
      rne=log10(den) 
      rte=log10(temp) 
      if (lpri.ge.1) then
        write (lun11,*)'in calt70:',rne,rte,nden,ntem,nxs,dmin,dmax, &
     &    tmin,tmax
        endif
!     test if nden>1
      if (nden.gt.1) then 
!       test if density > dmax          
        if (rne.gt.dmax) then 
          print*,'DENSITY TOO HIGH AT SUPREC' 
          print*,'z=',ic,' temp=',temp,' Ne=',den 
!          stop 
!         end of test if density > dmax          
          endif 
!       find intervals for density and temperature                            
!       in is the index for density interval for interpolation                
!       test if density < dmin
        if (rne.le.dmin) then 
            in=1 
          else 
            do i=1,nden-1 
              if (rne.ge.dens(i).and.rne.le.dens(i+1))                  &
     &          in=i                                                       
              enddo 
!         end of test if density < dmin
          endif 
          else
           in=1
          endif
        if (rte.lt.tmin.or.rte.gt.tmax) then 
          if (lpri.ne.0) then
            write (lun11,*)'TEMPERATURE OUT OF RANGE' 
            write (lun11,*)'z=',ic,' temp=',temp,' Tmin=',10.**tmin,    &
     &     ' Tmax=',10.**tmax                                           
            endif
!          stop 
          rte=min(0.999*tmax,max(1.001*tmin,rte))
          endif 
        do i=1,ntem-1 
          if (rte.ge.temps(i).and.rte.lt.temps(i+1))                    &
     &       it=i                                                       
          enddo 
! Interpolate in Temperature space
      kt=0
       x1=temps(it+kt)
       x2=temps(it+kt+1)
       y1=rcoefs(it+kt,in)
       y2=rcoefs(it+kt+1,in)
       call qcoefs(x1,x2,y1,y2,rte,y)
       rec1=y
       if (lpri.ge.1)                                                   &
     &  write (lun11,*)'interpolate in temperature',                    &
     &       it,kt,in,x1,x2,y1,y2,y
       if (in.eq.nden .or. in.le.1) then
        rec=10.**rec1
       if (lpri.gt.1)                                                   &
     &  write (lun11,*)'in=nden or in=1',rec,rec1
       else
        y1=rcoefs(it+kt,in+1)
        y2=rcoefs(it+kt+1,in+1)
        call qcoefs(x1,x2,y1,y2,rte,y)
        rec2=y
       if (lpri.gt.1)                                                   &
     &  write (lun11,*)'rec2=y',rec2,y1,y2
! interpolate in density
        dden=1./(dens(in+1)-dens(in))
        rm=(rec2-rec1)*dden
        rr=rec1+rm*(rne-dens(in))
        rec=10.**rr
       if (lpri.ge.1)                                                   &
     &  write (lun11,*)'interpolate in density',dden,rm,rr,rec
       endif
! if type70 recombination near zero flag error
      if (rec.lt.1.e-29) ierror=1
      if (lpri.ge.1)                                                    &
     & write (lun11,*)'in calt70:',rec,rr,rm,dden,rec2,y1,y2,x1,x2,  &
     &      it,kt,in
!                                                                       
! scale hydrogenic cross section
      i1=ntem*nden+ntem+nden 
      do i=1,nxs 
       xe8(i)=dtype70(i1+(i-1)*2+1) 
       xs8(i)=dtype70(i1+(i-1)*2+2) 
      enddo 
      call milne(temp8,nxs,xe8,xs8,eth8,al8,lun11,lpri) 
      scale=rec/sngl(al8)
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in calt70:',rec,al8,scale,xs8(1),nxs,eth8           
      rec8=rec 
      do i=1,nxs 
       xs8(i)=xs8(i)*scale 
      enddo 
      nx=nxs 
      return 
      END                                           
