      subroutine levwk(rniss,rnisse,bb,lpri,nlev,t,xee,xpx,lun11)    
!                                                                       
!     Name: levwk.f90  
!     Description:  
!           Calculates and collects quantities related to level populations
!           for one ion.  
!
!     List of Parameters:
!           Input:
!           t: temperature in 10^4K
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           nlev: number of levels for the ion
!           output:
!           rniss: lte level population
!           rnisse: lte level population relative to ground 
!                  with exponential removed
!           From Globaldata:
!           rlev(10,nd):  real data for levels of this ion
!           ilev(10,nd):  integer data for levels of this ion
!           nlpt(nd):
!           iltp(nd):
!

      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8) rniss(nd),rnisse(nd) 
      real(8) ergsev,bk,t,bktm,q2,rs,ethion,emltlv,                      &
     &     eexlv,ethsht,explev2,bb,expo                                 
      integer lpri,lprisv,nlev,lun11,ll 
!                                                                       
      real(8) xnx, xpx, xee, tm 
      integer mm 
                                                                        
      data ergsev/1.602197e-12/ 
      data bk/1.38062e-16/ 
!                                                                       
      lprisv=lpri 
!      lpri=0                                                           
      xnx=xpx*xee 
      bb=1. 
      tm=t*1.e4 
      bktm=bk*tm/ergsev 
      q2=2.07e-16*xnx*(tm**(-1.5)) 
      emltlv=leveltemp%rlev(2,nlev) 
      rs=q2/emltlv 
      ethion=leveltemp%rlev(1,nlev) 
      if (lpri.gt.1)                                                    &
     & write (lun11,9902)tm,bktm,q2,                                    &
     &    emltlv,rs,ethion,xnx                                          
 9902 format (1x,'in levwk',8(1pe11.3)) 
      rniss(nlev)=1. 
      rnisse(nlev)=1. 
      do ll=1,nlev-1 
        eexlv=leveltemp%rlev(1,ll) 
        emltlv=leveltemp%rlev(2,ll) 
        ethsht=(ethion-eexlv)/bktm 
        ethsht=max(ethsht,0.d0) 
        explev2=expo(-ethsht) 
        rniss(ll)=emltlv/(explev2/rs) 
        rnisse(ll)=emltlv*rs
        bb=bb+rniss(ll) 
        if (lpri.gt.1)                                                  &
     &   write (lun11,9901)ll,eexlv,emltlv,ethsht,explev2,              &
     &     rniss(ll),rs,bb,leveltemp%ilev(1,ll),leveltemp%iltp(ll),     &
     &     leveltemp%nlpt(ll),(leveltemp%klev(mm,ll),mm=1,8)         
 9901   format (1x,i4,7(1pe11.3),3i12,8a1) 
        enddo 
        if (lpri.gt.1) write (lun11,*)'rniss:'
        do ll=1,nlev 
          rniss(ll)=rniss(ll)/bb 
          if (lpri.gt.1) write (lun11,*)ll,rniss(ll)
          enddo 
       lpri=lprisv 
!                                                                       
      return 
      END                                           
