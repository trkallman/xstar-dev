      subroutine levwk(rniss,rnisse,bb,lpri,nlev,t,xee,xpx,lun11)    
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
 9901   format (1x,i4,7(1pe11.3),3i6,8a1) 
        enddo 
        do ll=1,nlev 
          rniss(ll)=rniss(ll)/bb 
          enddo 
       bb=1. 
       lpri=lprisv 
!                                                                       
      return 
      END                                           
