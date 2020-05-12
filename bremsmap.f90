      subroutine bremsmap(bremsa,bremsam,bremsint,epi,epim,ncn2,ncn2m,  &
     &          lpri,lun11)
!
!     this routine maps the continuum flux onto the small grid
!     for faster execution of calc_hmc_all
!     arguments:
!       bremsa:  flux stored at high resolution
!       bremsam: flux stored at low resolution
!       epi: high resolution energy grid
!       epim: low resolution energy grid
!       ncn2: number of high resolution bins
!       ncn2m: number of low resolution bins
!
      use globaldata
      implicit none 
!                                                                       
       real(8) bremsa(ncn),bremsam(ncn),bremsint(ncn)
       real(8) epi(ncn),epim(ncn),sumtmp
       integer ncn2,ncn2m
       integer mm,mmm,nbinc
       integer lun11,lpri,nskp,jk
!
       do mmm=1,ncn2m
         bremsam(mmm)=0.
         if (epim(mmm).le.1.e-39) stop 'epim error'
         enddo
       do mmm=1,ncn2m
          mm=nbinc(epim(mmm),epi,ncn2)
          bremsam(mmm)=bremsa(mm)
          enddo
       do mmm=1,ncn2m
          jk=ncn2m-mmm
          sumtmp=(bremsa(jk)+bremsa(jk+1))*(epi(jk+1)-epi(jk))/2. 
          bremsint(jk)=bremsint(jk+1)+sumtmp*(1.602197e-12) 
         enddo
!
       if (lpri.gt.1) then
         write (lun11,*)'in bremsmap:'
         nskp=max(1,int(ncn2/1000))
         do mm=1,ncn2,nskp
           write (lun11,*)mm,epi(mm),bremsa(mm)
           enddo
         do mm=1,ncn2m
           write (lun11,*)mm,epim(mm),bremsam(mm)
           enddo
         endif
!
       return
       end   
