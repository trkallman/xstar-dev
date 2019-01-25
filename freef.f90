!-----------------------------------------------------------------------
      subroutine freef(lpri,lun11,epi,ncn2,t,xpx,xee,opakc) 
!                                                                       
!     Name: freef.f90
!     Description:
!         this sub-routine computes the free-free opacity and               
!         include it into the total one (opakc)                             
!          author:  J. Garcia (July 2008)                                    
!     Parameters:
!         Input:
!         lpri=print switch
!         lun11=locial unit number for printing
!         epi(ncn): photon energy grid (ev)
!         ncn2: length of epi
!         t: temperature in 10^4K
!         xpx: H number density (cm^-3)
!         xee: electron fraction relative to H
!         Output:
!         opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!     Dependencies:  none (fbg)
!     Called by:  func.f90
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      real(8) opakc(ncn),epi(ncn) 
      real(8) t, opaff, ekt, t6, temp 
      real(8) xpx, xee, xnx, enz2, cc 
      real(8)  gau,  zz, gam
      integer numcon,lpri,lun11,ncn2,kk 
!                                                                       
!                                                                       
      data cc/2.614e-37/ 
!                                                                       
      if (lpri.gt.0) write (lun11,*)'in freef',t 
!                                                                       
      numcon=ncn2 
      xnx=xpx*xee 
      ekt = t*(0.861707) 
      t6 = t/100. 
      enz2=(1.4)*xnx 
      zz=1. 
      do kk=1,numcon 
         temp = epi(kk)/ekt 
         gam = zz*zz*(0.158)/t6 
         gau = 1. 
!         if ( temp.lt.100. ) gau = fbg(temp,gam)                       
!         if ( temp.lt.100. )                                           
!     1    gau = 10.**(0.2258*epi(kk)**(0.08)*(4.094-log10(epi(kk)))    
!     2      +log10(t6)*(0.133*(4.094-log10(epi(kk)))-0.2)-0.538)   !JG 
                                                                        
         opaff = cc*xnx*enz2*gau/sqrt(t)/epi(kk)**3.                    &
     &           *(1. - exp(-temp))                                     
                                                                        
         opakc(kk) = opakc(kk) + opaff 
!                                                                       
!!! THIS IS A TEST                                                      
!         if(opakc(kk).gt.6.65e-7)opakc(kk)=6.65e-7                     
!                                                                       
      enddo 
!                                                                       
      return 
      end                                           
