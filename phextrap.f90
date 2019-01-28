      subroutine phextrap(etmp,stmp,ntmp,ntmp2,ett,ncn2,lpri,lun11) 
!                                                                       
!     Name:  phextrap.f90
!     Description:
!       this routine does extrapolation of the ip photoionization cros
!       sections                                                          
!       Must be used with care because extrapolation can produce unphysical
!       values at high energies
!       author: T. Kallman                                                
!     Parameters:
!       etmp(ntmp)=energy grid (RY above threshold
!       stmp(ntmp)=cross section (cm^2)
!       ntmp=length of stmp
!       ett=threshold (eV)
!       ncn2=length of epi
!       lpri=print switch
!       lun11=logical unit number for print
!       Output:
!       ntmp2=length of stmp after extrapolation
!     Dependencies:  none
!     called by ucalc, for data types 53, 49.                           
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      integer ntmp,ntmp2 
!                                                                       
      real(8) etmp(ntmp2),stmp(ntmp2) 
      real(8) ett 
      integer ncn2,lpri,lun11 
      real(8) dele,dels,s1,e1,e2,s2 
      integer nadd 
!                                                                       
      dele=1.3 
      nadd=0 
      dels=dele**3 
      s1=stmp(ntmp-1) 
      e1=etmp(ntmp-1)*13.6+ett 
      if (lpri.gt.1) write (lun11,*)'in phextrap:',ntmp,e1,s1 
      do while ((s1.gt.1.e-27).and.(nadd+ntmp.lt.ncn2)                  &
     &         .and.(e1.lt.2.e+5))                                      
        e2=e1*dele 
        s2=s1/dels 
        nadd=nadd+1 
        stmp(nadd+ntmp-1)=s2 
        etmp(nadd+ntmp-1)=(e2-ett)/13.6 
        e1=e2 
        s1=s2 
        enddo 
      ntmp=nadd+ntmp-1 
!                                                                       
      return 
      END                                           
