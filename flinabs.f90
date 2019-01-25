      function flinabs(ptmp) 
!                                                                       
!     Name: flinabs.f90
!     Description:
!        this function will treat the absorption of incident               
!        continuum in the line. 
!        currently returns trivial answer                               
!     Parameters:
!        Input:
!        ptmp=escape probability
!        Output:
!        flinabs=fraction of line photons which escape local 
!                scattering/absorption
!     Dependencies: none
!     Called by:  ucalc
!                                                                       
      implicit none 
!                                                                       
      real(8) flinabs, ptmp 
!                                                                       
!      flinabs=0.                                                       
      flinabs=1. 
!                                                                       
      return 
      end                                           
