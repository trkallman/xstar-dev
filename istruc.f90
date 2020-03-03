      subroutine istruc(zeff,alpha,xitp,nnz,lpri,lun11) 
!                                                                       
!     Name: istruc.f90  
!     Description:  
!          computes ionization equilibrium                      
!          solves a system of ionization equations, attempting               
!          to avoid overflow problems.                                       
!     author:  T. Kallman                                               
!
!     List of Parameters
!         Input:
!         zeff(31): ionization rates
!         alpha(31): recombination rates
!         xitp(31): ion fractions
!         nnz: length of rate vectors
!         lpri:  print switch
!         lun11:  logical unit number for printing
!
!         Dependencies:
!         ioneqm.f90
!
!         Called by:
!         func.f90
!
      implicit none 
!                                                                       
      integer nnz,lpri,lun11 
!                                                                       
      real(8) zeff(30),alpha(30),xitp(31),xisum 
      real(8)  z8(30),a8(30),x8(31) 
      integer mm,nnzp1,ill 
!                                                                       
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'ion rates:',nnz                                  
      do mm=1,nnz 
          z8(mm)=dble(zeff(mm)) 
          a8(mm)=dble(alpha(mm)) 
          enddo 
!                                                                       
      nnzp1 = nnz + 1 
      ill=1 
      call ioneqm(z8,a8,x8,nnzp1,nnz,ill,lpri,lun11) 
!                                                                       
      xisum=0. 
      do mm=1,nnz 
          xitp(mm)=(x8(mm)) 
          xisum=xisum+xitp(mm) 
          if (lpri.ne.0)                                                &
     &     write (lun11,9901)mm,zeff(mm),alpha(mm),xitp(mm)
 9901     format (1x,i4,3(1pe11.3)) 
          enddo 
      xitp(nnz+1)=max(0.d0,1.-xisum) 
!                                                                       
      return 
      end                                           
