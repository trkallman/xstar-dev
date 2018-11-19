      subroutine ener(epi,ncn2) 
!                                                                       
!     sets up energy grin                                               
!     author: T. Kallman                                                
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      real(8) epi(ncn) 
      integer numcon,numcon2,numcon3,ncn2,ll,ll2 
      real(8) ebnd1,ebnd2,ebnd2o,dele 
!                                                                       
      numcon = ncn2 
      if (numcon.gt.ncn) stop 'ncn2 too large for arry dimension' 
      if (numcon.lt.4) write (6,*) 'in ener: numcon error' 
      numcon2=max(2,ncn2/50) 
      numcon3=numcon-numcon2 
      ebnd1=0.1 
!     nb changed energy grid for H only                                 
      ebnd2=4.e+5 
!      ebnd2=4.e+1                                                      
      ebnd2o=ebnd2 
      dele=(ebnd2/ebnd1)**(1./float(numcon3-1)) 
      epi(1)=ebnd1 
!      write (6,*)'in ener',ncn2,numcon,numcon2,numcon3                 
      do ll=2,numcon3 
        epi(ll)=epi(ll-1)*dele 
        enddo 
      ebnd2=1.e+6 
      ebnd1=ebnd2o 
      dele=(ebnd2/ebnd1)**(1./float(numcon2-1)) 
      do ll2=1,numcon2 
        ll=ll2+numcon3 
        epi(ll)=epi(ll-1)*dele 
        enddo 
!                                                                       
      return 
      end                                           
