      function nbinc(e,epi,ncn2) 
!                                                                       
!     Name:  nbinc.f90
!     Description:
!        this function bins the continuum                                  
!        lines between   epi(i) and   epi(i+1) are put in bin number i.    
!        energies between 0 and   epi(1) are put in bin number 50.         
!        author:  T. Kallman                                               
!     Parameters:
!        Input:
!        e=energy to be binned (ev)
!        epi(ncn)=energy grid (ev)
!        ncn2=length of epi
!        Output:
!        nbinc=bin number
!     Dependencies: huntf
!     Called by:  pprint.f90, ucalc.f90, ... many 
!
      use globaldata
!                                                                       
      real(8) e 
      integer jlo, lun11,nbinc, ncn2, numcon, numcon2,                  &
     &        numcon3                                                   
      real(8) epi(ncn) 
!                                                                       
      lun11=6 
      numcon=ncn2 
      numcon2=max0(2,ncn2/50) 
      numcon3=numcon-numcon2 
!      write (lun11,*)'in nbinc',e                                      
      call huntf(epi,numcon3,e,jlo,0,lun11) 
!      call hunt3(epi,numcon3,e,jlo,0,lun11)                            
!      if (abs(e-epi(jlo+1)).lt.abs(e-epi(jlo))) jlo=jlo+1              
      nbinc=jlo 
!                                                                       
      return 
      end                                           
