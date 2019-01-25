      subroutine ispec4(tp,xlum,epi,ncn2,zremsz,lpri,lun11) 
!                                                                       
!     Name: ispec4.f90
!     Description:
!       this subroutine generates the initial spectrum. 
!       Power law
!       brems stores the flux to be used    
!       author:  T. Kallman                   
!     Parameters:                            
!         Input:
!           tp: power law index (energy)
!           xlum: source luminosity integrated from 1-1000 Ry
!               in units of 10^38 erg/s
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           lpri: print switch
!           lun11: logical unit number for printing
!         Output:
!           zremsz:  input spectrum (erg s^-1 erg^-1 /10^38)
!     Dependencies: none
!     called by:  rread1
!                                                                       
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8) zremsz(ncn),epi(ncn) 
      integer ncn2,lpri,lun11 
      real(8) zremsi(ncn),ergsev,const,xlum 
      real(8) sum,ecut 
      integer i,numcon,nb1,nb2,nbinc,lprisv 
      real(8) tp 
!                                                                       
      data ergsev/1.602197e-12/ 
!                                                                       
      numcon=ncn2 
      ecut=0.01 
      sum=0. 
      lprisv=lpri 
      nb1=nbinc(13.6d0,epi,ncn2) 
      nb2=nbinc(1.36d+4,epi,ncn2) 
      if (lpri.ge.1) write (lun11,*)'in ispec4',tp,xlum,nb1,nb2 
      do i=1,numcon 
         zremsi(i)=1.e-24 
         if (epi(i).gt.ecut)                                            &
     &    zremsi(i)=epi(i)**tp                                          
         if (lpri.gt.1) write (lun11,*)i,epi(i),zremsi(i) 
         if ((i.ge.nb1).and.(i.le.nb2))                                 &
     &    sum=sum+(zremsi(i)+zremsi(i-1))*(epi(i)-epi(i-1))/2.          
         enddo 
!                                                                       
      const=xlum/sum/ergsev 
      do i=1,numcon 
         zremsz(i)=zremsz(i)+zremsi(i)*const 
         if (lpri.ge.1)                                                 &
     &        write (lun11,*)i,epi(i),zremsi(i),const,zremsz(i)         
         enddo 
      lpri=lprisv 
!                                                                       
      return 
      END                                           
