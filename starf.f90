      subroutine starf(tp,xlum,epi,ncn2,zremsz,lpri,lun11) 
!                                                                       
!                                                                       
!     Name: starf.f90
!     Description:
!       this subroutine generates the initial spectrum. 
!       blackbody
!       brems stores the flux to be used    
!       author:  T. Kallman                   
!     Parameters:                            
!         Input:
!           tp: radiation temperature in kev (for thermal spectrum) 
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
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8) epi(ncn),zremsz(ncn) 
      real(8) zremsi(ncn) 
      real(8) ergsev,del,q,xkt,sum,tempp,sum2,xlum,tp,const,expo 
      integer numcon,ncn2,lpri,lprisv,i,lun11 
!                                                                       
      data ergsev/1.602197e-12/ 
!                                                                       
      numcon=ncn2 
      del=1. 
      q=7.49e+08*del*xlum/(1.e-37+tp) 
      xkt=1.16e-03/(1.e-37+tp) 
      sum=0. 
      lprisv=lpri 
!      lpri=2                                                           
      if (lpri.gt.1) write (lun11,*)'in starf',tp,xlum,q,xkt 
      do i=1,numcon 
         tempp=epi(i)*xkt 
         zremsi(i)=0. 
!         zremsi(i)=(3.1415e+22)*epi(i)**3/exp(tempp)                   
         if (tempp.lt.1.e-3) then 
             zremsi(i)=epi(i)**3/tempp 
           else 
             if (tempp.lt.150.)                                         &
     &         zremsi(i)=epi(i)**3/(expo(tempp)-1.)                     
             if (tempp.gt.150.) zremsi(i)=xlum/epi(i)/ergsev/1.e+37 
           endif 
         zremsi(i)=(3.1415e+22)*zremsi(i) 
         if (lpri.gt.1) write (lun11,*)i,epi(i),zremsi(i) 
         if (.not.((epi(i).lt.13.6).or.(epi(i).gt.1.36e+4)              &
     &        .or.(i.le.1)))                                            &
     &    sum=sum+(zremsi(i)+zremsi(i-1))*(epi(i)-epi(i-1))/2.          
         enddo 
!                                                                       
      const=xlum/sum/ergsev 
      sum2=0. 
      do  i=1,numcon 
         zremsz(i)=zremsz(i)+zremsi(i)*const 
         if (i.gt.1)                                                    &
     &    sum2=sum2+(zremsz(i)+zremsz(i-1))*(epi(i)-epi(i-1))/2.        
         if (lpri.gt.1)                                                 &
     &        write (lun11,*)i,epi(i),zremsi(i),const,zremsz(i)         
         enddo 
      sum2=sum2*ergsev 
!      write (lun11,*)'normalization:',sum2                             
      lpri=lprisv 
!                                                                       
      return 
      end                                           
