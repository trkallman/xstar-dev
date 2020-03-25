      subroutine find53(stmpp,etmpp,ntmp,efnd,sg,jlo,lun11,lpri) 
!                                                                       
!     Name: find53.f90  
!     Description:  
!     this routine finds the continuum bin index for integrating over type 
!     53 data             
!     author T. Kallman                                                 
!     Parameters:
!          Input:
!          stmpp(ntmp)= type 53 photoionization cross section array (cm^2)
!          etmpp(ntmp)= type 53 energy array (eV)
!          ntmp=length of stmpp, etmpp
!          efnd=input energy(eV)
!          lpri= print switch
!          lun11= logical unit number
!          Output:
!          sg=cross section (cm^2)
!          jlo=bin in type 53 data
!     Dependencies:  hunt3
!     Called by:  phint53hunt
!                                                                       

      implicit none 
!                                                                       
      integer ntmp 
      real(8) stmpp(ntmp),etmpp(ntmp) 
      integer jlo,lun11,lpri 
      real(8) efnd,sg 
      integer ml2,mlp 
      real(8) del1,del2,alg1,alg2,algtmp 
!                                                                       
!      lpri=0                                                           
!                                                                       
!      if ((efnd.ge.etmpp(1)).and.(efnd.le.etmpp(ntmp))) then           
<<<<<<< HEAD
      if (lpri.gt.0) write (lun11,*)'in find53:',efnd,ntmp,          &
=======
      if (lpri.ne.0) write (lun11,*)'in find53:',efnd,ntmp,          &
>>>>>>> 2d75308c63b9789458ce092c697c7853fcdde44a
     &    etmpp(1),etmpp(ntmp),stmpp(1),stmpp(ntmp)                     
      if ((efnd.ge.0.).and.(efnd.le.etmpp(ntmp))) then 
        call hunt3(etmpp,ntmp,efnd,jlo,0,lun11) 
        ml2=max(jlo,1) 
        ml2=min(ml2,ntmp-1) 
        mlp=ml2+1 
        if (mlp.eq.ntmp) then 
            alg1=log(max(stmpp(mlp),1.d-26)/max(stmpp(ml2),1.d-26)) 
            alg2=log(max(etmpp(mlp),1.d-26)/max(etmpp(ml2),1.d-26)) 
            algtmp=alg1/alg2 
            sg=stmpp(ml2)*(efnd/etmpp(ml2))**algtmp 
          else 
            del1=(efnd-etmpp(ml2))/(etmpp(mlp)-etmpp(ml2)) 
            del2=(efnd-etmpp(mlp))/(etmpp(mlp)-etmpp(ml2)) 
            sg=-stmpp(ml2)*del2+stmpp(mlp)*del1 
          endif 
         sg=max(0.d0,sg) 
         if (lpri.gt.0)                                                 &
     &     write (lun11,*)sg,ml2,stmpp(ml2),stmpp(mlp),                 &
     &           del1,del2,efnd,etmpp(ml2)                              
         else 
              sg=0. 
         endif 
!                                                                       
      return 
      END                                           
