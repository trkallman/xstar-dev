      subroutine calt57(te,den2,e,ep,n,cion,crec,lun11,lpri) 
!                                                                       
!     Name: calt57.f90  
!     Description:  
!      this routine calculates rates for type 57 data, collisional ionizatio
!      author: M. Bautista                                              
!     Parameters:
!        Input:
!         temp=temperature in K                                              
!         den=electron density in cm-3                                      
!         e=level energy in eV (first real(8) in type 6)                   
!         ep=ionization potential in eV (forth real(8) in type 6)           
!         n= level's principal quatum number (first integer in type 6)     
!        lpri=print switch
!        lun11=logical unit number
!        Output:
!         cion=  ionization rate in s-1.cm+3                                   
!         crec=   3-body recombination rate in s-1.cm+6. THIS VALUE MUST BE     
!         MUTIPLIED BY (stat. weigth level/stat. weigth recombining     
!          level)                                                       
!     Dependencies: irc, eint
!     called by:  ucalc
!                                                                       
!                                                                       
!      author: M. Bautista                                              
!                                                                       
      use globaldata
       implicit none 
!                                                                       
       real(8) te,den,den2,e,ep,cion,crec 
       integer n,lun11,lpri 
       real(8) rk,cb,rio,rc,temp,tmin,rno,rno2,ciono,                    &
     &      beta,wte,wtm,e2,e3,ete,etm                                  
       real(8) rn 
!                                                                       
       data e3/0./,e2/0./ 
                                                                        
       cion=0. 
       crec=0. 
       rn=float(n) 
       rk=1.16058e+4 
       cb=13.605692*1.6021e-19/1.3805e-23 
       if (lpri.gt.1)                                                   &
     &  write (lun11,*)'entering calt57:',den2,te,                   &
     &   ep,e,n                                                         
       if (ep.lt.e) return 
       rio=(ep-e)/13.6 
       rc=sqrt(rio)*rn 
!       if (den.gt.1.e+18) then                                         
!        print*,'density too high for cics'                             
!        return                                                         
!       endif                                                           
       den=min(den2,1.d+18) 
       tmin=3.8e+4*rc*sqrt(rc) 
       temp=te 
       if (te.lt.tmin) then 
        temp=tmin 
       endif 
       rno=SQRT(1.8887E+8*rc/den**0.3333) 
       rno2=(1.814e26*(rc**6)/2./den)**0.13333 
       rno=min(rno,rno2) 
       if (int(rno).gt.n) then 
        call irc(n,temp,rc,rno,ciono,lpri,lun11) 
       if (lpri.gt.1)                                                   &
     &  write (lun11,*)'in calt57:',rno,den,rc,tmin,te,              &
     &   ciono                                                          
!                                                                       
! extrapolates to actual temperature below Tmin                         
!                                                                       
        cion=0. 
        crec=0. 
        if (te.lt.tmin) then 
         beta=.25*(sqrt((100.*rc+91.)/(4.*rc+3.))-5.) 
         wte=(log(1.+te/cb/rio))**(beta/(1.+te/cb*rio)) 
         wtm=(log(1.+tmin/cb/rio))**(beta/(1.+tmin/cb*rio)) 
         call eint(rio/te*cb,ete,e2,e3) 
         if (ete.lt.1.e-20) return 
         call eint(rio/tmin*cb,etm,e2,e3) 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*)'in calt57:',                               &
     &     n,temp,tmin,rc,rno,cion,te,ete,etm,wte,wtm,cion              
         cion=ciono*sqrt(tmin/te)*ete/(etm+1.e-30)*wte/(wtm+1.e-30) 
        else 
         cion=ciono 
        endif 
!                                                                       
        if (cion.le.1.e-24) return 
!                                                                       
!        se=log(cion)                                                   
!        sr=se-36.1136-1.5*log(te)+(13.6*                               
!     c                 rc*rc*(1./float(n*n)-1./rno/rno)*rk/te)         
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)te,n,rno,rk,te,cion                             
        cion=cion/float(n*n) 
!        crec=expo(sr)/float(n*n)                                       
        crec=cion*(2.0779e-16)*                                         &
     &    exp((ep-e)*rk/te)                                             &
     &              /te**(1.5)                                          
!     $    exp(13.6*rc*rc*(1./float(n*n)-1./rno/rno)*rk/te)             
       endif 
        return 
      END                                           
