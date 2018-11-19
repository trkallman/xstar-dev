      subroutine func2(jkk,kl,ilimh,lpriz,lun11,vturbi,                 &
     &       t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,                       &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       rniss,rnisse,nmat,nlev,                                    &
     &       ajisb,cjisb,indb,rrrtot,ipmat,nindb)
!                                                                       
!     this routine calculates rates affecting level populations         
!     author: T. Kallman                                                
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rcdum(2,ncn) 
!     continuum opacities                                               
      real(8) opdum(ncn)
      real(8) cedum(2,nnml),cadum(nnml),opadum(nnml) 
!     line emissivities                                                 
      real(8) rcedum(2,nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     element abundances                                                
!     state variables                                                   
      real(8) r,t,xpx,delr 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee 
      integer ncn2,lpri,lun11,lfpi,np2 
      integer nlsvn,ncsvn 
      character(49) kdesc2 
      real(8) tsq,ans1,ans2,xh1,xh0,cfrac,rrrtot 
      real(8) abund1,abund2,ptmp1,ptmp2,ans3,ans4,opakb1 
      integer idest1,idest2,idest3,idest4,ipmat,nindb 
!     continuum flux                                                    
      real(8) tauc(2,nnml) 
!                                                                       
      character(1) kblnk 
      real(8) rniss(nd),rnisse(nd) 
      real(8) ajisb(2,ndb),cjisb(ndb) 
      integer indb(2,ndb) 
      real(8) rrrt3,tau1,tau2,airtmp,rrrtot2,e1,e2,pescl,pescv,ptmp,eth 
      integer np1i,np1r,np1k 
      integer nlev,lpriz,mltype,ml,mllz,mlrdesc,lpriu,                  &
     &        llo,lup,ilimh,jkk,kl,nmat,ltyp,                           &
     &        lrtyp,lcon,nrdt,nidt,nkdt,kkkl,jkkl,mlpar,mlm
!                                                                       
      data kblnk/' '/ 
!                                                                       
!                                                                       
      lpri=lpriz 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in func2, inputs:',t,                           &
     &         xee,xpx                                                 
!                                                                       
!     lfpi value: photoionization and recombination, no opacities       
      lfpi=2 
!                                                                       
      tsq=sqrt(t) 
!     zero temporaries                                                  
      rrrt3=rrrtot 
      rrrtot=0. 
      rrrtot2=0. 
!                                                                       
!     now find the rates affecting this ion                             
      mltype=0 
      do while (mltype.lt.ntyp) 
        mltype=mltype+1 
        mlrdesc=mltype 
        if (((mlrdesc.lt.10).or.(mlrdesc.gt.13))                        &
     &    .and.(mlrdesc.gt.0)) then                                     
          ml=derivedpointers%npfi(mltype,jkk) 
          if (ml.gt.0) then 
            mllz=derivedpointers%npar(ml) 
!           step thru records of this type                              
            mlpar=derivedpointers%npar(ml) 
            do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
              if (masterdata%nptrs(3,ml).eq.mlrdesc) then 
                mlm=ml-1 
                call drd(ltyp,lrtyp,lcon,                               &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                    &
     &            0,lun11)                                        
!               calculate rates                                         
                lpriu=lpri 
                abund1=0. 
                abund2=0. 
                ptmp2=1. 
                ptmp1=0. 
                if ((lrtyp.eq.7).or.((lrtyp.eq.1).and.(ltyp.ne.53)))    &
     &            then                                                  
                  kkkl=derivedpointers%npconi2(ml) 
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)) then 
                    tau1=tauc(1,kkkl) 
                    tau2=tauc(2,kkkl) 
                    ptmp1=pescv(tau1)*(1.-cfrac) 
                    ptmp2=pescv(tau2)*(1.-cfrac)                        &
     &                  +2.*pescv(tau1+tau2)*cfrac                      
                    ptmp=(ptmp1+ptmp2) 
!                   note that radiative rates and emissivities will     
!                     have the escape probabilities in them             
!                     when calculated in ucalc                          
                    call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,           &
     &                nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,          &
     &                ans3,ans4,idest1,idest2,idest3,idest4,            &
     &                abund1,abund2,ptmp1,ptmp2,xpx,opakb1,             &
     &                opdum,opadum,rcdum,lpriu,kdesc2,                  &
     &                r,delr,t,trad,tsq,xee,xh1,xh0,                    &
     &                epi,ncn2,bremsa,bremsint,                         &
     &                rniss,rnisse,nlev,lfpi,lun11,                     &
     &                np2,ncsvn,nlsvn)               
!                   this Statement prevents double counting of PI from  
!                     excited levels.  All rate type 1 data should not h
!                     lower level associated with it, but some does.    
                    if ((lrtyp.eq.1).and.(idest1.ne.1)) ans1=0. 
                    if (lrtyp.eq.1) then 
                      ans2=0. 
                      ans4=0. 
                      endif 
                    llo=idest1 
                    lup=idest2 
                    if (kl.eq.ilimh) lup=min(nlev,lup) 
                    if ((idest1.gt.0).and.(idest2.gt.0).and.            &
     &                (llo.le.nd).and.(lup.le.nd)) then                 
                      eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
                      airtmp=ans2 
                      rrrt3=rrrt3-airtmp*max(0.,1.-ptmp) 
                      nindb=nindb+1 
                      if (nindb.gt.ndb) stop 'array indexing error' 
                      ajisb(1,nindb)=ans1 
                      ajisb(2,nindb)=airtmp 
                      cjisb(nindb)=0. 
                      indb(1,nindb)=lup+ipmat 
                      indb(2,nindb)=llo+ipmat 
                      if (lpri.gt.1)                                    &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                      nindb=nindb+1 
                      if (nindb.gt.ndb) stop 'array indexing error' 
                      ajisb(1,nindb)=airtmp 
                      ajisb(2,nindb)=ans1 
                      cjisb(nindb)=0. 
                      indb(1,nindb)=llo+ipmat 
                      indb(2,nindb)=lup+ipmat 
                      nindb=nindb+1 
                      if (lpri.gt.1)                                    &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                      if (nindb.gt.ndb) stop 'array indexing error' 
                      ajisb(1,nindb)=-ans1 
                      ajisb(2,nindb)=-ans1 
                      cjisb(nindb)=-ans3*xpx 
                      indb(1,nindb)=llo+ipmat 
                      indb(2,nindb)=llo+ipmat 
                      if (lpri.gt.1)                                    &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                      nindb=nindb+1 
                      if (nindb.gt.ndb) stop 'array indexing error' 
                      ajisb(1,nindb)=-airtmp 
                      ajisb(2,nindb)=-airtmp 
                      cjisb(nindb)=ans4*xpx 
                      indb(1,nindb)=lup+ipmat 
                      indb(2,nindb)=lup+ipmat 
                      if (lpri.gt.1)                                    &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                      nmat=max(nmat,llo) 
                      nmat=max(nmat,lup) 
                      rrrtot=rrrtot+airtmp 
                      if (idest2.le.nlev)                               &
     &                  rrrtot2=rrrtot2+airtmp                          
                      if ((lpri.ge.1))                                  &
     &                  write (lun11,9004)jkk,lrtyp,ltyp,idest1,        &
     &                    idest2,llo,lup,ml,ans1,ans2,ans3,ans4,        &
     &                    cedum(1,kkkl)+cedum(2,kkkl),opadum(kkkl),eth, &
     &                    kkkl,ptmp1,ptmp2                              
 9004                   format(1x,7i6,i8,' level',7(1pe10.3),i6,        &
     &                     7(1pe10.3))                                  
                      endif 
                     endif 
                   endif 
                 if ((lrtyp.eq.5).or.(lrtyp.eq.40)) then 
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,            &
     &               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,           &
     &               ans3,ans4,idest1,idest2,idest3,idest4,             &
     &               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,              &
     &               opdum,opadum,rcdum,lpriu,kdesc2,                   &
     &               r,delr,t,trad,tsq,xee,xh1,xh0,                     &
     &               epi,ncn2,bremsa,bremsint,                          &
     &               rniss,rnisse,nlev,lfpi,lun11,                      &
     &                np2,ncsvn,nlsvn)               
                   llo=idest1 
                   lup=idest2 
                   if (kl.eq.ilimh) lup=min(nlev,lup) 
                   if ((llo.ne.0).and.(lup.ne.0).and.                   &
     &               (llo.le.nd).and.(lup.le.nd)) then                  
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=ans1 
                     ajisb(2,nindb)=ans2 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=lup+ipmat 
                     indb(2,nindb)=llo+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=ans2 
                     ajisb(2,nindb)=ans1 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=lup+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=-ans1 
                     ajisb(2,nindb)=-ans1 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=llo+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=-ans2 
                     ajisb(2,nindb)=-ans2 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=lup+ipmat 
                     indb(2,nindb)=lup+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nmat=max(nmat,llo) 
                     nmat=max(nmat,lup) 
                     if (lpri.ge.1)                                     &
     &                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,idest2,  &
     &                   llo,lup,ml,ans1,ans2,ans3,ans4                 
                     endif 
                   endif 
                 if ((lrtyp.eq.3).or.(lrtyp.eq.23)) then 
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,            &
     &               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,           &
     &               ans3,ans4,idest1,idest2,idest3,idest4,             &
     &               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,              &
     &               opdum,opadum,rcdum,lpriu,kdesc2,                   &
     &               r,delr,t,trad,tsq,xee,xh1,xh0,                     &
     &               epi,ncn2,bremsa,bremsint,                          &
     &               rniss,rnisse,nlev,lfpi,lun11,                      &
     &                np2,ncsvn,nlsvn)               
                   if ((idest1.gt.0).and.(idest2.gt.0).and.             &
     &               (idest1.le.nlev).and.(idest2.le.nlev).and.         &
     &               (idest1.le.nd).and.(idest2.le.nd))                 &
     &               then                                               
                     e1=leveltemp%rlev(1,idest1) 
                     e2=leveltemp%rlev(1,idest2) 
                     if ((e1/(1.e-24+e2)-1.).lt.0.01) then 
                         lup=idest2 
                         llo=idest1 
                       else 
                         lup=idest1 
                         llo=idest2 
                       endif 
                     if (kl.eq.ilimh) lup=min(nlev,lup) 
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=ans2 
                     ajisb(2,nindb)=ans1 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=lup+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=ans1 
                     ajisb(2,nindb)=ans2 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=lup+ipmat 
                     indb(2,nindb)=llo+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=-ans1 
                     ajisb(2,nindb)=-ans1 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=llo+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=-ans2 
                     ajisb(2,nindb)=-ans2 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=lup+ipmat 
                     indb(2,nindb)=lup+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nmat=max(nmat,llo) 
                     nmat=max(nmat,lup) 
                     if (lpri.ge.1)                                     &
     &                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,         &
     &                   idest2,llo,lup,ml,ans1,ans2,ans3,ans4          
                     endif 
                   endif 
                 if ((lrtyp.eq.4).or.(lrtyp.eq.14).or.(lrtyp.eq.9))     &
     &             then                                                 
                   jkkl=derivedpointers%nplini(ml) 
                   if ((jkkl.lt.nnnl).and.(jkkl.gt.0)) then 
                     tau1=tau0(1,jkkl) 
                     tau2=tau0(2,jkkl) 
                     ptmp1=pescl(tau1)*(1.-cfrac) 
                     ptmp2=pescl(tau2)*(1.-cfrac)                       &
     &                      +2.*pescl(tau1+tau2)*cfrac                  
!                     ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac 
!                    note that radiative rates and emissivities will    
!                      have the escape probabilities in them            
!                      when calculated in ucalc                         
                     call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,          &
     &                 nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,         &
     &                 ans3,ans4,idest1,idest2,idest3,idest4,           &
     &                 abund1,abund2,ptmp1,ptmp2,xpx,opakb1,            &
     &                 opdum,opadum,rcdum,lpriu,kdesc2,                 &
     &                 r,delr,t,trad,tsq,xee,xh1,xh0,                   &
     &                 epi,ncn2,bremsa,bremsint,                        &
     &                 rniss,rnisse,nlev,lfpi,lun11,                    &
     &                 np2,ncsvn,nlsvn)               
                     if ((idest1.gt.0).and.(idest2.gt.0).and.           &
     &                 (idest1.le.nlev).and.(idest2.le.nlev).and.       &
     &                 (idest2.le.nd).and.(idest1.le.nd)) then          
                       e1=leveltemp%rlev(1,idest1) 
                       e2=leveltemp%rlev(1,idest2) 
                       if (e1.lt.e2) then 
                           lup=idest2 
                           llo=idest1 
                         else 
                           lup=idest1 
                           llo=idest2 
                         endif 
                       if (kl.eq.ilimh) lup=min(nlev,lup) 
!                      these terms represent radiative exctiation.      
!                       must turn them off to attain LTE.               
                       nindb=nindb+1 
                       if (nindb.gt.ndb) stop 'array indexing error' 
                       ajisb(1,nindb)=ans2 
                       ajisb(2,nindb)=ans1 
                       cjisb(nindb)=0. 
                       indb(1,nindb)=lup+ipmat 
                       indb(2,nindb)=llo+ipmat 
                       if (lpri.gt.1)                                   &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                               ajisb(1,nindb),ajisb(2,nindb)      
                       nindb=nindb+1 
                       if (nindb.gt.ndb) stop 'array indexing error' 
                       ajisb(1,nindb)=ans1 
                       ajisb(2,nindb)=ans2 
                       cjisb(nindb)=0. 
                       indb(1,nindb)=llo+ipmat 
                       indb(2,nindb)=lup+ipmat 
                       if (lpri.gt.1)                                   &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                               ajisb(1,nindb),ajisb(2,nindb)      
                       nindb=nindb+1 
                       if (nindb.gt.ndb) stop 'array indexing error' 
                       ajisb(1,nindb)=-ans2 
                       ajisb(2,nindb)=-ans2 
                       indb(1,nindb)=llo+ipmat 
                       indb(2,nindb)=llo+ipmat 
                       cjisb(nindb)=-ans3*xpx 
                       if (lpri.gt.1)                                   &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                       nindb=nindb+1 
                       if (nindb.gt.ndb) stop 'array indexing error' 
                       ajisb(1,nindb)=-ans1 
                       ajisb(2,nindb)=-ans1 
                       indb(1,nindb)=lup+ipmat 
                       indb(2,nindb)=lup+ipmat 
                       cjisb(nindb)=ans4*xpx 
                       if (lpri.gt.1)                                   &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                       nmat=max(nmat,llo) 
                       nmat=max(nmat,lup) 
                       if ((lpri.ge.1))                                 &
     &                 write (lun11,9009)jkk,lrtyp,ltyp,idest1,         &
     &                   idest2,llo,lup,ml,ans1,ans2,ans3,ans4,         &
     &                   tau1,tau2,jkkl,ptmp1,ptmp2                     
 9009                  format (1x,7i6,i8,' level',6(1pe10.3),           &
     &                       i6,2(1pe10.3))                             
                       endif 
                     endif 
                   endif 
                 if (lrtyp.eq.42) then 
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,            &
     &               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,           &
     &               ans3,ans4,idest1,idest2,idest3,idest4,             &
     &               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,              &
     &               opdum,opadum,rcdum,lpriu,kdesc2,                   &
     &               r,delr,t,trad,tsq,xee,xh1,xh0,                     &
     &               epi,ncn2,bremsa,bremsint,                          &
     &               rniss,rnisse,nlev,lfpi,lun11,                      &
     &                np2,ncsvn,nlsvn)               
                   if ((idest1.gt.0).and.(idest2.gt.0).and.             &
     &               (idest1.le.nlev).and.(idest2.le.nlev).and.         &
     &               (idest2.le.nd).and.(idest1.le.nd))                 &
     &               then                                               
                     e1=leveltemp%rlev(1,idest1) 
                     e2=leveltemp%rlev(1,idest2) 
                     if (e1.lt.e2) then 
                         lup=idest2 
                         llo=idest1 
                       else 
                         lup=idest1 
                         llo=idest2 
                       endif 
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=ans1 
                     ajisb(2,nindb)=0. 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=lup+ipmat 
                     indb(2,nindb)=llo+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=0. 
                     ajisb(2,nindb)=ans1 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=lup+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb)             
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=-ans1 
                     ajisb(2,nindb)=-ans1 
                     cjisb(nindb)=-ans3*xpx 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=llo+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                     nmat=max(nmat,llo) 
                     nmat=max(nmat,lup) 
                     if ((lpri.ge.1))                                   &
     &                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,idest2,  &
     &                 llo,lup,ml,ans1,ans2,ans3,ans4                   
                     endif 
                   endif 
                 if (lrtyp.eq.41) then 
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,            &
     &               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,           &
     &               ans3,ans4,idest1,idest2,idest3,idest4,             &
     &               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,              &
     &               opdum,opadum,rcdum,lpriu,kdesc2,                   &
     &               r,delr,t,trad,tsq,xee,xh1,xh0,                     &
     &               epi,ncn2,bremsa,bremsint,                          &
     &               rniss,rnisse,nlev,lfpi,lun11,                      &
     &                np2,ncsvn,nlsvn)               
                   if ((idest1.gt.0).and.(idest2.gt.0).and.             &
     &               (idest1.le.nlev).and.                              &
     &               (idest2.le.nd).and.(idest1.le.nd))                 &
     &               then                                               
                     lup=idest2 
                     llo=idest1 
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=ans1 
                     ajisb(2,nindb)=0. 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=lup+ipmat 
                     indb(2,nindb)=llo+ipmat 
                    if (lpri.gt.1)                                      &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                               ajisb(1,nindb),ajisb(2,nindb)      
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=0. 
                     ajisb(2,nindb)=ans1 
                     cjisb(nindb)=0. 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=lup+ipmat 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                               ajisb(1,nindb),ajisb(2,nindb)      
                     nindb=nindb+1 
                     if (nindb.gt.ndb) stop 'array indexing error' 
                     ajisb(1,nindb)=-ans1 
                     ajisb(2,nindb)=-ans1 
                     indb(1,nindb)=llo+ipmat 
                     indb(2,nindb)=llo+ipmat 
                     cjisb(nindb)=-ans3*xpx 
                     if (lpri.gt.1)                                     &
     &                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),&
     &                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                     nmat=max(nmat,llo) 
                     nmat=max(nmat,lup) 
                     if ((lpri.ge.1))                                   &
     &                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,idest2,  &
     &                   llo,lup,ml,ans1,ans2,ans3,ans4                 
                     endif 
                   endif 
                 endif 
               ml=derivedpointers%npnxt(ml) 
               if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
               enddo 
            endif 
          endif 
        enddo 
!                                                                       
      if (lpri.ne.0)                                                    &
     &   write (lun11,*)'rrrtt=',rrrtot,rrrtot2                         
!                                                                       
!                                                                       
      return 
      end                                           
