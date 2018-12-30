      subroutine func3(jkk,jkkl,lpri,lun11,vturbi,                      &
     &       t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,                       &
     &       epi,ncn2,bremsa,bremsint,tau0,tauc,                        &
     &       np2,ncsvn,nlsvn,                                           &
     &       rniss,rnisse,nlev,                                         &
     &       xeltp,rrcor,htt,cll,cllines,clcont,rrrt,                   &
     &       xilev,ipmat,ipmatsv,                                       &
     &       rcem,oplin,rccemis,opakc,opakcont,                         &
     &       cemab,cabab,opakab,fline,flinel)                  
!                                                                       
!     Name: func3.f90  
!     Description:  
!           Calculates rates affecting heating, cooling, emissivities
!           and opacities  for one ion.  
!           Puts results into emissivity and opacity arrays
!
!     List of Parameters:
!           Input:
!           jkk: index of ion in xstar scheme 1=H0, 432=Zn29+
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           vturbi: turbulent speed in km/s
!           t: temperature in 10^4K
!           trad: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!           r:  radius in nebula (cm)
!           delr: thickness of current spatial zone (cm)
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           xh1:  H+ number density (cm^-3)
!           xh0:  neutral H number density (cm^-3)
!           cfrac:  covering fraction (affects line and continuum 
!                forward-backward ratio
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!           tau0(nnnl):  line optical depths
!           tauc(nnml):  rrc optical depths
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           rniss:  lte level populations
!           rnisse:  lte level populations with exponential removed
!           nlev:  number of levels for the ion
!           xeltp:  abundance of current element relative to H=1
!           rrcor:  fudge factor for correcting recombination rates (not used)
!           also uses variables from globaldata
!           Output:
!           jkkl: index into line emissivity array of last line
!           htt:  total heating rate (erg cm^-3 s^-1)
!           cll:  total cooling rate (erg cm^-3 s^-1)
!           cllines:  cooling rate due to lines(erg cm^-3 s^-1)
!           clcont:  cooling rate due to rrcs (erg cm^-3 s^-1)
!           rrrt:  recombination rate coefficient (s^-1)
!           xilev:  level populations relative to element 
!           ipmat:  index of current ion into xilev
!           ipmatsv:  maximum index of current ion into xilev
!           rcem(2,nnnl):  line emissivities  (erg cm^-3 s^-1) /10^38
!                  inward and outward
!           oplin(nnnl):  line opacities  (cm^-1)
!           rccemis(2,ncn): continuum emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!                  inward and outward
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           opakcont(ncn):  continuum opacities lines excluded (cm^-1)
!           cemab(nnml):  rrc emissivities (erg cm^-3 s^-1) 
!           cabab(nnml):  total energy absorbed by rrc (erg cm^-3 s^-1) 
!           opakab(nnml):  rrc opacities (cm^-1)
!           fline(2,nnnl):  line emissivity (net radiative)
!              (erg cm^-3 s^-1) 
!           flinel(ncn):  line emissivity binned into continuum bins 
!              (erg cm^-3 s^-1 erg^-1)
!           
!        Dependencies:  Calls ucalc,drd
!        Called by: funcsyn
!
!     this routine calculates rates affecting emission and              
!        absorption                                                     
!     author: T. Kallman                                                
!                                                                       
!      note that the abundances are passed in in the array xilev        
!      this array is indexed for the element as a whole                 
!      and for each ion the offset is the index ipmat.                  
!      so that for each ion the levels begin at ipmat+1 ...             
!      the same does not go for rniss, the lte populations,             
!      which are numbered from 1 ...                                    
!     note that in this routine rniss indeces are relative to ground
!                                                                       
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     element abundances                                                
!     state variables                                                   
      real(8) r,t,xpx,delr 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) trad 
      real(8) vturbi,xee 
      integer ncn2,lpri,lun11,np2 
      integer nlsvn,ncsvn 
      character(49) kdesc2 
      real(8) fline(2,nnnl),flinel(ncn) 
!     level populations                                                 
      real(8) xilev(nd) 
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) tauc(2,nnml) 
      character(1) kblnk 
      real(8) tsq,ans1,ans2,xh1,xh0,cfrac 
      real(8) abund1,abund2,ptmp1,ptmp2,ans3,ans4,opakb1,                &
     &     xeltp,rrcor,cllines,clcont,htt,cll                           
      integer idest1,idest2,idest3,idest4 
      real(8) rniss(nd),rnisse(nd)
      real(8) abundtot,rrrt 
      real(8) tau1,tau2,e1,e2,pescl,pescv,                               &
     &     cemtmp1,cemtmp2,czzz,elin,ener,htsum,eth,opakbb,             &
     &     rcemm,rcsum,ergsev                                           
      integer nlev,mltype,ml,mllz,mlrdesc,lpriu,                        &
     &     llo,lup,jkk,ipmat,ltyp,                                      &
     &     lrtyp,lcon,nrdt,nidt,nkdt,kkkl,jkkl,ipmatsv,                 &
     &     lprisv,ml3,mm,nb1,nbinc,mlpar,mlm,lfpi
      integer np1i,np1r,np1k 
!                                                                       
      data kblnk/' '/ 
!                                                                       
      ergsev=1.602197e-12 
      lprisv=lpri 
!                                                                       
      if (lpri.gt.0)                                                    &
     &  write (lun11,*)'in func3, inputs:',t,xee,xpx,delr,ipmat,ipmatsv 
!                                                                       
                                                                        
!                                                                       
!     lfpi mode:  opacities only                                        
      lfpi=3 
!                                                                       
      abundtot=0. 
      rrrt=0. 
      tsq=sqrt(t) 
!                                                                       
      lprisv=lpri 
!                                                                       
      htt=0. 
      cll=0. 
      if (lpri.ne.0) then 
        write (lun11,*)'level populations:' 
        do mm=1,nlev 
          write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20),         &
     &      leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),                  &
     &      xilev(mm+ipmat),rniss(mm)                                   
 9022     format (i4,20a1,4(1pe10.3)) 
          enddo 
        endif 
!                                                                       
!                                                                       
!     now do other  rates                                               
      lpriu=lpri 
      mltype=9 
      mlrdesc=mltype 
      ml=derivedpointers%npfi(mltype,jkk) 
      if (ml.ne.0) then 
        mllz=derivedpointers%npar(ml) 
        mlpar=derivedpointers%npar(ml) 
        do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
          mlm=ml-1 
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          idest1=masterdata%idat1(np1i) 
          idest2=masterdata%idat1(np1i+1) 
          jkkl=derivedpointers%nplini(ml) 
          if ((masterdata%rdat1(np1r).gt.0.01).and.(jkkl.ne.0)          &
     &      .and.(idest1.gt.0).and.(idest2.gt.0).and.                   &
     &      (idest1.lt.nlev).and.(idest2.lt.nlev)) then                 
            e1=leveltemp%rlev(1,idest1) 
            e2=leveltemp%rlev(1,idest2) 
            eth=abs(e2-e1) 
            if (e1.lt.e2) then 
                lup=idest2+ipmat 
                llo=idest1+ipmat 
              else 
                lup=idest1+ipmat 
                llo=idest2+ipmat 
              endif 
            abund1=xilev(llo)*xpx*xeltp 
            abund2=xilev(lup)*xpx*xeltp 
            tau1=tau0(1,jkkl) 
            tau2=tau0(2,jkkl) 
            ptmp1=                                                      &
     &        pescl(tau1)*(1.-cfrac)                                    
            ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac 
!            ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac          
            call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,                  &
     &        ans3,ans4,idest1,idest2,idest3,idest4,                    &
     &        abund1,abund2,ptmp1,ptmp2,xpx,opakb1,                     &
     &        opakc,opakcont,rccemis,fline,lpriu,kdesc2,                &
     &        r,delr,t,trad,tsq,xee,xh1,xh0,                            &
     &        epi,ncn2,bremsa,bremsint,                                 &
     &        rniss,rnisse,nlev,lfpi,lun11,                             &
     &        np2,ncsvn,nlsvn)               
!                                                                       
            cemtmp1=abund2*ptmp1*ans4 
            cemtmp2=abund2*ptmp2*ans4 
            rcsum=ans4*abund2 
            cll=cll+rcsum 
            clcont=clcont+rcsum 
            if ((lpri.ge.1))                                            &
     &        write (lun11,9002)jkk,lrtyp,ltyp,idest1,idest2,           &
     &        llo,lup,ml,ans1,ans2,ans3,ans4,                           &
     &        cemtmp1+cemtmp2,opakb1,eth,                               &
     &        jkkl,cll,htt                                              
            endif 
          ml=derivedpointers%npnxt(ml) 
          mlpar=0 
          if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
          enddo 
        endif 
!                                                                       
                                                                        
!     now do other  rates                                               
      mltype=7 
      mlrdesc=mltype 
      ml=derivedpointers%npfi(mltype,jkk) 
      mllz=0 
      if (ml.ne.0) mllz=derivedpointers%npar(ml) 
      mlpar=mllz 
      do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
        mlm=ml-1 
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        idest1=masterdata%idat1(np1i+nidt-2) 
        idest2=nlev+masterdata%idat1(np1i-1+nidt-3)-1 
        kkkl=derivedpointers%npconi2(ml) 
        if ((kkkl.ne.0).and.(kkkl.le.ndat2)                             &
     &     .and.(idest1.gt.0)) then                                     
          llo=idest1+ipmat 
          lup=idest2+ipmat 
          eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
          abund1=xilev(llo)*xeltp 
          abund2=xilev(lup)*xeltp 
          if (lup.gt.ipmatsv) then 
            lup=min(lup,ipmatsv) 
            abund2=0. 
            endif 
          nb1=nbinc(eth,epi,ncn2) 
          if (lup.le.ipmatsv) then
            tau1=tauc(1,kkkl) 
            tau2=tauc(2,kkkl) 
            ptmp1=                                                      &
     &        pescv(tau1)*(1.-cfrac)                                    
            ptmp2=pescv(tau2)*(1.-cfrac)+2.*pescv(tau1+tau2)*cfrac 
            call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,                  &
     &        ans3,ans4,idest1,idest2,idest3,idest4,                    &
     &        abund1,abund2,ptmp1,ptmp2,xpx,opakab(kkkl),               &
     &        opakc,opakcont,rccemis,fline,lpriu,kdesc2,                &
     &        r,delr,t,trad,tsq,xee,xh1,xh0,                            &
     &        epi,ncn2,bremsa,bremsint,                                 &
     &        rniss,rnisse,nlev,lfpi,lun11,                             &
     &        np2,ncsvn,nlsvn)               
            cabab(kkkl)=ans3 
            cabab(kkkl)=cabab(kkkl)*abund1*xpx 
            cemab(1,kkkl)=ptmp1*ans4/(ptmp1+ptmp2) 
            cemab(2,kkkl)=ptmp2*ans4/(ptmp1+ptmp2) 
            cemab(1,kkkl)=cemab(1,kkkl)*abund2*xpx 
            cemab(2,kkkl)=cemab(2,kkkl)*abund2*xpx 
            if (lpri.ge.1)                                              &
     &        write (lun11,9002)jkk,lrtyp,ltyp,                         &
     &          idest1,idest2,llo,lup,ml,ans1,ans2,                     &
     &          ans3,ans4,cemab(1,kkkl)+cemab(2,kkkl),                  &
     &          opakab(kkkl),eth,                                       &
     &          kkkl,cll,htt,ptmp1,ptmp2                                            
 9002         format (1x,7i6,i8,' h-c ',                                &
     &          7(1pe10.3),i6,2(1pe10.3),4(1pe10.3))                    
            endif 
          endif 
        ml=derivedpointers%npnxt(ml) 
        mlpar=0 
        if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
        enddo 
!                                                                       
      mltype=42 
      mlrdesc=mltype 
      ml=derivedpointers%npfi(mltype,jkk) 
      if (ml.ne.0) then 
        mllz=derivedpointers%npar(ml) 
        mlpar=derivedpointers%npar(ml) 
        do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
          mlm=ml-1 
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          idest1=masterdata%idat1(np1i+nidt-2) 
          idest2=masterdata%idat1(np1i+nidt-3) 
          if (idest1.gt.0) then 
            llo=idest1+ipmat 
            lup=idest2+ipmat 
            eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
            abund1=xilev(llo)*xeltp 
            abund2=xilev(lup)*xeltp 
            if (lup.gt.ipmatsv) then 
              abund2=0. 
              lup=ipmatsv 
              endif 
            if ((lup.le.ipmatsv).and.                                   &
     &      (xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then        
              call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,                 &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,                &
     &          ans3,ans4,idest1,idest2,idest3,idest4,                  &
     &          abund1,abund2,ptmp1,ptmp2,xpx,opakbb,                   &
     &          opakc,opakcont,rccemis,fline,lpriu,kdesc2,              &
     &          r,delr,t,trad,tsq,xee,xh1,xh0,                          &
     &          epi,ncn2,bremsa,bremsint,                               &
     &          rniss,rnisse,nlev,lfpi,lun11,                           &
     &          np2,ncsvn,nlsvn)               
              htsum=ans3*xpx*abund1 
              htt=htt+htsum 
              czzz=0. 
              if (lpri.ge.1)                                            &
     &        write (lun11,9002)jkk,lrtyp,ltyp,                         &
     &          idest1,idest2,llo,lup,ml,ans1,ans2,                     &
     &          ans3,ans4,czzz,opakbb,eth,kkkl,cll,htt                  
              endif 
            endif 
          ml=derivedpointers%npnxt(ml) 
          if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
          enddo 
        endif 
!                                                                       
      mltype=1 
      mlrdesc=mltype 
      ml=derivedpointers%npfi(mltype,jkk) 
      mllz=0 
      if (ml.ne.0) mllz=derivedpointers%npar(ml) 
      mlpar=mllz 
      do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
        mlm=ml-1 
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        if (.not.((mlrdesc.eq.1).and.((ltyp.eq.93).or.(ltyp.eq.53)))    &
     &      .and.(.not.((mlrdesc.eq.7)                                  &
     &        .and.(masterdata%idat1(np1i+nidt-2).ne.1))))  then
           if (nidt.gt.3) then 
            idest1=masterdata%idat1(np1i+nidt-2) 
            idest2=nlev+masterdata%idat1(np1i+nidt-4)-1 
            kkkl=derivedpointers%npconi2(ml) 
!            if (lpri.ne.0) write (lun11,*)'kkkl=',kkkl,idest1,ltyp     
            if ((kkkl.ne.0).and.(kkkl.le.ndat2)                         &
     &        .and.(idest1.gt.0)) then                                  
              llo=idest1+ipmat 
              lup=idest2+ipmat 
!              eth=leveltemp%rlev(4,idest1)-rlev(1,idest1)                        
              eth=masterdata%rdat1(np1r) 
              abund1=xilev(llo)*xeltp 
              abund2=xilev(lup)*xeltp 
              if (lup.gt.ipmatsv) then 
                abund2=0. 
                lup=ipmatsv 
                endif 
!              if (lpri.ne.0) write (lun11,*)lup,ipmatsv,llo,xilev(llo),
!     $               xilev(1+ipmat)                                    
              if ((lup.le.ipmatsv).and.                                 &
     &          (xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then    
                tau1=tauc(1,kkkl) 
                tau2=tauc(2,kkkl) 
                ptmp1=                                                  &
     &          pescv(tau1)*(1.-cfrac)                                  
                ptmp2=pescv(tau2)*(1.-cfrac)+2.*pescv(tau1+tau2)*cfrac 
                call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,               &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,              &
     &            ans3,ans4,idest1,idest2,idest3,idest4,                &
     &            abund1,abund2,ptmp1,ptmp2,xpx,opakab(kkkl),           &
     &            opakc,opakcont,rccemis,fline,lpriu,kdesc2,            &
     &            r,delr,t,trad,tsq,xee,xh1,xh0,                        &
     &            epi,ncn2,bremsa,bremsint,                             &
     &            rniss,rnisse,nlev,lfpi,lun11,                         &
     &            np2,ncsvn,nlsvn)               
                htsum=ans3*xpx*abund1 
                rrrt=rrrt+xilev(lup)*(ptmp1+ptmp2)*ans2 
                abundtot=abundtot+xilev(lup) 
                htt=htt+htsum 
                cabab(kkkl)=ans3 
                cabab(kkkl)=cabab(kkkl)*abund1*xpx 
                cemab(1,kkkl)=ptmp1*ans4/(ptmp1+ptmp2) 
                cemab(2,kkkl)=ptmp2*ans4/(ptmp1+ptmp2) 
                cemab(1,kkkl)=cemab(1,kkkl)*abund2*xpx 
                cemab(2,kkkl)=cemab(2,kkkl)*abund2*xpx 
                if (lpri.ge.1)                                          &
     &            write (lun11,9002)jkk,lrtyp,ltyp,                     &
     &            idest1,idest2,llo,lup,ml,ans1,ans2,                   &
     &            ans3,ans4,cemab(1,kkkl)+cemab(2,kkkl),                &
     &            opakab(kkkl),eth,                                     &
     &            kkkl,cll,htt,ptmp1,ptmp2                                          
                endif 
              endif 
            endif 
          endif 
        ml=derivedpointers%npnxt(ml) 
        mlpar=0 
        if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
        enddo 
!                                                                       
      mltype=4 
      mlrdesc=mltype 
      ml=derivedpointers%npfi(mltype,jkk) 
      if (ml.ne.0) then 
        mllz=derivedpointers%npar(ml) 
        mlpar=derivedpointers%npar(ml) 
        do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
             mlm=ml-1 
             call drd(ltyp,lrtyp,lcon,                                  &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                    &
     &            0,lun11)                                        
             idest1=masterdata%idat1(np1i) 
             idest2=masterdata%idat1(np1i+1) 
             jkkl=derivedpointers%nplini(ml) 
             if ((masterdata%rdat1(np1r).gt.0.01).and.(jkkl.ne.0)       &
     &          .and.(idest1.gt.0).and.(idest2.gt.0).and.               &
     &          (idest1.lt.nlev).and.(idest2.lt.nlev)) then             
               e1=leveltemp%rlev(1,idest1) 
               e2=leveltemp%rlev(1,idest2) 
               if (e1.lt.e2) then 
                   lup=idest2+ipmat 
                   llo=idest1+ipmat 
                 else 
                   lup=idest1+ipmat 
                   llo=idest2+ipmat 
                 endif 
!               if (lpri.ne.0) write (lun11,*)idest1,idest2             
                nb1=nbinc(abs(e2-e1),epi,ncn2)
                mm=1
                 abund1=xilev(llo)*xpx*xeltp 
                 abund2=xilev(lup)*xpx*xeltp 
                 ml3=derivedpointers%nplin(jkkl) 
                 tau1=tau0(1,jkkl) 
                 tau2=tau0(2,jkkl) 
                 ptmp1=                                                 &
     &                pescl(tau1)*(1.-cfrac)                            
                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac 
!                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac     
                 lpriu=lpri 
!                we need to call ucalc again because rcem               
!                already has the abundance in from func3p               
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,              &
     &             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,             &
     &             ans3,ans4,idest1,idest2,idest3,idest4,               &
     &             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,                &
     &             opakc,opakcont,rccemis,fline,lpriu,kdesc2,           &
     &             r,delr,t,trad,tsq,xee,xh1,xh0,                       &
     &             epi,ncn2,bremsa,bremsint,                            &
     &             rniss,rnisse,nlev,lfpi,lun11,                        &
     &             np2,ncsvn,nlsvn)               
                 rcem(1,jkkl)=abund2*ans4*ptmp1/(ptmp1+ptmp2) 
                 rcem(2,jkkl)=abund2*ans4*ptmp2/(ptmp1+ptmp2) 
                 ml3=derivedpointers%nplin(jkkl) 
                 oplin(jkkl)=opakb1*abund1 
                 if (ml3.ne.0) then 
                   mlm=ml3-1 
                   call drd(ltyp,lrtyp,lcon,                            &
     &               nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                 &
     &               0,lun11)                                     
                   elin=abs(masterdata%rdat1(np1r)) 
                   ener=12398.41/elin 
                   nb1=nbinc(ener,epi,ncn2) 
                   nb1=max(2,min(ncn2-1,nb1)) 
!                   fline(1,jkkl)=(ans1*xilev(lup)-ans2*xilev(llo))     
!     $                *xpx*xeltp*ener*ergsev*ptmp1                     
!                   fline(2,jkkl)=(ans1*xilev(lup)-ans2*xilev(llo))     
!     $                *xpx*xeltp*ener*ergsev*ptmp2                     
!                   flinel(nb1)=flinel(nb1)+(fline(1,jkkl)+fline(2,jkkl)
!     $               *2./(epi(nb1+1)-epi(nb1-1))/ergsev                
                   cll=cll+rcem(1,jkkl)+rcem(2,jkkl) 
                   htt=htt+abund1*ans3 
                   cllines=cllines+rcem(1,jkkl)+rcem(2,jkkl) 
                   if ((lpri.ge.1))                                     &
     &               write (lun11,9002)jkk,lrtyp,                       &
     &                 ltyp,idest1,idest2,                              &
     &                 llo,lup,ml,ans1,ans2,ans3,ans4,                  &
     &                 rcem(1,jkkl)+rcem(2,jkkl),oplin(jkkl),           &
     &                 masterdata%rdat1(np1r),jkkl,cll,htt              &
     &                  ,ptmp1,ptmp2                                    
!     $                 ,rnrb,rnrb*xpx*xeltp*ener*ergsev                
                 endif 
               endif 
             ml=derivedpointers%npnxt(ml) 
             if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
             enddo 
        endif 
!                                                                       
           mltype=14 
           mlrdesc=mltype 
           ml=derivedpointers%npfi(mltype,jkk) 
           mllz=0 
           if (ml.ne.0)  mllz=derivedpointers%npar(ml) 
           mlpar=mllz 
           do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
             mlm=ml-1 
             call drd(ltyp,lrtyp,lcon,                                  &
     &         nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                       &
     &         0,lun11)                                           
             idest1=masterdata%idat1(np1i-1+nidt-3) 
             idest2=masterdata%idat1(np1i+nidt-3) 
             if ((idest1.gt.0).and.(idest2.gt.0).and.                   &
     &          (idest1.lt.nlev).and.(idest2.lt.nlev)) then             
               e1=leveltemp%rlev(1,idest1) 
               e2=leveltemp%rlev(1,idest2) 
               eth=abs(e2-e1) 
               if (e1.lt.e2) then 
                   lup=idest2+ipmat 
                   llo=idest1+ipmat 
                 else 
                   lup=idest1+ipmat 
                   llo=idest2+ipmat 
                 endif 
               if ((xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24).or.   &
     &             (xilev(lup)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then 
                 abund1=xilev(llo)*xpx*xeltp 
                 abund2=xilev(lup)*xpx*xeltp 
                 tau1=tau0(1,jkkl) 
                 tau2=tau0(2,jkkl) 
                 ptmp1=                                                 &
     &                pescl(tau1)*(1.-cfrac)                            
                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac 
                 lpriu=lpri 
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,              &
     &             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,             &
     &             ans3,ans4,idest1,idest2,idest3,idest4,               &
     &             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,                &
     &             opakc,opakcont,rccemis,fline,lpriu,kdesc2,           &
     &             r,delr,t,trad,tsq,xee,xh1,xh0,                       &
     &             epi,ncn2,bremsa,bremsint,                            &
     &             rniss,rnisse,nlev,lfpi,lun11,                        &
     &             np2,ncsvn,nlsvn)                
                 rcemm=abund2*ans4 
                 rccemis(2,3)=rccemis(2,3)+                             &
     &                  rcemm/(epi(4)-epi(3)+1.e-24)/ergsev/12.56       
!                if (lpri.ne.0) write (lun11,*)jkkl,tau0(1,jkkl),       
                 cll=cll+rcemm 
                 clcont=clcont+rcemm 
                 if ((lpri.ge.1))                                       &
     &              write (lun11,9002)jkk,lrtyp,ltyp,                   &
     &               idest1,idest2,llo,lup,ml,ans1,ans2,                &
     &               ans3,ans4,rcemm,opakb1,eth,                        &
     &               kkkl,cll,htt                                       
                 endif 
               endif 
             ml=derivedpointers%npnxt(ml) 
             mlpar=0 
             if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
             enddo 
!                                                                       
      rrrt=rrrt/max(1.e-24,abundtot) 
      lpri=lprisv 
!                                                                       
!                                                                       
      return 
      end                                           
