      subroutine heatt(lpri,lun11,                                      &
     &       t,r,cfrac,delr,xee,xpx,abel,                               &
     &       epi,ncn2,bremsa,                                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       zrems,zremso,elumab,elumabo,elum,elumo,                    &
     &       rcem,oplin,rccemis,opakc,opakcont,cemab,fline,flinel,      &
     &       brcems)
!                                                                       
!     this routine calculates heating and cooling.                      
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     state variables                                                   
      real(8) r,t,xpx,delr,delrl 
!     heating-cooling variables                                         
!     input parameters                                                  
      real(8) xee 
      integer ncn2,lpri,lun11 
      integer nlsvn,ncsvn 
      real(8) fline(2,nnnl),flinel(ncn) 
!     level populations                                                 
      real(8) abel(nl) 
      real(8) cemab(2,nnml) 
      character(1) kblnk 
      real(8) zrems(5,ncn)
      real(8) zremso(5,ncn)
      real(8) elum(3,nnnl),elumo(3,nnnl) 
      real(8) elumab(2,nnml),elumabo(2,nnml) 
      real(8) xeltp,cllines,clcont,cmp1,cmp2,cltot,                      &
     &     hmctot,htcomp,clcomp,clbrems,etst,ekt,epiio,                 &
     &     epii,fac,fpr2,optpp,optp2,tmpc1,tautmp,cfrac,                &
     &     tmpc2,tmpho,xnx,httot,tmpc,tmph,tmpco,hmctmp                 
      real(8) elin,ener,eth,r19,ergsev,hmctmpo 
      real(8) hpctot,hpctmp,hpctmpo 
      real(8) tmp2,tmp2o 
      integer lskp 
      integer idest1,jk,klel,kl,                                        &
     &     klion,mlleltp,mllel,mlel,mlion,mt2,nilin,nnz,numcon,         &
     &     np2
      integer np1i,np1r,np1k,np1ki 
      integer nlev,nlevmx,mltype,ml,mllz,jkk,ltyp,                      &
     &     lrtyp,lcon,nrdt,nidt,nkdt,lk,kkkl,                           &
     &     lprisv,mm,nbinc,mlpar,mlm                                    
!                                                                       
!                                                                       
      data kblnk/' '/ 
      data ergsev/1.602197e-12/ 
!                                                                       
!                                                                       
      lprisv=lpri 
!                                                                       
      xnx=xpx*xee 
      if (lpri.ge.1) lpri=2 
      fpr2=12.56*r19*r19 
      if (lpri.gt.1) write (lun11,*)'in heatt',httot,cltot,delr,r,fpr2
      if (lpri.gt.1) write (lun11,*)ncsvn,ncn2
      numcon=ncn2 
      r19=r*(1.e-19) 
!                                                                       

!     comment these out to implement scattering                         
      clbrems=0. 
      lskp=1 
      tmp2=0. 
      do kl=1,numcon 
        tmp2o=tmp2 
        tmp2 =brcems(kl) 
        if ( kl.ge.2 ) clbrems=clbrems+(tmp2+tmp2o)                     &
     &                *(epi(kl)-epi(kl-lskp))*ergsev/2.                 
        enddo 
!                                                                       
!      delrl=max(delr,1.e-24*r)                                         
        delrl=delr 
        httot=0. 
        cltot=0. 
        hmctot=0. 
        hpctot=0. 
        epii=epi(1) 
        hmctmp=0. 
        hpctmp=0. 
        tmpc=0. 
        tmph=0. 
        do kl=1,numcon 
          optpp=opakc(kl) 
          optp2=max(1.e-34,optpp) 
!         for outward only                                              
          epiio=epii 
          epii=epi(kl) 
          tmpho=tmph 
          tautmp=optp2*delrl 
          fac=1. 
          if (tautmp.gt.0.01)                                           &
     &      fac=(1.-exp(-tautmp))/tautmp                                
          tmph=bremsa(kl)*optp2 
          tmpco=tmpc 
!          tmpc1=rccemis(1,kl)                                          
!          tmpc2=rccemis(2,kl)                                          
          tmpc1=rccemis(1,kl)+brcems(kl)*(1.-cfrac)/2. 
          tmpc2=rccemis(2,kl)+brcems(kl)*(1.+cfrac)/2. 
          tmpc=(tmpc1+tmpc2)*12.56 
          hmctmpo=hmctmp 
          hpctmpo=hpctmp 
          hmctmp=(tmph-tmpc)*fac 
          hpctmp=(tmph+tmpc)*fac 
!         testing lte                                                   
!          zrems(1,kl)=12.56*(tmpc1+tmpc2)*fpr2/(1.e-34+optp2)          
!         this is the good expression                                   
          zrems(1,kl)=max(0.,zremso(1,kl)                               &
     &         -(tmph-12.56*(tmpc1+tmpc2)-flinel(kl))*fac*delrl*fpr2    &
     &        )                                                         
          zrems(2,kl)=zremso(2,kl)+12.56*tmpc1*fac*delrl*fpr2 
          zrems(3,kl)=zremso(3,kl)+12.56*tmpc2*fac*delrl*fpr2 
          if (kl.gt.1) then 
            hmctot=hmctot+(hmctmp+hmctmpo)                              &
     &       *(epii-epiio)*ergsev/2.                                    
            hpctot=hpctot+(hpctmp+hpctmpo)                              &
     &       *(epii-epiio)*ergsev/2.                                    
            endif 
          httot=(hmctot+hpctot)/2. 
          cltot=(-hmctot+hpctot)/2. 
          if (lpri.ge.1) write (lun11,9009)kl,epii,optpp,               &
     &     bremsa(kl),tmph,tmpc,flinel(kl),httot,cltot                  &
     &       ,rccemis(1,kl)+rccemis(2,kl),hmctot,tautmp,fac             &
     &       ,brcems(kl),zrems(2,kl),zremso(2,kl)
 9009     format (1x,i6,15(1pe12.5)) 
!         an afterthought:  continuum only
          optpp=opakcont(kl) 
          optp2=max(1.e-34,optpp) 
          tautmp=optp2*delrl 
          fac=1. 
          if (tautmp.gt.0.01)                                           &
     &      fac=(1.-exp(-tautmp))/tautmp                                
          zrems(4,kl)=zremso(4,kl)+12.56*tmpc1*fac*delrl*fpr2 
          zrems(5,kl)=zremso(5,kl)+12.56*tmpc2*fac*delrl*fpr2 
          enddo 
        if (lpri.ge.1)                                                  &
     &   write (lun11,*)'continuum heating, cooling:',                  &
     &       hmctot,hpctot,(hmctot+hpctot)/2.,-(hmctot-hpctot)/2.       
        clcont=cltot 
        do jkk=1,nlsvn 
          jk=jkk 
          ml=derivedpointers%nplin(jk) 
          mlm=ml-1 
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          elin=abs(masterdata%rdat1(np1r)) 
          if ((elin.lt.1.e+8).and.(elin.gt.1.)                          &
     &      .and.(ml.ne.0).and.(lrtyp.eq.4)) then                       
            nilin=derivedpointers%npar(ml) 
            ener=ergsev*(12398.41)/max(elin,1.e-24) 
            tmph=elumo(1,jk)*optp2/fpr2 
            tmpc1=rcem(1,jk) 
            tmpc2=rcem(2,jk) 
            tmpc=tmpc1+tmpc2 
            etst=ener/ergsev 
!           this statement is in v221boxx                               
            optp2=0. 
!            nblin=nbinc(etst,epi,ncn2) 
!            optp2=opakc(nblin)                                         
            tautmp=optp2*delrl 
            fac=1. 
!            if (tautmp.gt.1.e-4)                                       
!     $        fac=(1.-exp(-tautmp))/tautmp                             
!           must remove these statements when radiative excitation is in
            hmctot=hmctot-tmpc*fac 
            hpctot=hpctot+tmpc*fac 
            cltot=cltot+tmpc 
            elum(1,jk)=max(0.,elumo(1,jk)+(-tmph+tmpc1)*fac*delrl*fpr2) 
            tmph=elumo(2,jk)*optp2/fpr2 
            elum(2,jk)=max(0.,elumo(2,jk)+(-tmph+tmpc2)*fac*delrl*fpr2) 
            if ((elum(1,jk).gt.1.e-49).and.(lpri.ge.1))                 &
     &       write (lun11,9019)jk,nilin,etst,                           &
     &        rcem(1,jk),rcem(2,jk),delrl,fpr2,cltot,                   &
     &        elumo(2,jk),elum(2,jk),optp2,tmph                         
!            write (lun11,*)'line sum',jk,nilin,nblin,etst,             
!     $        rcem(2,jk),optp2,fac,delrl,fpr2,tmph,tmpc2,              
!     $        elumo(2,jk),elum(2,jk),(elum(2,jk)-elumo(2,jk))/delrl/fpr
 9019     format (1x,2i6,14(1pe12.4)) 
            endif 
          enddo 
!                                                                       
!       calculate rrc luminosities                                      
!       First look for element data (jk is element index)               
          if (lpri.ge.2)                                                &
     &     write (lun11,*)'rrc print:'                                  
          klel=11 
          mlel=derivedpointers%npfirst(klel) 
          jk=0 
          jkk=0 
          do while (mlel.ne.0) 
            jk=jk+1 
            mt2=mlel-1 
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                        &
     &        0,lun11)                                            
            if (nidt.gt.0) then 
              mllel=masterdata%idat1(np1i+nidt-1) 
              xeltp=masterdata%rdat1(np1r) 
              xeltp=abel(mllel) 
              nnz=masterdata%idat1(np1i) 
              if (lpri.ge.2)                                            &
     &        write (lun11,*)'element:',jk,mlel,mllel,nnz,              &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)                    
!           ignore if the abundance is small                            
            if (xeltp.lt.1.e-10) then 
                jkk=jkk+nnz 
              else 
!               now step thru ions (jkk is ion index)                   
                klion=12 
                mlion=derivedpointers%npfirst(klion) 
                jkk=0 
                kl=0 
                do while ((mlion.ne.0).and.(kl.lt.nnz)) 
                  jkk=jkk+1 
!                 retrieve ion name from kdati                          
                  mlm=mlion-1 
                  call drd(ltyp,lrtyp,lcon,                             &
     &              nrdt,np1r,nidt,np1i,nkdt,np1ki,mlm,                 &
     &              0,lun11)                                      
!                 if not accessing the same element, skip to the next el
                  mlleltp=masterdata%idat1(np1i+nidt-2) 
                  if (mlleltp.eq.mllel) then 
                    kl=kl+1 
                    if (lpri.ge.2)                                      &
     &              write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,       &
     &                (masterdata%kdat1(np1ki+mm-1),mm=1,nkdt)           
!                   now find level data                                 
!                   step thru types                                     
                    nlevmx=0 
                    mltype=13 
                    ml=derivedpointers%npfi(mltype,jkk) 
                    mllz=0 
                    if (ml.ne.0) mllz=derivedpointers%npar(ml) 
!                   step thru records of this type                      
                    mlpar=0 
                    if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                    do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
                      mlm=ml-1 
                      call drd(ltyp,lrtyp,lcon,                         &
     &                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,              &
     &                  0,lun11)                                  
                      nlev=masterdata%idat1(np1i+nidt-2) 
                      nlevmx=max(nlevmx,nlev) 
                      if ((nlev.gt.0).and.(nlev.le.nd)) then 
                        do  lk=1,nrdt 
                          leveltemp%rlev(lk,nlev)                       &
     &                        =masterdata%rdat1(np1r+lk-1) 
                          enddo 
                        do lk=1,nidt 
                          leveltemp%ilev(lk,nlev)                       &
     &                       =masterdata%idat1(np1i+lk-1) 
                          enddo 
                        do lk=1,nkdt 
                          leveltemp%klev(lk,nlev)                       &
     &                       =masterdata%kdat1(np1k+lk-1) 
                          enddo 
                        do lk=nkdt+1,20 
                          leveltemp%klev(lk,nlev)=kblnk 
                          enddo 
                        endif 
                      ml=derivedpointers%npnxt(ml) 
                      if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                      enddo 
                    nlev=nlevmx 
                    mltype=7 
                    ml=derivedpointers%npfi(mltype,jkk) 
                    mllz=0 
                    if (ml.ne.0) mllz=derivedpointers%npar(ml) 
                    mlpar=0 
                    if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                    do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
!                     step thru records of this type                    
                      mlm=ml-1 
                      call drd(ltyp,lrtyp,lcon,                         &
     &                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,              &
     &                  0,lun11)                                  
                      kkkl=derivedpointers%npconi2(ml) 
                      idest1=masterdata%idat1(np1i+nidt-2) 
                      if ((kkkl.gt.0).and.(kkkl.le.ndat2)               &
     &                .and.((cemab(1,kkkl).gt.1.e-36)                   &
     &                .or.(cemab(2,kkkl).gt.1.e-36))) then              
                        eth=leveltemp%rlev(4,idest1)                    &
     &                      -leveltemp%rlev(1,idest1) 
                        tmpc=(cemab(1,kkkl)+cemab(2,kkkl)) 
                        optp2=0. 
                        fac=delrl 
                        tmph=elumabo(1,kkkl)*optp2/fpr2 
                        elumab(1,kkkl)=max(0.,elumabo(1,kkkl)           &
     &                            +(-tmph+tmpc)*fac*fpr2/2.)            
                        elumab(2,kkkl)=max(0.,elumabo(2,kkkl)           &
     &                            +(-tmph+tmpc)*fac*fpr2/2.)            
                        if (lpri.ge.2)                                  &
     &                  write (lun11,981)kkkl,eth,idest1,               &
     &                    cemab(1,kkkl),cemab(2,kkkl),                  &
     &                    elumabo(1,kkkl),elumab(1,kkkl)                
  981                   format (1x,i6,1pe11.3,i6,6(1pe11.3)) 
                        endif 
                      ml=derivedpointers%npnxt(ml) 
                      if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                      enddo 
                    endif 
!                 Go to next ion                                        
                  mlion=derivedpointers%npnxt(mlion) 
                  enddo 
                endif 
              endif 
            mlel=derivedpointers%npnxt(mlel) 
!           Go to next element                                          
            enddo 
        cllines=cltot-clcont 
        if (lpri.ge.1)                                                  &
     &   write (lun11,*)'line cooling',cltot,cllines                    
!                                                                       
!                                                                       
      fac=delrl 
      ekt = t*(0.861707) 
      htcomp = cmp1*xnx*ergsev 
      clcomp = ekt*cmp2*xnx*ergsev 
!      these statements needed to implement scattering                  
!      httot=htcomp*fac*fpr2+httot                                      
!      cltot=clcomp*fac*fpr2+cltot                                      
      hmctot=hmctot+(htcomp-clcomp) 
      hpctot=hpctot+(htcomp+clcomp) 
!      hmctot=hmctot+(htcomp-clcomp-clbrems)                            
!      hpctot=hpctot+(htcomp+clcomp+clbrems)                            
      httot=httot+htcomp 
!      cltot=cltot+clcomp+clbrems                                       
      cltot=cltot+clcomp 
      if (lpri.ge.1) write (lun11,9953)htcomp,clcomp,cmp1,cmp2,         &
     &   clbrems,(hmctot+hpctot)/2.,-(hmctot-hpctot)/2.                 
      hmctot=2.*hmctot/(1.e-37+hpctot) 
      if (lpri.ge.1) write (lun11,*)hmctot 
 9953 format (1h , ' compton heating, cooling=',8e12.4) 
      lpri=lprisv 
!                                                                       
                                                                        
      return 
      end                                           
