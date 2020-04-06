      subroutine phint53(stmpp,etmpp,ntmp,ethi,pirt,rrrt,piht,rrcl,  &
     & piht2,rrcl2,                                                     &
     & abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,                      &
     & opakc,opakcont,rccemis,lpri,epi,ncn2,bremsa,t,swrat,xnx,         &
     & lfast,lun11)                                                     
!                                                                       
!                                                                       
!     Name:  phint53.f90
!     Description:
!       this routine does the integration over the spectrum to calculate 
!       photoionization rates, milne rates, heating and cooling rates
!       and return opacity and emissivity arrays.
!       uses power law piecewise analytic integrals, 
!       author T. Kallman
!     Parameters:
!       stmpp(ntmp)=cross section (Mb)
!       etmpp(ntmp)=energy (Ry above threshold)
!       ntmp=length of stmpp
!       ethi=threshold energy (eV)
!       abund1=lower level population
!       abund2=upper level population
!       lpri=print switch
!       ptmp1=backward rrc escape probability
!       ptmp2=forward rrc escape probability
!       xpx=H nucleus density (cm^-3)
!       rnist=LTE level population ratio
!       epi(ncn2)=continuum bins (eV)
!       ncn2=length of epi
!       bremsa(ncn)=radiation flux (erg s^-1 cm^-2 erg^-1)
!       t=temperature in 10^4 K
!       swrat=statistical weight ratio relative to continuum
!       xnx=electron number sensity (cm^-3)
!       lfast=fast switch, >=2 --> include milne integral
!       lun11=logical unit number
!       Output:
!       pirt=photoionization rate (s^-1)
!       rrrt=recombination rate (s^-1)
!       piht=photoionization heating rate (erg s^-1)
!       rrcl=recombination cooling rate (erg s^-1)
!       opakab=opacity at threshold (cm^-1)
!       opakc(ncn)=opacity (cm^-1)
!       opakcont(ncn)=continuum only opacity (cm^-1)
!       rccemis(2,ncn)=recombination emissivity in energy bins 
!                       (erg cm^-1 erg^-1)
!     dependencies:  none
!     called by:  ucalc
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      integer ntmp 
!                                                                       
!                                                                       
      real(8)  bremsa(ncn),epi(ncn),stmpp(ntmp),etmpp(ntmp) 
      real(8) rccemis(2,ncn),opakc(ncn),opakcont(ncn) 
      real(8) sgbar(ncn)
      integer lpri,ncn2,lfast,lun11 
      real(8) ethi,pirt,rrrt,piht,rrcl,abund1,abund2,ptmp1,ptmp2,xpx,   &
     &     opakab,t,swrat,xnx,piht2,rrcl2
      real(8) eth,ergsev,bk,tm,bktm,ener,sgtmp,epii,sgtp,optmp,         &
     &     sumr,sumh,sumh2,sumi,sumc,sumc2,tempi,                       &
     &     bremtmp,tempr,exptst,                                        &
     &     atmp2,rnist,ethsht,optmp2,                                   &
     &     wwir,wwih,bbnurjp,e2t,e2,e1,e1o,bremtmpp,                    &
     &     epiip,e2o,tempc,tempcp,exptmpp,rctmp2,rctmp1,                &
     &     s2,s2to,s2o,s2t,sgtpp,sum,tempip,temphp,temprp,temph,        &
     &     temphp2,temph2,tempcp2,tempc2,                               &
     &     wwirp,wwicp,wwic,wwii,wwihp,wwiip,enermx,exptsto,expo        
      integer lprisv,numcon2,nphint,nb1,kl,jk,nbn,klmax,nbinc 
!                                                                       
      logical done 
!                                                                       
      data ergsev/1.602197e-12/ 
      save ergsev
      data bk/1.38062e-16/ 
      save bk
!                                                                       
!                                                                       
!     internal to this routine we use the threshold binned              
!      eth=ethi                                                         
      nb1=nbinc(ethi,epi,ncn2) 
      eth=epi(nb1) 
      lprisv=lpri 
!                                                                       
      tm=t*1.e4 
      bktm=bk*tm/ergsev 
!      tsq = sqrt(t) 
!      rnist=(2.61e-21)*swrat/t/tsq                                     
!                                                                       
      ethsht=eth/bktm 
      ethsht=max(ethsht,0.d0) 
!                                                                       
      numcon2=max(2,ncn2/50) 
      nphint=ncn2-numcon2 
!                                                                       
      sumr = 0. 
      sumh = 0. 
      sumh2 = 0. 
      sumc = 0. 
      sumc2 = 0. 
      sumi=0. 
      ener=eth+etmpp(1)*(13.605692) 
      nb1=nbinc(ener,epi,ncn2) 
      do while ((epi(nb1).lt.ener).and.(nb1.lt.nphint)) 
        nb1=nb1+1 
        enddo 
      nb1=nb1-1 
      nb1=max(nb1,1) 
      if (lpri.ge.1) then
        write (lun11,*)'entering phint53:'         
        write (lun11,*)'  eth=',eth
        write (lun11,*)'  xnx=',xnx
        write (lun11,*)'  swrat=',swrat
        write (lun11,*)'  t=',t
        write (lun11,*)'  sigmath=',etmpp(1),stmpp(1),ntmp
        write (lun11,*)'  lfast=',lfast
        write (lun11,*)'  ptmp1,ptmp2=',ptmp1,ptmp2
        write (lun11,*)'  abund1,abund2=',abund1,abund2
        write (lun11,*)'  rnist=',rnist
        write (lun11,*)'  nb1,ener,nphint=',nb1,ener,nphint 
        endif
      if (nb1.ge.nphint) return 
!      if (epi(kl+1).lt.ethi) return                                    
      tempr=0. 
                                                                        
!     step through cross section, map onto kl grid                      
      jk=1 
      enermx=eth+etmpp(ntmp)*(13.605692) 
      nbn=nbinc(enermx,epi,ncn2) 
      nbn=max(nbn,min(nb1+1,ncn2-1)) 
      ener=eth+etmpp(jk)*(13.605692) 
      sgtmp=stmpp(1) 
      sgbar(max(1,nb1-1))=0. 
      kl=nb1 
      jk=1 
      e1=epi(kl) 
      e2=eth+etmpp(jk)*(13.605692) 
      s2=stmpp(jk) 
      if (e1.lt.e2) then 
        kl=kl+1 
        e1=epi(kl) 
        endif 
      e1o=e2 
      sum=0. 
      done=.false. 
      do while (.not.done) 
!       step through where jk grid is finer                             
!        if (lpri.gt.0) write (lun11,*)'mapping loop:',jk,kl,e2,e1      
        do while ((e2.lt.e1).and.(jk.lt.(ntmp-1))) 
          jk=jk+1 
          e2o=e2 
          s2o=s2 
          e2=eth+etmpp(jk)*(13.605692) 
          s2=stmpp(jk) 
          sum=sum+(s2+s2o)*(e2-e2o)/2. 
!          if (lpri.gt.0) write (lun11,*)'jk loop',jk,e2,s2,e2o,s2o,sum 
          enddo 
!       kl bin exceeds jk bin, subtract off extra                       
        sum=sum-(s2+s2o)*(e2-e2o)/2. 
!       now interpolate to find value at kl bin                         
        e2t=e1 
        if (e2-e2o.gt.1.e-8) then 
             s2t=s2o+(s2-s2o)*(e2t-e2o)/(e2-e2o+1.e-24) 
           else 
             s2t=s2o 
           endif 
        s2to=s2t 
!       now update sum at kl bin                                        
        sum=sum+(s2t+s2o)*(e2t-e2o)/2. 
!       save                                                            
        if (abs(e1-e1o).gt.1.d-36) then 
            sgbar(kl)=sum/(e1-e1o) 
          else 
            sgbar(kl)=0. 
          endif 
!        if (lpri.gt.0) write (lun11,*)'saving:',kl,e1,sgbar(kl),sum,s2t
        e1o=e1 
!       increment kl                                                    
        kl=kl+1 
        e1=epi(kl) 
!       step through where kl grid is finer                             
        do while ((e1.lt.e2).and.(kl.lt.ncn2)) 
          e2t=e1 
          if (e2-e2o.gt.1.e-8) then 
              s2t=s2o+(s2-s2o)*(e2t-e2o)/(e2-e2o) 
            else 
              s2t=s2o 
            endif 
          s2to=s2t 
          sum=(s2t+s2to)*(e1-e1o)/2. 
          sgbar(kl)=sum/(e1-e1o) 
!          if (lpri.gt.0) write (lun11,*)'kl loop',kl,e1,               
!     $                                   sgbar(kl),sum,s2t             
          e1o=e1 
          kl=kl+1 
          e1=epi(kl) 
          enddo 
!       update sum for remaining bit                                    
        sum=(s2+s2t)*(e2-e2t)/2. 
!        if (lpri.gt.0) write (lun11,*)'testing for done:',kl,nphint,   
!     $                                                    jk,ntmp      
        if ((kl.gt.nphint-1).or.(jk.ge.ntmp-1))                         &
     &       done=.true.                                                
        enddo 
      klmax=kl-1 
!                                                                       
!                                                                       
!     preliminary setup                                                 
      sgtpp=sgbar(nb1) 
      bremtmpp=bremsa(nb1)/(12.56) 
      epiip=epi(nb1) 
      temprp=(12.56)*sgtpp*bremtmpp/epiip 
      temphp=temprp*epiip 
      temphp2=temprp*(epiip-eth) 
      exptst=(epiip-eth)/bktm 
      exptmpp=expo(-exptst) 
      bbnurjp=(min(2.d+4,epiip))**3*(1.571e+22)*2. 
                                                                        
      tempip=rnist*(bremtmpp+bbnurjp)                                   &
     &  *sgtpp*exptmpp/epiip*(ptmp1+ptmp2)                              
      if (lpri.gt.0) write (lun11,*)tempip,rnist,bremtmpp,bbnurjp,      &
     & sgtpp,exptmpp,ptmp1,ptmp2,exptst                                 
      tempcp=tempip*epiip 
      tempcp2=tempip*(epiip-eth)
!                                                                       
      kl=nb1 
      epii=epi(kl) 
!      if (lpri.gt.0) write (lun11,*)'kl=',kl,klmax,sumh
      rctmp1=0. 
      rctmp2=0. 
      if (lpri.gt.0)                                                    &
     &  write (lun11,*)'  jk,kl,epi(kl),sgtp,bremtmp,tempr,sumr,exptsto,&
     & tempi,sumi,tempip,wwir,,opakab,optmp,optmp2rctmp1,rccemis(1,kl)'
      do while (kl.lt.klmax) 
!                                                                       
!       the basics                                                      
        sgtmp=max(0.d0,sgbar(kl)) 
        sgtp=sgtmp 
        sgtpp=sgbar(kl+1) 
        bremtmp=bremsa(kl)/(12.56) 
        bremtmpp=bremsa(kl+1)/(12.56) 
        epii=epi(kl) 
        epiip=epi(kl+1) 
!                                                                       
!       pi rate                                                         
        tempr=temprp 
        temprp=(12.56)*sgtpp*bremtmpp/epiip 
        wwir=(epiip-epii)/2. 
        wwirp=wwir 
        sumr = sumr + (tempr*wwir+temprp*wwirp) 
!                                                                       
!       heat                                                            
        temph=temphp 
        temph2=temphp2 
        temphp=temprp*epiip
        temphp2=temprp*(epiip-eth)
        wwih=wwir 
        wwihp=wwih 
        sumh = sumh + (temph*wwih+temphp*wwihp) 
        sumh2 = sumh2 + (temph2*wwih+temphp2*wwihp) 
!                                                                       
!       rec                                                             
        exptsto=exptst 
        exptst=(epiip-eth)/bktm 
        exptmpp=0. 
        if ((exptsto.lt.200.).and.(lfast.ge.2)) then 
          bremtmpp=bremsa(kl+1)/(12.56) 
          exptmpp=expo(-exptst) 
          bbnurjp=(min(2.d+4,epiip))**3*(1.571e+22)*2. 
          tempi=tempip 
          tempip=rnist*(bremtmpp+bbnurjp)                               &
     &        *sgtpp*exptmpp*12.56/epiip                                
          atmp2=tempip*epiip 
          tempip=tempip*(ptmp1+ptmp2) 
          if (lpri.ge.1)                                                &
     &       write (lun11,*)'tempip',rnist,bremtmpp,bbnurjp,sgtpp,      &
     &          exptmpp,epiip,ptmp1,ptmp2
          wwii=wwir 
          wwiip=wwir 
          sumi = sumi + (tempi*wwii+tempip*wwiip) 
!                                                                       
!         cool                                                          
          tempc=tempcp 
          tempc2=tempcp2 
          tempcp=tempip*epiip 
          tempcp2=tempip*(epiip-eth)
          wwic=wwir 
          wwicp=wwir 
          sumc = sumc+tempc*wwic+tempcp*wwicp 
          sumc2 = sumc2+tempc2*wwic+tempcp2*wwicp 
!                                                                       
          rctmp1=abund2*atmp2*ptmp1*xpx/12.56 
          rccemis(1,kl)=rccemis(1,kl)+rctmp1 
          rctmp2=abund2*atmp2*ptmp2*xpx/12.56 
          rccemis(2,kl)=rccemis(2,kl)+rctmp2 
!                                                                       
          endif 
!                                                                       
!       emiss and opac                                                  
!       the emission must be fudged to get the right cooling with a     
!         trapezoid integration.                                        
        optmp=abund1*xpx*sgtp 
        opakc(kl)=opakc(kl)+optmp 
        opakcont(kl)=opakcont(kl)+optmp 
!                                                                       
        if (kl.eq.(nb1+2)) then 
          optmp2=rnist*exptmpp*sgtp*abund2*(ptmp1+ptmp2)*xpx 
          opakab=optmp-optmp2 
!          if (lpri.ge.1)                                                &
!     &       write (lun11,*)'otmp,optmp2',optmp,optmp2,exptmpp,         &
!     &       sgtp,abund2,ptmp1,ptmp2,xpx,rnist,exptst                   
          endif 
!                                                                       
!       print                                                           
        if (lpri.ge.1) then 
!          write (lun11,*)jk,kl,epi(kl),kl,nphint 
          write (lun11,901)jk,kl,epi(kl),sgtp,bremtmp,tempr,sumr        &
     &         ,exptsto,tempi,sumi,tempip,wwir,opakab,optmp,optmp2,     &
     &         rctmp1,rccemis(1,kl)       
  901     format(1x,2i6,15(1pe11.3)) 
          endif 
!     $      write (lun11,901)kl,epii,sgtp,bremtmp,                     
!     $         tempr,temprp,wwir,wwirp,sumr,                           
!     $         temph,temphp,wwih,wwihp,sumh,                           
!     $         tempi,tempip,wwii,wwiip,sumi,                           
!     $         tempc,tempcp,wwic,wwicp,sumc                            
!     $         ,rctmp1,rctmp2,exptst                                   
! 901       format(1x,'found something',i6,25(1pe11.3))                 
        if (sgtp.lt.0.) then 
          write (6,*) 'phint error' 
          return 
          endif 
!                                                                       
        kl=kl+1 
        enddo 
!                                                                       
      pirt = pirt + sumr 
      rrrt = rrrt + sumi 
      piht = piht + sumh*ergsev 
      piht2 = piht2 + sumh2*ergsev 
      rrcl = rrcl + sumc*ergsev 
      rrcl2 = rrcl2 + sumc2*ergsev 
                                                                        
      if (lpri.ge.1) write (lun11,*)'in phint53:',eth,pirt,rrrt      &
     &         ,piht,rrcl                                               
      lpri=lprisv 
!                                                                       
      do kl=nb1,nbn 
        sgbar(kl)=0. 
        enddo 
!                                                                       
!                                                                       
      return 
      END                                           
