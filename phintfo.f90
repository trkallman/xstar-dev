      subroutine phintfo(sigc,ethi,pirt,rrrt,piht,rrcl,piht2,rrcl2,  &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lpri,epi,ncn2,bremsa,t,swrat,xnx,lfast,lun11)  
!                                                                       
!     Name:  phintfo.f90
!     Description:
!       this routine does the integration over the spectrum as required by
!       type 59 data                                                      
!       author T. Kallman
!     Parameters:
!       sigc(ncn)=cross section (Mb)
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
!     dependencies:  none
!     called by:  ucalc
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
      real(8)  bremsa(ncn),epi(ncn) 
      real(8) opakc(ncn),opakcont(ncn) 
      real(8) sigc(ncn) 
      integer lpri,ncn2,lfast,lun11 
      real(8) ethi,pirt,rrrt,piht,rrcl,abund1,abund2,xpx,               &
     &     opakab,t,swrat,xnx,piht2,rrcl2                                           
      real(8) eth,ergsev,bk,tm,bktm,ener,sgtmp,epii,ansar1,optmp,       &
     &     sumr,sumh,sumi,sumc,tempi,enero,sumh2,sumc2,                 &
     &     bremtmp,tempr,tempro,exptst,                                 &
     &     tempi1,tempi2,tempio,atmp2,rnist,tsq,ethsht,                 &
     &     atmp2o,bbnurj,exptmp,deld,tsti,atmp22,atmp22o,               &
     &     bbnu,htmpo,htmp,tmpp,sigth,htsum,eps,tst,expo                
      integer lprisv,nphint,nb1,kl,                                     &
     &     lprie,nskp,nbinc,lrcalc,numcon2                              
!                                                                       
!                                                                       
      data ergsev/1.602197e-12/ 
      data bk/1.38062e-16/ 
!                                                                       
      data eps/1.e-2/ 
!                                                                       
      nb1=nbinc(ethi,epi,ncn2) 
      eth=ethi 
!                                                                       
      tm=t*1.e4 
      bktm=bk*tm/ergsev 
      tsq = sqrt(t) 
      rnist=(5.216e-21)*swrat/t/tsq 
!                                                                       
!     from levwk                                                        
      ethsht=eth/bktm 
      ethsht=max(ethsht,0.d0) 
!                                                                       
!     initialize.                                                       
      pirt =0. 
      rrrt =0. 
!                                                                       
      lprisv=lpri 
!                                                                       
!     continuum                                                         
      tempro = 0. 
      tempr = 0. 
      sumr = 0. 
      sumh = 0. 
      sumh2 = 0. 
      sumc = 0. 
      sumc2 = 0. 
      tempio = 0. 
      tempi=0. 
      sumi=0. 
      htsum=0. 
      htmp=0. 
      sgtmp = 0. 
      numcon2=max(2,ncn2/50) 
      nphint=ncn2-numcon2 
      nphint=max(nphint,nb1+1) 
      if (lpri.ge.1) write (lun11,*)'in phintf:',                       &
     &      eth,rnist,xnx,swrat,t,abund1,abund2,                        &
     &      nphint,nb1                                                  
      ener = epi(nb1) 
      sgtmp=0. 
      bbnu=0. 
      kl=nb1 
      atmp2=0. 
      atmp22=0. 
      tst=1.e+10 
      do while ((kl.le.nphint).and.                                     &
     &            ((lfast.le.2).or.(tst.gt.eps)))                       
            enero=ener 
            ener = epi(kl) 
            epii=ener 
            sgtmp = sigc(kl) 
            if (kl.eq.nb1) sigth=sgtmp 
            bremtmp=bremsa(kl)/(25.3) 
            tempro=tempr 
            tempr=(25.3)*sgtmp*bremtmp/epii 
            deld = ener - enero 
            tst=(tempr+tempro)*deld/2. 
            ansar1=sgtmp 
            sumr = sumr + tst 
            sumh=sumh+(tempr*ener+tempro*enero)*deld*ergsev/2. 
            sumh2=sumh2+(tempr*(ener-eth)+tempro*(enero-eth))           &
     &               *deld*ergsev/2. 
            exptst=max(1.d-36,(epii-eth)/bktm) 
            exptmp=expo(-exptst) 
            bbnurj=(min(epii,2.d+4))**3*(1.571e+22) 
            tempi1=rnist*bbnurj*exptmp*sgtmp/epii 
            tempi2=rnist*bremtmp*exptmp*sgtmp/epii 
            tempi=tempi1+tempi2 
            atmp2o=atmp2 
            atmp2=tempi1*epii 
            atmp22o=atmp22 
            atmp22=tempi1*(epii-eth)
            tsti = (tempi+tempio)*deld/2. 
            sumi = sumi + tsti 
            sumc = sumc+(atmp2+atmp2o)*deld*ergsev/2. 
            sumc2 = sumc2+(atmp22+atmp22o)*deld*ergsev/2. 
!            rctmp1=abund2*ansar2*ptmp1*xpx*xnx                         
!            rctmp2=abund2*ansar2*ptmp2*xpx*xnx                         
!            rccemis(1,kl)=rccemis(1,kl)+rctmp1                         
!            rccemis(2,kl)=rccemis(2,kl)+rctmp2                         
            optmp=abund1*ansar1*xpx 
            if (kl.le.nb1+1) opakab=optmp 
            htmpo=htmp 
            tmpp=optmp 
            htmp=bremsa(kl)*tmpp 
            htsum=htsum+(htmp+htmpo)                                    &
     &               *(epi(kl)-epi(kl-1))*(1.602197e-12)/2.             
            opakc(kl)=opakc(kl)+optmp 
            opakcont(kl)=opakcont(kl)+optmp 
            if (lpri.gt.1) then 
              write (lun11,*)kl,ener,bremtmp,bbnu,bbnurj 
              write (lun11,*) tempr,tempi,rnist,exptmp,bbnurj 
              write (lun11,*) sumr,sumi,sumh,sumc,sgtmp 
              endif 
            tempro = tempr 
            tempio = tempi 
            enero = ener 
            tst=1. 
            lprie=0 
            call enxt(ethi,nb1,lprie,epi,ncn2,t,lfast,lun11,         &
     &                  kl,nskp,nphint,lrcalc)                          
            kl=kl+nskp 
            enddo 
!                                                                       
         pirt = pirt + sumr 
         rrrt = rrrt + xnx*sumi 
         piht = piht + sumh 
         rrcl = rrcl + xnx*sumc 
         piht2 = piht2 + sumh2 
         rrcl2 = rrcl2 + xnx*sumc2 
                                                                        
         if (lpri.ge.1)                                                 &
     &    write (lun11,*)'in phintf:',eth,sigth,                        &
     &     pirt,rrrt,piht,rrcl                                          
         lpri=lprisv 
!                                                                       
      return 
      END                                           
