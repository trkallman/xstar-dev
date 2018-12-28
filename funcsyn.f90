      subroutine funcsyn(lpri,lun11,vturbi,critf,                       &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xiin,                                                      &
     &       xilevt,bilevt,rnist,                                       &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                    
                                                                        
!                                                                       
!     Name: funcsyn.f90  
!     Description:  
!           Master routine which steps throuch elements and ions, 
!           calls routines which calculate rates,
!           calculates  emissivities and opacities
!           this routine is called after heating=cooling solution
!
!     List of Parameters:
!           Input:
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           vturbi: turbulent speed in km/s
!           critf: threshold value for ion fraction to be included in 
!                   level population calculation
!           t: temperature in 10^4K
!           trad: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!           r:  radius in nebula (cm)
!           delr: thickness of current spatial zone (cm)
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           abel(nl):  element abundances relative to H=1
!           cfrac:  covering fraction (affects line and continuum 
!                forward-backward ratio
!           p:  pressure in dynes/cm^2
!           lcdd: constant pressure switch, 1=constant pressure 
!                      0=constant density
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
!           Output:
!           xiin(nni):  ion fractions, xiin(1)=H, xiin(2)=He0, xiin(3)=He+ etc
!           xilevt(nnml):  level populations (relative to parent element)
!           bilevt(nnml):  departure coefficients for levels
!           rnist(nnml): lte level populations
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
!           also uses variables from globaldata
!           
!        Dependencies:  Calls func2i,func2l,func3
!
!     calculates opacities and emissivities and does transfer           
!     level populations and integrates continuum emissivities           
!     and opacities are assumed as input                                
!     note that in this routine rniss indeces are relative to ground
!                                                                       
!     author: T. Kallman                                                
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
      real(8) fline(2,nnnl),flinel(ncn) 
!     level populations                                                 
      real(8) xilevt(nnml),bilevt(nnml),rnist(nnml) 
!     ion abundances                                                    
      real(8) xiin(nni) 
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) tauc(2,nnml) 
!     element abundances                                                
      real(8) abel(nl) 
!     the saved rates                                                   
!                                                                       
      real(8) rniss(nd),rnisse(nd) 
      real(8) rrcor(nni),xin(31)
      real(8) bmat(nd),bmatl(nd) 
      integer lpri,lun11,lcdd,ncn2,np2,ncsvn,nlsvn 
      real(8) vturbi,critf,t,trad,r,delr,xee,xpx,cfrac,p
      real(8) xh1,xh0
      real(8) httmp,cltmp,cllines,clcont,                                &
     &     rtdm,xeltp,ximax,                                            &
     &     xisum,xipp
      integer nlev,                                                     &
     &     jkk,ipmat,ltyp,ldir,llp,imax,ilimh,                          &
     &     lrtyp,lcon,nrdt,nidt,nkdt,ll,jkkl,ipmatsv,                   &
     &     iliml,jk,kl1,mm,kl,klion,klel,klp,llm,lp,lm,                 &
     &     lprif,lsum,lsumt,ndtmp,mlel,                                 &
     &     mmt,mllel,mlion,mleltp,mmtmp,                                &
     &     nlevm,nnz,mlm,np1i,np1r,np1k,lprisv                 
!                                                                       
      lprisv=lpri 
      lprif=lpri 
      if (lprif.ne.0)                                                   &
     &  write (lun11,*)'in funcsyn, inputs:',t,                         &
     &         xee,xpx,lcdd,p,abel(1),delr                              
       if (lcdd.ne.1)                                                   &
     &   xpx = p/1.38e-12/max(t,1.e-24)                                 
!                                                                       

      xh0=xpx*xiin(1)*abel(1) 
      xh1=xpx*(1.-xiin(1))*abel(1) 
!                                                                       
!      zero emissivitiesd and opacities                                 
!      note that here variables                                         
!      on the level grid  are                                           
!      already calculated in func3p                                     
       do ll=1,nnml 
         cemab(1,ll)=0.
         cemab(2,ll)=0.
         cabab(ll)=0.
         opakab(ll)=0.
         enddo
       do ll=1,nnnl 
         fline(1,ll)=0. 
         fline(2,ll)=0. 
         rcem(1,ll)=0. 
         rcem(2,ll)=0. 
         oplin(ll)=0. 
         enddo 
       do ll=1,ncn2 
         rccemis(1,ll)=0. 
         rccemis(2,ll)=0. 
         opakc(ll)=0. 
         opakcont(ll)=0. 
         enddo 
       do ll=1,ncn2 
         flinel(ll)=0. 
         enddo 
!                                                                       
!                                                                       
!      now calculate.  first step thru elements                         
       jkk=0 
       jkkl=0 
       klel=11 
       mlel=derivedpointers%npfirst(klel) 
       jk=0 
       do while (mlel.ne.0) 
         mlm=mlel-1 
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         mllel=0 
         if (nidt.gt.0) then 
           if (lprif.ne.0) write (lun11,9339)                           &
     &         (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)           
 9339      format (1x, ' element:',12a1) 
           mllel=masterdata%idat1(np1i) 
           jk=mllel 
           nnz=masterdata%idat1(np1i) 
           xeltp=0. 
           if (jk.gt.0) xeltp=abel(jk) 
           if (xeltp.gt.1.e-24) then 
!                                                                       
!            find ion indeces                                           
             if (lprif.ne.0) write (lun11,*)' finding ion indeces' 
             klion=12 
             mlion=derivedpointers%npfirst(klion) 
             jkk=0 
             kl=0 
             do while ((mlion.ne.0).and.(kl.lt.nnz)) 
               jkk=jkk+1 
               mlm=mlion-1 
               call drd(ltyp,lrtyp,lcon,                                &
     &           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                     &
     &           0,lun11)                                         
               mleltp=derivedpointers%npar(mlion) 
               if (mleltp.eq.mlel) then 
                 kl=kl+1 
                 endif 
               mlion=derivedpointers%npnxt(mlion) 
               enddo 
!                                                                       
!            unpack ion fractions                                       
             xisum=0. 
             ldir=+1 
             do mm=1,kl 
                 kl1=mm 
                 xin(kl1)=xiin(jkk-kl+kl1) 
                 xisum=xisum+xin(kl1) 
                 enddo 
             xipp=1.-xisum 
             klp=kl+1 
             xin(klp)=xipp 
!                                                                       
!            find iliml, ilimh                                          
             iliml=0 
             ilimh=0 
             imax=0 
             ximax=0. 
             do mm=1,klp 
               if (xin(mm).gt.ximax) then 
                 ximax=xin(mm) 
                 imax=mm 
                 endif 
               enddo 
             imax=max(min(nnz,imax),1) 
             llp=imax 
             llm=imax 
             iliml=imax 
             ilimh=imax 
             lp=0 
             lm=0 
             if (imax.ne.klp) then 
                 lsumt=derivedpointers%nlevs(jkk-kl+imax) 
               else 
                 lsumt=0 
               endif 
             ndtmp=nd 
             mmt=0 
             do while ((lsumt.lt.ndtmp).and.(ldir.ne.0)) 
               mmt=mmt+1 
               lsum=lsumt 
               iliml=min(iliml,llm) 
               ilimh=max(ilimh,llp) 
               if ((llp.ge.klp).or.(xin(llp).lt.critf)) then 
                 ldir=-1 
                 lp=1 
                 endif 
               if ((llm.le.1).or.(xin(llm).lt.critf)) then 
                 ldir=+1 
                 lm=1 
                 endif 
               if ((lm.ne.1).and.(lp.ne.1)) then 
                 if (xin(llp+1).gt.xin(llm-1)) then 
                     ldir=+1 
                   else 
                     ldir=-1 
                   endif 
                 endif 
               if ((lp.eq.1).and.(lm.eq.1)) ldir=0 
               if (ldir.eq.+1) then 
                   llp=llp+1 
                   if (llp.ne.klp) then 
                       lsumt=lsum+derivedpointers%nlevs(jkk-kl+llp) 
                     else 
                       lsumt=lsum 
                     endif 
                   endif 
               if (ldir.eq.-1) then 
                   llm=llm-1 
                   lsumt=lsum+derivedpointers%nlevs(jkk-kl+llm) 
                   endif 
               ilimh=max(ilimh-1,iliml+1) 
               enddo 
             if (lpri.ne.0) write (lun11,*)'iliml,ilimh:',iliml,ilimh 
!                                                                       
!            step thru ions                                             
             if (lprif.ne.0) write (lun11,*)' third pass',ipmat 
             mlion=derivedpointers%npfirst(klion) 
             ipmat=ipmat+1 
             ipmat=0 
             jkk=jkk-nnz 
             kl=0 
             do while ((mlion.ne.0).and.(kl.lt.nnz)) 
               ltyp=klion 
               mlm=mlion-1 
               call drd(ltyp,lrtyp,lcon,                                &
     &           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                     &
     &           0,lun11)                                         
               mleltp=derivedpointers%npar(mlion) 
               if (mleltp.eq.mlel) then 
                 kl=kl+1 
                 jkk=jkk+1 
                 call func2i(jkk,                                       &
     &             nlev)            
                 if ((kl.ge.iliml).and.(kl.le.ilimh)) then 
                   if (lprif.ne.0)write (lun11,9338)                    &
     &              (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)      
 9338              format (1x, ' ion:',8a1) 
                   call func2l(jkk,lpri,lun11,t,xee,xpx,                &
     &              rniss,rnisse,nlev)
                   nlevm=nlev-1 
!                  retrieve saved abundances                            
                   do mm=1,nlevm 
                     mmtmp=derivedpointers%npilev(mm,jkk) 
                     if (mmtmp.gt.0) then 
                       if (mmtmp.gt.nnml) stop 'mmtmp error' 
                       bmatl(mm+ipmat)=xilevt(mmtmp) 
                       if (lpri.gt.1)                                   &
     &                  write (lun11,*)mm,mmtmp,xilevt(mmtmp)           
                       endif 
!                                                                       
!                    nb this code makes H and He fully ionized          
!                     if (jkk.le.3) then                                
!                       bmatl(mm+ipmat)=0.                              
!                       endif                                           
!                                                                       
                     enddo 
                   mmtmp=derivedpointers%npilev(nlev,jkk) 
                   if (mmtmp.gt.nnml) stop 'mmtmp error' 
                   bmatl(ipmat+nlev)=xilevt(mmtmp) 
!                                                                       
!                  nb this code makes H and He fully ionized            
!                   if ((jkk.eq.1).or.(jkk.eq.3))                       
!     $               bmatl(nlev+ipmat)=1.                              
!                                                                       
                   rrcor(jkk)=1. 
                   ipmatsv=ipmat+nlev 
                   call func3(jkk,jkkl,lpri,lun11,vturbi,               &
     &                 t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,             &
     &                 epi,ncn2,bremsa,bremsint,tau0,tauc,              &
     &                 np2,ncsvn,nlsvn,                                 &
     &                 rniss,rnisse,nlev,                               &
     &                 xeltp,rrcor(jkk),httmp,cltmp,cllines,clcont,rtdm,&
     &                 bmatl,ipmat,ipmatsv,                             &
     &                 rcem,oplin,rccemis,opakc,opakcont,               &
     &                 cemab,cabab,opakab,fline,flinel)      
!                                                                       
                   ipmat=ipmat+nlev-1 
                   endif 
                 endif 
!                                                                       
               mlion=derivedpointers%npnxt(mlion) 
               enddo 
!                                                                       
             endif 
                                                                        
           endif 
!                                                                       
         if  (mlel.ne.0) mlel=derivedpointers%npnxt(mlel) 
!                                                                       
         enddo 
!                                                                       
!                                                                       
      if (lprif.ne.0) write (lun11,*)'leaving funcsyn' 
!                                                                       
      lprisv=lpri 
!                                                                       
      return 
      end                                           
