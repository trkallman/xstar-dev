      subroutine func(lpri,lun11,vturbi,critf,                          &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xiin,rrrts,pirts,htt,cll,httot,cltot,hmctot,elcter,        &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       xilevt,bilevt,rnist)
                                                                        
!                                                                       
!     Name: func.f90  
!     Description:  
!           Master routine which steps throuch elements and ions, 
!           calls routines which calculate rates,
!           calculate ion fractions, level populations, heating cooling
!           this is called from within the heating=cooling loop
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
!           rrrts(nni): total recombination rates for each ion (s^-1)
!           pirts(nni): total photoionization rates for each ion(s^-1)
!           htt(nni): total heating rate for each ion (approximate) 
!                       (erg s^-1 cm^-3)
!           cll(nni): total cooling rate for each ion (approximate) 
!           httot: total heating rate (erg s^-1 cm^-3) 
!           cltot: total cooling rate (erg s^-1 cm^-3) 
!           hmctot:  (httot-cltot)*2./(httot+cltot)
!           elcter:  charge conservation error (relative to H)
!           cllines:  total cooling rate due to lines (erg s^-1 cm^-3) 
!           clcont:  total cooling rate due to continuum (erg s^-1 cm^-3) 
!           cllines:  total cooling rate due to lines (erg s^-1 cm^-3) 
!           htcomp:  compton heating rate (erg s^-1 cm^-3) 
!           clcomp:  compton cooling rate (erg s^-1 cm^-3) 
!           clbrems:  bremsstrahlung cooling rate (erg s^-1 cm^-3) 
!           xilevt(nnml):  level populations (relative to parent element)
!           bilevt(nnml):  departure coefficients for levels
!           rnist(nnml): lte level populations
!           also uses variables from globaldata
!           
!        Dependencies:  Calls func1,func2,func2i,func2l,istruc,
!                   msolvelucy,chisq,comp2,bremem,heatf
!
!     this routine steps through data and calculates                    
!     new version attempts to avoid rates for unabundant ions           
!     author: T. Kallman                                                
!                                                                       
!     with data structures designed for Lucy's iterative method         
!       nsup is a pointer from level n to superlevel N                  
!                                                                       
!     no longer calls full func3 in main loop.                          
!     func3 calls moved to funcsyn                                      
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     level populations                                                 
      real(8) xilevt(nnml),bilevt(nnml),rnist(nnml)
!     ion abundances                                                    
      real(8) xiin(nni) 
      real(8) rrrts(nni),pirts(nni) 
      real(8) tauc(2,nnml) 
      real(8) htt(nni),cll(nni) 
!     element abundances                                                
      real(8) abel(nl) 
      real(8) brcems(ncn)
      real(8) opakc(ncn)
!                                                                       
      real(8) rrrt(31),pirt(31),xin(31),xitmp(31) 
      real(8) ajisb(2,ndb),cjisb(ndb) 
      integer indb(2,ndb) 
      real(8) xilev(nd),rniss(nd),rnisse(nd) 
      real(8) rrcor(nni),pirtt(31) 
      real(8) bmat(nd),bmatl(nd) 
      real(8) rnisl(nd) 
      real(8) x(nd) 
      integer ipsv(31),nsup(nd) 
      integer lpri,lun11,lcdd,ncn2,np2,ncsvn,nmat,nlsvn 
      real(8) vturbi,critf,t,trad,r,delr,xee,xpx,cfrac,p,                &
     &     hmctot,elcter,cllines,clcont,htcomp,clcomp,clbrems           
      real(8) xh1,xh0,httot,cltot 
      real(8) rnisum,cltmp,cmp1,cmp2,                                    &
     &     enelec,httmp,pirtsum,rniss2,rrrtt,                           &
     &     rtdm,tt1,tt2,xintp,xeltp,ximax,xintp2,                       &
     &     xisum,xipp,cl,ht                                             
      integer nlev,nindb,                                               &
     &     jkk,ipmat,ltyp,ldir,llp,imax,ilimh,                          &
     &     lrtyp,lcon,nrdt,nidt,nkdt,ll,jkkl,ipmatsv,                   &
     &     iliml,jk,kl1,mm,kl,kl2,klion,klel,klp,llm,lp,lm,             &
     &     lprim,lprif,lpril,lpritp,lsum,lsumt,ndtmp,mlel,              &
     &     ml1,mmt,mllel,mlion,mleltp,mmtmp,nit,nit2,nit3,nitmx,        &
     &     nitmx2,nlevm,nnz,nnzp,nsp,mlm,np1i,np1r,np1k,lprisv                 
!                                                                       
!                                                                       
      lprisv=lpri 
      lprif=lpri 
      lpritp=0 
!      if (lpri.ge.1) lpritp=2                                          
      if (lprif.ne.0)                                                   &
     &  write (lun11,*)'in func, inputs:',t,                            &
     &         xee,xpx,lcdd,p,abel(1),delr                              
       if (lcdd.ne.1)                                                   &
     &   xpx = p/1.38e-12/max(t,1.e-24)                                 
!                                                                       

      xh0=xpx*xiin(1)*abel(1) 
      xh1=xpx*(1.-xiin(1))*abel(1) 
!                                                                       
!      zero emissivitiesd and opacities                                 
       do ll=1,nni 
         htt(ll)=0. 
         cll(ll)=0. 
         xiin(ll)=0. 
         enddo 
       do ll=1,29 
         rrrt(ll)=0. 
         pirt(ll)=0. 
         enddo 
       elcter=0. 
       httot=0. 
       cltot=0. 
       clcont=0. 
       cllines=0. 
!                                                                       
!                                                                       
!      now calculate.  first step thru elements                         
       nmat=nds
       jkk=0 
       jkkl=0 
       klel=11 
       mlel=derivedpointers%npfirst(klel) 
       jk=0 
       httot=0. 
       cltot=0. 
       do while (mlel.ne.0) 
         mlm=mlel-1 
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         mllel=0 
         if (nidt.gt.0) then 
           mllel=masterdata%idat1(np1i) 
           jk=mllel 
           nnz=masterdata%idat1(np1i) 
           nnzp=nnz+1 
           xeltp=0. 
           if (jk.gt.0) xeltp=abel(jk) 
           if (lprif.ne.0) then 
              write (lun11,9339)                                        &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,12),mllel,xeltp      
             endif 
 9339      format (1x, ' element:',12a1,i4,1pe11.3) 
           if (xeltp.gt.1.e-24) then 
!
!            now step thru ions first pass: func1                       
             if (lprif.ne.0) write (lun11,*)' first pass' 
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
                 call func2l(jkk,lpritp,lun11,t,xee,xpx,                &
     &              rniss,rnisse,nlev)
                 if (lprif.ne.0) write (lun11,9338)                     &
     &             (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)        
 9338            format (1x, ' ion:',8a1) 
                 if (lprif.ne.0) write (lun11,9328) 
                 call func1(jkk,kl,nnz,lpri,lun11,vturbi,               &
     &              t,trad,r,delr,xee,xpx,xh1,xh0,                      &
     &              epi,ncn2,bremsa,bremsint,                           &
     &              np2,ncsvn,nlsvn,                                    &
     &              rniss,rnisse,nlev,pirtt,rrrtt) 
                 rrcor(jkk)=1. 
                 pirtsum=0. 
                 do mm=kl,nnzp 
                   pirtsum=pirtsum+pirtt(mm) 
                   enddo 
                 pirt(kl)=pirtsum 
                 rrrt(kl)=rrrtt 
                 endif 
               mlion=derivedpointers%npnxt(mlion) 
               enddo 
!                                                                       
!            do ion balance                                             
             kl2=kl 
             call istruc(pirt,rrrt,xin,kl2,lpritp,lun11) 
             xisum=0. 
             ldir=+1 
             do mm=1,kl 
                 kl1=mm 
                 xiin(jkk-kl+kl1)=xin(kl1) 
                 xisum=xisum+xin(kl1) 
                 enddo 
             xipp=1.-xisum 
             klp=kl+1 
             xin(klp)=xipp 
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
             if (lprif.ne.0) then 
                 write (lun11,*)'ion fractions:',iliml,ilimh,lsum 
                 write (lun11,*)'ion, pi rate,    rec rate,   fraction' 
                 do mm=1,kl 
                   write (lun11,9023)mm,pirt(mm),rrrt(mm),xin(mm) 
 9023              format (1x,i4,3(1pe10.2)) 
                   enddo 
                 endif 
!                                                                       
!            now step thru ions for second pass                         
             if (lprif.ne.0) write (lun11,*)' second pass' 
 9328        format ('     ion    process     d1    d2 ',               &
     &         '    rec use    ans1      ans2  ',                       &
     &     '   ionization  recombination')                              
 9329        format ('     ion    process     d1    d2 ',               &
     &  '             rec use     ans1      ans2     ans3       ans4  ',&
     &         '  aji(lo,up) aji(up,lo)',                               &
     &         ' aji(lo,lo) aji(up,up)')                                
 9330        format ('     ion    process     d1    d2 ',               &
     &   '             rec use     ans1      ans2     ans3       ans4 ',&
     &      '  emiss.    opac.     energy   rec   heat      cool    ')  
!                                                                       
             if (lprif.ne.0) write (lun11,*)'zeroing:',nmat 
             do ml1=1,nmat 
               bmat(ml1)=0. 
               nsup(ml1)=0 
               enddo 
!            a new treatment
             nmat=0                                                                       
             klion=12 
             mlion=derivedpointers%npfirst(klion) 
             jkk=0 
             kl=0 
             nindb=0 
             ipmat=0 
             nsp=1 
             do while ((mlion.ne.0).and.(kl.lt.nnz)) 
               jkk=jkk+1 
               mlm=mlion-1 
               call drd(ltyp,lrtyp,lcon,                                &
     &           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                     &
     &           0,lun11)                                         
               mleltp=derivedpointers%npar(mlion) 
               if (mleltp.eq.mlel) then 
                 kl=kl+1 
                 ipsv(kl)=-1 
                 if (lprif.ne.0)write (lun11,9338)                      &
     &             (masterdata%kdat1(np1k-1+mm),mm=1,nkdt)       
                 if (lprif.ne.0) write (lun11,9329) 
                 if ((kl.ge.iliml).and.(kl.le.ilimh)) then 
                   call func2l(jkk,lpritp,lun11,t,xee,xpx,              &
     &              rniss,rnisse,nlev)
                   if (lprif.ne.0) write (lun11,*)'ipmat=',ipmat 
                   if (lprif.ne.0) write (lun11,*)'before func2',nindb 
                   ipsv(kl)=ipmat 
                   rrrtt=rrrt(kl) 
                   call func2(jkk,kl,ilimh,lpri,lun11,vturbi,           &
     &                   t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,           &
     &                   epi,ncn2,bremsa,bremsint,                      &
     &                   tau0,tauc,                                     &
     &                   np2,ncsvn,nlsvn,                               &
     &                   rniss,rnisse,nmat,nlev,                        &
     &                   ajisb,cjisb,indb,rrrtt,ipmat,nindb)
!                   call func2a(jkk,kl,ilimh,                           
!       $                 lpri,lun11,lfpi,vturbi,                       
!       $                 t,trad,r,xee,xpx,xh1,xh0,                     
!       $                 epi,ncn2,bremsa,bremsint,                     
!       $                 np2,                  
!       $                 npar,npnxt,npfi,npfirst,                      
!       $                 nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi
!       $                 npconi2,ncsvn,                                
!       $                 rlev,ilev,                                    
!       $                 nlpt,iltp,nlev,klev,                          
!       $                 ajis,cjis,ipmat,ipsv)                         
!                                                                       
!                  condense the reaction matrix to omit levels without t
!                    some definition:                                   
!                       nlev is the number of levels in the ion (includi
!                       nmat is the maximum index of levels linked by th
!                       ipmat is the base index into the master arrays f
!                          i.e. the array elements for  this ion start a
!                          and go to ipmat+nmat                         
!                                                                       
                   if (lpri.ge.2) write (lun11,*)'level mapping:',      &
     &                ipmat,nmat                                        
                   do mm=1,nmat 
                     if (mm.le.nlev) then 
                       mmtmp=derivedpointers%npilev(mm,jkk) 
                       if (mmtmp.gt.0) then 
                         x(mm+ipmat)=xilevt(mmtmp) 
                         if (lpri.ge.2)                                 &
     &                    write (lun11,*)mm,mmtmp,mm+ipmat,xilevt(mmtmp)
                      endif 
                       endif 
                     enddo 
!                                                                       
!                  set up superlevel pointers                           
                   nsup(1+ipmat)=nsp 
                   nsp=nsp+1 
                   do mm=2,nlev-1 
                     nsup(mm+ipmat)=nsp 
                     enddo 
                   nsp=nsp+1 
                   ipmat=ipmat+nlev-1 
                   if (ipmat.gt.nd) stop 'ipmat too large.              &
     &                                    Increase critf value.'        
                   endif 
                 endif 
               mlion=derivedpointers%npnxt(mlion) 
               enddo 
!                                                                       
!                                                                       
!                                                                       
             nsup(ipmat+1)=nsp 
             x(ipmat+1)=0. 
             ipmat=ipmat+1 
             nmat=ipmat 
!                                                                       
!            now calculate populations:  second pass, full list         
             call remtms(tt1) 
             nitmx=400 
             nitmx2=400 
             lprim=0 
!             if (lpri.ne.0) lprim=4                                    
             if (lpri.ne.0)                                             &
     &         write (lun11,*)'before msolvelucy',ipmat,nindb                 
             call msolvelucy(ajisb,cjisb,indb,nindb,nsup,nsp,ipmat,     &
     &          bmat,x,ht,cl,nit,nit2,nit3,nitmx,nitmx2,lun11,lprim)    
             if (lpri.ge.1)                                             &
     &         call chisq(ajisb,cjisb,indb,nindb,                       &
     &                    ipmat,x,lun11,lpri)                           
             call remtms(tt2) 
             if (lpri.gt.0)                                             &
     &        write (lun11,981)abs(tt2-tt1),nit,nit2,nit3,ht,cl         
  981        format (1x,'after msolvelucy',(1pe11.3),3i4,2(1pe11.3)) 
             do mm=1,ipmat 
               bmatl(mm)=x(mm) 
               enddo 
             cltot=cltot+cl*xeltp 
             httot=httot+ht*xeltp 
             htt(jk)=ht*xeltp 
             cll(jk)=cl*xeltp 
!                                                                       
             xisum=0. 
             do kl1=1,kl 
               xisum=xisum+xin(kl1) 
               enelec=float(kl1-1) 
               elcter=elcter+xin(kl1)*enelec*xeltp 
               enddo 
             enelec=float(kl) 
             elcter=elcter+max(0.,1.-xisum)*enelec*xeltp 
!                                                                       
             endif 
                                                                        
           endif 
!                                                                       
         if  (mlel.ne.0) mlel=derivedpointers%npnxt(mlel) 
!                                                                       
         enddo 
!                                                                       
      lpril=lpri
      call comp2(lpril,lun11,epi,ncn2,bremsa,t,cmp1,cmp2)                               
!     nonrelativistic compton                                           
!     call comp(lpri,lun11,epi,ncn2,bremsa,cmp1,cmp2)                 
!     ferland compton                                                   
!      call comp3(lpri,lun11,epi,ncn2,bremsa,cmp1,cmp2)               
!      call freef(lpri,lun11,epi,ncn2,t,xpx,xee,opakc)                  
      call bremem(lpril,lun11,xee,xpx,t,epi,ncn2,brcems,opakc) 
      call heatf(jkk,lpri,lun11,                                        &
     &       t,r,cfrac,delr,xee,xpx,abel,                               &
     &       epi,ncn2,bremsa,                                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       nlev,                                                      &
     &       brcems,cmp1,cmp2,httot,cltot,hmctot,                       &
     &             cllines,clcont,htcomp,clcomp,clbrems)                
!                                                                       
       elcter=xee-elcter 
!                                                                       
      if (lprif.ne.0) write (lun11,*)'leaving func' 
!                                                                       
      lprisv=lpri 
!                                                                       
!                                                                       
      return 
      end                                           
