      subroutine dsec(lnerr,nlim,                                       &
     &       lpri,lppri,lun11,tinf,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,zeta,              &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,elcter, &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       httot2,cltot2,                                             &
     &       xilev,bilev,rnist,                                         &
     &       rcem,oplin,brcems,opakc,cemab,                             &
     &       cabab,opakab)                         
!                                                                       
!     Name: dsec.f90  
!     Description:  
!           this routine solves for thermal equilibrium and charge 
!           conservation  by the   
!          double secant method                                              
!          author:  T. Kallman                                               
!     List of Parameters:
!           Input:
!           nlim=maximum interations allowed
!           lppri=print switch for iteration information
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           vturbi: turbulent speed in km/s
!           critf: threshold value for ion fraction to be included in 
!                   level population calculation
!           trad: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!           r:  radius in nebula (cm)
!           delr: thickness of current spatial zone (cm)
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
!           xee= electron fraction relative to H
!           t= temperature in 10^4K
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
!           lnerr=error flag
!           ntotit=total number of iterations 
!           also uses variables from globaldata
!           
!           
!        Dependencies:  Calls calc_ion_rates,calc_rates_level,
!                   calc_num_level,calc_rates_level_lte,istruc,
!                   msolvelucy,chisq,comp2,bremem,heatf
!        Called by: xstar, dsec
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,nd) 
        integer:: ilev(10,nd),nlpt(nd),iltp(nd) 
        character(1) :: klev(100,nd) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
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
      real(8) brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn)
!     level populations                                                 
      real(8) xilev(nnml),bilev(nnml),rnist(nnml)
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating and cooling                                               
      real(8) htt(nni),cll(nni) 
      real(8) htt2(nni),cll2(nni) 
      real(8) rrrt(nni),pirt(nni) 
!     element abundances                                                
      real(8) abel(nl) 
!     limits on ion indeces vs element
      integer mml(nl),mmu(nl)
!     state variables                                                   
      real(8) p,r,t,xpx,delr 
!     heating-cooling variables                                         
      real(8) httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,          &
     &     clcont,hmctot,httot2,cltot2
!     input parameters                                                  
      real(8) trad,tinf 
      real(8) cfrac,critf,vturbi,xee,zeta
      integer lcdd,ncn2,lpri,lun11,np2,nlim 
!     variables associated with thermal equilibrium solution            
      integer ntotit 
!     temporary for xwrite                                              
      character(133) tmpst 
      integer nlsvn,ncsvn 
!                                                                       
!     local variables                                                   
      integer nnt,nntt,lnerr,lppri0,lppri,nlimt,nlimx,nnxx,             &
     &        nlimtt,nlimxx,iht,ilt,iuht,iult,ihx,ilx,nnx               
      real(8) crite,crith,critt,fact,facx,epst,epsx,epstt,to,            &
     &     tl,th,xeel,xeeh,elctrl,elctrh,hmctth,hmcttl,tst,             &
     &     testt                                                        
!                                                                       
!                                                                       
      if (lpri.gt.0) write (lun11,*)'in dsec' 
!                                                                       
      crite=1.e-03 
!      crite=1.e-06                                                     
      crith=1.e-02 
!      crith=5.e-03                                                     
!      crith=1.e-05                                                     
      critt=2.e-09 
!                                                                       
      ntotit=0 
      nnt = 0 
      nntt=0 
      lnerr = 0 
      lppri0 = lppri 
      nlimt =max(nlim,0) 
      nlimx=abs(nlim) 
      nlimtt=max0(nlimt,1) 
      nlimxx=max0(nlimx,1) 
      if (lpri.gt.0)                                                    &
     & write (lun11,*)'nlimtt,nlimxx,lppri--',nlimtt,nlimxx,lppri       
      fact = 1.2 
      facx = 1.2 
      epst = crith 
      epsx = crite 
      epstt = critt 
      to = 1.e+30 
      tl = 0. 
      th = 0. 
      xeel = 0. 
      xeeh = 1. 
      elctrl = 1. 
      elctrh = -1. 
      hmctth = 0. 
      hmcttl = 0. 
!                                                                       
      iht = 0 
      ilt = 0 
      iuht = 0 
      iult = 0 
!                                                                       
  100 nnx = 0 
      t=max(t,tinf) 
      if (t.lt.tinf*1.01) then 
          nlimt=0 
          nlimtt=0 
          nlimx=0 
          nlimxx=0 
        else 
          nlimt =max(nlim,0) 
          nlimx=abs(nlim) 
          nlimxx=nlimx 
        endif 
!      if (t.lt.tinf) return                                            
      nnxx=0 
      ihx = 0 
      ilx = 0 
  200 continue 
      if ( lppri.ne.0 ) then 
        write (lun11,99001)                                             &
     &   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,                        &
     &   nnt,t,tl,th,hmctot,hmcttl,hmctth                               
        write (tmpst,99001)                                             &
     &   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,                        &
     &   nnt,t,tl,th,hmctot,hmcttl,hmctth                               
        call xwrite(tmpst,10) 
        endif 
      call calc_hmc_all(lpri,lun11,vturbi,critf,                        &
     &       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,zeta,              &
     &       mml,mmu,                                                   &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       leveltemp,                                                 &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,elcter, &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       httot2,cltot2,                                             &
     &       xilev,bilev,rnist)
      if ( lppri.ne.0 ) then 
        write (lun11,99001)                                             &
     &   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,                        &
     &   nnt,t,tl,th,hmctot,hmcttl,hmctth                               
        write (tmpst,99001)                                             &
     &   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,                        &
     &   nnt,t,tl,th,hmctot,hmcttl,hmctth                               
        call xwrite(tmpst,10) 
99001   format (' in dsec -- ',i4,6(1pe9.2),i4,6(1pe9.2)) 
        endif 
      ntotit=ntotit+1 
      nnx = nnx + 1 
      nnxx=nnxx+1 
      if (nnxx.ge.nlimxx) go to 300 
      tst=abs(elcter)/max(1.d-48,xee) 
      if (tst.lt.epsx) go to 300 
      if ( elcter.lt.0 ) then 
            ihx = 1 
            xeeh = xee 
            elctrh = elcter 
            if ( ilx.ne.1 ) then 
               xee = xee*facx 
               goto 200 
               endif 
         else 
            ilx = 1 
            xeel = xee 
            elctrl = elcter 
            if ( ihx.ne.1 ) then 
               xee = xee/facx 
               goto 200 
               endif 
         endif 
         xee = (xeel*elctrh-xeeh*elctrl)/(elctrh-elctrl) 
         goto 200 
!                                                                       
!                                                                       
  300 continue 
      nntt=nntt+1 
      nnt = nnt + 1 
      if ( abs(hmctot).le.epst )  goto 500 
      if (nntt.ge.nlimtt) go to 500 
      if ( nnt.lt.nlimt ) then 
         if ( hmctot.lt.0 ) then 
            iht = 1 
            th = t 
            hmctth = hmctot 
            iuht = 1 
            if ( iult.eq.0 ) hmcttl = hmcttl/2. 
            iult = 0 
            if ( ilt.ne.1 ) then 
               t = t/fact
!              doubling the step for far from equilibrium
               if (abs(hmctot).gt.0.9) t=t/fact 
               goto 100 
            endif 
         else 
            ilt = 1 
            tl = t 
            hmcttl = hmctot 
            iult = 1 
            if ( iuht.eq.0 ) hmctth = hmctth/2. 
            iuht = 0 
            if ( iht.ne.1 ) then 
               t = t*fact 
!              doubling the step for far from equilibrium
               if (abs(hmctot).gt.0.9) t=t*fact 
               goto 100 
            endif 
         endif 
         testt = abs(1.-t/to) 
         if ( testt.lt.epstt ) then 
            lnerr = -2 
            if ( lppri.ne.0 ) then 
               write (lun11,99004) 
               write (lun11,99006) nnt,t,tl,th,hmctot,hmcttl,           &
     &                         hmctth                                   
            endif 
            goto 500 
         else 
            to = t 
            t = (tl*hmctth-th*hmcttl)/(hmctth-hmcttl) 
            goto 100 
         endif 
      endif 
!                                                                       
      lnerr = 2 
      write (lun11,99002) 
      write (lun11,99006) nnt,t,tl,th,hmctot,hmcttl,hmctth 
!                                                                       
  500 if ( lppri.ne.0 ) write (lun11,99007) testt,epst,hmctot 
      lppri = lppri0 
!                                                                       
      return 
99002 format (' ','**** note: in dsec --  too many iterations **** ') 
99004 format (' ',' warrning -- dsec not converging ') 
99006 format (' ',' temperature ',i4,6(1pe16.8)) 
99007 format (' ',' finishing dsec -- test,epst,hmctot',3(1pe16.8)) 
      end                                           
