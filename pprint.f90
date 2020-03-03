                                                                        
!-----------------------------------------------------------------------
      subroutine pprint(jj,jkstep,                                   &
     & tp,xlum,lwri,lpri,r,t,xpx,p,lcdd,numrec,npass,                   &
     & nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,zremsz,epi,ncn2,            &
     & abel,cfrac,emult,taumax,xeemin,spectype,specfile,specunit,       &
     & kmodelname,nloopctl,nparms,parname,partype,parms,parcomm,        &
     & atcredate,lun11,tinf,xcol,vturbi,critf,radexp,                   &
     & delr,rdel,enlum,xee,ababs,                                       &
     & bremsa,bremsint,tau0,dpthc,tauc,                                 &
     & ncsvn,nlsvn,                                                     &
     & ntotit,lnerrd,                                                   &
     & xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,              &
     & cllines,clcont,htcomp,clcomp,clbrems,                            &
     &       httot2,cltot2,                                             &
     & xilev,bilev,rnist,                                               &
     & rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,           &
     & cabab,elumab,elum,zrems,                                         &
     & zrtmp,zrtmpcol,zrtmph,zrtmpc)
!                                                                       
!                                                                       
!     Name: pprint.f90  
!     Description:  
!           prints out all quantities of interest
!
!     List of Parameters:
!           Input:
!           jj:  index which selects print type
!           jkstep:  spatial step index 
!           tp: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!           xlum: source luminosity integrated from 1-1000 Ry
!               in units of 10^38 erg/s
!           lunlog: logical unit number for printing
!           lwri: write switch
!           lpri: print switch
!           r:  radius in nebula (cm)
!           t: temperature in 10^4K
!           xpx: H number density (cm^-3)
!           p:  pressure in dynes/cm^2
!           lcdd: constant pressure switch, 1=constant pressure 
!                      0=constant density
!           numrec:  total number of spatial zones
!           npass:  number of global passes
!           nnmax:  obsolete
!           nlimd:  maximum number of temperature iterations
!           rmax:  maxumum radius (cm)
!           xpxcol:  column density (cm^-2)
!           xi:  ionization parameter (erg cm s^-1)
!           zeta:  log10(xi)
!           lfix:  step size switch.  obsolete.
!           zremsz:  input spectrum (erg s^-1 erg^-1 /10^38)
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           abel(nl):  element abundances relative to xstar 
!                standard values
!           cfrac:  covering fraction (affects line and continuum 
!                forward-backward ratio
!           emult:  courant condition step multiplier
!           taumax:  maximum optical depth for inclusion in courant condition
!           xeemin:  minimum electron fraction
!           spectrype:  spectrum type
!           specfile:  file name for input spectrum
!           specunit:  units of input spectrum 0=energy, 1=photons 2=log(energy)
!           kmodelname:  model name (char*72)
!           nloopctl:  loop control parameter (for xstar2xspec)
!           nparms: number of input parameters
!           parname(nparms): names of input parameters
!           partype(nparms): types of input parameters
!           parval(nparms): values of input parameters
!           parcomm(nparms): comments of input parameters
!           atcredate:  atomic data file creation date (string length 63)
!           lun11: logical unit number for printing
!           tinf:  temperature lower limit
!           xcol:  maximum allowed column density (cm^-2)
!           vturbi:  turbulent speed (km/s)
!           critf: threshold value for ion fraction to be included in 
!                   level population calculation
!           radexp:  power law for radial dependence of density
!           delr: thickness of current spatial zone (cm)
!           rdel:  distance from illuminated face of cloud (cm)
!           xee: electron fraction relative to H
!           ababs(nl):  element abundances relative to H=1
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!           tau0(2,nnnl):  line optical depths
!           tauc(2,nnml):  rrc optical depths
!           ncsvn:  number of rrcs
!           nlsvn:  number of lines
!           xii(nni):  ion fractions, xiin(1)=H, xiin(2)=He0, xiin(3)=He+ etc
!           rrrt(nni): total recombination rates for each ion (s^-1)
!           pirt(nni): total photoionization rates for each ion(s^-1)
!           htt(nni): total heating rate for each ion (approximate) 
!                       (erg s^-1 cm^-3)
!           cll(nni): total cooling rate for each ion (approximate) 
!           httot: total heating rate (erg s^-1 cm^-3) 
!           cltot: total cooling rate (erg s^-1 cm^-3) 
!           hmctot:  2*(httot-cltot)/(httot+cltot)
!           cllines:  total cooling rate due to lines (erg s^-1 cm^-3) 
!           clcont:  total cooling rate due to continuum (erg s^-1 cm^-3) 
!           htcomp:  compton heating rate (erg s^-1 cm^-3) 
!           clcomp:  compton cooling rate (erg s^-1 cm^-3) 
!           clbrems:  bremsstrahlung cooling rate (erg s^-1 cm^-3) 
!           xilev(nnml):  level populations (relative to all elements)
!           bilev(nnml):  level departure coefficients
!           rnist(nnml):  lte level populations (relative to parent element)
!           rcem(2,nnnl):  line emissivities  (erg cm^-3 s^-1) /10^38
!                  inward and outward
!           oplin(nnnl):  line opacities  (cm^-1)
!           rccemis(2,ncn): continuum emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!                  inward and outward
!           brcems(ncn):  bremsstrahlung emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           opakcont(ncn):  continuum opacities lines excluded (cm^-1)
!           cemab(nnml):  rrc emissivities (erg cm^-3 s^-1) 
!           cabab(nnml):  total energy absorbed by rrc (erg cm^-3 s^-1) 
!           opakab(nnml):  rrc opacities (cm^-1)
!           dpthc(2,ncn):  continuum optical depths in continuum bins
!           dpthcont(2,ncn):  continuum optical depths in continuum bins 
!                          without lines
!           elumab(2,nnml):  rrc luminosities (erg s^-1)/10^38 
!           elumabo(2,nnml):  old rrc luminosities (erg s^-1)/10^38 
!           elum(2,nnnl):  line luminosities (erg/s/10^38)
!           elum(2,nnnl):  old line luminosities (erg/s/10^38)
!           zrems(5,ncn):  radiation field in continuum bins 
!                          (erg/s/erg)/10^38
!     Output:  none
!     Dependencies: drd, calc_rates_level_lte, xwrite
!     Called by:  xstar, xstarsetup
!
!
!     this routine prints                                               
!     author:  T. Kallman                                               
!                                                                       
!     variable categories:                                              
!     step indeces:                                                     
!      jj,jkstep,                                                       
!     input parameters:                                                 
!      tp,xlum,lwri,lpri,xpx,p,lcdd,numrec,npass,                       
!      nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,                            
!      abel,cfrac,emult,taumax,xeemin,                                  
!      tinf,xcol,vturbi,critf,ababs,                                    
!     other input                                                       
!      spectype,specfile,specunit,kmodelname,nloopctl,                  
!      nparms,parname,partype,parms,parcomm,lun11,                      
!     input spectrum                                                    
!      zremsz,epi,ncn2,bremsa,                                          
!     database quantities                                               
!      np2,                                                  
!      npar,npnxt,npfi,npfirst,nplin,nplini,                            
!      npcon,npconi,npilev,npilevi,npconi2,                             
!      nlsvn,ncsvn,                                                     
!     step diagnostics                                                  
!      ntotit,lnerrd,                                                   
!     state variables                                                   
!      r,t,xee,xii,xilev,bilev,rnist,
!     derived state variables                                           
!      delr,rdel,enlum,                                                 
!     rates                                                             
!      rrrt,pirt,htt,cll,httot,cltot,hmctot,                            
!      cllines,clcont,htcomp,clcomp,clbrems,                            
!     emissivities and opacities                                        
!      rcem,oplin,rccemis,brcems,opakc,cemab,opakab,                    
!     optical depths and luminosities                                   
!      tau0,dpthc,tauc,elumab,elum,zrems                                
!     note that in this routine rniss indeces are relative to ground
!     the lte abundances with global indeces are in rnist                                                                  
!                                                                       
!     A plethora of printing options...                                 
!                                                                       
!        jj                                                             
!         1 - 500 strongest emissions lines, sorted by strength         
!         19 - RRC luminosities                                         
!         23 - 500 strongest absorption lines, sorted by strength       
!         24 - absorption edges                                         
!         2 - print input parameter list!                               
!         4 - continuum opacity and emissivity                          
!         5 - energy sums                                               
!         6 - continuum luminosities and depths                         
!         8 - line list                                                 
!         10 - ion abundances and thermal rates (erg/sec)               
!        11 - Write FITS file with summary of ion abundances            
!        12 - append abundance values to the data array for xout_abund1.
!             Doesn't actually write the file, just accumulates values. 
!        13 - blank space                                               
!        14 - line opacities and emissivities                           
!        15 - line luminosities                                         
!        18 - line wavelengths and levels                               
!        20 - line finding list                                         
!        21 - level opacities and emissivities                          
!         7 - level populations                                         
!        17 - print column headings for this pass                       
!         9 - print short summary line of the radial zone               
!        16 - times                                                     
!        22 - ionization parameter etc.                                 
!        25 - outputting to common block                                
!        26 - ferland print                                             
!                                                                       
!     Modifications:                                                    
!        1998/12/17, WTB: Fix FITS keyword format problem for           
!                       writeascii routine.  Removed dependence on      
!                       writeimage routine.                             
!        1999/01/04, WTB: Added model name to writeascii parameter list 
!        1999/01/05, WTB: Removed log(Xi)& log(U1) columns from calls #1
!                       & #12                                           
!                       Inserted '_' in spaces for ion names in call #11
!                                                                       
      use globaldata
      use times
      implicit none 
!                                                                       
!                                                                       
      integer nptmpdim 
      parameter (nptmpdim=400000) 
!                                                                       
!     line luminosities                                                 
      real(8) elum(3,nnnl) 
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum lum                                                     
      real(8) zrems(5,ncn),                                              &
     &          zremsz(ncn)                                             
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
!     continuum flux                                                    
      real(8) bremsa(ncn) 
      real(8) bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     level populations                                                 
      real(8) xilev(nnml),bilev(nnml),rnist(nnml)
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) elumab(2,nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating/cooling                                                   
      real(8) htt(nni),cll(nni) 
      real(8) htt2(nni),cll2(nni) 
!     the atomic data creation date                                     
      character(63) atcredate 
      real(8) rrrt(nni),pirt(nni) 
      real(8) abel(nl),ababs(nl) 
      real(8) xcoltmp(nni) 
      integer kltmp(5000) 
      real(8) zrtmp(999,3999)
      real(8)  zrtmpcol(999,1)
      real(8)  zrtmpc(999,3999)
      real(8)  zrtmph(999,3999)                              

!                                                                       
!      common /ewout/newout,lnewo(nnnl),kdewo(8,nnnl),                  
!     $  kdewol(20,nnnl),kdewou(20,nnnl),aijewo(nnnl),flinewo(nnnl),    
!     $  ggloewo(nnnl),ggupewo(nnnl),                                   
!     $  elewo(nnnl),tau0ewo(nnnl),elout(2,nnnl),zrtmp,epi2,zrems2      
!     A feature added November 2007 is output of the strongest lines,   
!     sorted by element and ion into a common block called 'ewout'      
!     The contents of the common block are:                             
!       newout:  number of lines in the list.                           
!       lnewo:   array conatining line indexes.                         
!       kdewo:   character array containing the name of the ion         
!       kdewol:  character array containing the name of the lower level 
!       kdewou:  character array containing the name of the upper level 
!       aijewo:  array containing A values for the lines                
!       flinewo: array containing f values for the lines                
!       ggloewo: array containing statistical weights for the lower leve
!       ggupewo: array containing statistical weights for the upper leve
!       elewo:   array containing the line wavelengths                  
!       tau0ewo: array containing the line center depths                
!       elout:   array containing line luminosities in xstar units (erg/
!       zrtmp: a 2d array 999x3999, containing a zone-by-zone summary of
!                state of the gas.  The second index is the zone number 
!                the first index is as follows:                         
!                1: radius (cm)                                         
!                2: log(xi)                                             
!                3: electron fraction (relative to nuclei)              
!                4: nucleus number density                              
!                5: pressure (dynes/cm^2)                               
!                6: tempeatre/10^4 K                                    
!                7:  heating-cooling/(heating+cooling)                  
!                8-..: ion fractions for all the ions in the model      
!                    a model with non-zero abundance for all elements   
!                    will have 168 ions (excluding bare nuclei)         
!                    numbered such that 1=H0, 2=He0, 3=He+, 4=C0,       
!                    ... 168=Ni27+.  In this case the upper limit       
!                    for these columns will be 168+8=176                
!       epi: energy grid in eV, length=99999                            
!       zrems: spectrum in erg/s/erg/10**38.  This is a 2d array 3x99999
!           where column 1=transmitted outward flux, including diffuse  
!           emission, column 2=diffuse inward emission (no direct flux) 
!           column 3=diffuse outward emission (no direct flux)          
!                                                                       
!                                                                       
!       character(1) kdewo,kdewol,kdewou                                
!       character(1) klevl(20),klevu(20)                                
!       real(8)  flinewo, aijewo,                                        
!      real(8) ggloewo, ggupewo, elewo, tau0ewo                          
!      real(8) elout                                                     
!      integer nilino, jkktmp, lup, lnewo, lupfnd, llofnd               
!      integer newout                                                   
!                                                                       
      character(1), dimension(:), allocatable :: kdat
      real(8), dimension(:), allocatable ::  elsv
      integer, dimension(:), allocatable :: jpnt
      character(16), dimension(:), allocatable :: klabs,kunits,kform
      character(16) knam,ktmp,kblnk16                                             
      character(20) parname(55) 
      character(10) partype(55) 
      real(8) parms(55) 
      character(30) parcomm(55),kmodelname 
      character(8) spectype, specfile 
      character(1) kdtmp(100),kblnk,klablo(20),klabup(20) 
      integer nparms, specunit, nloopctl 
      character(8) kabstring(30) 
      character(8) ktmp8 
      character(20) ktmp20,klevu,klevl,kblnk20 
      character(9) kinam1 
      character(133) tmpst 
      integer klen, unit, status 
! jg                                                                    
      real(8) xnx, xpx, xee, eliml, elimh 
      real(8) elmmtpp, elcomp, xeltp
      real(8)  cfrac, t, p, tp, xlum 
      real(8) xpxcol, zeta, taumax, xeemin, critf, radexp 
      real(8) vturbi, opsum, tstar, fstr, rsum1 
      real(8) rsum2, sgtmp, tmp, crayj, fstro 
      real(8) emult, ekkr, delte, rssmn, elsum, ergsev 
      real(8) sumtmp1, sumtmp2, r19, tmp1, tmp2, tmp1o 
      real(8) tmp2o, err, sum1, sum2, sum3, sum4 
      real(8) r, ener, etst, aij, gglo, ggup, flin 
      real(8) httot, cltot, htcomp, clcont, cllines,                    &
     &       httot2,cltot2
      real(8) clcomp, clbrems, uu1, enlum, alguu1 
      real(8) skse, ecc, ekt, sksec, zetac, enn0 
      real(8) egam, rdel, hmctot,  expo 
      real(8) elmtp, elmtpb,  terr 
      real(8) ett, optpp, optppo, tmp2c, xcol,fpr2 
      real(8) tmp3, tmp4, tmp5, tmp6, tmp7, tmp8 
      real(8) ttot, enlumx 
      real(8) uux, alguux, eth, abundel 
      real(8)  xi, delr, elin, flux1, flux2 
      real(8) rmax, tinf, rdum, dep, rss, bbe, rocc 
      real(8) flux,ab12
      real(8) vtherm,t1,deleth,aasmall,delea,a,delearad
      real(8) ttmpi
                                                                        
      integer lfnd
      integer jj, lun11, lpril, kltmpo, nlplmx, lm 
      integer nlpl, lnn, nlsvn, ln, ml, ltyp, lrtyp 
      integer lcon, nrdt, nidt, nkdt, nilin 
      integer lmm, kl2, k 
      integer kltmpn, mm, kk, j,  klel, mlel 
      integer jkk, nnz, klion,mlion, jk, mt2, mllel 
      integer kl, nkdti, mlleltp,  mltype 
      integer mlpar, lk, nilin2, mllz, nlev, kkkl, idest1 
      integer  nlevp, idest2, jkk3, ndtmp 
      integer  iltmp, ltyp2, lrtyp2, lcon2, nrdt2 
      integer nidt2, nkdt2, lcdd, numrec, nlimd, lwri, lpri 
      integer lfix, npass, numcon, ncn2, i, jlk
      integer ktt,  nell, lkk, nelin, jkl, mmlv, nidti 
      integer nlyc, nry, jkstep, ilevup, ilevlo
      integer niter, jjj, jp1
      integer    mlcu, mm2, mmtmp, ll 
      integer ntotit, nb1, nb10 
      integer lnerrd, nbinc 
      integer nnmax, ncsvn,mlm,nnzel
      integer np1i,np1r,np1k,np1i2,np1r2,np1k2,np1ki 
!                                                                       
      logical done 
!                                                                       
!     Not used                                                          
      real(8) javir 
      integer javi 
      character(20) javik                                              
!                                                                       
      data kblnk20/'                    '/ 
      data kblnk/' '/,kblnk16/'                '/ 
      data kabstring/'H abund=','Heabund=','Liabund=',                  &
     &               'Beabund=','B abund=','C abund=','N abund=',       &
     &               'O abund=','F abund=','Neabund=','Naabund=',       &
     &               'Mgabund=','Alabund=','Siabund=','P abund=',       &
     &               'S abund=','Clabund=','Arabund=','K abund=',       &
     &               'Caabund=','Scabund=','Tiabund=','V abund=',       &
     &               'Crabund=','Mnabund=','Feabund=','Coabund=',       &
     &               'Niabund=','Cuabund=','Znabund='/                  
!                                                                       
      allocate(kdat(nptmpdim))
      allocate(elsv(nnnl))
      allocate(jpnt(nnnl))
      allocate(klabs(3999))
      allocate(kunits(3999))
      allocate(kform(3999))
!                                                                       
!                                                                       
      javi=nnmax 
      javir=rmax 
      javi=nparms 
      javik=parname(1)                                                 
      javik=partype(1)                                                 
      javir=tinf 
      javir=bremsa(1) 
      bremsa(1)=javir 
      javir=bremsint(1) 
      bremsint(1)=javir 
      javi=derivedpointers%nplini(1) 
      javi=derivedpointers%npilevi(1) 
      javi=derivedpointers%npconi(1) 
      javi=ncsvn 
      javir=bilev(1) 
      javir=xi 
      javi=lnerrd 
!      lnerrd=javi                                                      
!                                                                       
                                                                        
      xnx=xpx*xee 
!                                                                       
      if ((jj.ne.9).and.(jj.ne.12))                                     &
     & write (lun11,9211)jj                                             
 9211 format (1x, 'print option:',i2) 
      if ((jj.le.0).or.(jj.gt.27)) go to 9000
!                                                                       
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,   &
     &23,24,25,26,27),                                                  &
     &  jj                                                              
!                                                                       
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
!                                                                       
    1 continue 
!                                                                       
      lpril=0 
!     500 strongest emission lines, sorted by strength                  
      write (lun11,*)'emission line luminosities (erg/sec/10**38))' 
      kltmpo=0 
      nlplmx=500 
      eliml=0.1 
      elimh=1.0e10 
!     find the strongest lines.                                         
      do  lm=1,nlplmx 
        kltmp(lm)=0 
        enddo 
!                                                                       
!     step through lines                                                
      nlpl=1 
      do lnn=1,nlsvn 
!                                                                       
!       get line data                                                   
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
!       exclude rate type 14                                            
        if ((lrtyp.ne.14).and.(lrtyp.ne.9)) then 
!                                                                       
!         get ion data                                                  
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          nilin2=masterdata%idat1(np1i+nidt-1) 
!                                                                       
!         get lum and test for strength, wavelength                     
          elmmtpp=(elum(2,ln)+elum(1,ln))/2. 
          if (lpril.ne.0)                                               &
     &       write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml                
          if ((ln.gt.0).and.(ln.lt.nnnl)                                &
     &       .and.(elin.ge.eliml).and.(elin.le.elimh)                   &
     &       .and.(elin.le.8.9e+6)                                      &
     &       .and.(elmmtpp.gt.1.d-49)                                   &
     &       .and.(nilin2.gt.0).and.(nilin2.le.nni))                    &
     &        then                                                      
!                                                                       
!           insertion sort                                              
            lmm=0 
            elcomp=1.e+10 
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp)) 
              lmm=lmm+1 
              kl2=kltmp(lmm) 
              elcomp=0. 
              if (kl2.gt.0)                                             &
     &          elcomp=(elum(2,kl2)+elum(1,kl2))/2.                     
              enddo 
            if (lpril.ne.0)                                             &
     &       write (lun11,8516)ln,elin,elmmtpp,lmm,nlpl,kl2,elcomp      
 8516       format (1h ,i4,2e12.4,3i4,e12.4) 
            kltmpo=ln 
            do  k=lmm,min(nlplmx,nlpl) 
              if ((lpril.ne.0).and.(kltmp(k).ne.0))                     &
     &          write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo          
              kltmpn=kltmp(k) 
              kltmp(k)=kltmpo 
              kltmpo=kltmpn 
              enddo 
           nlpl=min(nlplmx,nlpl+1) 
           if (lpril.ne.0)                                              &
     &       write (lun11,*)'done with 557 loop',lm                     
            endif 
!           end of insertion                                            
!                                                                       
          endif 
!                                                                       
        enddo 
!                                                                       
      if (nlpl.gt.0) kltmp(nlpl)=kltmpo 
!                                                                       
!     printing loop                                                     
      write (lun11,959) 
!                                                                       
!     step through lines                                                
      do  kk=1,nlpl 
        if (lpril.ne.0)                                                 &
     &    write (lun11,*)'kk=',kk                                       
        ln=kltmp(kk) 
        if (ln.ne.0) then 
!                                                                       
!         get line data                                                 
          ml=derivedpointers%nplin(ln) 
          if (ml.ne.0) then 
            if (lpril.ne.0)                                             &
     &      write (lun11,*)'   ',ln,ml                                  
            mlm=ml
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
!           get ion data                                                
            nilin=derivedpointers%npar(ml) 
            mlm=nilin
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            do mm=1,nkdt 
              kdtmp(mm)=masterdata%kdat1(np1k-1+mm) 
              enddo 
            do mm=nkdt+1,9 
              kdtmp(mm)=kblnk 
              enddo 
!            nilin=masterdata%idat1(np1i+2)                                        
            if (lpril.ne.0)                                             &
     &      write (lun11,*)ml,nilin,derivedpointers%npar(ml)                   
!                                                                       
!           print                                                       
            write (lun11,9955)kk,ln,(kdtmp(mm),mm=1,9),elin,            &
     &      elum(1,ln),elum(2,ln)                                       
 9955       format (1x,2i8,1x,9a1,3(1pe13.5)) 
!                                                                       
            endif 
!                                                                       
          endif 
!                                                                       
        enddo 
!                                                                       
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   19 continue 
!                                                                       
!     print 500 strongest recombination continua                        
      lpril=0 
      write (lun11,*)'recombination continuum luminosities',            &
     &  '(erg/sec/10**38))'                                             
      write (lun11,*)'index, ion, level, energy (eV), RRC luminosity ' 
!                                                                       
!     lpril is flag for printing debug information                      
      if (lpril.ne.0) then 
        write (lun11,*)'raw data' 
        do j=1,nnml 
          if (xilev(j).gt.1.e-37)                                       &
     &     write (lun11,*)j,xilev(j),elumab(1,j)                        
          enddo 
        endif 
!                                                                       
!     First look for element data (jk is element index)                 
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      jk=0 
      kk=0 
      jkk=0 
!                                                                       
!     step through elements                                             
      do while (mlel.ne.0) 
!                                                                       
!       get element data                                                
        jk=jk+1 
        mt2=mlel
        call drd(ltyp,lrtyp,lcon,                                       &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                        &
     &        0,lun11)                                            
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=masterdata%rdat1(np1r) 
        xeltp=abel(mllel) 
        nnz=masterdata%idat1(np1i) 
        if (lpril.ge.1) then
              write (lun11,902)jk,mlel,nnz,                             &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,min(8,nkdt))
902           format (1x,'  element:',3(i12,1x),8(1a1))
           endif
!                                                                       
!       ignore if the abundance is small                                
        if (xeltp.lt.1.e-10) then 
            jkk=jkk+nnz 
          else 
!                                                                       
!           now step thru ions (jkk is ion index)                       
            klion=12 
            mlion=derivedpointers%npfirst(klion) 
            jkk=0 
            kl=0 
            do while ((mlion.ne.0).and.(kl.lt.nnz)) 
              jkk=jkk+1 
!                                                                       
!             retrieve ion name from kdati                              
              mlm=mlion
              call drd(ltyp,lrtyp,lcon,                                 &
     &            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,                 &
     &            0,lun11)                                        
!                                                                       
!             if not accessing the same element, skip to the next elemen
              mlleltp=masterdata%idat1(np1i+nidti-2) 
              if (mlleltp.eq.mllel) then 
!                                                                       
                kl=kl+1 
                if (lpril.ge.1)                                         &
     &            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,         &
     &                (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)            
!                                                                       
!               now find level data                                     
                jkk=masterdata%idat1(nidt+np1i-1)
                call calc_rates_level_lte(jkk,lpril,lun11,t,xee,xpx,    &
     &              nlev)
!                                                                       
!               now step through rate type 7 data                       
                mltype=7 
                ml=derivedpointers%npfi(mltype,jkk) 
                mllz=0 
                if (ml.ne.0) mllz=derivedpointers%npar(ml) 
                mlpar=0 
                if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
!                                                                       
!                 get rrc data                                          
                  kkkl=derivedpointers%npconi2(ml) 
                  if (lpril.ne.0) write (lun11,*)kkkl,ml,idest1,        &
     &                    elumab(1,kkkl),elumab(2,kkkl)                 
!                                                                       
!                 test for non-zero rrc data                            
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)                   &
     &                .and.((elumab(1,kkkl).gt.1.d-49)                  &
     &                .or.(elumab(2,kkkl).gt.1.d-49))) then             
!                                                                       
!                   get rrc  data                                       
                    mlm=ml
                    call drd(ltyp,lrtyp,lcon,                           &
     &                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                &
     &                0,lun11)                                    
                    idest1=masterdata%idat1(np1i+nidt-2) 
                    nlevp=nlev 
                    idest2=nlevp+masterdata%idat1(np1i-1+nidt-3)-1 
!                                                                       
!                   label for lower level                               
                    do lk=1,20 
                      write (ktmp20(lk:lk),'(a1)')                      &
     &                    leveltemp%klev(lk,idest1) 
                      enddo 
                    klevl=ktmp20 
!                                                                       
!                   label for upper level                               
                    write (ktmp20(1:20),'(a20)')'continuum           ' 
                    klevu=ktmp20 
!                                                                       
!                   ion label                                           
                    do lk=1,nkdti 
                       write (ktmp8(lk:lk),'(a1)')                      &
     &                       masterdata%kdat1(np1ki+lk-1) 
                      enddo 
                    do lk=nkdti+1,8 
                      write (ktmp8(lk:lk),'(a1)')kblnk 
                      enddo 
!                                                                       
                    eth=leveltemp%rlev(4,idest1)                        &
     &                  -leveltemp%rlev(1,idest1) 
                    ett=eth 
!                                                                       
!                   get upper level data                                
                    go to 9089
                    if (idest2.gt.nlevp) then 
                      jkk3=jkk+1 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*)jkk3,ndtmp,nlevp,idest2          
                      ndtmp=derivedpointers%npfi(13,jkk3) 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*)jkk3,ndtmp,nlevp,idest2          
                      if (ndtmp.le.0) stop 'ndtmp error' 
                      mllz=derivedpointers%npar(ndtmp) 
                      iltmp=0 
                      do while ((ndtmp.ne.0).and.                       &
     &                    (iltmp.ne.(idest2-nlevp+1)).and.              &
     &                    (derivedpointers%npar(ndtmp).eq.mllz))               
                        mlm=ndtmp
                        call drd(ltyp2,lrtyp2,lcon2,                    &
     &                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,      &
     &                    0,lun11)                                
                        iltmp=masterdata%idat1(np1i2+nidt2-2) 
                        if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
                        ndtmp=derivedpointers%npnxt(ndtmp) 
                        if (ndtmp.le.0) stop 'ndtmp error' 
                        enddo 
!                     NB fix to excited level PI and rec                
                      ett=ett+masterdata%rdat1(np1r2) 
                      eth=ett 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*) ndtmp,iltmp,idest2,ett          
!                     label for lower level                             
                      ktmp20=kblnk20 
                      do lk=1,nkdt2 
                       write (ktmp8(lk:lk),'(a1)')                      &
     &                       masterdata%kdat1(np1ki+lk-1) 
                        enddo 
                      klevu=ktmp20 
                      endif 
9089                continue
!                                                                       
!                   other data                                          
                    mmlv=derivedpointers%npilev(idest1,jkk) 
                    write (lun11,9293)kkkl,mmlv,ktmp8,idest1,idest2,    &
     &                  klevl,klevu,eth,elumab(1,kkkl),elumab(2,kkkl)   
!                                                                       
!                   done with this rr!                                  
                    endif 
!                                                                       
!                 end of loop over rrcs                                 
                  ml=derivedpointers%npnxt(ml) 
                  if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                  enddo 
!                                                                       
!               end of test for element                                 
                endif 
!                                                                       
!             Go to next ion                                            
              mlion=derivedpointers%npnxt(mlion) 
              enddo 
!                                                                       
!         end of test for non-zero element abund                        
          endif 
!                                                                       
        mlel=derivedpointers%npnxt(mlel) 
!       Go to next element                                              
        enddo 
!                                                                       
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   23 continue 
!                                                                       
      lpril=0 
!     print 500 strongest absoprtion lines                              
      write (lun11,*)'line depths' 
      kltmpo=0 
      nlplmx=500 
      eliml=0.1 
      elimh=1.0e10 
!     find the strongest lines.                                         
      do  lm=1,nlplmx 
        kltmp(lm)=0 
        enddo 
!                                                                       
!     step through lines                                                
      nlpl=1 
      do lnn=1,nlsvn 
!                                                                       
!       get line data                                                   
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
!       exclude rate type 14                                            
        if ((lrtyp.ne.14).and.(lrtyp.ne.9)) then 
!                                                                       
!         get ion data                                                  
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          nilin2=masterdata%idat1(np1i+nidt-1) 
!                                                                       
!         get lum and test for strength, wavelength                     
          elmmtpp=tau0(1,ln) 
          if (lpril.ne.0)                                               &
     &       write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml                
          if ((ln.gt.0).and.(ln.lt.nnnl)                                &
     &       .and.(elin.ge.eliml).and.(elin.le.elimh)                   &
     &       .and.(elin.le.8.9e+6)                                      &
     &       .and.(elmmtpp.gt.1.d-49)                                   &
     &       .and.(nilin2.gt.0).and.(nilin2.le.nni))                    &
     &        then                                                      
!                                                                       
!           insertion sort                                              
            lmm=0 
            elcomp=1.e+10 
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp)) 
              lmm=lmm+1 
              kl2=kltmp(lmm) 
              elcomp=0. 
              if (kl2.gt.0)                                             &
     &          elcomp=tau0(1,kl2)                                      
              enddo 
            if (lpril.ne.0)                                             &
     &       write (lun11,8516)ln,elin,elmmtpp,lmm,nlpl,kl2,elcomp      
            kltmpo=ln 
            do  k=lmm,min(nlplmx,nlpl) 
              if ((lpril.ne.0).and.(kltmp(k).ne.0))                     &
     &          write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo          
              kltmpn=kltmp(k) 
              kltmp(k)=kltmpo 
              kltmpo=kltmpn 
              enddo 
           nlpl=min(nlplmx,nlpl+1) 
           if (lpril.ne.0)                                              &
     &       write (lun11,*)'done with 557 loop',lm                     
            endif 
!           end of insertion                                            
!                                                                       
          endif 
!                                                                       
        enddo 
!                                                                       
      if (nlpl.gt.0) kltmp(nlpl)=kltmpo 
!                                                                       
!     printing loop                                                     
      write (lun11,959) 
  959 format (1x,'index, ion, wavelength, reflected, transmitted') 
!                                                                       
!     step through lines                                                
      do  kk=1,nlpl 
        if (lpril.ne.0)                                                 &
     &    write (lun11,*)'kk=',kk                                       
        ln=kltmp(kk) 
        if (ln.ne.0) then 
!                                                                       
!         get line data                                                 
          ml=derivedpointers%nplin(ln) 
          if (ml.ne.0) then 
            if (lpril.ne.0)                                             &
     &      write (lun11,*)'   ',ln,ml                                  
            mlm=ml
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
!           get ion data                                                
            nilin=derivedpointers%npar(ml) 
            mlm=nilin
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            do mm=1,nkdt 
              kdtmp(mm)=masterdata%kdat1(np1k-1+mm) 
              enddo 
            do mm=nkdt+1,9 
              kdtmp(mm)=kblnk 
              enddo 
!            nilin=masterdata%idat1(np1i+2)     
            if (lpril.ne.0)                                             &
     &      write (lun11,*)ml,nilin,derivedpointers%npar(ml)                   
!                                                                       
!           print                                                       
            write (lun11,9955)kk,ln,(kdtmp(mm),mm=1,9),elin,            &
     &      tau0(1,ln),tau0(2,ln)                                       
!                                                                       
            endif 
!                                                                       
          endif 
!                                                                       
        enddo 
!                                                                       
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   24 continue 
!                                                                       
      lpril=0 
!     print 500 strongest absorption edges                              
      write (lun11,*)'absorption edge depths' 
      write (lun11,*)'index, ion, level, energy (eV), depth ' 
!                                                                       
!     lpril is flag for printing debug information                      
      if (lpril.ne.0) then 
        write (lun11,*)'raw data' 
        do j=1,nnml 
          if (xilev(j).gt.1.e-37)                                       &
     &     write (lun11,*)j,xilev(j),tauc(1,j)                          
          enddo 
        endif 
!                                                                       
!     First look for element data (jk is element index)                 
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      jk=0 
      kk=0 
      jkk=0 
!                                                                       
!     step through elements                                             
      do while (mlel.ne.0) 
!                                                                       
!       get element data                                                
        jk=jk+1 
        mt2=mlel
        call drd(ltyp,lrtyp,lcon,                                       &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                        &
     &        0,lun11)                                            
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=masterdata%rdat1(np1r) 
        xeltp=abel(mllel) 
        nnz=masterdata%idat1(np1i) 
        if (lpril.ge.1) then
             write (lun11,902)jk,mlel,nnz,                              &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,nkdt),xeltp                    
          endif
!                                                                       
!       ignore if the abundance is small                                
        if (xeltp.lt.1.e-10) then 
            jkk=jkk+nnz 
          else 
!                                                                       
!           now step thru ions (jkk is ion index)                       
            klion=12 
            mlion=derivedpointers%npfirst(klion) 
            jkk=0 
            kl=0 
            do while ((mlion.ne.0).and.(kl.lt.nnz)) 
              jkk=jkk+1 
!                                                                       
!             retrieve ion name from kdati                              
              mlm=mlion
              call drd(ltyp,lrtyp,lcon,                                 &
     &            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,                 &
     &            0,lun11)                                        
!                                                                       
!             if not accessing the same element, skip to the next elemen
              mlleltp=masterdata%idat1(np1i+nidti-2) 
              if (mlleltp.eq.mllel) then 
!                                                                       
                kl=kl+1 
                if (lpril.ge.1)                                         &
     &            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,         &
     &               (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)            
!                                                                       
!               now find level data                                     
                jkk=masterdata%idat1(nidt+np1i-1)
                call calc_rates_level_lte(jkk,lpril,lun11,t,xee,xpx,    &
     &              nlev)
!                                                                       
!               now step through rate type 7 data                       
                mltype=7 
                ml=derivedpointers%npfi(mltype,jkk) 
                mllz=0 
                if (ml.ne.0) mllz=derivedpointers%npar(ml) 
                mlpar=0 
                if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
!                                                                       
!                 get rrc data                                          
                  kkkl=derivedpointers%npconi2(ml) 
                  if (lpril.ne.0) write (lun11,*)kkkl,ml,idest1,        &
     &                    elumab(1,kkkl),elumab(2,kkkl)                 
!                                                                       
!                 test for non-zero rrc data                            
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)                   &
     &                .and.((tauc(1,kkkl).gt.1.d-49)                    &
     &                .or.(tauc(2,kkkl).gt.1.d-49))) then             
!                                                                       
!                   get rrc  data                                       
                    mlm=ml
                    call drd(ltyp,lrtyp,lcon,                           &
     &                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                &
     &                0,lun11)                                    
                    idest1=masterdata%idat1(np1i+nidt-2) 
                    nlevp=nlev 
                    idest2=nlevp+masterdata%idat1(np1i-1+nidt-3)-1 
!                                                                       
!                   label for lower level                               
                    do lk=1,20 
                      write (ktmp20(lk:lk),'(a1)')                      &
     &                    leveltemp%klev(lk,idest1) 
                      enddo 
                    klevl=ktmp20 
!                                                                       
!                   label for upper level                               
                    write (ktmp20(1:20),'(a20)')'continuum           ' 
                    klevu=ktmp20 
!                                                                       
!                   ion label                                           
                    do lk=1,nkdti 
                       write (ktmp8(lk:lk),'(a1)')                      &
     &                       masterdata%kdat1(np1ki+lk-1) 
                      enddo 
                    do lk=nkdti+1,8 
                      write (ktmp8(lk:lk),'(a1)')kblnk 
                      enddo 
!                                                                       
                    eth=leveltemp%rlev(4,idest1)                        &
     &                 -leveltemp%rlev(1,idest1) 
                    ett=eth 
!                                                                       
!                   get upper level data                                
                    go to 9088
                    if (idest2.gt.nlevp) then 
                      jkk3=jkk+1 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*)jkk3,ndtmp,nlevp,idest2          
                      ndtmp=derivedpointers%npfi(13,jkk3) 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*)jkk3,ndtmp,nlevp,idest2          
                      if (ndtmp.le.0) stop 'ndtmp error' 
                      mllz=derivedpointers%npar(ndtmp) 
                      iltmp=0 
                      do while ((ndtmp.ne.0).and.                       &
     &                    (iltmp.ne.(idest2-nlevp+1)).and.              &
     &                    (derivedpointers%npar(ndtmp).eq.mllz))               
                        mlm=ndtmp
                        call drd(ltyp2,lrtyp2,lcon2,                    &
     &                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,      &
     &                    0,lun11)                                
                        iltmp=masterdata%idat1(np1i2+nidt2-2) 
                        if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
                        ndtmp=derivedpointers%npnxt(ndtmp) 
                        if (ndtmp.le.0) stop 'ndtmp error' 
                        enddo 
!                     NB fix to excited level PI and rec                
                      ett=ett+masterdata%rdat1(np1r2) 
                      eth=ett 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*) ndtmp,iltmp,idest2,ett          
!                     label for lower level                             
                      ktmp20=kblnk20 
                      do lk=1,nkdt2 
                         write (ktmp20(lk:lk),'(a1)')                   &
     &                         masterdata%kdat1(np1k2+lk-1) 
                        enddo 
                      klevu=ktmp20 
                      endif 
9088                  continue
!                                                                       
!                   other data                                          
                    mmlv=derivedpointers%npilev(idest1,jkk) 
                    write (lun11,9293)kkkl,mmlv,ktmp8,idest1,idest2,    &
     &                  klevl,klevu,eth,tauc(1,kkkl),tauc(2,kkkl)       
 9293               format(1x,2i6,1x,a8,2i6,1x,2(a20,1x),3(1pe11.3)) 
!                                                                       
!                   done with this rrc                                  
                    endif 
!                                                                       
!                 end of loop over rrcs                                 
                  ml=derivedpointers%npnxt(ml) 
                  if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                  enddo 
!                                                                       
!               end of test for element                                 
                endif 
!                                                                       
!             Go to next ion                                            
              mlion=derivedpointers%npnxt(mlion) 
              enddo 
!                                                                       
!         end of test for non-zero element abund                        
          endif 
!                                                                       
        mlel=derivedpointers%npnxt(mlel) 
!       Go to next element                                              
        enddo 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
    2 continue 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
!                                                                       
!     Print list of input parameters                                    
!                                                                       
!                                                                       
!                                                                       
      write (lun11,'(1x)') 
      write (lun11,*)'input parameters:' 
      write (lun11,'(a22,1pe11.3)')'covering fraction=    ',cfrac 
      write (lun11,'(a22,1pe11.3)')'temperature (/10**4K)=',t 
      write (lun11,*)'constant pressure switch (1=yes, 0=no)=',1-lcdd 
      write (lun11,'(a22,1pe11.3)')'pressure (dyne/cm**2)=',p 
      write (lun11,'(a22,1pe11.3)')'density (cm**-3)=     ',xpx 
      write (lun11,*)'spectrum type=',spectype 
      write (lun11,*)'spectrum file=',specfile 
      write (lun11,*)'spectrum units? (0=energy, 1=photons)',specunit 
      write (lun11,'(a31,1pe11.3)')'radiation temperature or alpha=',tp 
      write (lun11,'(a27,1pe11.3)')'luminosity (/10**38 erg/s)=',xlum 
      write (lun11,'(a22,1pe11.3)')'column density (cm**-2)=',xpxcol 
      write (lun11,'(a26,1pe11.3)')'log(ionization parameter)=',zeta 
      r19=r/1.e+19
      flux=xlum/12.56/r19/r19
      write (lun11,'(a22,1pe11.3)')'flux=                 ',flux
      write (lun11,*)'abundances:'
      write (lun11,*)'element,   rel.to cosmic,     rel. to H,     H=12'
      do j=1,nl 
         ab12=12.+log10(max(1.d-12,ababs(j))/max(1.d-48,ababs(1)))
         write (lun11,9322)kabstring(j),abel(j),ababs(j),ab12 
9322     format (1x,a8,3(1pe11.3))
         enddo 
      write (lun11,*)'model name=',kmodelname 
      write (lun11,*)'number of steps=',numrec 
      write (lun11,*)'number of iterations=',nlimd 
      write (lun11,*)'write switch (1=yes, 0=no)=',lwri 
      write (lun11,*)'print switch (1=yes, 0=no)=',lpri 
      write (lun11,*)'step size choice switch=',lfix 
      write (lun11,*)'loop control (0=standalone)=',nloopctl 
      write (lun11,*)'number of passes=',npass 
      write (lun11,*)'emult=',emult 
      write (lun11,*)'taumax=',taumax 
      write (lun11,*)'xeemin=',xeemin 
      write (lun11,*)'critf=',critf 
      write (lun11,*)'vturbi=',vturbi 
      write (lun11,*)'ncn2=',ncn2 
      write (lun11,*)'radexp=',radexp 
      write (lun11,'(1x)') 
      write (lun11,993) 
      go to 9000
                                                                        
!      r19=r*(1.d-97)                                                   
!      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10                             
!      alguu1=log10(max(1.e-24,uu1))                                    
!      skse=xlum/(xpx*r19*r19)                                          
!      zeta=log10(max(1.e-24,skse))                                     
!      ecc=2.998e+10                                                    
!      ekt=t*(0.861707)*ergsev                                          
!      sksec=skse/12.56/((1.+xee)*ekt*ecc)                              
!      zetac=log10(max(1.e-24,sksec))                                   
!      enn0=xpx                                                         
!      nlyc=nbinc(13.7,epi,ncn2)                                        
!      nry=nlyc+1                                                       
!      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24)              
!      nry=nbinc(13.6,epi,ncn2)+1                                       
!      write (lun11,993)                                                
!      write (lun11,*)'input parameters:'                               
!      write (lun11,*)'continuum luminosity=',xlum                      
!      write (lun11,*)'pressure or density=',xpx                        
!      write (lun11,*)'radius=',r                                       
!      write (lun11,*)'ionization parameter=',skse                      
!      write (lun11,993)                                                
!      write (tmpst,993)                                                
!      call xwrite(tmpst,10)                                            
!      write (tmpst,*)'input parameters:'                               
!      call xwrite(tmpst,10)                                            
!      write (tmpst,*)'continuum luminosity=',xlum                      
!      call xwrite(tmpst,10)                                            
!      write (tmpst,*)'pressure or density=',xpx                        
!      call xwrite(tmpst,10)                                            
!      write (tmpst,*)'radius=',r                                       
!      call xwrite(tmpst,10)                                            
!      write (tmpst,*)'ionization parameter=',skse                      
!      call xwrite(tmpst,10)                                            
!      write (tmpst,993)                                                
!      call xwrite(tmpst,10)                                            
                                                                        
                                                                        
    3 continue 
!                                                                       
!     Write the parameter list to the FITS file                         
!     When changing this list, make sure nparms, parname, partype,      
!     and parcomm are also properly updated.  Watch out for the         
!     model name which is currently parcomm(37)                         
      parms(1)=cfrac 
      parms(2)=t 
      parms(3)=1-lcdd 
      parms(4)=p 
      parms(5)=xpx 
      parms(6)=0.0 
      parcomm(6)=spectype 
      parms(7)=0.0 
      parcomm(7)=specfile 
      parms(8)=specunit 
      parms(9)=tp 
      parms(10)=xlum 
      parms(11)=xpxcol 
      parms(12)=zeta 
      parms(13)=numrec 
      parms(14)=nlimd 
      parms(15)=lwri 
      parms(16)=lpri 
      parms(17)=lfix 
      do j=1,nl 
         parms(17+j)=abel(j) 
         enddo 
      parms(17+nl+1)=emult 
      parms(17+nl+2)=taumax 
      parms(17+nl+3)=xeemin 
      parms(17+nl+4)=critf 
      parms(17+nl+5)=vturbi 
      parms(17+nl+6)=npass 
      parms(17+nl+7)=0.0 
      parcomm(17+nl+7)=kmodelname 
      parms(17+nl+8)=nloopctl 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
    4 continue 
!                                                                       
!     print continuum opacity and emissivity                            
      write (lun11,*)                                                   &
     & 'continuum opacity and emissivities (/cm**3/sec/10**38)'         
      write (lun11,*)'channel, energy,      opacity,    sigma*e**3,'    &
     &, 'scattered,  rec. in,   rec. out,  brem. em., source, bbe,',    &
     &'photon occ'                                                      
      opsum=0. 
      tstar=t 
      ekkr=xnx*(6.65e-25) 
      ekkr=max(1.d-20,ekkr) 
      optpp=max(opakc(1),ekkr) 
      fstr=0. 
      rsum1=0. 
      rsum2=0. 
      numcon=ncn2 
!                                                                       
!     step thru continuum bins                                          
      do 135 kl=2,numcon 
!                                                                       
!        sigma*e**3                                                     
         sgtmp=(opakc(kl)*(epi(kl)/1000.)**3)/max(1.d-24,xpx) 
!                                                                       
!        calculate sum and rosseland mean                               
         if ((kl.gt.1).and.(epi(kl).gt.100.))                           &
     &    opsum=opsum+(opakc(kl)+opakc(kl-1))*(epi(kl)-epi(kl-1))/2.    
         i=kl 
         tmp = epi(i)*1.16/tstar 
         crayj = 1./tmp 
         fstro = fstr 
         if ( tmp.le.50. ) then 
            if ( tmp.gt.1.e-4 ) crayj = 1./(expo(tmp)-1.) 
            crayj=crayj*crayj 
!            fstr= cconst*tmp*crayj*epi(i)**3/tstar                     
            fstr= tmp*crayj*epi(i)**3/tstar 
         endif 
         optppo=optpp 
         optpp=max(opakc(kl),ekkr) 
         delte=epi(kl)-epi(kl-1) 
         rsum1=min(1.d+20,rsum1+(fstr/optpp+fstro/optppo)*delte/2.) 
         rsum2=min(1.d+20,rsum2+(fstr+fstro)*delte/2.) 
!                                                                       
!        source function                                                
         rss=(rccemis(1,kl)+rccemis(2,kl)+brcems(kl)/12.56)/            &
     &           (1.d-36+opakc(kl))                                     
!         rss=(rccemis(1,kl)+rccemis(2,kl))/(1.d-36+opakc(kl))          
!                                                                       
!        planck function                                                
         bbe=2.*min(2.d+4,epi(kl))**3*(1.5642e+22)                      &
     &       /(expo(epi(kl)/(0.861707*t))-1.+1.d-36)                    
!                                                                       
!        photon occupation number                                       
         rocc=rss/(bbe+1.d-36) 
!                                                                       
!        print                                                          
         write (lun11,967)kl,epi(kl),opakc(kl),sgtmp,opakcont(kl),      &
     &            rccemis(1,kl),rccemis(2,kl), brcems(kl),rss,bbe,rocc  
  967    format (1h ,i6,11(1pe13.5)) 
!                                                                       
  135    continue 
!                                                                       
!     print summed opacities                                            
      write (lun11,*)'opsum cont=',opsum 
      rssmn=rsum2/rsum1 
!      ens1=rssmn/(t*1.e+4)**(-3.5)/xpx/xpx/1.66e-24                    
      write (lun11,*)'rosseland mean opacity=',t,rssmn 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
    5 continue 
!                                                                       
!     print energy sums                                                 
      elsum=0. 
      ergsev=1.602197e-12 
      do jlk=1,nlsvn 
         ln=jlk 
         ml=derivedpointers%nplin(ln) 
         mlm=ml
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         elin=abs(masterdata%rdat1(np1r)) 
         if ((elin.lt.1.e+8).and.(elin.gt.1.)) then 
           elsum=elsum+elum(1,ln)+elum(2,ln) 
           endif 
         enddo 
      sumtmp1=0. 
      sumtmp2=0. 
      ergsev=1.602197e-12 
      r19=r*1.e-19
      tmp1=zremsz(1)*(1.-exp(-dpthc(1,1))) 
      tmp2=zrems(3,1)+zrems(2,1) 
      do jk=2,ncn2 
         tmp1o=tmp1 
         tmp1=zremsz(jk)*(1.-exp(-dpthc(1,jk))) 
         sumtmp1=sumtmp1+(tmp1+tmp1o)*(epi(jk)-epi(jk-1))*ergsev/2. 
         tmp2o=tmp2 
         tmp2=zrems(2,jk)+zrems(3,jk) 
         sumtmp2=sumtmp2+(tmp2+tmp2o)*(epi(jk)-epi(jk-1))*ergsev/2. 
         enddo 
      err=(sumtmp1-sumtmp2-elsum)/(sumtmp1+1.e-24) 
      write (lun11,9981)sumtmp1,sumtmp2,elsum,err 
 9981 format (1x,'energy sums: abs, cont, line, err:',4(1pe13.5)) 
      write (lun11,993) 
!      write (lun11,*),httot,cltot,fpr2dr,httot/(sumtmp1+1.e-24),       
!     $  cltot/(elsum+sumtmp2+1.e-24)                                   
!                                                                       
      go to 9000
!                                                                       
!                                                                       
6     continue
!c
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c
!     print continuum luminosities and depths
      write (lun11,*)'continuum luminosities (/sec/10**38) and depths'
      write (lun11,*)'real quantities are as follows:'
      write (lun11,*)                                                   &
     &'1 ) photon energy in eV'
      write (lun11,*)                                                   &
     &'2) Incident radiation field'
      write (lun11,*)                                                   &
     &'3) Radiation field used by xstar for internal calculations of rat&
     &es, etc.'
      write (lun11,*)                                                   &
     &'   The quantity defined in equation (12) in the xstar manual, cha&
     &pter 8, the Physics of xstar.'
      write (lun11,*)                                                   &
     &'4) and 5) the quantities in equations (15) and (16), with lines a&
     &dded'
      write (lun11,*)                                                   &
     &'6) and 7) the quantities in equations (15) and (16), without line&
     &s'
      write (lun11,*)                                                   &
     &'8) inward depth'
      write (lun11,*)                                                   &
     &'9) outward depth'
      write (lun11,*)                                                   &
     &'10) Planck function at local gas temperature'
      write (lun11,*)                                                   &
     &'11) ratio of the xstar internal radiation field (in flux units) t&
     &o Planck function at local gas temperature.'
!
      numcon=ncn2
      sum1=0.
      sum2=0.
      sum3=0.
      sum4=0.
      ergsev=1.602197e-12
      r19=r*1.e-19
      fpr2=12.56*r19*r19
!
!     step thru continuum bins
      do kl=1,numcon
!
!       planck function
        bbe=2.*min(2.d+4,epi(kl))**3*(1.5642e+22)                       &
     &       /(expo(epi(kl)/(0.861707*t))-1.+1.d-36)
!
!       photon occupation number
        rocc=zrems(1,kl)/(bbe+1.d-36)/fpr2/12.56
!
        write (lun11,968)kl,epi(kl),zremsz(kl),                         &
     &    zrems(1,kl),zrems(2,kl),zrems(3,kl),zrems(4,kl),              &
     &    zrems(5,kl),dpthc(1,kl),dpthc(2,kl),bbe,rocc
!
!       sums
        if (kl.gt.1) then
           sum1=sum1+(zremsz(kl)+zremsz(kl-1))                          &
     &         *(epi(kl)-epi(kl-1))*ergsev/2.
           sum2=sum2+(zrems(1,kl)+zrems(1,kl-1))                        &
     &         *(epi(kl)-epi(kl-1))*ergsev/2.
           sum3=sum3+(zrems(2,kl)+zrems(2,kl-1))                        &
     &         *(epi(kl)-epi(kl-1))*ergsev/2.
           sum4=sum4+(zrems(3,kl)+zrems(3,kl-1))                        &
     &         *(epi(kl)-epi(kl-1))*ergsev/2.
           endif
968      format (1h ,i6,11(1pe13.5))
!
         enddo
!
      write (lun11,*)'norms:'
      write (lun11,9698)sum1,sum2,sum3,sum4
 9698 format (20x,4(1pe13.5))
      write (lun11,993)
!
!                                                                       
      go to 9000
!                                                                       
!                                                                       
    8 continue 
!                                                                       
                                                                        
      lpril=0 
      write (lun11,*)'line list' 
      write (lun11,9943) 
 9943 format (1x,'           wave (A)   ion        Aul        fij       &
     &   glo         gup        level_lo             level_up           &
     &   DE_aug      DE_th       a_damp')
!                                                                       
!     step through lines                                                
      nlpl=1 
      do lnn=1,nlsvn 
!                                                                       
!       get line data                                                   
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
!       get ion data                                                  
        nilin=derivedpointers%npar(ml) 
        mlm=nilin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
        do ktt=1,min(8,nkdt) 
          write (kinam1(ktt:ktt),'(a1)')masterdata%kdat1(np1k-1+ktt) 
          enddo 
        do ktt=nkdt+1,9 
          write (kinam1(ktt:ktt),'(a1)')kblnk 
          enddo 
!                                                                       
!          if (lpri.ge.1)                                               
!     $      write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,              
!     $          (kdat1(np1ki+mm-1),mm=1,nkdti)                         
!                                                                       
        nelin=derivedpointers%npar(nilin) 
        mlm=nelin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                       &
     &      0,lun11)                                                  
        nnzel=masterdata%idat1(np1i2)
        a=masterdata%rdat1(np1r2+1) 
!
!       exclude rate type 14                                            
        if ((lrtyp.ne.14).and.(abs(elin).gt.0.1).and.(lrtyp.ne.9)       &
     &       .and.(abs(elin).lt.9.e+9)) then
!                                                                       
!         get line data                                                   
          ln=lnn 
          ml=derivedpointers%nplin(ln) 
          mlm=ml
          call drd(ltyp,lrtyp,lcon,                                     &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
          elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
          ergsev=1.602197e-12 
          ener=ergsev*(12398.41)/max(elin,1.d-24) 
          etst=ener/ergsev 
          idest1=masterdata%idat1(np1i) 
          idest2=masterdata%idat1(np1i+1) 
          aij=masterdata%rdat1(np1r+2) 
          if (ltyp.eq.82) aij=masterdata%rdat1(np1r+3) 
!                                                                       
!                                                                       
!         now find level data                                           
          jkk=masterdata%idat1(np1i+nidt-1) 
          call calc_rates_level_lte(jkk,lpril,lun11,t,xee,xpx,          &
     &              nlev)
          lpri=0
          call deleafnd(jkk,idest2,delea,lfnd,lpri,lun11)           
          delearad=6.6262e-27*aij/6.28/1.602197e-12
          t1=1.        
          vtherm=1.29e+6/sqrt(a/t1)
          deleth=etst*vtherm/3.e+10 
          aasmall=(delea+delearad)/(1.E-24+deleth)/12.56 
                                                                        
          ggup=leveltemp%rlev(2,idest2) 
          gglo=leveltemp%rlev(2,idest1) 
          do lk=1,20 
            klablo(lk)=leveltemp%klev(lk,idest1) 
            klabup(lk)=leveltemp%klev(lk,idest2) 
            enddo 
          flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo) 
          write (lun11,9944)lnn,elin,kinam1,aij,flin,gglo,ggup,         &
     &             klablo,klabup,delea,deleth,aasmall 
 9944     format (1h ,i9,e12.4,1x,a9,4(1pe12.4),1x,20a1,1x,20a1,        &
     &        3(1pe12.4)) 
!                                                                       
          endif 
        enddo 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   10 continue 
!                                                                       
      write (lun11,*)'ion abundances and  rates (/sec)' 
      write (lun11,947) 
  947 format (1x,'index, ion, abundance, recombination, ionization,') 
!                                                                       
!     step thru ions                                                    
      klion=12 
      mlion=derivedpointers%npfirst(klion) 
      lk=0 
      do while (mlion.ne.0) 
!                                                                       
!        get ion data                                                   
         lk=lk+1 
         ltyp=klion 
         mlm=mlion
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         do mm=1,nkdt 
           kdtmp(mm)=masterdata%kdat1(np1k-1+mm) 
           enddo 
         do mm=nkdt+1,9 
           kdtmp(mm)=kblnk 
           enddo 
!                                                                       
!        get element data                                               
         nell=derivedpointers%npar(mlion) 
         mlm=nell
         call drd(ltyp,lrtyp,lcon,                                      &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
!         write (lun11,*)mlion,lk,np1i,nidt,np1i+nidt-1,                
!     $      idat1(np1i+nidt-1),ababs(idat1(np1i+nidt-1)),mlm           
         if ((masterdata%idat1(np1i+nidt-1).gt.0)                       &
     &     .and.(masterdata%idat1(np1i+nidt-1).le.nl)) then             
           abundel=ababs(masterdata%idat1(np1i+nidt-1)) 
!                                                                       
!          print out                                                    
           if (abundel.gt.1.e-15)                                       &
     &      write (lun11,9046)lk,(kdtmp(mm),mm=1,9),                    &
     &      xii(lk),rrrt(lk),pirt(lk)                                   
 9046      format (1x,i4,1x,9a1,5(1pe16.8)) 
!                                                                       
           endif 
!                                                                       
         mlion=derivedpointers%npnxt(mlion) 
         enddo 
!                                                                       
      write (lun11,*)'heating and cooling rates (erg/sec)' 
      write (lun11,9947) 
 9947 format (1x,'index, element, heating, cooling: ') 
!                                                                       
!     step thru elements                                                
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      lk=0 
      do while (mlel.ne.0) 
!                                                                       
!        get element data                                               
         ltyp=klel 
         mlm=mlel
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         do mm=1,nkdt 
           kdtmp(mm)=masterdata%kdat1(np1k-1+mm) 
           enddo 
         do mm=nkdt+1,9 
           kdtmp(mm)=kblnk 
           enddo 
         lk=masterdata%idat1(np1i) 
!                                                                       
!                                                                       
         if ((masterdata%idat1(np1i+nidt-1).gt.0)                       &
     &     .and.(masterdata%idat1(np1i+nidt-1).le.nl)) then                
!                                                                       
           abundel=ababs(masterdata%idat1(np1i+nidt-1))                  
           if (abundel.gt.1.d-36)                                       &
     &      write (lun11,9046)lk,(kdtmp(mm),mm=1,9),htt(lk),cll(lk),    &
     &          htt2(lk),cll2(lk)
!                                                                       
           endif 
!                                                                       
         mlel=derivedpointers%npnxt(mlel) 
         enddo 
!                                                                       
      write (lun11,*)'total heating, cooling:',                         &
     &            httot,cltot                                           
      write (lun11,*)'partial heating rates: photo,compton',            &
     &            httot-htcomp,htcomp                                   
      write (lun11,*)'partial cooling rates: rec,lines,brems,compton',  &
     &            clcont,cllines,clcomp,clbrems                         
      write (lun11,*)'total heating, cooling, electron point of view:', &
     &            httot2,cltot2                                           
      write (lun11,993) 
!                                                                       
!                                                                       
      go to 9000
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!        
!                                                                       
!     Write FITS file with summary of ion abundances                    
!     by radial zone                                                    
!                                                                       
   11 continue 
!                                                                       
      knam='xout_abund1.fits' 
      call fheader(unit,knam,atcredate,kmodelname,status) 
                                                                        
!      write (lun11,*)'in pprint, 11: numrec=',numrec                   
!      call writeimage(knam)                                            
      do mm=1,999 
        kunits(mm)=kblnk16 
        klabs(mm)=kblnk16 
        kform(mm)=kblnk16 
        enddo 
      klabs(1)='radius          ' 
      kform(1)='E11.3' 
      kunits(1)='cm' 
      klabs(2)='delta_r         ' 
      kform(2)='E11.3' 
      kunits(2)='cm' 
      klabs(3)='ion_parameter   ' 
      kform(3)='E11.3' 
      kunits(3)='erg*cm/s' 
      klabs(4)='x_e             ' 
      kform(4)='E11.3' 
      kunits(4)=' ' 
      klabs(5)='n_p             ' 
      kform(5)='E11.3' 
      kunits(6)='cm**(-3)' 
      klabs(6)='pressure        ' 
      kform(6)='E11.3' 
      kunits(6)='dynes/cm**2' 
      klabs(7)='temperature     ' 
      kform(7)='E11.3' 
      kunits(7)='10**4 K' 
      klabs(8)='frac_heat_error' 
      kform(8)='E11.3' 
      kunits(8)=' ' 
!     Search for the ion names in the database                          
      klion=12 
      mlion=derivedpointers%npfirst(klion) 
      do lkk=1,nni 
        xcoltmp(lkk)=0. 
        enddo 
      lk=0 
      do while (mlion.ne.0) 
           lk=lk+1 
           ltyp=klion 
           mlm=mlion
           call drd(ltyp,lrtyp,lcon,                                    &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
!                                                                       
!          get element abundance                                        
           nelin=derivedpointers%npar(mlion) 
           ml=nelin 
           mlm=ml
           call drd(ltyp,lrtyp,lcon,                                    &
     &       nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                         &
     &       0,lun11)                                             
           mllel=masterdata%idat1(np1i+nidt-1) 
           xeltp=ababs(mllel) 
!                                                                       
!          go back to ion data                                          
           mlm=mlion
           call drd(ltyp,lrtyp,lcon,                                    &
     &       nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                         &
     &       0,lun11)                                             
                                                                        
!          Compute string length from character array by searching backw
!          for first non-blank character                                
           klen=0 
           do mm=1,nkdt 
             kdat(mm)=masterdata%kdat1(np1k-1+mm) 
             enddo 
           do mm=nkdt+1,9 
             kdat(mm)=kblnk 
             enddo 
           do mm=1,9 
             if(kdat(10-mm).ne.' '.and.klen.eq.0) then 
                klen=10-mm 
             endif 
             enddo 
!           write (lun11,*)'kdat:',(kdat(mm),mm=1,9)                    
!           write (lun11,*)'klen:',klen                                 
!          Replace ' ' in ion names to '_' to match FITS standard       
           do mm=1,9 
             if(kdat(mm).eq.' '.and.mm.lt.klen) then 
                write (ktmp(mm:mm),'(a1)')'_' 
             else 
                write (ktmp(mm:mm),'(a1)')kdat(mm) 
             endif 
             enddo 
           do mm=10,16 
             write (ktmp(mm:mm),'(a1)')' ' 
             enddo 
           do jkl=2,numrec 
!             write (lun11,*)jkl,lk,zrtmp(2,jkl),zrtmp(8+lk,jkl),       
!     $                             zrtmp(5,jkl),xeltp,xcoltmp(lk)      
             xcoltmp(lk)=xcoltmp(lk)                                    &
     &         +(zrtmp(8+lk,jkl)*zrtmp(5,jkl)                           &
     &             +zrtmp(8+lk,jkl-1)*zrtmp(5,jkl-1))                   &
     &         *(zrtmp(2,jkl)-zrtmp(2,jkl-1))*xeltp/2.                  
             enddo 
           klabs(8+lk)=ktmp 
           kform(8+lk)='E11.3' 
           kunits(8+lk)=' ' 
           mlion=derivedpointers%npnxt(mlion) 
           enddo 
!                                                                       
      call fwrtascii(unit,'ABUNDANCES',zrtmp,8+lk,                   &
     &                  numrec,klabs,kform,kunits,lun11)                
!                                                                       
!     calculate columns                                                 
      numrec=1 
      do lkk=1,lk 
        zrtmpcol(8+lkk,1)=xcoltmp(lkk) 
        enddo 
      do lkk=1,8 
        zrtmpcol(lkk,1)=0. 
        enddo 
!                                                                       
      call fwrtascii(unit,'COLUMNS                                 ',&
     & zrtmpcol,8+lk,1,klabs,kform,kunits,lun11)                
!                                                                       
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      lk=0 
      do while (mlel.ne.0) 
           lk=lk+1 
           ltyp=klel 
           mlm=mlel
           call drd(ltyp,lrtyp,lcon,                                    &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
!                                                                       
!          Compute string length from character array by searching backw
!          for first non-blank character                                
           klen=0 
           do mm=1,nkdt 
             kdat(mm)=masterdata%kdat1(np1k-1+mm) 
             enddo 
           do mm=nkdt+1,9 
             kdat(mm)=kblnk 
             enddo 
           do mm=1,9 
             if(kdat(10-mm).ne.' '.and.klen.eq.0) then 
                klen=10-mm 
             endif 
             enddo 
           do mm=1,9 
             if(kdat(mm).eq.' '.and.mm.lt.klen) then 
                write (ktmp(mm:mm),'(a1)')'_' 
             else 
                write (ktmp(mm:mm),'(a1)')kdat(mm) 
             endif 
             enddo 
           do mm=10,16 
             write (ktmp(mm:mm),'(a1)')' ' 
             enddo 
           klabs(8+lk)=ktmp 
           kform(8+lk)='E11.3' 
           kunits(8+lk)=' ' 
           mlel=derivedpointers%npnxt(mlel) 
           enddo 
      klabs(8+lk+1)='compton' 
      kform(8+lk+1)='E11.3' 
      kunits(8+lk+1)=' ' 
      klabs(8+lk+2)='total' 
      kform(8+lk+2)='E11.3' 
      kunits(8+lk+2)=' ' 
      call fwrtascii(unit,'HEATING                                ', &
     & zrtmph,8+lk+2,numrec,klabs,kform,kunits,lun11)                   
!                                                                       
      klabs(8+lk+1)='compton' 
      kform(8+lk+1)='E11.3' 
      kunits(8+lk+1)=' ' 
      klabs(8+lk+2)='brems' 
      kform(8+lk+2)='E11.3' 
      kunits(8+lk+2)=' ' 
      klabs(8+lk+3)='total' 
      kform(8+lk+3)='E11.3' 
      kunits(8+lk+3)=' ' 
      call fwrtascii(unit,'COOLING                                 ',&
     & zrtmpc,8+lk+3,numrec,klabs,kform,kunits,lun11)                   
!                                                                       
      call fitsclose(lun11,unit,status) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   12 continue 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!        
!                                                                       
!     Add ionic abundances info in this radial zone to array for        
!     eventual inclusion in xout_abund1.fits                            
!     Modifies zrtmp                                                    
!                                                                       
      ergsev=1.602197e-12 
      r19=r*(1.e-19)
!      write (lun11,*)enlum,xpx,r,xlum,jkstep                           
      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10 
!      write (lun11,*)uu1                                               
      alguu1=log10(max(1.d-24,uu1)) 
      skse=xlum/(xpx*r19*r19) 
      zeta=log10(max(1.d-24,skse)) 
      ecc=2.998e+10 
      ekt=t*(0.861707)*ergsev 
!      sksec=skse/(12.56*((1.+xee)*ekt+pradl/(1.e-24+xpx))*ecc)         
      sksec=skse/12.56/((1.+xee)*ekt*ecc) 
      zetac=log10(max(1.d-24,sksec)) 
      enn0=xpx 
      nlyc=nbinc(13.7d0,epi,ncn2) 
      nry=nlyc+1 
      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24) 
      nry=nbinc(13.6d0,epi,ncn2)+1 
!     Copy the values for radial zone jkstep                            
      if (jkstep.gt.3999) go to 9000
      zrtmp(1,jkstep)=r 
      zrtmp(2,jkstep)=rdel 
      zrtmp(3,jkstep)=zeta 
      zrtmp(4,jkstep)=xee 
      zrtmp(5,jkstep)=xpx 
      zrtmp(6,jkstep)=p 
      zrtmp(7,jkstep)=t 
      zrtmp(8,jkstep)=hmctot 
      do lk=1,8 
        zrtmpc(lk,jkstep)=zrtmp(lk,jkstep) 
        zrtmph(lk,jkstep)=zrtmp(lk,jkstep) 
        enddo 
      klion=12 
      mlion=derivedpointers%npfirst(klion) 
      lk=0 
      do while (mlion.ne.0) 
        lk=lk+1 
        zrtmp(8+lk,jkstep)=xii(lk) 
        mlion=derivedpointers%npnxt(mlion) 
        enddo 
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      lk=0 
      do while (mlel.ne.0) 
        lk=lk+1 
        zrtmpc(8+lk,jkstep)=cll(lk) 
        zrtmph(8+lk,jkstep)=htt(lk) 
        mlel=derivedpointers%npnxt(mlel) 
        enddo 
      zrtmph(8+lk+1,jkstep)=htcomp 
      zrtmph(8+lk+2,jkstep)=httot 
      zrtmpc(8+lk+1,jkstep)=clcomp 
      zrtmpc(8+lk+2,jkstep)=clbrems 
      zrtmpc(8+lk+3,jkstep)=cltot 
!                                                                       
      go to 9000
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!      
!                                                                       
!                                                                       
   13 continue 
!                                                                       
  993 format (1h ) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   14 continue 
!                                                                       
      write (lun11,900) 
  900 format ('line opacities and emissivities',                        &
     & ' (erg/cm**3/sec/10**38)')                                       
      write (lun11,915) 
  915 format (1x,'index,wavelength,energy,ion,opacity,rec. em.,',       &
     &'coll. em.,fl. em.,di. em.,cx. em.')                              
!                                                                       
!     step through lines                                                
      nlpl=1 
      write (lun11,*)nlsvn 
      do lnn=1,nlsvn 
!                                                                       
!       get line data                                                   
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        elin=abs(masterdata%rdat1(np1r)) 

!       get ion data                                                  
        nilin=derivedpointers%npar(ml) 
        mlm=nilin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
        nilin2=masterdata%idat1(np1i+nidt-1) 

!       get element data                                                  
        nelin=derivedpointers%npar(nilin) 
        mlm=nelin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=abel(mllel) 

        if ((lrtyp.ne.14).and.(lrtyp.ne.9).and.(abs(elin).gt.0.1)       &
     &       .and.(abs(elin).lt.9.e+9).and.(xeltp.gt.1.d-36)) then             
                                                                        
          ergsev=1.602197e-12 
          ener=ergsev*(12398.41)/max(elin,1.d-24) 
          etst=ener/ergsev 
          idest1=masterdata%idat1(np1i) 
          idest2=masterdata%idat1(np1i+1) 
          aij=masterdata%rdat1(np1r+2) 
!                                                                       
!         get ion data                                                  
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
           do ktt=1,min(8,nkdt) 
            write (kinam1(ktt:ktt),'(a1)')masterdata%kdat1(np1k-1+ktt) 
            enddo 
          do ktt=nkdt+1,9 
            write (kinam1(ktt:ktt),'(a1)')kblnk 
            enddo 
!                                                                       
          j=ln 
          write (lun11,904)j,elin,etst,kinam1,oplin(j),rcem(1,j),       &
     &                      rcem(2,j)                                   
  904     format (1h ,i9,2(1pe13.5),1x,a9,6(1pe13.5)) 
!                                                                       
          endif 
        enddo 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   15 continue 
!                                                                       
      write (lun11,*)'line luminosities (erg/sec/10**38) and depths' 
      write (lun11,9923) 
 9923 format (1x,' line, wavelength, ion, ref. lum.,trn. lum.,',        &
     &'backward depth, forward depth')                                  
!     step through lines                                                
      nlpl=1 
      do lnn=1,nlsvn 
!                                                                       
!       get line data                                                   
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt2,np1k2,mlm,                          &
     &    0,lun11)                                                
!        call dprinto(ltyp,lrtyp2,lcon,                                 
!     $    nrdt,np1r,nidt,np1i,nkdt2,np1k2,rdat1,idat1,kdat1,lun11)     
        elin=abs(masterdata%rdat1(np1r)) 

!       get ion data                                                  
        nilin=derivedpointers%npar(ml) 
        mlm=nilin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
        nilin2=masterdata%idat1(np1i+nidt-1) 

!       get element data                                                  
        nelin=derivedpointers%npar(nilin) 
        mlm=nelin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=abel(mllel) 

        if ((lrtyp.ne.14).and.(lrtyp.ne.9).and.(abs(elin).gt.0.1)       &
     &       .and.(abs(elin).lt.9.e+9).and.(xeltp.gt.1.d-36)) then            
                                                                        
          ergsev=1.602197e-12 
          ener=ergsev*(12398.41)/max(elin,1.d-24) 
          etst=ener/ergsev 
          idest1=masterdata%idat1(np1i) 
          idest2=masterdata%idat1(np1i+1) 
          aij=masterdata%rdat1(np1r+2) 
!                                                                       
!         get ion data                                                  
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          do ktt=1,min(8,nkdt) 
            write (kinam1(ktt:ktt),'(a1)')masterdata%kdat1(np1k-1+ktt) 
            enddo 
          do ktt=nkdt+1,9 
            write (kinam1(ktt:ktt),'(a1)')kblnk 
            enddo 
!                                                                       
          j=ln 
          elmtp=elum(1,j) 
          elmtpb=elum(2,j) 
          write (lun11,9924)j,elin,kinam1,                              &
     &       elmtp,elmtpb,tau0(1,j), tau0(2,j),                         &
     &      (masterdata%kdat1(mm),mm=np1k2,np1k2+nkdt2-1)                       
 9924     format (1h ,i9,1pe13.5,1x,a9,1x,4(1pe13.5),20a1) 
!                                                                       
          endif 
        enddo 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   18 continue 
!                                                                       
      lpril=0 
      write (lun11,*)'line wavelengths and levels' 
!     step through lines                                                
      nlpl=1 
      do lnn=1,nlsvn 
!                                                                       
!       get line data                                                   
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
!       get ion data                                                  
        nilin=derivedpointers%npar(ml) 
        mlm=nilin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
        nilin2=masterdata%idat1(np1i+nidt-1) 

!       get element data                                                  
        nelin=derivedpointers%npar(nilin) 
        mlm=nelin
        call drd(ltyp,lrtyp,lcon,                                       &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=abel(mllel) 

!       exclude rate type 14                                            
        if ((lrtyp.ne.14).and.(lrtyp.ne.9).and.(abs(elin).gt.0.1)       &
     &       .and.(abs(elin).lt.9.e+9).and.(xeltp.gt.1.d-36)) then             
!                                                                       
!         get line data                                                   
          ln=lnn 
          ml=derivedpointers%nplin(ln) 
          mlm=ml
          call drd(ltyp,lrtyp,lcon,                                     &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
          elin=abs(masterdata%rdat1(np1r)) 
!                                                                       
          ergsev=1.602197e-12 
          ener=ergsev*(12398.41)/max(elin,1.d-24) 
          etst=ener/ergsev 
          idest1=masterdata%idat1(np1i) 
          idest2=masterdata%idat1(np1i+1) 
          aij=masterdata%rdat1(np1r+2) 
!                                                                       
!         get ion data                                                  
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          do ktt=1,min(8,nkdt) 
            write (kinam1(ktt:ktt),'(a1)')masterdata%kdat1(np1k-1+ktt) 
            enddo 
          do ktt=nkdt+1,9 
            write (kinam1(ktt:ktt),'(a1)')kblnk 
            enddo 
!                                                                       
          if (lpril.ge.1)                                               &
     &      write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,               &
     &         (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)                       
!                                                                       
!         now find level data                                           
          jkk=masterdata%idat1(np1i+nidt-1) 
          call calc_rates_level_lte(jkk,lpril,lun11,t,xee,xpx,          &
     &              nlev)
                                                                        
          ggup=leveltemp%rlev(2,idest1) 
          gglo=leveltemp%rlev(2,idest2) 
          do lk=1,20 
            klablo(lk)=leveltemp%klev(lk,idest1) 
            klabup(lk)=leveltemp%klev(lk,idest2) 
            enddo 
          flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo) 
          ilevlo=idest1 
          ilevup=idest2 
!                                                                       
          j=ln 
            if (lpri.ne.0)                                              &
     &        write (lun11,9929)j,elin,kinam1,                          &
     &        (leveltemp%klev(mm,ilevlo),mm=1,20),                      &
     &        (leveltemp%klev(mm,ilevup),mm=1,20),                      &
     &        (leveltemp%rlev(mm,ilevlo),leveltemp%rlev(mm,ilevup),     &
     &         mm=1,3),                                                 &
     &        (leveltemp%ilev(mm,ilevlo),leveltemp%ilev(mm,ilevup),     &
     &         mm=1,3)
 9929       format (1h ,i9,1pe13.5,1x,a9,1x,2(20a1,1x),6(1pe13.5),      &
     &          6i6)                                                    
!                                                                       
          endif 
        enddo 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   20 continue 
!                                                                       
      lpril=0 
      write (lun11,*)'line finding list' 
      do jlk=1,nlsvn 
         j=jlk 
         ml=derivedpointers%nplin(j) 
         mlm=ml
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         elin=abs(masterdata%rdat1(np1r)) 
         jpnt(j)=j 
         elsv(j)=abs(elin) 
         enddo 
!     sort                                                              
      done=.false. 
      niter=0 
      do while (.not.done) 
        niter=niter+1 
        done=.true. 
!        do jjj=1,100                                                   
        do jjj=1,nlsvn-1 
          j=jpnt(jjj) 
          jp1=jpnt(jjj+1) 
!          write (lun11,*)jjj,j,jp1,                                    
!     $           elsv(jp1),elsv(j)                                     
          if (elsv(jp1).lt.elsv(j)) then 
            jpnt(jjj)=jp1 
            jpnt(jjj+1)=j 
            done=.false. 
            endif 
          enddo 
        enddo 
!                                                                       
!     print out sorted list                                             
      do jlk=1,nlsvn 
        j=jpnt(jlk) 
        lnn=j 
!                                                                       
!       get line data                                                   
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
!                                                                       
!       exclude rate type 14                                            
        if ((lrtyp.ne.14).and.(lrtyp.ne.9).and.(abs(elin).gt.0.1)       &
     &       .and.(abs(elin).lt.9.e+9)) then                            
!                                                                       
          elin=abs(masterdata%rdat1(np1r)) 
          ergsev=1.602197e-12 
          ener=ergsev*(12398.41)/max(elin,1.d-24) 
          etst=ener/ergsev 
          idest1=masterdata%idat1(np1i) 
          idest2=masterdata%idat1(np1i+1) 
          aij=masterdata%rdat1(np1r+2) 
!                                                                       
!         get ion data                                                  
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          do ktt=1,min(8,nkdt) 
            write (kinam1(ktt:ktt),'(a1)')masterdata%kdat1(np1k-1+ktt) 
            enddo 
          do ktt=nkdt+1,9 
            write (kinam1(ktt:ktt),'(a1)')kblnk 
            enddo 
!                                                                       
          if (lpril.ge.1)                                               &
     &      write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,               &
     &        (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)                          
!                                                                       
!         now find level data                                           
          jkk=masterdata%idat1(np1i+nidt-1) 
          call calc_rates_level_lte(jkk,lpril,lun11,t,xee,xpx,          &
     &              nlev)
                                                                        
          ggup=leveltemp%rlev(2,idest1) 
          gglo=leveltemp%rlev(2,idest2) 
          do lk=1,20 
            klablo(lk)=leveltemp%klev(lk,idest1) 
            klabup(lk)=leveltemp%klev(lk,idest2) 
            enddo 
          flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo) 
          ilevlo=idest1 
          ilevup=idest2 
!                                                                       
            if (lpri.ne.0)                                              &
     &        write (lun11,9929)j,elin,kinam1,                          &
     &        (leveltemp%klev(mm,ilevlo),mm=1,20),                      &
     &        (leveltemp%klev(mm,ilevup),mm=1,20),                      &
     &        (leveltemp%rlev(mm,ilevlo),leveltemp%rlev(mm,ilevup),     &
     &         mm=1,3),                                                 &
     &        (leveltemp%ilev(mm,ilevlo),leveltemp%ilev(mm,ilevup),     &
     &         mm=1,3)
           endif 
         enddo 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
   21 continue 
!                                                                       
      lpril=0 
      write (lun11,*)' level opacities and emissivities' 
      write (lun11,*)'index,energy,ion,level,index,emiss in,emiss out,th&
     &reshold opacity,absorbed energy,depth in, depth out'              
!                                                                       
!     First look for element data (jk is element index)                 
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      jk=0 
      kk=0 
      jkk=0 
!                                                                       
!     step through elements                                             
      do while (mlel.ne.0) 
!                                                                       
!       get element data                                                
        jk=jk+1 
        mt2=mlel
        call drd(ltyp,lrtyp,lcon,                                       &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                        &
     &        0,lun11)                                            
        mllel=masterdata%idat1(np1i+nidt-1) 
        xeltp=masterdata%rdat1(np1r) 
        xeltp=abel(mllel) 
        nnz=masterdata%idat1(np1i) 
        if (lpril.ge.1) then
             write (lun11,902)jk,mlel,nnz,                              &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,nkdt),xeltp                    
          endif
!                                                                       
!       ignore if the abundance is small                                
        if (xeltp.lt.1.e-10) then 
            jkk=jkk+nnz 
          else 
!                                                                       
!           now step thru ions (jkk is ion index)                       
            klion=12 
            mlion=derivedpointers%npfirst(klion) 
            jkk=0 
            kl=0 
            do while ((mlion.ne.0).and.(kl.lt.nnz)) 
              jkk=jkk+1 
!                                                                       
!             retrieve ion name from kdati                              
              mlm=mlion
              call drd(ltyp,lrtyp,lcon,                                 &
     &            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,                 &
     &            0,lun11)                                        
!                                                                       
!             if not accessing the same element, skip to the next elemen
              mlleltp=masterdata%idat1(np1i+nidti-2) 
              if (mlleltp.eq.mllel) then 
!                                                                       
                kl=kl+1 
                if (lpril.ge.1)                                         &
     &            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,         &
     &               (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)            
!                                                                       
!               now find level data                                     
                jkk=masterdata%idat1(nidt+np1i-1)
                call calc_rates_level_lte(jkk,lpril,lun11,t,xee,xpx,    &
     &              nlev)
!                                                                       
!               now step through rate type 7 data                       
                mltype=7 
                ml=derivedpointers%npfi(mltype,jkk) 
                mllz=0 
                if (ml.ne.0) mllz=derivedpointers%npar(ml) 
                mlpar=0 
                if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                do while ((ml.ne.0).and.(mlpar.eq.mllz)) 
!                                                                       
!                 get rrc data                                          
                  kkkl=derivedpointers%npconi2(ml) 
                  if (lpril.ne.0) write (lun11,*)kkkl,ml,idest1,        &
     &                    elumab(1,kkkl),elumab(2,kkkl)                 
!                                                                       
!                 test for non-zero rrc data                            
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)                   &
     &                .and.((opakab(kkkl).gt.1.d-49)                    &
     &                .or.(cabab(kkkl).gt.1.d-49)                       &
     &                .or.(cemab(1,kkkl).gt.1.d-49)                     &
     &                .or.(cemab(2,kkkl).gt.1.d-49))) then             
!                                                                       
!                   get rrc  data                                       
                    mlm=ml
                    call drd(ltyp,lrtyp,lcon,                           &
     &                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                &
     &                0,lun11)                                    
                    idest1=masterdata%idat1(np1i+nidt-2) 
                    nlevp=nlev 
                    idest2=nlevp+masterdata%idat1(np1i-1+nidt-3)-1 
!                                                                       
!                   label for lower level                               
                    do lk=1,20 
                      write (ktmp20(lk:lk),'(a1)')                      &
     &                    leveltemp%klev(lk,idest1) 
                      enddo 
                    klevl=ktmp20 
!                                                                       
!                   label for upper level                               
                    write (ktmp20(1:20),'(a20)')'continuum           ' 
                    klevu=ktmp20 
!                                                                       
!                   ion label                                           
                    do lk=1,nkdti 
                      write (ktmp8(lk:lk),'(a1)')                       &
     &                      masterdata%kdat1(np1ki+lk-1) 
                      enddo 
                    do lk=nkdti+1,8 
                      write (ktmp8(lk:lk),'(a1)')kblnk 
                      enddo 
!                                                                       
                    eth=leveltemp%rlev(4,idest1)                        &
     &                  -leveltemp%rlev(1,idest1) 
                    ett=eth 
!                                                                       
!                   get upper level data                                
                    go to 9087
                    if (idest2.gt.nlevp) then 
                      jkk3=jkk+1 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*)jkk3,ndtmp,nlevp,idest2          
                      ndtmp=derivedpointers%npfi(13,jkk3) 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*)jkk3,ndtmp,nlevp,idest2          
                      if (ndtmp.le.0) stop 'ndtmp error' 
                      mllz=derivedpointers%npar(ndtmp) 
                      iltmp=0 
                      do while ((ndtmp.ne.0).and.                       &
     &                    (iltmp.ne.(idest2-nlevp+1)).and.              &
     &                    (derivedpointers%npar(ndtmp).eq.mllz))               
                        mlm=ndtmp
                        call drd(ltyp2,lrtyp2,lcon2,                    &
     &                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,      &
     &                    0,lun11)                                
                        iltmp=masterdata%idat1(np1i2+nidt2-2) 
                        if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
                        ndtmp=derivedpointers%npnxt(ndtmp) 
                        if (ndtmp.le.0) stop 'ndtmp error' 
                        enddo 
!                     NB fix to excited level PI and rec                
                      ett=ett+masterdata%rdat1(np1r2) 
                      eth=ett 
                      if (lpril.gt.1)                                   &
     &                  write (lun11,*) ndtmp,iltmp,idest2,ett          
!                     label for lower level                             
                      ktmp20=kblnk20 
                      do lk=1,nkdt2 
                         write (ktmp20(lk:lk),'(a1)')                    &
     &                         masterdata%kdat1(np1k2+lk-1) 
                        enddo 
                      klevu=ktmp20 
                      endif 
9087                  continue
!                                                                       
!                   other data                                          
                    mmlv=derivedpointers%npilev(idest1,jkk) 
                    mlcu=kkkl 
                    write (lun11,969)kkkl,mmlv,ktmp8,idest1,idest2,     &
     &                  klevl,klevu,eth,                                &
     &                  cemab(1,mlcu),cemab(2,mlcu),opakab(mlcu),       &
     &                  cabab(mlcu),tauc(1,mlcu),tauc(2,mlcu)           
  969               format (1x,2i6,1x,a8,1x,2i6,1x,2(a20,1x),           &
     &                  8(1pe13.5),2i6)                                 
!                                                                       
!                   done with this rrc                                  
                    endif 
!                                                                       
!                 end of loop over rrcs                                 
                  ml=derivedpointers%npnxt(ml) 
                  if (ml.ne.0) mlpar=derivedpointers%npar(ml) 
                  enddo 
!                                                                       
!               end of test for element                                 
                endif 
!                                                                       
!             Go to next ion                                            
              mlion=derivedpointers%npnxt(mlion) 
              enddo 
!                                                                       
!         end of test for non-zero element abund                        
          endif 
!                                                                       
        mlel=derivedpointers%npnxt(mlel) 
!       Go to next element                                              
        enddo 
!                                                                       
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
    7 continue 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!       
!                                                                       
!     Write level populations                                           
!                                                                       
      write (lun11,9985) 
 9985 format (1x,' level populations ') 
      write (lun11,9986) 
 9986 format (1x,' ion                      level              '        &
     &,' e_exc population')                                             
                                                                        
!     lpril is flag for printing debug information                      
      lpril=0
      if (lpril.ne.0) then 
        write (lun11,*)'raw data' 
        do j=1,nnml 
!          if (xilev(j).gt.1.e-37)                                       &
          write (lun11,*)j,xilev(j),bilev(j),rnist(j),elumab(1,j)                        
          enddo 
        endif 
!                                                                       
!     First look for element data (jk is element index)                 
      klel=11 
      mlel=derivedpointers%npfirst(klel) 
      jk=0 
!                                                                       
!     step through elements                                             
      do while (mlel.ne.0) 
!                                                                       
!       get element data                                                
        jk=jk+1 
        mt2=mlel
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                            &
     &    0,lun11)                                                
        mllel=masterdata%idat1(np1i+nidt-1) 
        nnz=masterdata%idat1(np1i) 
        xeltp=masterdata%rdat1(np1r) 
        xeltp=abel(mllel) 
        if (lpril.ne.0)                                                 &
     &        write (lun11,*)'element:',jk,mlel,mllel,nnz,              &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,nkdt),xeltp            
!                                                                       
!       ignore if the abundance is small                                
        if (xeltp.lt.1.e-10) then 
            jkk=jkk+nnz 
          else 
!                                                                       
!           now step thru ions (jkk is ion index)                       
            klion=12 
            mlion=derivedpointers%npfirst(klion) 
            jkk=0 
            kl=0 
            do while ((mlion.ne.0).and.(kl.lt.nnz)) 
!                                                                       
              jkk=jkk+1 
!             retrieve ion name from kdati                              
              mlm=mlion
              call drd(ltyp,lrtyp,lcon,                                 &
     &            nrdt,np1r,nidt,np1i,nkdti,np1ki,mlm,                  &
     &            0,lun11)                                        
!                                                                       
!             if not accessing the same element, skip to the next elemen
              mlleltp=masterdata%idat1(np1i+nidt-2) 
              if (mlleltp.eq.mllel) then 
!                                                                       
                kl=kl+1 
                if (lpril.ne.0)                                         &
     &            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,         &
     &               (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)            
                do ktt=1,min(8,nkdti) 
                   write (kinam1(ktt:ktt),'(a1)')                       &
     &                   masterdata%kdat1(np1ki-1+ktt) 
                  enddo 
                do ktt=nkdti+1,9 
                  write (kinam1(ktt:ktt),'(a1)')kblnk 
                  enddo 
!                                                                       
!               get level data                                          
                jkk=masterdata%idat1(nidt+np1i-1)
                call calc_rates_level_lte(jkk,lpril,lun11,t,xee,xpx,    &
     &              nlev)
!                                                                       
!               step thru levels                                        
                do mm2=1,nlev 
!                                                                       
!                 get level pointer                                     
                  mmtmp=derivedpointers%npilev(mm2,jkk) 
                  if (mmtmp.ne.0) then 
                    kkkl=mmtmp 
                    mmlv=mmtmp 
!                                                                       
!                   test for level pop                                  
                    if (xilev(kkkl).gt.1.d-64) then 
!                                                                       
!                     get data                                          
                      eth=leveltemp%rlev(1,mm2) 
                      dep=xilev(kkkl)/(rnist(kkkl)+1.d-36) 
                      write (lun11,9296)kkkl,kinam1,                    &
     &                   (leveltemp%klev(lk,mm2),lk=1,20),              &
     &                   eth,xilev(kkkl),rnist(kkkl),dep,bilev(kkkl)                                
 9296                 format (1x,i6,1x,a8,1x,(20a1),7(1pe13.5)) 
!                                                                       
                                                                        
!                     end of test for level pop                         
                      endif 
                                                                        
!                   end of test for level pointer                       
                    endif 
!                                                                       
!                 end of step thru levels                               
                  enddo 
!                                                                       
!               end of test for element                                 
                endif 
!                                                                       
!             Go to next ion                                            
              mlion=derivedpointers%npnxt(mlion) 
              enddo 
!                                                                       
!           end of test for abundance                                   
            endif 
!                                                                       
!       Go to next element                                              
        if (mlel.ne.0) mlel=derivedpointers%npnxt(mlel) 
        enddo 
!                                                                       
      write (lun11,*)'done with 7' 
      write (lun11,993) 
!                                                                       
!                                                                       
      go to 9000
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!   
!                                                                       
!     Print a short summary line of the radial calculation              
!                                                                       
    9 continue 
!                                                                       
      elsum=0. 
      ergsev=1.602197e-12 
      do jlk=1,nlsvn 
         ln=jlk 
         ml=derivedpointers%nplin(ln) 
         mlm=ml 
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         elin=abs(masterdata%rdat1(np1r)) 
         if ((elin.lt.1.e+8).and.(elin.gt.1.)) then 
           elsum=elsum+elum(1,ln)+elum(2,ln) 
           endif 
         enddo 
      sumtmp1=0. 
      sumtmp2=0. 
      ergsev=1.602197e-12 
      r19=r*1.e-19
      tmp1=zremsz(1) 
      tmp2=zrems(1,1) 
      do jk=2,ncn2 
         tmp1o=tmp1 
         tmp1=zremsz(jk) 
         sumtmp1=sumtmp1+(tmp1+tmp1o)*(epi(jk)-epi(jk-1))*ergsev/2. 
         tmp2o=tmp2 
         tmp2=zrems(1,jk) 
         sumtmp2=sumtmp2+(tmp2+tmp2o)*(epi(jk)-epi(jk-1))*ergsev/2. 
         enddo 
!      terr=(sumtmp1-sumtmp2-elsum)/(sumtmp1+1.e-24)                    
      terr=(sumtmp1-sumtmp2)/(sumtmp1+1.e-24) 
      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10 
      alguu1=log10(max(1.d-24,uu1)) 
      skse=xlum/(xpx*r19*r19) 
      zeta=log10(max(1.d-24,skse)) 
      ecc=2.998e+10 
      ekt=t*(0.861707)*ergsev 
      sksec=skse/12.56/((1.+xee)*ekt*ecc) 
      zetac=log10(max(1.d-24,sksec)) 
      enn0=xpx 
      nlyc=nbinc(13.7d0,epi,ncn2) 
      nry=nlyc+1 
      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24) 
      nry=nbinc(13.6d0,epi,ncn2)+1 
!      write (lun11,9969)r,rdel,zeta,xee,xpx,t,hmctot,                  
!     $ dpthc(1,nry),dpthc(2,nry),ntotit,lnerrd                         
!9969  format (1x,9(1pe10.3),2i3)                                       
      tmp1=log10(r) 
      tmp2=log10(max(1.d-36,min(99.d0,rdel/r))) 
      tmp2c=log10(max(xcol,1.d-10)) 
      tmp3=log10(xpx) 
      tmp4=log10(t)+4. 
      tmp5=log10(max(dpthc(1,nry),1.d-10)) 
      tmp6=log10(max(dpthc(2,nry),1.d-10)) 
      tmp7=min(99.99d0,max(-99.99d0,hmctot*100.)) 
      tmp8=min(99.99d0,max(-99.99d0,terr*100.)) 
      write (tmpst,9889)tmp1,tmp2,tmp2c,zeta,xee,tmp3,tmp4,tmp7,        &
     & tmp8,tmp5,tmp6,ntotit                                            
      write (lun11,9889)tmp1,tmp2,tmp2c,zeta,xee,tmp3,tmp4,tmp7,        &
     & tmp8,tmp5,tmp6,ntotit                                            
 9889  format (1x,11(1x,f6.2),2i3) 
      call xwrite(tmpst,10) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   16 continue 
!                                                                       
!     times                                                             
      write (lun11,*)'times:',tread,tloop,tfunc,trates1,thcor,trates2,  &
     &          theat                                                  
      ttot=0. 
      do ll=1,ntyp 
        ttmpi=tucalc(ll)/max(1,ncall(ll))   !jg                        
        ttot=ttot+tucalc(ll)   !jg                                     
        write (lun11,9892)ll,ncall(ll),tucalc(ll),ttmpi   !jg          
 9892   format (1x,2i8,2(1pe11.3))                                     
        enddo 
      write (lun11,*)'total ucalc=',ttot 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   17 continue 
!                                                                       
!     column headings                                                   
      klabs(1)='log(r)' 
      klabs(2)='delr/r' 
      klabs(3)='log(N)' 
      klabs(4)='log(xi)' 
      klabs(5)=' x_e  ' 
      klabs(6)='log(n)' 
      klabs(7)='log(t)' 
      klabs(8)='h-c(%)' 
      klabs(9)='h-c(%)' 
      klabs(10)='log(tau)' 
!      klabs(10)='ntotit'                                               
      write (lun11,9979)(klabs(mm),mm=1,10) 
      write (tmpst,9979)(klabs(mm),mm=1,10) 
 9979 format (2x,3(1x,a6),1x,a7,a6,4(1x,a6),(1x,a9)) 
      call xwrite(tmpst,10) 
      klabs(1)='      ' 
      klabs(2)='      ' 
      klabs(3)='      ' 
      klabs(4)='      ' 
      klabs(5)='      ' 
      klabs(6)='      ' 
      klabs(7)='      ' 
      klabs(8)='      ' 
      klabs(9)='      ' 
      klabs(10)='fwd   ' 
      klabs(11)='rev   ' 
      write (lun11,9989)(klabs(mm),mm=1,11) 
      write (tmpst,9989)(klabs(mm),mm=1,11) 
      call xwrite(tmpst,10) 
 9989 format (3x,11a7) 
!                                                                       
      go to 9000
!                                                                       
!                                                                       
   22 continue 
!                                                                       
!     ionization parameter etc.                                         
      rdum=delr 
      delr=rdum 
      ergsev=1.602197e-12 
      r19=r*(1.e-19)
!      write (lun11,*)enlum,xpx,r,xlum                                  
      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10 
!      write (lun11,*)uu1                                               
      enlumx=0. 
      nb1=nbinc(100.d0,epi,ncn2) 
      nb10=nbinc(10000.d0,epi,ncn2) 
!      write (lun11,*)'nb1=',nb1,nb10                                   
      do kl=nb1,nb10 
!        write (lun11,*)kl,epi(kl),zremsz(kl),enlumx                    
        enlumx=enlumx+(zremsz(kl)/epi(kl)+zremsz(kl-1)/epi(kl-1))       &
     &                *(epi(kl)-epi(kl-1))/2.                           
        enddo 
      uux=enlumx/(12.56*xpx*r19*r19)/3.e+10 
      alguux=log10(max(1.d-24,uux)) 
      alguu1=log10(max(1.d-24,uu1)) 
      skse=xlum/(xpx*r19*r19) 
      zeta=log10(max(1.d-24,skse)) 
      ecc=2.998e+10 
      ekt=t*(0.861707)*ergsev 
!      sksec=skse/(12.56*((1.+xee)*ekt+pradl/(1.e-24+xpx))*ecc)         
      sksec=skse/12.56/((1.+xee)*ekt*ecc) 
      zetac=log10(max(1.d-24,sksec)) 
      enn0=xpx 
      nlyc=nbinc(13.7d0,epi,ncn2) 
      nry=nlyc+1 
      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24) 
      nry=nbinc(13.6d0,epi,ncn2)+1 
 9968 format (1x,' log(Xi)=',1pe11.3, ' log(u1)=',1pe11.3,              &
     & ' log(ux)=',1pe11.3,' gamma=',1pe11.3, ' rdel=',1pe11.3)         
 9965 format (1x,' r=',1pe11.3,' t=',1pe11.3,' log(xi)=',1pe11.3,       &
     & ' n_e=',1pe11.3,' n_p=',1pe11.3)                                 
 9966 format (1x,'httot=',1pe11.3,' cltot=',1pe11.3,                    &
     &      'taulc=',1pe11.3,'taulcb=',1pe11.3)                         
      write(lun11,9965)r,t,zeta,xnx,xpx 
      write(lun11,9966)httot,cltot,dpthc(1,nry),dpthc(2,nry) 
      write(lun11,9968)zetac,alguu1,alguux,egam,rdel 
      write (lun11,993) 
!                                                                       
      go to 9000
!                                                                       
   25 continue 
!                                                                       
!      write (lun11,*)'outputting to the common block',nlsvn            
!      do mm=1,ncn2                                                     
!        epi2(mm)=epi(mm)                                               
!        do ll=1,3                                                      
!          zrems2(ll,mm)=zrems(ll,mm)                                   
!          enddo                                                        
!        enddo                                                          
!      lpril=1                                                          
!      nilino=0                                                         
!      jkktmp=0                                                         
!      do j=1,nlsvn                                                     
!          kk=j                                                         
!          ln=nplin(j)                                                  
!          ml=ln                                                        
!          if (ml.ne.0) then                                            
!!            write (lun11,*)'   ',j,ml                                 
!            mlm=ml-1                                                   
!            call drd(ltyp,lrtyp,lcon,                                  
!     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                       
!     $                     0,lun11)                         
!            elin=rdat1(np1r)                                           
!            llo=idat1(np1i)                                            
!            lup=idat1(np1i+1)                                          
!            elin=rdat1(np1r)                                           
!            aij=rdat1(np1r+2)                                          
!            nilin=npar(ml)                                             
!            if ((nilin.gt.0).and.(nilin.lt.ndat2)) then                
!                if (nilin.ne.nilino) jkktmp=jkktmp+1                   
!                nilino=nilin                                           
!                mlm=nilin-1                                            
!                call drd(ltyp,lrtyp,lcon,                              
!     $            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                   
!     $                         0,lun11)                     
!                do mm=1,nkdt                                           
!                  kdtmp(mm)=kdat1(np1k-1+mm)                           
!                  enddo                                                
!                do mm=nkdt+1,9                                         
!                  kdtmp(mm)=kblnk                                      
!                  enddo                                                
!                nilin=idat1(np1i+2)                                    
!!                write (lun11,*)ml,nilin,npar(ml)                      
!                newout=newout+1                                        
!                newout=min(newout,nnnl)                                
!                lnewo(newout)=j                                        
!                ml=npfi(13,jkktmp)                                     
!                mllz=npar(ml)                                          
!                lupfnd=0                                               
!                llofnd=0                                               
!                mlpar=npar(ml)                                         
!                do while ((ml.ne.0).and.(mlpar.eq.mllz)                
!     $            .and.((llofnd.ne.1).or.(lupfnd.ne.1)))               
!!                    write (lun11,*)ml,masterdata(2,ml),mltype,jkk          
!                  mlm=ml-1                                             
!                  call drd(ltyp,lrtyp,lcon,                            
!     $              nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                 
!     $                           0,lun11)                   
!                  nlevmx=nlevmx+1                                      
!                  nlev=idat1(np1i+nidt-2)                              
!!                  write (lun11,*)ml,nlev,llo,lup,(kdat1(np1k-1+mm),mm=
!                  if (nlev.eq.llo) then                                
!                    do mm=1,20                                         
!                      if (mm.le.nkdt) then                             
!                          klevl(mm)=kdat1(np1k-1+mm)                   
!                        else                                           
!                          klevl(mm)=kblnk                              
!                        endif                                          
!                      enddo                                            
!!                   write (lun11,*)kk,ktmp2                            
!                    llofnd=1                                           
!                    gglo=rdat1(np1r+1)                                 
!                    endif                                              
!                  if (nlev.eq.lup) then                                
!                    do mm=1,20                                         
!                      if (mm.le.nkdt) then                             
!                          klevu(mm)=kdat1(np1k-1+mm)                   
!                        else                                           
!                          klevu(mm)=kblnk                              
!                        endif                                          
!                      enddo                                            
!                    lupfnd=1                                           
!                    ggup=rdat1(np1r+1)                                 
!                    endif                                              
!                  ml=npnxt(ml)                                         
!                  if (ml.ne.0) mlpar=npar(ml)                          
!                  enddo                                                
!                if ((llofnd.eq.1).and.(lupfnd.eq.1)) then              
!                  flinewo(newout)=(1.e-16)*aij*ggup*elin*elin          
!     $                             /((0.667274)*gglo)                  
!                  aijewo(newout)=aij                                   
!                  ggloewo(newout)=gglo                                 
!                  ggupewo(newout)=ggup                                 
!                  do mm=1,8                                            
!                    kdewo(mm,newout)=kdtmp(mm)                         
!                    enddo                                              
!                  do mm=1,20                                           
!                    kdewol(mm,newout)=klevl(mm)                        
!                    enddo                                              
!                  do mm=1,20                                           
!                    kdewou(mm,newout)=klevu(mm)                        
!                    enddo                                              
!                  elewo(newout)=elin                                   
!                  tau0ewo(newout)=tau0(1,j)                            
!                  elout(1,newout)=elum(1,j)                            
!                  elout(2,newout)=elum(2,j)                            
!!                  write (lun11,*)kk,ln,j,(kdtmp(mm),mm=1,8),elin,     
!!     $             tau0(1,j),elum(1,j),elum(2,j),newout               
!c9955             format (1x,2i8,1x,8a1,3(1pe11.3))                    
!                  endif                                                
!              endif                                                    
!            endif                                                      
!          enddo                                                        
!      call commonprint(lun11)                                          
!                                                                       
      go to 9000
!                                                                       
   26 continue 
!                                                                       
      go to 9000
!                                                                       
!     ferland print                                                     
      lpril=0 
!     print 500 strongest emission lines                                
      write (lun11,*)'log(emission line fluxes (erg/sec/cm^2))' 
      kltmpo=0 
      nlplmx=500 
      eliml=0.1 
      elimh=1.0e10 
!     find the strongest lines.                                         
      do  lm=1,nlplmx 
        kltmp(lm)=0 
        enddo 
      nlpl=1 
      do lnn=1,nlsvn 
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        mlm=ml 
        call drd(ltyp,lrtyp,lcon,                                       &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                            &
     &    0,lun11)                                                
        elin=abs(masterdata%rdat1(np1r)) 
        if ((lrtyp.ne.14).and.(lrtyp.ne.9)) then 
          nilin=derivedpointers%npar(ml) 
          mlm=nilin
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          nilin2=masterdata%idat1(np1i+nidt-1) 
          elmmtpp=(elum(2,ln)+elum(1,ln))/2. 
          if (lpril.ne.0)                                               &
     &       write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml                
          if ((ln.gt.0).and.(ln.lt.nnnl)                                &
     &       .and.(elin.ge.eliml).and.(elin.le.elimh)                   &
     &       .and.(elin.le.8.9e+6)                                      &
     &       .and.(elmmtpp.gt.1.d-36)                                   &
     &       .and.(nilin2.gt.0).and.(nilin2.le.nni))                    &
     &        then                                                      
            lmm=0 
            elcomp=1.e+10 
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp)) 
              lmm=lmm+1 
              kl2=kltmp(lmm) 
              elcomp=0. 
              if (kl2.gt.0)                                             &
     &          elcomp=(elum(2,kl2)+elum(1,kl2))/2.                     
              enddo 
            if (lpril.ne.0)                                             &
     &       write (lun11,8516)ln,elin,elmmtpp,lmm,nlpl,kl2,elcomp      
! 8516       format (1h ,i4,2e12.4,3i4,e12.4)                           
            kltmpo=ln 
            do  k=lmm,min(nlplmx,nlpl) 
              if ((lpril.ne.0).and.(kltmp(k).ne.0))                     &
     &          write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo          
              kltmpn=kltmp(k) 
              kltmp(k)=kltmpo 
              kltmpo=kltmpn 
              enddo 
           nlpl=min(nlplmx,nlpl+1) 
           if (lpril.ne.0)                                              &
     &       write (lun11,*)'done with 557 loop',lm                     
            endif 
          endif 
        enddo 
       if (nlpl.gt.0) kltmp(nlpl)=kltmpo 
!      nlpl=nlpl-1                                                      
      write (lun11,9599) 
 9599 format (1x,'index, ion, wavelength, reflected, transmitted,total') 
      r19=r*1.e-19
      do  kk=1,nlpl 
        if (lpril.ne.0)                                                 &
     &    write (lun11,*)'kk=',kk                                       
        ln=kltmp(kk) 
        if (ln.ne.0) then 
          ml=derivedpointers%nplin(ln) 
          if (ml.ne.0) then 
            if (lpril.ne.0)                                             &
     &      write (lun11,*)'   ',ln,ml                                  
            mlm=ml 
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            elin=abs(masterdata%rdat1(np1r)) 
            nilin=derivedpointers%npar(ml) 
            mlm=nilin
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            do mm=1,nkdt 
              kdtmp(mm)=masterdata%kdat1(np1k-1+mm) 
              enddo 
            do mm=nkdt+1,9 
              kdtmp(mm)=kblnk 
              enddo 
             flux1=elum(1,ln)/12.56/r19/r19 
             flux2=elum(2,ln)/12.56/r19/r19 
!            nilin=idat1(np1i+2)                                        
            if (lpril.ne.0)                                             &
     &      write (lun11,*)ml,nilin,derivedpointers%npar(ml)                   
            write (lun11,9956)kk,ln,(kdtmp(mm),mm=1,9),elin,            &
     &      log10(flux1),log10(flux2),log10(flux1+flux2)                
 9956       format (1x,2i8,1x,9a1,4(1pe13.5)) 
            endif 
          endif 
        enddo 
      write (lun11,993) 
      go to 9000
!                                                                       
                                                                        
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!        
!                                                                       
!     Write  ion column densities                                       
!     requires that zrtmp be filled by calling pprint(12)               
!                                                                       
   27 continue 
!                                                                       
      write (lun11,*)'ion column densities' 
      write (lun11,9447) 
 9447  format (1x,'index, ion, column density') 
!                                                                       
      do lk=1,nni 
        xcoltmp(lk)=0. 
        enddo 
!                                                                       
!     step thru ions                                                    
      klion=12 
      mlion=derivedpointers%npfirst(klion) 
      lk=0 
      do while (mlion.ne.0) 
!                                                                       
!        get ion data                                                   
         lk=lk+1 
         ltyp=klion 
         mlm=mlion
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         do mm=1,nkdt 
           kdtmp(mm)=masterdata%kdat1(np1k-1+mm) 
           enddo 
         do mm=nkdt+1,9 
           kdtmp(mm)=kblnk 
           enddo 
!                                                                       
!        get element data                                               
         nell=derivedpointers%npar(mlion) 
         mlm=nell
         call drd(ltyp,lrtyp,lcon,                                      &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
!         write (lun11,*)mlion,lk,np1i,nidt,np1i+nidt-1,                
!     $      idat1(np1i+nidt-1),ababs(idat1(np1i+nidt-1)),mlm           
         if ((masterdata%idat1(np1i+nidt-1).gt.0)                       &
     &     .and.(masterdata%idat1(np1i+nidt-1).le.nl)) then            
           abundel=ababs(masterdata%idat1(np1i+nidt-1)) 
           xeltp=abundel 
!                                                                       
           do jkl=2,numrec 
!             write (lun11,*)jkl,lk,zrtmp(2,jkl),zrtmp(8+lk,jkl),       
!     $                             zrtmp(5,jkl),xeltp,xcoltmp(lk)      
             xcoltmp(lk)=xcoltmp(lk)                                    &
     &         +(zrtmp(8+lk,jkl)*zrtmp(5,jkl)                           &
     &             +zrtmp(8+lk,jkl-1)*zrtmp(5,jkl-1))                   &
     &         *(zrtmp(2,jkl)-zrtmp(2,jkl-1))*xeltp/2.                  
             enddo 
!                                                                       
!          print out                                                    
           if (xcoltmp(lk).gt.1.e-15)                                   &
     &      write (lun11,9446)lk,(kdtmp(mm),mm=1,9),                    &
     &      xcoltmp(lk)                                                 
 9446      format (1x,i4,1x,9a1,1pe16.8) 
!                                                                       
           endif 
!                                                                       
         mlion=derivedpointers%npnxt(mlion) 
         enddo 
!       
9000  continue
      deallocate(kdat)
      deallocate(elsv)
      deallocate(jpnt)
      deallocate(klabs)
      deallocate(kunits)
      deallocate(kform)
!                                                                       
!                                                                       
      return 
      END                                           
