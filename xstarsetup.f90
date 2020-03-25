      subroutine xstarsetup(lnerrd,nlimd,                               &
     &       lpri,lprid,lunlog,tinf,critf,                              &
     &       t,tp,r,delr,xee,xpx,ababs,abel,cfrac,xlum,p,lcdd,          &
     &       epi,ncn2,bremsa,bremsint,atcredate,                        &
     &       zrems,zremso,zremsz,                                       &
     &       tau0,dpthc,dpthcont,tauc,                                  &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,hmctot,elcter, &
     &       httot2,cltot2,                                             &
     &       xilev,bilev,rnist,elum,                                    &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,nlin,elin)                                    
!                                                                       
!     Name: xstarsetup.f90  
!     Description:  
!           performs initialization needed for xstar run. 
!           initializes all physical variables, reads in atomic data,
!           compton heating tables,   sets up pointers.
!
!     List of Parameters:
!           Input: 
!           lnerrd: thermal equilibrium error switch
!           nlimd: thermal equilibrium solver iteration limit
!           lpri: print switch
!           lprid: thermal equilibrium solver print switch
!           lunlog: logical unit number for printing
!           tinf:  temperature lower limit
!           critf: threshold value for ion fraction to be included in 
!                   level population calculation
!           t: temperature in 10^4K
!           tp: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!           r:  radius in nebula (cm)
!           delr: thickness of current spatial zone (cm)
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           ababs(nl):  element abundances relative to H=1
!           abel(nl):  element abundances relative to xstar 
!                standard values
!           cfrac:  covering fraction (affects line and continuum 
!                forward-backward ratio
!           xlum: source luminosity integrated from 1-1000 Ry
!               in units of 10^38 erg/s
!           p:  pressure in dynes/cm^2
!           lcdd: constant pressure switch, 1=constant pressure 
!                      0=constant density
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!       Output:
!           bremsa(ncn):  Ionizing flux (erg/s/cm^2/erg)
!           bremsint(ncn):  Integral of bremsa from each bin to epi(ncn2)
!               (erg/s/cm^2)
!           atcredate:  atomic data file creation date (string length 63)
!           zrems(4,ncn):  master spectrum array.  (erg/s/erg/10^38)
!           zremso(4,ncn):  master spectrum array previous step.  
!                     (erg/s/erg/10^38)
!           zremsz(ncn):  spectrum array, incident spectrum
!                     (erg/s/erg/10^38)
!           tau0(2,nnnl):  line optical depths
!           dpthc(2,ncn):  continuum optical depths in continuum bins
!           dpthcont(2,ncn):  continuum optical depths in continuum bins 
!                          without lines
!           tauc(2,nnml):  rrc optical depths
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           ntotit:  number of iterations for thermal equilibrium
!           xii(nni):  ion fractions, xiin(1)=H, xiin(2)=He0, xiin(3)=He+ etc
!           rrrt(nni): total recombination rates for each ion (s^-1)
!           pirt(nni): total photoionization rates for each ion(s^-1)
!           htt(nni): total heating rate for each ion (approximate) 
!                       (erg s^-1 cm^-3)
!           cll(nni): total cooling rate for each ion (approximate) 
!           httot: total heating rate (erg s^-1 cm^-3) 
!           cltot: total cooling rate (erg s^-1 cm^-3) 
!           cllines:  total cooling rate due to lines (erg s^-1 cm^-3) 
!           hmctot: 2*(heating-cooling)/(heating+cooling)
!           elcter:  charge conservation error
!           xilev(nnml):  level populations (relative to parent element)
!           bilev(nnml):  level departure coefficient
!           rnist(nnml):  lte level populations (relative to parent element)
!           elum(nnnl):  line luminosities (erg/s/10^38)
!           rcem(2,nnnl):  line emissivities  (erg cm^-3 s^-1) /10^38
!                  inward and outward
!           oplin(nnnl):  line opacities  (cm^-1)
!           rccemis(2,ncn): continuum emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!                  inward and outward
!           brcems(ncn):  bremsstrahlung emissivity (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           opakcont(ncn):  continuum opacities lines excluded (cm^-1)
!           cemab(nnml):  rrc emissivities (erg cm^-3 s^-1) 
!           cabab(nnml):  total energy absorbed by rrc (erg cm^-3 s^-1) 
!           opakab(nnml):  rrc opacities (cm^-1)
!           nlin:  number of lines in special line array
!           elin:  special line array, line wavelengths
!     Dependencies: getenv (system), xwrite, readtbl, getlunx, setptrs, 
!          init, trnfrc, drd
!     Called by:  xstar
!
!      this routine does many of the setup chores: read in atomic       
!        data, set up pointers, zeroing variables.                      
!      NB: no input parameters are affected                             
!                                                                       
      use globaldata
      use times
      implicit none 
!                                                                       
!                                                                       
!     global xstar data

!     line luminosities                                                 
      real(8) elum(3,nnnl),elumo(3,nnnl) 
!     line emissivities                                                 
      real(8) rcem(2,nnnl) 
!     line opacities                                                    
      real(8) oplin(nnnl) 
      real(8) fline(2,nnnl),flinel(ncn) 
!     line optical depths                                               
      real(8) tau0(2,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!      continuum lum                                                    
      real(8) zrems(5,ncn),zremsz(ncn)                 
      real(8) zremso(5,ncn)
!     continuum optical depths                                          
      real(8) dpthc(2,ncn),dpthcont(2,ncn)
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     level populations                                                 
      real(8) xilev(nnml),bilev(nnml),rnist(nnml)
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) elumab(2,nnml),elumabo(2,nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating and cooling                                               
      real(8) htt(nni),cll(nni) 
      real(8) htt2(nni),cll2(nni) 
      real(8) rrrt(nni),pirt(nni) 
!     element abundances                                                
      real(8) abel(nl),abcosmic(30),ababs(nl) 
      integer nlin(nnnl) 
      real(8) elin(nnnl) 
!     the atomic data creation date                                     
      character(63) atcredate 
!                                                                       
!     local variables                                                   
!     state variables                                                   
      real(8) p,r,t,xpx,delr,tp 
!     heating-cooling variables                                         
      real(8) httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,          &
     &     clcont,hmctot,httot2,cltot2
!     input parameters                                                  
      real(8) xlum,xpxcol 
      real(8) cfrac,critf,xee 
      integer lcdd,ncn2 
!     variables associated with thermal equilibrium solution            
      integer ntotit,lnerrd 
!     switches                                                          
      integer lprid,lpri 
      integer  nlimd,lunlog,lun11,lun25 
!     temporary for xwrite                                              
      character(133) tmpst 
!     strings for atomic data read                                      
      character(256) datafil3,datafil4,datafile 
      logical ex3,ex4 
!     temporary integers                                                
      integer ll,mm,ldir,jk,mlm 
      integer nlsvn,ncsvn, lenact 
!     times                                                             
      real(8) tinf,t1s 
      integer np1r,np1i,np1k,np2 
      integer jlk,j,ml,ltyp,lrtyp,lcon,kl,                              &
     &        nrdt,nidt,nkdt                                            
!     storing info for parameters                                       
      character(20) parname(55) 
      real(8) abeltmp(nl)
!                                                                       
!     Not used                                                          
      integer javi 
      real(8) javir 
!                                                                       
!     these are the anders and grevesse abundances from the xspec    
!       manual                                                          
!      data abcosmic/'1.00d+00 ','9.77d-02 ','1.45d-11 ','1.41d-11 ',   
!     $ '3.98d-10 ','3.63d-04 ','1.12d-04 ','8.51d-04 ','3.63d-08 ',    
!     $ '1.23d-04 ','2.14d-06 ','3.80d-05 ','2.95d-06 ','3.55d-05 ',    
!     $ '2.82d-07 ','1.62d-05 ','1.88d-07 ','3.63d-06 ','1.32d-07 ',    
!     $ '2.29d-06 ','1.26d-09 ','9.77d-08 ','1.00d-08 ','4.84d-07 ',    
!     $ '2.45d-07 ','4.68d-05 ','8.60d-08 ','1.78d-06 ','1.62d-08 ',    
!     $ '3.98d-08 '/                                                    
!     old values                                                        
      data abcosmic/1.00E+000, 1.00D-001,                               &
     & 1.d-10, 1.d-10, 1.d-10, 3.70D-004, 1.10D-004,                    &
     & 6.80D-004, 3.98D-008, 2.80D-005, 1.78D-006, 3.50D-005,           &
     & 2.45D-006, 3.50D-005, 3.31D-007, 1.60D-005, 3.98D-007,           &
     & 4.50D-006, 8.91D-008, 2.10D-006, 1.66D-009, 1.35D-007,           &
     & 2.51D-008, 7.08D-007, 2.51D-007, 2.50d-005, 1.26D-007,           &
     & 2.00D-006, 3.16D-008, 1.58D-008/                                 
!                                                                       
!                                                                       
!     Parameter Names                                                   
!                                                                       
      data parname/'cfrac','temperature',                               &
     &   'lcpres','pressure','density','spectrum',                      &
     &   'spectrum_file','spectun','trad',                              &
     &   'rlrad38','column','rlogxi',                                   &
     &   'nsteps','niter','lwrite',                                     &
     &   'lprint','lstep',                                              &
     &   'habund','heabund',                                            &
     &   'liabund','beabund','babund','cabund',                         &
     &   'nabund','oabund','fabund','neabund',                          &
     &   'naabund','mgabund','alabund',                                 &
     &   'siabund','pabund','sabund',                                   &
     &   'clabund','arabund','kabund',                                  &
     &   'caabund','scabund','tiabund',                                 &
     &   'vabund','crabund','mnabund ',                                 &
     &   'feabund','coabund','niabund',                                 &
     &   'cuabund','znabund','emult','taumax','xeemin',                 &
     &   'critf','vturbi','npass','modelname',                          &
     &   'loopcontrol'/                                                 
                                                                        
!                                                                       
      javir=t 
      javir=tp 
      javi=nlimd 
      javi=lcdd 
      lcdd=javi 
      javir=delr 
      javir=cfrac 
      javir=xlum 
      javir=p 
      p=javir 
      javir=critf 
!                                                                       
      call remtms(t1s) 
!                                                                       
!     opening message                                                   
<<<<<<< HEAD
      write (lunlog,*)'xstar version 2.56c' 
=======
      write (lunlog,*)'xstar version 2.56a' 
>>>>>>> 2d75308c63b9789458ce092c697c7853fcdde44a
!                                                                       
!     Test if atomic database files are available.  Abort if not.       
      call getenv('LHEA_DATA', datafile) 
      datafil4 = datafile(1:lenact(datafile))//'/atdb.fits' 
      datafil3 = datafile(1:lenact(datafile))//'/coheat.dat' 
      inquire(file=datafil3,exist=ex3) 
      inquire(file=datafil4,exist=ex4) 
      if (.not.(ex3 .and. ex4 )) then 
         write(tmpst,*)'xstar: One or more of the Atomic Database files' 
         write(lunlog,*)tmpst 
         call xwrite(tmpst,10) 
         write(tmpst,*)'xstar: are missing.' 
         write(lunlog,*)tmpst 
         call xwrite(tmpst,10) 
         write(tmpst,*)'xstar: ',datafil4(1:lenact(datafil4)) 
         write(lunlog,*)tmpst 
         call xwrite(tmpst,10) 
         write(tmpst,*)'xstar: ',datafil3(1:lenact(datafil3)) 
         write(lunlog,*)tmpst 
         call xwrite(tmpst,10) 
         write(tmpst,*)'Program aborting...' 
         write(lunlog,*)tmpst 
         call xwrite(tmpst,10) 
         close(lunlog) 
         return 
      endif 
!                                                                       
!                                                                       
!                                                                       
      tread=0.                                                         
      trates1=0.                                                       
      thcor=0.                                                         
      trates2=0.                                                       
      theat=0.                                                         
      do kl=1,ntyp                                                     
         tucalc(kl)=0.                                                 
         ncall(kl)=0                                                   
         enddo                                                         
!                                                                       
!                                                                       
!     read in                                                           
      write (lunlog,*)'Loading Atomic Database...' 
      write (tmpst,*)'Loading Atomic Database...' 
      call xwrite(tmpst,10) 
                                                                        
      call readtbl(np1r,np1i,np1k,                                   &
     &       np2,                                                       &
     &      datafil4,atcredate,lpri,lunlog)                             
!                                                                       
!                                                                       
      call getlunx(lun25) 
      open(unit=lun25,file=datafil3,status='unknown') 
      rewind(lun25) 
      read (lun25,901) 
  901 format (1x) 
      do mm=1,ncomp 
        do ll=1,ncomp 
          read (lun25,902)sxcomp(mm),                                   &
     &            ecomp(ll),decomp(mm,ll)                               
  902     format (9x,e12.4,12x,2(e12.4)) 
          enddo 
        enddo 
      close(lun25) 
!                                                                       
!     Initialize the database                                           
      write (lunlog,*)'initializng database...' 
      write (tmpst,*)'initializng database...' 
      call xwrite(tmpst,10) 
!     need to fool setptrs into using all the data
      do mm=1,nl
        abeltmp(mm)=1.
        enddo
      call setptrs(lunlog,lpri,                                         &
     &       np2,ncsvn,nlsvn,                                           &
     &       abcosmic,abeltmp)                         
!
!                                                                       
!     read in parameter values                                          
!      write(lunlog,*)'Atomic Abundances'                               
!      write(lunlog,*)'Element      Solar    Hydrogen'                  
      do ll=1,nl 
        ababs(ll)=abel(ll)*abcosmic(ll) 
        write(lunlog,9990)parname(17+ll),abel(ll),ababs(ll) 
        enddo 
      write(lunlog,*)' ' 
 9990 format(A10,2(E12.4)) 
!                                                                       
!                                                                       
!     set up and initialize                                             
      tinf=0.31 
      call init(lunlog,abel,bremsa,bremsint,tau0,dpthc,dpthcont,tauc,   &
     &   xii,rrrt,pirt,htt,cll,htt2,cll2,httot,cltot,httot2,cltot2,     &
     &   cllines,clcont,htcomp,clcomp,clbrems,                          &
     &   xilev,rcem,oplin,rccemis,brcems,opakc,opakcont,                &
     &   cemab,cabab,opakab,elumab,elumabo,elum,elumo,                  &
     &   zrems,zremso,fline,flinel)                                     
!                                                                       
      do jk=1,ncn 
        zrems(1,jk)=0. 
        zrems(2,jk)=0. 
        dpthc(1,jk)=0. 
        dpthc(2,jk)=0. 
        dpthcont(1,jk)=0.
        dpthcont(2,jk)=0.
        enddo 
      do jk=1,nnnl 
        elum(1,jk)=0. 
        elum(2,jk)=0. 
        tau0(1,jk)=0. 
        tau0(2,jk)=0. 
        enddo 
      do jk=1,nnml 
        bilev(jk)=0. 
        rnist(jk)=0. 
        tauc(1,jk)=0. 
        tauc(2,jk)=0. 
        enddo 
!                                                                       
      lnerrd=0 
      ntotit=0 
      lprid=0 
      xee=1.21 
      elcter=0. 
      hmctot=0. 
!                                                                       
          ldir=1 
          lun11=lunlog 
          call trnfrc(lpri,lun11,ldir,                                  &
     &      r,xpxcol,xpx,                                               &
     &      epi,ncn2,zremsz,dpthc,opakc,                                &
     &      zrems,bremsa,bremsint)                               
!                                                                       
!                                                                       
      do jlk=1,nlsvn 
         j=jlk 
         ml=derivedpointers%nplin(j) 
         mlm=ml 
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         elin(jlk)=0. 
         nlin(jlk)=0 
         if ((lrtyp.ne.14).and.(lrtyp.ne.9)) then 
           elin(jlk)=masterdata%rdat1(np1r) 
           nlin(jlk)=masterdata%idat1(np1i+nidt-1) 
           endif 
         enddo 
!                                                                       
      return 
      END                                           
