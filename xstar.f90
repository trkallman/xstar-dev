      module globaldata
      implicit none 
!                                                                       
      include './PARAM' 
!                                                                       
!     global xstar data
!     master data
      TYPE :: master_data
        integer :: idat1(nidat1) ! integer data
        real(8) :: rdat1(nrdat1)  ! real data
        integer :: nptrs(nptt,ndat2) ! pointer data
        character(1) :: kdat1(nkdat1) ! character data
      END TYPE master_data
      TYPE(master_data) :: masterdata
      TYPE :: derived_pointers
        integer :: npar(ndat2)   !    pointers to master data
        integer :: npnxt(ndat2)  !    pointers to master data
        integer :: npfirst(ntyp) !    pointers to master data
        integer :: npfi(ntyp,nni)!    pointers to master data
        integer :: nplin(nnnl)   ! pointers to line data
        integer :: nplini(ndat2) ! pointers to line data
        integer :: npcon(nnml)
        integer :: npconi2(ndat2) 
        integer :: npconi(ndat2)
        integer :: npilev(nd,nni)
        integer :: npilevi(nnml)
        integer :: nlevs(nni) 
      END TYPE derived_pointers
      TYPE(derived_pointers) :: derivedpointers
      TYPE :: level_temp
        real(8) :: rlev(10,nd) 
        integer:: ilev(10,nd),nlpt(nd),iltp(nd) 
        character(1) :: klev(100,nd) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
!     compton heating data                                              
      real(8) decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp) 
      end module globaldata

      program x116n5
      use globaldata
      implicit none

!      this version resurrects the continuum escape probabilities
!      also allows ncn=1e6
!      also allows density > 1.e18 (with reservations)

!      common /lpass/lunu

!      lunu=20

      open (unit=5,file='xin116no.dat',status='unknown')
      call xstar

      stop
      end program
      subroutine xstar 
!                                                                       
!      based on attenuate
!
      use globaldata
      implicit none 
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
      real(8) zrems(5,ncn),zremso(5,ncn),zremsz(ncn)
!     continuum optical depths                                          
      real(8) dpthc(2,ncn),dpthcont(2,ncn)
!     continuum flux                                                    
      real(8) bremsa(ncn),bremsint(ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn),brcems(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn),opakcont(ncn)
!     level populations                                                 
      real(8) xilev(nnml),bilev(nnml),rniss(nnml) 
      real(8) cemab(2,nnml),cabab(nnml),opakab(nnml) 
      real(8) elumab(2,nnml),elumabo(2,nnml) 
      real(8) tauc(2,nnml) 
!     ion abundances                                                    
      real(8) xii(nni) 
!     heating and cooling                                               
      real(8) htt(nni),cll(nni) 
      real(8) rrrt(nni),pirt(nni) 
!     element abundances                                                
      real(8) abel(nl),ababs(nl) 
!     the atomic data creation date                                     
      character(63) atcredate 
!                                                                       
!     local variables                                                   
!     state variables                                                   
      real(8) p,rdel,r19,xi,xcol,zeta,rdelo,r,t,xpx,delr,xpx0,r0 
      real(8) rnew,dennew 
!     heating-cooling variables                                         
      real(8) httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,          &
     &     clcont,hmctot
!     input parameters                                                  
      character(16) knam,knam2,knam3,knam4 
      character(80) kmodelname,specfile,spectype 
      real(8) enlum,emult,taumax,xeemin,xlum,rmax,xpxcol,trad 
      real(8) cfrac,critf,vturbi,xlum2,xee,radexp 
      integer lcdd,ncn2 
!     variables associated with thermal equilibrium solution            
      integer ntotit,lnerrd 
!     switches                                                          
      integer lprid,lpril,lpri,lwri,lpriu,nlimdt,lpris 
      integer lprisv,lpri2,lpri3
      integer  nnmax,nlimd,lunlog,nloopctl,numrec,npass,lfix 
!     temporary for xwrite                                              
      character(133) tmpst 
      real(8) ectt 
!     temporary for spectrum                                            
      real(8) eptmp(ncn),zrtmp(ncn) 
      integer nlprnt(19),nlnprnt 
!     temporary integers                                                
      integer ll,mm,kk,ldir,jk,jkp
      integer nlsvn,ncsvn 
      real(8) eliml,elimh 
      integer istatus, iunit,iunit2,numcon2,iunit3,iunit4 
      integer iunito,iunit2o,iunit3o,iunit4o,ierr,lmax
!     times                                                             
      real(8) tinf,ttot,t1s,t2s 
      integer ntmp,np2,lun40 
      integer lpril2 
!     storing info for parameters                                       
      character(20) parname(55) 
      character(10) partype(55) 
      real(8) parms(55) 
      character(30) parcomm(55) 
      integer nparms, specunit
!     temporary line pointers                                           
      integer nlin(nnnl) 
      real(8) elin(nnnl)
!
      integer status
!                                                                       
!     local definitions                                                 
!                                                                       
!                                                                       
!     Warning!!  #11 & #7 must be run before #5 since #5 changes        
!                the variable zrtmp                                     
!                Also make sure that if you change the number of        
!                entries in nlprnt, that you also update nlnprnt        
!                and the real(4) statement for nlprnt                    
      data nlnprnt/10/,nlprnt/22,11,1,23,24,5,10,16,15,19,14,4,6,       &
     &            21,7,18,27,8,26/                                      
!     $            21,7,26,2*0/                                         
!                                                                       
!     Parameter Names                                                   
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
      data partype/'real','real',                                       &
     &    'integer','real','real','string',                             &
     &    'string','integer','real',                                    &
     &    'real','real','real',                                         &
     &    'integer','integer','integer',                                &
     &    'integer','integer',                                          &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real',                                                &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real','real',                                         &
     &    'real','real','real','integer','string',                      &
     &    'integer'/                                                    
      data parcomm/' ','Units of 10**4 K',                              &
     &     '1=yes, 0=no','dynes/cm**2','cm**(-3)',' ',                  &
     &     ' ','0=energy, 1=photons','or alpha',                        &
     &     '/10**38 erg/sec','cm**(-2)',' ',                            &
     &     ' ',' ','1=yes, 0=no',                                       &
     &     '1=yes, 0=no',' ',                                           &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',                                                     &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',' ',                                                 &
     &     ' ',' ',                                                     &
     &     ' '/                                                         
!
!                                                                       
      nparms=55 
!                                                                       
!     Allocate temporary and logging files                              
      call getlunx(lunlog) 
      open(unit=lunlog,file='xout_step.log',status='unknown') 
!                                                                       
      call remtms(t1s) 
!                                                                       
!     opening message                                                   
      write (tmpst,*)'xstar version 2.54' 
      call xwrite(tmpst,10) 
!                                                                       
!     default parameter values                                          
!                                                                       
      trad=1. 
      xlum=1. 
      lpri=0 
      lwri=0 
      r=0. 
      r0=r 
      t=1. 
      xpx=1000. 
      xpx0=xpx 
      p=0.03 
      lcdd=1 
      numrec=2 
      npass=0 
      nnmax=1 
      nlimd=0 
      rmax=1. 
      xpxcol=1.e+16 
      zeta=0. 
      xi=10.**zeta 
      lfix=0 
      ncn2=min(ncn,999) 
      call ener(epi,ncn2) 
      do mm=1,ncn2 
        zremsz(mm)=0. 
        enddo 
      cfrac=0. 
      emult=0.75 
      taumax=5. 
      xeemin=0.1 
      nloopctl=0 
      critf=1.d-8 
      vturbi=1. 
      xcol=0. 
      radexp=0. 
!                                                                       
!     read in parameter values                                          
      call rread1(trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,numrec,npass,      &
     & nlimd,rmax,xpxcol,xi,zeta,lfix,                                  &
     & lunlog,abel,cfrac,emult,taumax,xeemin,spectype,specfile,         &
     & specunit,kmodelname,nloopctl,critf,vturbi,eptmp,zrtmp,numcon2,   &
     &     ncn2,radexp)                                                 
!            
      if (radexp.lt.-99.) then 
        call getlunx(lun40) 
        open(lun40,iostat=ierr,file='density.dat') 
        if (ierr.ne.0) stop 'missing density file' 
        read (lun40,*,iostat=ierr)rnew,dennew 
        r=rnew 
        xpx=dennew 
        endif 
!                                                                       
      lprisv=lpri 
      lpri=0 
      call ener(epi,ncn2) 
      do mm=1,ncn2 
        zremsz(mm)=0. 
        enddo 
!                                                                       
      write (lunlog,*)'main loop' 
      write (tmpst,*)'main loop' 
      call xwrite(tmpst,10) 
!                                                                       
      call xstarsetup(lnerrd,nlimd,                                     &
     &       lpri,lprid,lunlog,tinf,critf,                              &
     &       t,trad,r,delr,xee,xpx,ababs,abel,cfrac,xlum,p,lcdd,        &
     &       epi,ncn2,bremsa,bremsint,atcredate,                        &
     &       zrems,zremso,zremsz,                                       &
     &       tau0,dpthc,dpthcont,tauc,                                  &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &       xilev,bilev,rniss,elum,                                    &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,nlin,elin)                                    
!                                                                       
      xee=1.21 
      xi=10.**zeta 
!                                                                       
!     Write the parameter list to the FITS file                         
!     When changing this list, make sure nparms, parname, partype,      
!     and parcomm are also properly updated.  Watch out for the         
!     model name which is currently parcomm(37)                         
      call pprint(3,1,trad,xlum,lwri,lprisv,r,t,xpx,p,lcdd,             &
     &        numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,        &
     &        zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,           &
     &        spectype,specfile,specunit,kmodelname,nloopctl,           &
     &        nparms,parname,partype,parms,parcomm,atcredate,           &
     &        lunlog,tinf,xcol,vturbi,critf,radexp,                     &
     &        delr,rdel,enlum,xee,ababs,                                &
     &        bremsa,bremsint,tau0,dpthc,tauc,                          &
     &        ncsvn,nlsvn,                                              &
     &        ntotit,lnerrd,                                            &
     &        xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                 &
     &        cllines,clcont,htcomp,clcomp,clbrems,                     &
     &        xilev,bilev,rniss,                                        &
     &        rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,    &
     &        cabab,elumab,elum,zrems)                                  
!                                                                       
!     set up and initialize                                             
      rdel = 0. 
      call pprint(2,1,trad,xlum,lwri,lprisv,r,t,xpx,p,lcdd,             &
     &        numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,        &
     &        zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,           &
     &        spectype,specfile,specunit,kmodelname,nloopctl,           &
     &        nparms,parname,partype,parms,parcomm,atcredate,           &
     &        lunlog,tinf,xcol,vturbi,critf,radexp,                     &
     &        delr,rdel,enlum,xee,ababs,                                &
     &        bremsa,bremsint,tau0,dpthc,tauc,                          &
     &        ncsvn,nlsvn,                                              &
     &        ntotit,lnerrd,                                            &
     &        xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                 &
     &        cllines,clcont,htcomp,clcomp,clbrems,                     &
     &        xilev,bilev,rniss,                                        &
     &        rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,    &
     &        cabab,elumab,elum,zrems)                                  
!                                                                     
      write (lunlog,*) 'running ...' 
      write (tmpst,*) 'running ...' 
!                                                                       
!     Entering the main processing loop                                 
      if (numrec.le.0) npass=1 
      do kk=1,npass 
!                                                                       
        if ((lwri.gt.0).or.(npass.gt.1)) then 
!         Open and initialize FITS file for detailed ionic data         
          knam='xout_detail.fits' 
          call fnappend(knam,kk) 
          iunito=iunit 
          call fheader(iunit,knam,atcredate,kmodelname,istatus) 
          knam2='xout_detal2.fits' 
          call fnappend(knam2,kk) 
          iunit2o=iunit2 
          call fheader(iunit2,knam2,atcredate,kmodelname,istatus) 
          if (istatus .gt. 0)call printerror(lunlog,istatus) 
          knam3='xout_detal3.fits' 
          call fnappend(knam3,kk) 
          iunit3o=iunit3 
          call fheader(iunit3,knam3,atcredate,kmodelname,istatus) 
          knam4='xout_detal4.fits' 
          call fnappend(knam4,kk) 
          iunit4o=iunit4 
          call fheader(iunit4,knam4,atcredate,kmodelname,istatus) 
          call fparmlist(iunit,1,kmodelname,nparms,parname,partype,     &
     &               parms,parcomm,nloopctl,istatus,lunlog)             
          call fparmlist(iunit2,1,kmodelname,nparms,parname,partype,    &
     &               parms,parcomm,nloopctl,istatus,lunlog)             
          call fparmlist(iunit3,1,kmodelname,nparms,parname,partype,    &
     &               parms,parcomm,nloopctl,istatus,lunlog)             
          call fparmlist(iunit4,1,kmodelname,nparms,parname,partype,    &
     &               parms,parcomm,nloopctl,istatus,lunlog)             
          endif 
!                                                                       
        if (kk.eq.1) then 
        if (spectype.eq.'file    ')                                     &
     &    call ispecg(eptmp,zrtmp,numcon2,epi,ncn2,zremsz,xlum,         &
     &                lpri,lunlog)                                    
        xlum2=1. 
        if (spectype.eq.'bbody   ') then 
          call starf(trad,xlum2,epi,ncn2,zremsz,lpri,lunlog) 
          endif 
        if (spectype.eq.'pow     ')                                     &
     &    call ispec4(trad,xlum,epi,ncn2,zremsz,lpri,lunlog)          
        if (spectype.eq.'brems   ')                                     &
     &    call ispec(trad,xlum,epi,ncn2,zremsz,lpri,lunlog)           
        call ispecgg(xlum,epi,ncn2,zremsz,                              &
     &             lpri,lunlog)                                       
        call ispcg2(zremsz,epi,ncn2,enlum,lpri,lunlog) 
          endif 
!                                                                     
      call init(lunlog,bremsa,bremsint,tau0,dpthc,dpthcont,tauc,        &
     &  xii,rrrt,pirt,htt,cll,httot,cltot,                              &
     &  cllines,clcont,htcomp,clcomp,clbrems,                           &
     &  xilev,rcem,oplin,rccemis,brcems,opakc,opakcont,                 &
     &  cemab,cabab,opakab,elumab,elumabo,elum,elumo,                   &
     &  zrems,zremso,fline,flinel)                                    
      lpri=0 
!                                                                     
!                                                                       
!       nb suppressing line escape                                      
!        write (lunlog,*) ' nb suppressing line escape'                 
!        write (tmpst,*) ' nb suppressing line escape'                  
!        call xwrite(tmpst,10)                                          
!        do mm=1,nnnl                                                   
!          tau0(1,mm)=1.e+10                                            
!          tau0(2,mm)=1.e+10                                            
!          enddo                                                        
!                                                                       
        ldir=(-1)**kk 
        write (lunlog,*) ' ' 
        write (tmpst,*) ' ' 
        call xwrite(tmpst,10) 
        write (lunlog,*) 'pass number=',kk,ldir 
        write (tmpst,*) 'pass number=',kk,ldir 
        call xwrite(tmpst,10) 
!                                                                       
!     Print a heading for this pass                                   
      call pprint(17,1,trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,              &
     &      numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,          &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
!                                                                     
      do jk=1,ncn2 
        zrems(1,jk)=zremsz(jk) 
        zremso(1,jk)=zremsz(jk) 
        enddo 
!
!       main loop, step thru radius zones                               
        jkp=0 
        rdelo=0. 
        rdel=0. 
        ierr=0 
        do while ((((kk.eq.1).and.(xcol.lt.xpxcol).and.(xee.gt.xeemin)  &
     &         .and.(t.gt.tinf*(0.99)).and.(numrec.gt.0))               &
     &         .or.((kk.gt.1).and.(jkp.lt.numrec))).and.(ierr.eq.0))    
                                                                        
!                                                                       
!         save for variable density                                     
          if ((kk.eq.1).and.(lcdd.eq.1).and.(jkp.eq.0)) then 
            xpx0=xpx 
            r0=r 
            endif 
                                                                        
          jkp=jkp+1 
          if (jkp.gt.3999) stop 'too many steps: buffer filled' 
!                                                                       
!         step forward                                                  
          if ((kk.gt.1).and.(kk.le.npass)) then 
             jk=numrec+1-jkp 
             lpriu=0 
             rdelo=rdel 
             if (ldir.gt.0) then 
                 nlimdt=0 
               else 
                 nlimdt=nlimd 
               endif 
             lpriu=0 
             call unsavd(jk+2,ldir,                                     &
     &       lpri,iunito,iunit2o,iunit3o,iunit4o,                       &
     &       t,p,r,rdel,delr,xcol,xee,xpx,zeta,                         &
     &       xilev,rniss,                                               &
     &       rcem,oplin,tau0,                                           &
     &       cemab,cabab,opakab,tauc,                                   &
     &       epi,ncn2,zrems,dpthc,opakc,rccemis,                        &
     &       lunlog,status)                                             
             endif 
          if (kk.eq.1) then 
             nlimdt=nlimd 
             delr=0. 
             ectt=1. 
             jk=jkp 
             lpris=0 
             if (jk.gt.1)                                               &
     &         call step(ectt,emult,epi,ncn2,opakc,rccemis,fline,       &
     &           zrems,lpris,delr,dpthc,r,                              &
     &           xpxcol,xcol,xpx,taumax,numrec,lunlog)                  
            endif 
!
!          calculate flux                                                
          r19=r*(1.e-19) 
          xi=xlum/r19/r19/xpx 
          zeta=log10(xi) 
           lpri3=0
           call trnfrc(lpri3,lunlog,ldir,                               &
     &      r,xpxcol,xpx,                                               &
     &      epi,ncn2,zremsz,dpthc,opakc,                                &
     &      zrems,bremsa,bremsint)                            
!                                                                     
          if (t.lt.tinf*(1.02)) then 
            t=tinf*(1.01) 
            endif 
!                                                                       
!          calculate temperature, ionization, etc.                       
           lpri2=0
           lprid=0
           call xstarcalc(lpri2,lnerrd,nlimd,                           &
     &       lpri3,lprid,lunlog,                                        &
     &       tinf,vturbi,critf,t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,&
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       ntotit,                                                    &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       xilev,bilev,rniss,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                     
!                                                                     
!
!          do tranfer.  assumes comp2 and brems have been called             
!          already                                                           
           call heatt(lpri,lunlog,                                      &
     &       t,r,cfrac,delr,xee,xpx,abel,                               &
     &       epi,ncn2,bremsa,                                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       zrems,zremso,elumab,elumabo,elum,elumo,                    &
     &       rcem,oplin,rccemis,opakc,opakcont,cemab,fline,flinel,      &
     &       brcems)
!                                                                       
!
!          printout for each step:                                           
           call pprint(9,jkp,trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,        &
     &      numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,          &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
!                                                                     
!         If this is the final pass over all zones                      
          if (kk.eq.npass) then 
!           Add the abundances information to the output array zrtmp    
            call pprint(12,jkp,trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,      &
     &      numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,          &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
            endif 
!         Write radial profile FITS extension in here...                
          if ((lwri.gt.0).or.(npass.gt.1))                              &
     &     call savd(jkp+1,ldir,                                        &
     &       lpri,iunit,iunit2,iunit3,iunit4,                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       t,p,r,rdel,delr,xcol,xee,xpx,zeta,abel,                    &
     &       xilev,rniss,                                               &
     &       rcem,oplin,tau0,                                           &
     &       cemab,cabab,opakab,tauc,                                   &
     &       epi,ncn2,zrems,dpthc,opakc,rccemis,                        &
     &       lunlog,status)                                             
!                                                                       
!         new position etc                                              
          if (radexp.lt.-99.) then 
              read (lun40,*,iostat=ierr)rnew,dennew 
              delr=rnew-r 
              if (delr.lt.0.) stop 'radius error' 
              r=rnew 
              xpx=dennew 
!              write (lunlog,*)'rnew,dennew:',rnew,dennew,ierr          
            else 
              r=r+delr 
              if (lcdd.eq.1)                                            &
     &          xpx=xpx0*(r/r0)**radexp                                 
            endif 
          rdel=rdel+delr 
          xcol=xcol+xpx*delr 
!                                                                       
!         transfer                                                      
          eliml=1. 
          elimh=1.0e6 
          lpril2=0 
          call stpcut(ldir,lpril2,lunlog,vturbi,                        &
     &         np2,ncsvn,nlsvn,                                         &
     &       epi,ncn2,opakc,opakcont,oplin,opakab,delr,t,               &
     &       dpthc,dpthcont,tau0,tauc,eliml,elimh) 
          call trnfrn(lpri,lunlog,                                      &
     &       nlsvn,ncsvn,ncn2,                                          &
     &       zrems,zremso,elumab,elumabo,elum,elumo)                    
!
!          All done looping over the radial zones                        
           enddo 
!

          if (kk.eq.1) numrec=jkp 
!                                                                       
!       another printout to get the last step                           
        call pprint(9,jkp,trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,           &
     &      numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,          &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
!       If this is the final pass over all zones                        
        if (kk.eq.npass) then 
!         Add the abundances information to the output array zrtmp      
          call pprint(12,jkp,trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,        &
     &      numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,          &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
            endif 
          if ((lwri.gt.0).or.(npass.gt.1))                              &
     &     call savd(jkp+1,ldir,                                        &
     &       lpri,iunit,iunit2,iunit3,iunit4,                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       t,p,r,rdel,delr,xcol,xee,xpx,zeta,abel,                    &
     &       xilev,rniss,                                               &
     &       rcem,oplin,tau0,                                           &
     &       cemab,cabab,opakab,tauc,                                   &
     &       epi,ncn2,zrems,dpthc,opakc,rccemis,                        &
     &       lunlog,status)                                             
        if (kk.gt.1) then 
          call fitsclose(lunlog,iunito,istatus) 
          call fitsclose(lunlog,iunit2o,istatus) 
          call fitsclose(lunlog,iunit3o,istatus) 
          call fitsclose(lunlog,iunit4o,istatus) 
          endif 
                                                                        
        lpriu=0 
!                                                                       
!       End of the main iteration loop                                  
        enddo 

!                                                                       
!                                                                       
      write (lunlog,*)' ' 
      write (lunlog,*)' final print:',lprisv 
      write (tmpst,*)' final print:',lprisv 
      lpri=lprisv 
      call xwrite(tmpst,10) 
      lpril=0
      if (lpri.ge.2) lpril=1 
      delr=1.e-15
      call func(lpril,lunlog,vturbi,critf,                              &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,           &
     &       cllines,clcont,htcomp,clcomp,clbrems,                      &
     &       xilev,bilev,rniss,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                            
      call funcsyn(lpril,lunlog,vturbi,critf,                           &
     &       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,                  &
     &       epi,ncn2,bremsa,bremsint,                                  &
     &       tau0,tauc,                                                 &
     &       np2,ncsvn,nlsvn,                                           &
     &       xii,                                                       &
     &       xilev,bilev,rniss,                                         &
     &       rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,            &
     &       cabab,opakab,fline,flinel)                           
      call heatt(lpril,lunlog,                                          &
     &       t,r,cfrac,delr,xee,xpx,abel,                               &
     &       epi,ncn2,bremsa,                                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       zrems,zremso,elumab,elumabo,elum,elumo,                    &
     &       rcem,oplin,rccemis,opakc,opakcont,cemab,fline,flinel,      &
     &       brcems)
      lpril2=0 
      call stpcut(ldir,lpril2,lunlog,vturbi,                            &
     &         np2,ncsvn,nlsvn,                                         &
     &       epi,ncn2,opakc,opakcont,oplin,opakab,delr,t,               &
     &       dpthc,dpthcont,tau0,tauc,eliml,elimh) 
      if (lpri.gt.0) then 
        write (lunlog,9902)t,httot,cltot,hmctot 
 9902   format (4(1pe16.8)) 
        write (lunlog,*)' ' 
        do ll=11,17
          ntmp=nlprnt(ll) 
          call pprint(ntmp,numrec,trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,   &
     &      numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,          &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
          enddo 
        endif 
!                                                                       
!                                                                       
!     final print                                                       
!     Iterate through the list of 'final' reports                       
      call pprint(22,numrec,trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,         &
     &      numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,          &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
      if (lpri.ge.0) then
        do ll=2,nlnprnt 
          call pprint(nlprnt(ll),numrec,trad,xlum,lwri,lpri,r,t,xpx,p,  &
     &      lcdd,numrec,npass,nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,     &
     &      zremsz,epi,ncn2,abel,cfrac,emult,taumax,xeemin,             &
     &      spectype,specfile,specunit,kmodelname,nloopctl,             &
     &      nparms,parname,partype,parms,parcomm,atcredate,             &
     &      lunlog,tinf,xcol,vturbi,critf,radexp,                       &
     &      delr,rdel,enlum,xee,ababs,                                  &
     &      bremsa,bremsint,tau0,dpthc,tauc,                            &
     &      ncsvn,nlsvn,                                                &
     &      ntotit,lnerrd,                                              &
     &      xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,                   &
     &      cllines,clcont,htcomp,clcomp,clbrems,                       &
     &      xilev,bilev,rniss,                                          &
     &      rcem,oplin,rccemis,brcems,opakc,opakcont,cemab,opakab,      &
     &      cabab,elumab,elum,zrems)                                  
            enddo 
          endif
!                                                                       
!                                                                       
!     Write spectral data file xout_spect1.fits                         
      if (lwri.ge.0) then 
      write(6,*)'xstar: Prepping to write spectral data ' 
      lpril=0
      call writespectra(lunlog,lpril,nparms,parname,partype,parms,      &
     &        parcomm,atcredate,t,vturbi,epi,ncn2,dpthc,                &
     &        np2,nlsvn,                                                &
     &        elum,zrems,zremsz,kmodelname,nloopctl)            
      write (lunlog,*)'after writespectra' 
      lpril=0
      call writespectra2(lunlog,lpril,nparms,parname,partype,parms,     &
     &        parcomm,atcredate,epi,ncn2,dpthc,                         &
     &        np2,nlsvn,                                                &
     &        elum,tau0,kmodelname,nloopctl)                    
      write (lunlog,*)'after writespectra2' 
      lpril=0
      call writespectra3(lunlog,lpril,nparms,parname,partype,parms,     &
     &        parcomm,atcredate,epi,ncn2,dpthc,dpthcont,                &
     &        np2,                                                      &
     &        elum,zrems,zremsz,kmodelname,nloopctl)            
      write (lunlog,*)'after writespectra3' 
      call writespectra4(lunlog,lpril,nparms,parname,partype,parms,     &
     &        parcomm,atcredate,epi,ncn2,dpthc,ababs,                   &
     &        np2,ncsvn,                                                &
     &        elumab,tauc,kmodelname,nloopctl)                  
      write(6,*)'xstar: Done writing spectral data' 
      endif 
!                                                                       
!                                                                       
      call remtms(t2s) 
      ttot=abs(t2s-t1s) 
      write (lunlog,*)'total time',ttot 
      write (tmpst,*)'total time',ttot 
      call xwrite(tmpst,10) 
!                                                                       
      if ((lwri.gt.0).or.(npass.gt.1)) then 
        call fitsclose(lunlog,iunit,istatus) 
        call fitsclose(lunlog,iunit2,istatus) 
        call fitsclose(lunlog,iunit3,istatus) 
        call fitsclose(lunlog,iunit4,istatus) 
        endif 
      if (npass.gt.1) then 
        call fitsclose(lunlog,iunito,istatus) 
        call fitsclose(lunlog,iunit2o,istatus) 
        call fitsclose(lunlog,iunit3o,istatus) 
        call fitsclose(lunlog,iunit4o,istatus) 
        endif 
                                                                        
      close(unit=lunlog) 
!
      return 
      end                                           
