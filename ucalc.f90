      subroutine ucalc(ndesc,nrdesc,ml,lcon,jkion,vturbi,cfrac,         &
     &   nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,                       &
     &   ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,               &
     &   abund1,abund2,ptmp1,ptmp2,xpx,opakab,                          &
     &   opakc,opakcont,rccemis,lpriu,kdesc2,                           &
     &   rr,delr,t,trad,tsq,xee,xh1,xh0,                                &
     &   epi,ncn2,bremsa,bremsint,                                      &
     &   leveltemp,                                                     &
     &   nlev,lfast,lun11,                                              &
     &   np2,ncsvn,nlsvn)                
!                                                                       
!     Name: ucalc.f90  
!     Description:  
!       this routine calculates rates for all atomic processes            
!       author:  T. Kallman                                               
!     List of Parameters:
!       Input:
!       ndesc:  data type
!       nrdesc=rate type
!       ml=database record index
!       lcon=continuation flag (not used)
!       jkion=ion index
!       vturbi:  turbulent speed (km/s)
!       nrdt=number of reals
!       np1r=real pointer
!       nidt=number of integers
!       np1i=integer pointer
!       nkdt=number of chars
!       np1k=char pointer
!       idest1=level index of initial level
!       idest2=level index of final level
!       idest3=ion index of initial level
!       idest4=ion index of final level
!       abund1=population of initial levle
!       abund2=population of final level
!       ptmp1=escape probability in reverse direction
!       ptmp2=escape probability in the forward direction
!       xpx= H number density (cm^-3)
!       opakab(nnml):  rrc opacities (cm^-1)
!       opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!       opakcont(ncn):  continuum opacities lines excluded (cm^-1)
!       rccemis(2,ncn): continuum emissivities (erg cm^-3 s^-1 erg^-1) 
!                   /10^38
!                  inward and outward
!       lpriu=print switch
!       kdesc2=description string
!       r:  radius in nebula (cm)
!       delr=step size (cm)
!       t: temperature in 10^4K
!       trad: radiation temperature (for thermal spectrum) 
!               or energy index (for power law).
!       tsq=sqrt(t)
!       xee=electron fractional abundance relative to neutral H
!       xh1=H+ density
!       xh0=H0 density
!       epi(ncn)=energy grid (eV)
!       ncn2=length of epi
!       bremsa=flux (erg cm^-2 s^-1 erg^-1)
!       bremsint=integrated flux
!       nlev:  number of levels for the ion
!       lfast=fast switch
!       lun11=logical unit for printing
!       np2=number of atomic database records
!       ncsvn=number of rrcs
!       nlsvn=number of lines
!       Output:
!       ans1=forward rate
!       ans2=reverse rate
!       ans3=heating rate
!       ans4=cooling rate
!     Dependencies: many
!     Called by:  calc_ion_rates,calc_rates_level,calc_emis_ion
!
      use globaldata
      use times
      implicit none 
!                                                                       
      integer nptmpdim 
      parameter (nptmpdim=max(10000,ncn)) 
!
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,ndl) 
        integer:: ilev(10,ndl),nlpt(ndl),iltp(ndl) 
        character(1) :: klev(100,ndl) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
      character(49) kdesc(ntyp),kdesc2 
      character(29) krdesc(ntyp) 
      real(8) epi(ncn) 
      real(8) bremsa(ncn),bremsint(ncn) 
      real(8) rccemis(2,ncn),opakc(ncn),opakcont(ncn)
      real(8) aa(11),aaa(11,10),bbb(10),sg(ncn) 
      real(8) rstorey(5),dcfe(8),defe(8),alhe(2),alh(2) 
      real(8) etmpp(nptmpdim),stmpp(nptmpdim),ttmp(400),xsec(100) 
      real(8)  zc,eion,far,gam,scal,etmp8(ncn),stmp8(ncn) 
      real(8) scal2 
      real(8) a, aa1,aarec, aax, abund1, abund2, adi, aij,              &
     &     airt, al, algt, alm, alp, alph                               
      real(8) alpha, alppp, ans1, ans1o, ans2, ans2d, ans3,             &
     &     ans3d, ans4, ans4s, ansar2, ansar2o, ap, arad, atan, atmp, b,&
     &     bb, ans5, ans6
      real(8) bbb2, bbrec, bbx, bdi, beta, bethe, c, beth,              &
     &     cai, ccrec, ccx, ch, ch2, ch3, chi, chir, chitmp,            &
     &     cii                                                          
      real(8) cij, cijpp, cion, citmp1, citmp2, cji, clu, cn, cno,crate,&
     &     crec, crit53, csum, csum2, cul, d, ddd, ddx, delea           
      real(8) del1, del2, dele, delev, delt, den, dirt, dirtemp,        &
     &     e, e0, e1, eai, ecm, ediff,  ee1expo, eelo, eeup      
      real(8) eex, eexc, efnd, eij, eijry, ekt, elammu, elin, elo, em1, &
     &     em2ph, emax, enelec, ener, enn, ep, epii, erel, ergsev,eijkev
      real(8) eta, eth, etkh, etmp, ett, ett2, ettry, eup, exp10,       &
     &     expo,exptmp, f1, f2, fchi, ff, ff2, fh2lke, fi, flin,        &
     &     enerm
      real(8) float, fudge, gamma, gflin, ggl, gglo, ggu, ggup, hecxrt, &
     &     hij, opakab, texp, flinabs, opakb1,                          &
     &     p1, p2, p3, p4, p5, phi, phi1, phi2                          
      real(8) pi, pp, ppp, psi, ptmp1, ptmp2, q2, qq, r19, rate,        &
     &     rcemsum, rctmp1, rctmp2, rec,                                &
     &     rinf,rcem1,rcem2,rnist, rnissel, rnisseu, emltlv, ethion,    &
     &     bktm, bk                                       
      real(8) rm, rr, rrrt, rs, s0, scale, sd, se, sg0, delr,           &
     &     sgth, sigma, sigvtherm, sqrt, sscal, sth, sum                
      real(8) swrat, t, t0, t1, t3s2, t6, tbig, temp, tfnd,             &
     &     time1, time2, tk, tm, tmr, trad                          
      real(8) tsq, ttz, tz,                                             &
     &     upsil, upsiln, vth, vtherm, vturb, vturbi, wav, xee,   &
     &      xh0, xh1, xhe1, upsilon, cfrac
      real(8) xkt, xnx, xpx, xx, y, ya, ypow, yw, ywsq, yy, z1,         &
     &     zap, zeff, zz, zzz, y0,y1,yyqq                               
      real(8) dc,dt4,t2,term1,term2,term3,optst,opcrit,hcxrt 
!      real(8) er,ee1,ee2,ee3,er1,co,z2s,qij,sig,cr,crp,cr1
      real(8)  tmin,tmax,alphamilne,amilnerr                                
      real(8) tstr(100),cstr(100),rdattmp(100) 
      real(8) tt,e3,e2,rho,ee,term4,bbnurjp,bremtmpp,epiip,anstmp
!      real(8) min       !jg                                             
      integer nspline 
      integer indonly
      integer i57, ic, idest1, idest2, idest3, idest4,                  &
     &     ierr, ik, il, iltmp, int, iq, ist,                           &
     &     itmp, iz                                                     
      integer jj, jkion, jkk, jkk2, jkk3, jkkl, jlo, kdim, kl, l2, lcon,&
     &     lcon2, lf, lfast, lfastl, lfasto, lff,  lforce, li, li1,lii  
      integer lk, ll, lm, lorb, lpri, lprib, lpric, lpril, lprim,       &
     &     lprisv, lprit, lpriu, lrcalc, lrtyp, lrtyp2, lskp, ltyp,     &
     &     ltyp2, lun11, luse8, nkdti                                   
      integer lz, m, ml, ml2, ml3, mlion, mllz, mlp,                    &
     &     mm, mm5, mml, n, na, nb1, nbinc, nbmx,mlm                    
      integer ncn2, ncsvn, ndesc, ndtmp, nelin, nf, ni,                 &
     &     nidt, nidt2, nidti, nilin, nind, nistage, njj, nptmp         
      integer nkdt, nkdt2, nlev, nlevp, nll, nlsvn, nmin, nmx,          &
     &     nn,  nnz, nphint, npr, nprn, ndtmpo                          
      integer nq, nrdt, nrdti, nrdesc, nrdt2, nsh, nskp, ntcs,          &
     &     ntmp, ntmp2, nu, numcon2, nzel, nterm, np2                        
      integer lunsv,lfnd 
      integer np1r,np1i,np1k,np1r2,np1i2,np1k2 
      integer lctype
!     integer ncase,npts 
      integer ntem 
!                                                                       
!     Not used                                                          
      real(8) javir 
      integer javi 
!      character(80) javik                                              
!                                                                       
!      save aa,bb,ddd,ett,ggup,gglo,hij,opcrit,                          &
!     &         swrat,elin,pi,c,ergsev,etmp8,luse8                       
!                                                                       
      data bk/1.38062e-16/ 
!      data opcrit/1.d-39/                                              
      data opcrit/1.d-26/ 
      data ergsev/1.602197e-12/ 
      data pi/3.1415927/,c/2.997925e10/,luse8/0/ 
      data krdesc(1)/'ground state ionization      '/ 
      data krdesc(2)/'level ionization/recombinatio'/ 
      data krdesc(3)/'bound-bound collision        '/ 
      data krdesc(4)/'bound-bound radiative        '/ 
      data krdesc(5)/'bound-free collision (level) '/ 
      data krdesc(6)/'total recombination          '/ 
      data krdesc(8)/'total recombination          '/ 
      data krdesc(7)/'bound-free radiative (level) '/ 
      data krdesc(9)/'2 photon decay               '/ 
      data krdesc(11)/'element data                 '/ 
      data krdesc(12)/'ion data                     '/ 
      data krdesc(13)/'level data                   '/ 
      data krdesc(23)/'collisional superlevel->spect'/ 
      data krdesc(14)/'radiative superlevel->spect  '/ 
      data krdesc(15)/'CI total rate                '/ 
      data krdesc(40)/'CI from superlevels          '/ 
      data krdesc(41)/'non-radiative auger transtion'/ 
      data krdesc(42)/'Inner shell photoabsorption  '/ 
      data kdesc(1)/'radiative recombination:  aldrovandi and pequign '/ 
      data kdesc(2)/'charge exch. h0: Kingdon and Ferland             '/ 
      data kdesc(3)/'autoionization: hamilton, sarazin chevalier      '/ 
      data kdesc(4)/'line data radiative: mendosa; raymond and smith  '/ 
      data kdesc(5)/'2 photon transition collisional                  '/ 
      data kdesc(6)/'level data                                       '/ 
      data kdesc(7)/'dielectronic recombination: aldrovandi and pequi '/ 
      data kdesc(8)/'dielectronic recombination: arnaud and raymond   '/ 
      data kdesc(9)/'charge exch. H0 Kingdon and Ferland              '/ 
      data kdesc(10)/'charge exchange H+ Kingdon and Ferland          '/ 
      data kdesc(11)/'2 photon radiative                              '/ 
      data kdesc(12)/'photoionization, excited levels: hydrogenic     '/ 
      data kdesc(13)/'element data:                                   '/ 
      data kdesc(14)/'ion data:                                       '/ 
      data kdesc(15)/'photoionization: barfield koontz and huebner    '/ 
      data kdesc(16)/'arnaud and raymond ci                           '/ 
      data kdesc(17)/'collisional excitation hydrogenic: cota         '/ 
      data kdesc(18)/'radiative recombination hydrogenic: cota        '/ 
      data kdesc(19)/'photoionization: hullac                         '/ 
      data kdesc(20)/'charge exchange H+ Kingdon and Ferland          '/ 
      data kdesc(21)/'pixc bkh continued 3                            '/ 
      data kdesc(22)/'dielectronic recombination: storey              '/ 
      data kdesc(23)/'photoionization, excited levels: clark          '/ 
      data kdesc(24)/'pi xc clark continued                           '/ 
      data kdesc(25)/'collisional ionization: raymond and smith       '/ 
      data kdesc(26)/'collisional ionization hydrogenic: cota         '/ 
      data kdesc(27)/'photoionization: hydrogenic                     '/ 
      data kdesc(28)/'line data collisional: mendosa; raymond and smi '/ 
      data kdesc(29)/'collisional ionization data: scaled hydrogenic  '/ 
      data kdesc(30)/'radiative recombination hydrogenic: gould and t '/ 
      data kdesc(31)/'line data no levels                             '/ 
      data kdesc(32)/'collisional ionization: cota                    '/ 
      data kdesc(33)/'line data collisional: hullac                   '/ 
      data kdesc(34)/'line data radiative: mendosa; raymond and smitha'/ 
      data kdesc(35)/'photoionization: table (from bkh)               '/ 
      data kdesc(36)/'photoionization, excited levels:hydrogenic(no l)'/ 
      data kdesc(37)/'iron 3pq dr data from badnell                   '/ 
      data kdesc(38)/'total rr  from badnell amdpp.phys.strath.ac.uk  '/ 
      data kdesc(39)/'total dr  from badnell amdpp.phys.strath.ac.uk  '/ 
      data kdesc(40)/'                                                '/ 
      data kdesc(41)/'                                                '/ 
      data kdesc(42)/'                                                '/ 
      data kdesc(43)/'total photoionization cross sections tabulated  '/ 
      data kdesc(44)/'                                                '/ 
      data kdesc(45)/'                                                '/ 
      data kdesc(46)/'                                                '/ 
      data kdesc(47)/'                                                '/ 
      data kdesc(48)/'                                                '/ 
      data kdesc(49)/'op pi xsections for inner shells                '/ 
      data kdesc(50)/'op line rad. rates                              '/ 
      data kdesc(51)/'op and chianti line coll rates                  '/ 
      data kdesc(52)/'same as 59 but rate type 7                      '/ 
      data kdesc(53)/'op pi xsections                                 '/ 
      data kdesc(54)/'h-like cij, bautista (hlike ion)                '/ 
      data kdesc(55)/'hydrogenic pi xsections, bautista format        '/ 
      data kdesc(56)/'tabulated collision strength, bautista          '/ 
      data kdesc(57)/'effective charge to be used in coll. ion.       '/ 
      data kdesc(58)/'hlike rec rates, bautista                       '/ 
      data kdesc(59)/'verner pi x!                                    '/ 
      data kdesc(60)/'calloway h-like coll. strength                  '/ 
      data kdesc(62)/'calloway h-like coll. strength                  '/ 
      data kdesc(61)/'h-like cij, bautista (non-hlike ion)            '/ 
      data kdesc(63)/'h-like cij, bautista (hlike ion)                '/ 
      data kdesc(64)/'hydrogenic pi xsections, bautista format        '/ 
      data kdesc(65)/'effective charge to be used in coll. ion.       '/ 
      data kdesc(66)/'Like type 69 but, data in fine structure.       '/ 
      data kdesc(67)/'Effective collision strengths from Keenan et al.'/ 
      data kdesc(68)/'coll. strength He-like ions by Zhang & Sampason '/ 
      data kdesc(69)/'Kato & Nakazaki (1996) fit to Helike coll. strgt'/ 
      data kdesc(70)/'Coefficients for phot x-section of suplevels    '/ 
      data kdesc(71)/'Transition rates from superlevel to spect. lvls '/ 
      data kdesc(72)/'Autoinization rates (in s^-1) for satellite lvls'/ 
      data kdesc(73)/'Fit to coll. strengths satellite lvls Helike ion'/ 
      data kdesc(74)/'Delta functions to add to phot. x-sections  DR  '/ 
      data kdesc(75)/'autoionization data for Fe XXiV satellites      '/ 
      data kdesc(76)/'2 photon decay                                  '/ 
      data kdesc(77)/'coll rates from 71                              '/ 
      data kdesc(78)/'Auger level data                                '/ 
      data kdesc(79)/'fluorescence line data                          '/ 
      data kdesc(80)/' Collisional ionization rates gnd of Fe and Ni  '/ 
      data kdesc(81)/' Bhatia Fe XIX collision strengths              '/ 
      data kdesc(82)/' Fe UTA rad rates                               '/ 
      data kdesc(83)/' Fe UTA level data                              '/ 
      data kdesc(84)/' Iron K Pi xsections, spectator Auger binned    '/ 
      data kdesc(85)/' Iron K Pi xsections, spectator Auger summed    '/ 
      data kdesc(86)/' Iron K Auger data from Patrick                 '/ 
      data kdesc(88)/' Iron inner shell resonance excitation (Patrick)'/ 
      data kdesc(89)/' saf line wavelengths same as 50                '/ 
      data kdesc(91)/' aped line wavelengths same as 50               '/ 
      data kdesc(92)/' aped collision strengths                       '/ 
      data kdesc(93)/' OP PI xsections?                               '/ 
      data kdesc(94)/' OP PI xsections?                               '/ 
      data kdesc(95)/' Bryans CI rates                                '/ 
      data kdesc(96)/' Fe XXiV satellites from safranova              '/ 
      data kdesc(97)/' CI rates from inner shells from palmeri 2016   '/ 
      data kdesc(98)/' chianti2016 collisional rates                  '/ 
      data kdesc(99)/' new type 70                                    '/ 

      save kdesc,ergsev,krdesc,pi,c,luse8,bk,opcrit
                                                                        
      javir=trad 
!      trad=javir                                                       
      javi=nlsvn 
!      nlsvn=javi                                                       
      javi=derivedpointers%npcon(1) 
      javi=derivedpointers%npconi(1) 
      javi=derivedpointers%npilev(1,1) 
      javi=derivedpointers%npilevi(1) 
      javi=derivedpointers%npconi2(1) 
      javi=ncsvn 
!      javik=krdesc(1)                                                  
!                                                                       
      call remtms(time1) 
!                                                                       
      xnx=xpx*xee 
!                                                                       
      lpri=lpriu 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc:',ndesc,lcon,nrdt,nidt,nkdt,        &
     &  ml,(masterdata%rdat1(np1r+mm-1),mm=1,nrdt),                     &
     &  (masterdata%idat1(np1i+mm-1),mm=1,nidt),                        &
     &  (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)                        
       if (lpri.gt.1) write (lun11,*)'in ucalc, inputs:',            &
     &   t,xee,xpx,xnx                                                  
!
!                                                                       
      vturb=vturbi 
!                                                                       
      kdesc2=kdesc(ndesc) 
!                                                                       
      if (luse8.eq.0) then 
        luse8=1 
        do mm=1,ncn2 
          etmp8(mm)=dble(epi(mm)/13.605692) 
          enddo 
        endif 
      nlevp=nlev 
      ans1=0. 
      ans2=0. 
      ans3=0. 
      ans4=0. 
      ans5=0. 
      ans6=0. 
      idest1=0 
      idest2=0 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=masterdata%idat1(np1i+nidt-1)+1 
      opakab=0. 
      lforce=1 
!
!     nb indonly will be needed later
      indonly=0
!
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,                    &
     &  17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,       &
     &  36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,    &
     &  56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,    &
     &  76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,    &
     &  96,97,98,99),                                                   &
     &  ndesc                                                           
!                                                                       
!                                                                       
!     rr, a&p formula                                                   
    1 continue 
!      write (lun11,*)'in ucalc, ndesc=1'                               
      arad=masterdata%rdat1(np1r) 
      eta=masterdata%rdat1(np1r+1) 
      rrrt=arad/t**eta 
      ans1=rrrt*xnx 
      idest1=1 
      idest2=0 
!      write (lun11,*)'in ucalc, ndesc=1'                               
      go to 9000 
!                                                                       
!     h charge exchange recombination                                   
    2 continue 
      if (t.gt.5.) go to 9000 
      aax=masterdata%rdat1(np1r) 
      bbx=masterdata%rdat1(np1r+1) 
      ccx=masterdata%rdat1(np1r+2) 
      ddx=masterdata%rdat1(np1r+3) 
      rate=aax*expo(log(t)*bbx)*max(0.d0,(1.+ccx*expo(ddx*t)))    &
     & *(1.d-9) 
      ans1=rate*xh0 
!      if (lpri.ge.1) write (lun11,*)'note turning off charge exchange' 
!      ans1=0.                                                          
      ans2=0. 
      if (nrdesc.eq.5) then 
        ans2=rate*xh0 
        ans1=0. 
        endif 
      idest1=1 
      idest2=nlevp 
      if (lpri.gt.1) write (lun11,*)'type 2 data',aax,bbx,ccx,ddx,rate, &
     &                               xh0,ans1,idest2                    
      go to 9000 
!      beth=rdat1(np1r)                                                 
!      alh(1)=rdat1(np1r+1)                                             
!      alh(2)=rdat1(np1r+2)                                             
!      ntcs=2                                                           
!      if (t.lt.1.) ntcs=1                                              
!      hcxrt = beth*t**alh(ntcs)                                        
!      xh1 = xiin(1)*xpx                                                
!      xh1=0.                                                           
!      xh2 =max(0.,(1.-xiin(1)))*xpx                                    
!      ans1=hcxrt*xh1                                                   
!      idest1=1                                                         
!      idest2=0                                                         
!      go to 9000                                                       
!                                                                       
    3 continue 
!     autoionization rates                                              
      ekt = t*(0.861707) 
      cai=masterdata%rdat1(np1r) 
      eai=masterdata%rdat1(np1r+1) 
      airt = cai*expo(-eai/ekt)/tsq 
      ans1=airt*xnx 
      idest1=1 
      idest2=1 
      go to 9000 
!                                                                       
    4 continue 
!     line rates, coll and rad                                          
!           write (lun11,*)'level data'                                 
!           do 1906 ll=1,nlev                                           
!             write (lun11,*)ll,(rlev(mm,ll),mm=1,3),                   
!     $          (ilev(mm,ll),mm=1,3),(klev(mm,ll),mm=1,3)              
! 1906        continue                                                  
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      elin=abs(masterdata%rdat1(np1r)) 
      flin=masterdata%rdat1(np1r+1) 
!      if (flin.le.1.d-10) flin=1.                                      
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         endif 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      a=masterdata%rdat1(np1r+4) 
      hij=elin*1.d-8 
      elammu=elin*1.d-4 
      aij=(6.67e+7)*gglo*flin/ggup/elammu/elammu 
!     this is a fudge to avoid badnumerics from fine structure.         
      if (flin.le.1.01e-12) aij=1.e+5 
      if (elin.ge.1.e+9) aij=1.e+5 
      ans1=aij*(ptmp1+ptmp2) 
      ans4=aij*(ptmp1+ptmp2)*ergsev*12398.41/abs(elin) 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      sigma=(0.02655)*flin*elin*(1.d-8)/vtherm 
      sigvtherm=sigma 
      ener=12398.41/abs(elin) 
      nb1=nbinc(ener,epi,ncn2) 
      ans2=0. 
!      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10                         
      if (elin.gt.0.99e+9) then 
         ans2=0. 
         sigvtherm=0. 
         endif 
      ans1=ans1+ans2*ggup/(1.d-48+gglo) 
      opakab=sigvtherm 
      ans3=ans2*ener*ergsev 
      ans4=ans1*ener*ergsev 
!      write (lun11,*)'ltyp=4',idest1,idest2,elin,flin,ggup,gglo        
      go to 9000 
!                                                                       
    5 continue 
!     2photon rates, col                                                
      idest1=masterdata%idat1(np1i+1) 
      idest2=masterdata%idat1(np1i) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      ans1=0. 
      ans2=0. 
      ggup=leveltemp%rlev(2,masterdata%idat1(np1i+1)) 
      gglo=leveltemp%rlev(2,masterdata%idat1(np1i)) 
      hij=elin*1.d-8 
      ekt=t*(0.861707) 
      eex=abs(leveltemp%rlev(1,masterdata%idat1(np1i))                  &
     &      -leveltemp%rlev(1,masterdata%idat1(np1i+1))) 
      ans2=(8.629e-8)*masterdata%rdat1(np1r+4)                          &
     &      *t**masterdata%rdat1(np1r+5)/ggup 
      exptmp=expo(-eex/ekt) 
      ans1=(8.629e-8)*masterdata%rdat1(np1r+4)                          &
     &      *t**masterdata%rdat1(np1r+5)*exptmp/gglo 
!      write (lun11,*)'ltyp=5',idest1,idest2,elin,flin,ggup,gglo,       
!     $       ans2,eex,ekt,exptmp,ans1                                  
      go to 9000 
!      write (lun11,*)'in ucalc, ltyp=17',idat1(np1i),idat1(np1i+1),ggup
!     $          gglo,eex,ans2,ans1                                     
!                                                                       
!     level quantities, partition function                              
    6 continue 
      go to 9000 
!                                                                       
!     dr, a&p formula                                                   
    7 continue 
      adi=masterdata%rdat1(np1r) 
      bdi=masterdata%rdat1(np1r+1) 
      t0=masterdata%rdat1(np1r+2) 
      t1=masterdata%rdat1(np1r+3) 
      ap=1. 
      dirt=adi*ap*(1.d-06)*expo(-t0/t)                                  &
     &  *(1.+bdi*expo(-t1/t))/(t*sqrt(t))                               
      ans1=dirt*xnx 
      if (lpri.gt.1) write (lun11,*)'type 7 data:',                     &
     &  adi,bdi,t0,t1,dirt,xnx,ans1                                    
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
!     dr, arnaud and raymond                                            
    8 continue 
      dirt=0. 
      ekt=0.861707*t 
      t3s2=t**(-1.5) 
      tmr = 1.d-6*t3s2 
      do 820 n = 1,4 
        dcfe(n)=masterdata%rdat1(np1r+n-1) 
        defe(n)=masterdata%rdat1(np1r-1+n+4) 
        dirt = dirt + dcfe(n)*expo(-defe(n)/ekt) 
  820    continue 
      dirt = dirt*tmr 
      ans1=dirt*xnx 
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
!     he charge exchange                                                
    9 continue 
      go to 9000
      aax=masterdata%rdat1(np1r) 
      bbx=masterdata%rdat1(np1r+1) 
      ccx=masterdata%rdat1(np1r+2) 
      ddx=masterdata%rdat1(np1r-1+4) 
      texp=min(t,1000.d0)**bbx 
      rate=aax*texp*(1.+ccx*expo(ddx*t))*(1.d-9) 
      ans2=rate*xh0 
      ans1=0. 
      idest1=1 
      idest2=nlevp 
      if (nidt.gt.1) then 
        idest1=masterdata%idat1(np1i) 
        idest2=nlevp+masterdata%idat1(np1i+1)-1 
        ans2=ans2/6. 
        endif 
!      if (jkion.eq.18) idest1=3                                        
      go to 9000 
      bethe=masterdata%rdat1(np1r) 
      alhe(1)=masterdata%rdat1(np1r+1) 
      alhe(2)=masterdata%rdat1(np1r+2) 
      ntcs=2 
      if (t.lt.1.) ntcs=1 
      hecxrt = bethe*t**alhe(ntcs) 
!      xhe1 = xiin(2)*xpx                                               
      xhe1=0. 
      ans1=hecxrt*xhe1 
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
!                                                                       
   10 continue 
      go to 9000
!     charge transfer ionzation as used in calc_rates_level, for level rates
      aax=masterdata%rdat1(np1r) 
      bbx=masterdata%rdat1(np1r+1) 
      ccx=masterdata%rdat1(np1r+2) 
      ddx=masterdata%rdat1(np1r-1+4) 
      eex=masterdata%rdat1(np1r-1+7) 
      rate=aax*t**bbx*(1.+ccx*expo(ddx*t))                           &
     &             *expo(-eex/t)*(1.d-9)                                
      ans1=rate*xh1 
      ans2=0. 
!     this is tricky: calc_ion_rates only uses the rate type 15 rate 
!       if idest1=1 
!     for O I we have idest1=1,2,3, so it's OK                          
      idest1=masterdata%idat1(np1i) 
      idest2=nlevp 
      go to 9000 
!                                                                       
   11 continue 
      idest2=masterdata%idat1(np1i) 
      idest1=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      dele=abs(leveltemp%rlev(1,idest2)-leveltemp%rlev(1,idest1)) 
      ggl=masterdata%rdat1(np1r+2) 
      ggu=masterdata%rdat1(np1r-1+4) 
      ans1=(6.669e+15)*masterdata%rdat1(np1r+1)*ggl                     &
     &   /(ggu*masterdata%rdat1(np1r)*masterdata%rdat1(np1r)) 
      ans4=ans1*dele*ergsev 
      go to 9000 
!                                                                       
   12 continue 
      go to 36 
!                                                                       
   13 continue 
      go to 9000 
!                                                                       
   14 continue 
      go to 9000 
!                                                                       
   15 continue 
      lprisv=lpri 
      if (lpri.gt.1) write(lun11,*)'ltyp=15',ml,derivedpointers%npar(ml) 
!      if (lpri.gt.1) write (lun11,*)(rdat1(np1r-1+jj),jj=1,nrdt)       
!      if (lpri.gt.0) write (lun11,*)(idat1(np1i-1+jj),jj=1,nidt)       
!      if (lpri.gt.0) write (lun11,*)(kdat1(np1k-1+jj),jj=1,nkdt)       
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      ntmp=nrdt/2 
      do ml2=1,ntmp+1 
        etmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2-1) 
        stmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2)*1.d-18 
        enddo 
      mlm=nilin
      call drd(ltyp2,lrtyp2,lcon2,                                      &
     &  nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                        &
     &  0,lun11)                                                  
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=masterdata%idat1(np1i+nidt-3)-masterdata%idat1(np1i+nidt-1) 
      if (indonly.eq.1) return
      if (lpri.gt.1)                                                    &
     & write (lun11,*)ml,nilin,masterdata%rdat1(np1r),idest1                 
      ett=masterdata%rdat1(np1r2) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'ett=',ett                                        
      nb1=nbinc(ett,epi,ncn2) 
      gglo=leveltemp%rlev(2,1) 
      ggup=leveltemp%rlev(2,nlevp) 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      d=masterdata%rdat1(np1r+1) 
      do mm=1,11 
        aa(mm)=masterdata%rdat1(np1r-1+3+mm) 
        enddo 
!      aa(1)=min(max(aa(1),-6.),6.)                                     
      ekt=t*(0.861707) 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc, ind=15:',lcon,                     &
     &               nrdt,nidt,nkdt                                     
      if (lpri.gt.1)                                                    &
     & write (lun11,891)(masterdata%rdat1(np1r-1+mm),mm=1,nrdt)        
  891 format (1x,10(1pe10.3)) 
      if (lpri.gt.1)                                                    &
     & write (lun11,892)(masterdata%idat1(np1i-1+mm),mm=1,nidt)              
  892 format (1x,10(i6)) 
      if (lpri.gt.1)                                                    &
     & write (lun11,893)(masterdata%kdat1(np1k-1+mm),mm=1,nkdt)             
  893 format (1x,100a1) 
      na=masterdata%idat1(np1i-1+nidt-5) 
      nsh=masterdata%idat1(np1i-1+nidt-4) 
      do lk=1,na 
        ll=masterdata%idat1(np1i-1+9*lk-5) 
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)'ll=',ll,lk                                     
        lz=15*(lk-1) 
        ett=masterdata%rdat1(np1r-1+1+lz) 
        ddd=masterdata%rdat1(np1r-1+2+lz) 
        bb=masterdata%rdat1(np1r-1+3+lz) 
        aa(1)=masterdata%rdat1(np1r-1+4+lz) 
        aa(2)=masterdata%rdat1(np1r-1+5+lz) 
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)'ltest=2',ett,ddd,bb,aa(1)                      
        do 1011 mml=1,5 
          aa(2+mml)=masterdata%rdat1(np1r-1+mml+5+lz) 
 1011     continue 
        do 1012 mml=1,4 
          aa(7+mml)=masterdata%rdat1(np1r-1+mml+10+lz) 
 1012     continue 
        bbb(lk)=bb 
        do 1013 mml=1,11 
          aaa(mml,lk)=aa(mml) 
          if (lpri.gt.1)                                                &
     &     write (lun11,*)mml,aa(mml)                                   
 1013     continue 
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)'ltest=0',ll,na,nsh,aa(7)                       
        enddo 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'calling bkhsgo:',ett,t,                       &
     & ddd,(bbb(mm),mm=1,3),na,                                         &
     & aaa(1,1),aaa(7,1)                                                
      lfastl=1 
      call bkhsgo(sg,ett,ddd,bbb,na,                                 &
     &         aaa,epi,ncn2,t,lprib,lfastl,lun11)                       
      lprib=0 
!      if (lpri.gt.1) lprib=lpri                                        
      gglo=leveltemp%rlev(2,1) 
      ggup=leveltemp%rlev(2,nlevp) 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      if (lpri.ge.1) then 
        npr=nbinc(ett,epi,ncn2) 
        write (lun11,*)'bkh threshold xsection:',                       &
     &         npr,ett,sg(npr)                                          
        endif 
      lpri=lprisv 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      go to 9000 
!                                                                       
   16 continue 
      ekt=t*(0.861707) 
      njj=int(nrdt/5) 
      lpriu=0 
      if (lpriu.ne.0)                                                   &
     & write (lun11,*)'ltyp=16:',masterdata%idat1(np1i+nidt-1)          
      csum=0. 
      csum2=0. 
      do mm=1,njj 
        mm5=5*(mm-1) 
        eth=masterdata%rdat1(np1r-1+mm5+1) 
        a=masterdata%rdat1(np1r-1+mm5+2) 
        b=masterdata%rdat1(np1r-1+mm5+3) 
        c=masterdata%rdat1(np1r-1+mm5+4) 
        d=masterdata%rdat1(np1r-1+mm5+5) 
        xx=eth/ekt 
        if (lpriu.ne.0)                                                 &
     &   write (lun11,*)'xx=',xx,eth,ekt,a,b,c,d,np1r,mm5               
        em1=ee1expo(xx) 
        f1=em1/xx 
        if (lpriu.ne.0)                                                 &
     &   write (lun11,*)'before ff2:',f1,em1,xx                         
        f2=ff2(xx,lpriu,lun11) 
        if (lpriu.ne.0)                                                 &
     &   write (lun11,*)xx,a,b,c,d,em1,f1,f2                            
        term1=a*(1.-xx*f1) 
        term2=b*(1.+xx-xx*(2.+xx)*f1) 
        term3=c*f1 
        term4=d*xx*f2 
        fi=term1+term2+term3+term4 
        fi=max(fi,0.d0) 
        csum=csum+fi*expo(-xx)/xx 
        csum2=csum2+fi/xx 
        if (lpriu.ne.0)                                                 &
     &   write (lun11,*)term1,term2,term3,term4                         
        if (lpriu.ne.0)                                                 &
     &   write (lun11,*)mm,mm5,a,b,c,d,xx,f1,fi,csum                    
        enddo 
      citmp1=csum*(6.69e-7)/ekt**(1.5) 
      ans1=citmp1*xnx 
      citmp2=csum2*(6.69e-7)/ekt**(1.5) 
      idest1=1 
      idest2=1 
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,1) 
!     note that rinf has exponential removed                            
      rinf=(2.08e-22)*gglo/ggup/t/tsq 
      ans2=citmp2*xnx 
      ans2=ans2*rinf*xnx 
      if (nrdesc.eq.5) then 
         idest2=nlevp 
        else 
          idest2=1 
        endif 
      if (lpriu.ne.0)                                                   &
     &   write (lun11,*)csum,citmp1,citmp2,ans1,                        &
     &     ggup,gglo,rinf,ans2,idest1,idest2                            
      go to 9000 
!                                                                       
   17 continue 
!     line rates, col                                                   
      ans1=0. 
      ans2=0. 
      hij=elin*1.d-8 
!      write (lun11,*)'ltyp=4',idest1,idest2,elin,flin,ggup,gglo        
      ekt=t*(0.861707) 
      idest1=masterdata%idat1(np1i+1) 
      idest2=masterdata%idat1(np1i) 
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      if (eelo.gt.eeup) then 
        idest1=masterdata%idat1(np1i) 
        idest2=masterdata%idat1(np1i+1) 
        eeup=leveltemp%rlev(1,idest2) 
        eelo=leveltemp%rlev(1,idest1) 
        endif 
      if (indonly.eq.1) return
      eex=eeup-eelo 
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      ans2=(8.629e-8)*masterdata%rdat1(np1r)                            &
     &      *t**masterdata%rdat1(np1r+1)/ggup 
      ans1=0. 
      exptmp=expo(-eex/ekt) 
      exptmp=1. 
      if (ekt.gt.eex/20.)                                               &
     & ans1=ans2*ggup*exptmp/gglo                                       
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in ucalc, ltyp=17',masterdata%idat1(np1i),    &
     &   masterdata%idat1(np1i+1),ggup,                                 &
     &   gglo,eex,ans2,ans1                                             
      go to 9000 
!                                                                       
   18 continue 
      aarec=masterdata%rdat1(np1r) 
      bbrec=masterdata%rdat1(np1r+1) 
      ccrec=masterdata%rdat1(np1r+2) 
      ttz=masterdata%rdat1(np1r-1+4) 
      algt=log10(t/(1.d-48+ttz))+4. 
      algt=max(algt,3.5d0) 
      algt=min(algt,7.5d0) 
      idest1=masterdata%idat1(np1i) 
      ans1=exp10(aarec+bbrec*(algt-ccrec)**2)/t/1.e+4 
      ans1=ans1*xnx 
!      ans1=0.                                                          
      ans2=0. 
      idest2=0 
      go to 9000 
!                                                                       
   19 continue 
      etkh=masterdata%rdat1(np1r-1+5) 
      enelec=1. 
      eth=etkh 
      nb1=nbinc(eth,epi,ncn2) 
      idest1=masterdata%idat1(np1i) 
      idest2=nlevp 
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,idest1) 
      swrat=gglo/ggup 
      ekt=t*(0.861707) 
      lm=nb1 
      do while (lm.le.nphint) 
         bbb2=epi(lm)/max(etkh,1.d-48) 
         etmp=log(bbb2) 
         alppp=masterdata%rdat1(np1r)+etmp*(masterdata%rdat1(np1r+1)+   &
     &         etmp*(etmp*masterdata%rdat1(np1r+2)                      &
     &        +etmp*masterdata%rdat1(np1r-1+4)))                    
         ppp=expo(alppp) 
         sg(lm)=(1.d-18)*enelec*ppp*(13.606)/etkh 
         call enxt(eth,nb1,lpri,epi,ncn2,t,lfastl,lun11,             &
     &                  lm,nskp,nphint,lrcalc)                          
        lm=lm+nskp 
        enddo 
      call phintfo(sg,eth,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      go to 9000 
!                                                                       
   20 continue 
!     charge transfer ionzation as used in calc_ion_rates, for total rate        
      aax=masterdata%rdat1(np1r) 
      bbx=masterdata%rdat1(np1r+1) 
      ccx=masterdata%rdat1(np1r+2) 
      ddx=masterdata%rdat1(np1r-1+4) 
      eex=masterdata%rdat1(np1r-1+5) 
      rate=aax*t**bbx*(1.+ccx*expo(ddx*t))                           &
     &             *expo(-eex/t)*(1.d-9)                                
      ans1=rate*xh1 
      ans2=0. 
      idest1=1 
      idest2=nlevp 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
      go to 9000 
!                                                                       
   21 continue 
!     h charge exchange dalgarno and butler                             
      beth=masterdata%rdat1(np1r) 
      alh(1)=masterdata%rdat1(np1r+1) 
      alh(2)=masterdata%rdat1(np1r+2) 
      ntcs=2 
      if (t.lt.1.) ntcs=1 
      hcxrt = beth*t**alh(ntcs) 
      ans1=hcxrt*xh1 
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
!     dr storey                                                         
   22 continue 
      ans1=0. 
      ans2=0. 
      idest1=1 
      idest2=0 
      if (t.gt.6.) go to 9000 
      do 221 kl=1,5 
        rstorey(kl)=masterdata%rdat1(np1r-1+kl) 
  221   continue 
!      if (rstorey(5).lt.0.) go to 9000                                 
      t3s2=t**(-1.5) 
      dirtemp=                                                          &
     &   (1.d-12)*(rstorey(1)/t+rstorey(2)                              &
     &   +t*(rstorey(3)+t*rstorey(4)))*t3s2                             &
     &   *expo(-rstorey(5)/t)                                           
      dirtemp=max(dirtemp,0.d0) 
      if (lpri.gt.1) write (lun11,*)'in ucalc, ltyp=22:',            &
     &   ndesc,lcon,nrdt,nidt,nkdt,                                     &
     &  ml,(masterdata%rdat1(np1r-1+mm),mm=1,nrdt),                     &
     &  (masterdata%idat1(np1i-1+mm),mm=1,nidt),                        &
     &  (masterdata%kdat1(np1k-1+mm),mm=1,nkdt),dirtemp,xnx                    
      ans1=dirtemp*xnx 
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
   23 continue 
      lfastl=1 
      lprisv=lpri 
!      lpri=2                                                           
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in ucalc, 23:',masterdata%rdat1(np1r),        &
     &       masterdata%rdat1(np1r+1),                                  &
     &  masterdata%rdat1(np1r+2),masterdata%rdat1(np1r-1+4),            &
     &  masterdata%rdat1(np1r-1+5),swrat,masterdata%idat1(np1i),        &
     &  masterdata%idat1(np1i+1),masterdata%idat1(np1i+2),              &
     &  masterdata%idat1(np1i-1+4),masterdata%idat1(np1i-1+5)
      ett=leveltemp%rlev(1,masterdata%idat1(np1i)) 
      eth=ett 
      nb1=nbinc(eth,epi,ncn2) 
      gglo=leveltemp%rlev(2,masterdata%idat1(np1i)) 
      ggup=leveltemp%rlev(2,nlevp) 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      ekt=t*(0.861707) 
      jkk2=masterdata%idat1(np1i+nidt-1) 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      jkk3=0 
      jkk=0 
      ndtmp=derivedpointers%npfirst(12) 
      do while ((jkk.ne.jkk2).and.(ndtmp.ne.0)) 
        jkk3=jkk3+1 
        mlm=ndtmp
        call drd(ltyp2,lrtyp2,lcon2,                                    &
     &    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                      &
     &    0,lun11)                                                
        jkk=masterdata%idat1(np1i2+nidt2-1) 
        ndtmp=derivedpointers%npnxt(ndtmp) 
        enddo 
      if (ndtmp.le.0) go to 9000 
      zzz=float(masterdata%idat1(np1i2)) 
      ndtmp=derivedpointers%npfi(13,jkk3) 
      mllz=derivedpointers%npar(ndtmp) 
      if (lpri.gt.1) write (lun11,*)jkk,jkk2,jkk3,zzz,ndtmp 
      iltmp=1 
      mlm=ndtmp
      call drd(ltyp2,lrtyp2,lcon2,                                      &
     &    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                      &
     &    0,lun11)                                                
      ndtmpo=ndtmp 
      ndtmp=derivedpointers%npnxt(ndtmp) 
      do while ((ndtmp.ne.0).and.(iltmp.ne.masterdata%idat1(np1i))      &
     &      .and.(derivedpointers%npar(ndtmp).eq.mllz))               
        mlm=ndtmp
        call drd(ltyp2,lrtyp2,lcon2,                                    &
     &    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                      &
     &    0,lun11)                                                
        iltmp=masterdata%idat1(np1i2+nidt2-2) 
        ndtmpo=ndtmp 
        ndtmp=derivedpointers%npnxt(ndtmp) 
        enddo 
      ndtmp=ndtmpo 
      nprn=masterdata%idat1(np1i) 
      enn=float(nprn) 
      if ((enn.le.1.d-24).or.(zzz.le.1.d-24).or.(ett.le.1.d-6))         &
     &  go to 9000                                                      
      sg0=6.3e-18*enn/zzz/zzz 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'ind=23:',ml,nilin,zzz,jkk,nprn,sg0,ett           
      ll=nb1 
      do while (ll.le.nphint) 
        epii=epi(ll) 
        sg(ll)=sg0*(epii/ett)**(-3) 
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,              &
     &                  ll,nskp,nphint,lrcalc)                          
        ll=ll+nskp 
        enddo 
      lprib=0 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
      lpri=lprisv 
!      ans2=ans2*xnx                                                    
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=nlevp 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      go to 9000 
!                                                                       
   24 continue 
      go to 9000 
!                                                                       
   25 continue 
      idest4=masterdata%idat1(np1i+nidt-1)+1 
      idest3=masterdata%idat1(np1i+nidt-1) 
      if (nrdesc.eq.5) then 
         idest1=masterdata%idat1(np1i+nidt-2) 
         idest2=nlevp 
        else 
         idest2=1 
         idest1=1 
        endif 
      if (indonly.eq.1) return
      e=masterdata%rdat1(np1r) 
      a=masterdata%rdat1(np1r+1) 
      b=masterdata%rdat1(np1r+2) 
      c=masterdata%rdat1(np1r-1+4) 
      d=masterdata%rdat1(np1r-1+5) 
      cion = 0. 
      chir = (t*1.e+4)/(11590.*e) 
      citmp1=cion 
      ans1=citmp1*xnx 
      ans2=0. 
!      idest2=1                                                         
      if ( chir.le..0115 ) go to 9000 
      chi = max(chir,0.1d0) 
      ch2 = chi*chi 
      ch3 = ch2*chi 
      alpha = (.001193+.9764*chi+.6604*ch2+.02590*ch3)                  &
     &        /(1.0+1.488*chi+.2972*ch2+.004925*ch3)                    
      beta = (-.0005725+.01345*chi+.8691*ch2+.03404*ch3)                &
     &       /(1.0+2.197*chi+.2457*ch2+.002503*ch3)                     
      ch = 1./chi 
      fchi = 0.3*ch*(a+b*(1.+ch)+(c-(a+b*(2.+ch))*ch)*alpha+d*beta*ch) 
      chitmp=expo(-1./chir) 
      cion = 2.2e-6*sqrt(chir)*fchi/(e*sqrt(e)) 
      citmp1=cion 
      ans1=citmp1*xnx 
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,nidt-1) 
!     note that rinf has exponential removed                            
      rinf=(2.08e-22)*gglo/ggup/t/tsq 
      ans2=ans1*rinf*xnx 
      ans1=ans1*chitmp 
      eth=e
      ans6=ans1*eth*ergsev
      ans5=ans2*eth*ergsev
      go to 9000 
!                                                                       
   26 continue 
      go to 9000 
!      ekt=t*(0.861707)                                                 
!      idest1=idat1(np1i)                                               
!      gglo=leveltemp%rlev(2,idest1)                                              
!      edelt=abs(rlev(1,idest1)-rlev(1,nlev))                           
!      exptmp=expo(-edelt/ekt)                                          
!      ans1=(4.1416e-9)*rdat1(np1r)*t**rdat1(np1r+1)/gglo               
!      ggup=leveltemp%rlev(2,nlev)                                                
!      rinf=(2.08e-22)*gglo/ggup/t/tsq                                  
!      ans2=ans1*rinf                                                   
!      ans1=ans1*exptmp                                                 
!      write (lun11,*)'ltyp=26',idest1,gglo,ggup,                       
!     $   edelt,rdat1(np1r),rdat1(np1r+1),ans1                          
!      idest2=nlev                                                      
!      go to 9000                                                       
!                                                                       
   27 continue 
      lprisv=lpri 
!      ett=masterdata%rdat1(np1r+1)                                           
      lfastl=1 
      idest1=1 
      if (nrdesc.eq.1) then 
          idest2=0 
        else 
          idest2=nlevp 
        endif 
      ett=abs(leveltemp%rlev(1,nlev)-leveltemp%rlev(1,1)) 
!      if (lpri.gt.1)                                                   
!      write (lun11,*)'in ucalc, ind=27:',rlev(1,nlev),                 
!     $     rlev(1,1),nlev,ett                                          
      if (ett.le.1.d-5) go to 9000 
      eth=ett 
      nb1=nbinc(eth,epi,ncn2) 
      gglo=leveltemp%rlev(2,masterdata%idat1(np1i)) 
      swrat=gglo 
      ekt=t*(0.861707) 
      ll=nb1 
      do while (ll.le.nphint) 
        epii=epi(ll) 
        e=epii 
        eth=ett 
        zap = e/eth - 1. 
        y = e/eth 
        yy=sqrt(zap) 
        yy=max(yy,1.d-04) 
        fh2lke=((6.3e-18)/masterdata%rdat1(np1r)**2                     &
     &   *y**(-4)*expo(4.-4.*atan(yy)/yy)                            &
     &   /(1.-expo(-6.2832/yy)))                                         
!        fh2lke=((6.3e-18)/rdat1(np1r)/rdat1(np1r))*y**(-3)             
        sg(ll)=fh2lke 
        if (lpri.ge.2) write (lun11,*)ll,epii,zap,y,yy,fh2lke 
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,              &
     &                  ll,nskp,nphint,lrcalc)                          
        ll=ll+nskp 
        enddo 
      lprib=0 
      if (lpri.gt.1) lprib=1 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      lpri=lprisv 
      go to 9000 
!                                                                       
   28 continue 
!     line rates, col                                                   
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1)).lt.                &
     &         leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      nind=5 
      if (nrdt.ge.12) then 
        do ll=1,4 
          ttmp(ll)=masterdata%rdat1(np1r-1+nrdt-4+ll) 
          enddo 
        jlo=0 
        call hunt3(ttmp,4,t,jlo,0,lun11) 
        nind=nrdt-8+jlo 
        endif 
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      elin=abs(masterdata%rdat1(np1r)) 
      hij=elin*1.d-8 
      if (elin.le.1.d-24) go to 9000 
!      nind=nrdt-2                                                      
      cijpp=masterdata%rdat1(np1r-1+nind) 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      cji=(8.626e-8)*cijpp/tsq/ggup 
      cij=0. 
      exptmp=expo(-delt) 
      cij=cji*ggup*exptmp/tsq/gglo 
      ans1=cij*xnx 
      ans2=cji*xnx 
      if (lpri.gt.1) then 
        write (lun11,*)'ltyp=28',idest1,idest2,elin,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,(masterdata%rdat1(np1r-1+mm)        &
     &      ,mm=1,8),nind,jlo 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp 
        endif 
      elin=0. 
      go  to 9000 
!                                                                       
   29 continue 
      go to 9000 
!      anstmp=rdat1(np1r+1)*(8.626e-8)/tsq                              
!      ans2=anstmp*(2.08e-22)*(rdat1(np1r+2)/rdat1(np1r-1+4))/t/tsq     
!      ans1=0.                                                          
!      delt=rdat1(np1r)/t                                               
!      if (delt.lt.50.) then                                            
!         exptmp=1.                                                     
!         exptmp=expo(-delt)                                            
!         ans1=anstmp*exptmp                                            
!         endif                                                         
!      idest1=idat1(np1i+1)                                             
!      idest2=nlev                                                      
!      write (lun11,*)'ltyp=29',ans1,ans2,(rdat1(np1r-1+ii),ii=1,4),anst
!      go to 9000                                                       
!                                                                       
   30 continue 
!      write (lun11,*)'ltyp=30',idat1(np1i)                             
        nmx=masterdata%idat1(np1i) 
        t6=t/100. 
        zeff=float(nmx) 
        beta=zeff*zeff/(6.34*t6) 
        yy=beta 
        vth=(3.10782e+7)*sqrt(t) 
!       fudge factor makes the 2 expressions join smoothly              
        ypow=min(1.d0,(0.06376)/yy/yy) 
        fudge=0.9*(1.-ypow)+(1./1.5)*ypow 
        phi1=(1.735+log(yy)+1./6./yy)*fudge/2. 
        phi2=yy*(-1.202*log(yy)-0.298) 
        phi=phi1 
        if (yy.lt.0.2525) phi=phi2 
        rrrt=2.*(2.105e-22)*vth*yy*phi 
        ans1=rrrt*xnx 
        ans2=0. 
        idest1=1 
        idest2=0 
      go to 9000 
!                                                                       
   31 continue 
!     line rates, coll and rad                                          
!           write (lun11,*)'level data'                                 
!           do 1906 ll=1,nlev                                           
!             write (lun11,*)ll,(rlev(mm,ll),mm=1,3),                   
!     $          (ilev(mm,ll),mm=1,3),(klev(mm,ll),mm=1,3)              
! 1906        continue                                                  
      idest1=masterdata%idat1(np1i+1) 
      idest2=masterdata%idat1(np1i) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      elin=abs(masterdata%rdat1(np1r)) 
      flin=masterdata%rdat1(np1r+1) 
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
!      a=rdat1(np1r-1+5)                                                
      ans1=0. 
      ans2=0. 
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000 
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      a=masterdata%rdat1(np1r2+1) 
      hij=elin*1.d-8 
      aij=(0.02655)*flin*8.*pi/hij/hij*gglo/(1.d-24+ggup) 
      ans1=aij*(ptmp1+ptmp2) 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      sigma=(0.02655)*flin*elin*(1.d-8)/vtherm 
      sigvtherm=sigma 
      ener=12398.41/abs(elin) 
      nb1=nbinc(ener,epi,ncn2) 
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10 
      if (elin.gt.0.99e+9) then 
         ans2=0. 
         sigvtherm=0. 
         endif 
      ans1=ans1+ans2*ggup/(1.d-48+gglo) 
!     notice that opakab does not have abundance in                     
      opakab=sigvtherm 
      ans3=ans2*ener*ergsev 
      ans4=ans1*ener*ergsev 
      delea=0. 
      lfasto=4 
      if (lfasto.ge.4) ans2=0. 
!      if (opakab.gt.1.d-34)                                            
!     $  call linopac(lpri,lun11,opakab,rcem1,rcem2,elin,vturb,t,a,     
!     $               delea,epi,ncn2,opakc,rccemis,
!     $               lfasto)                                           
      ans4=ans1*ener*ergsev 
!      write (lun11,*)'ltyp=31',idest1,idest2,elin,flin,ggup,gglo       
      go to 9000 
!                                                                       
   32 continue 
      idest1=masterdata%idat1(np1i) 
      if (indonly.eq.1) return
      gglo=masterdata%rdat1(np1r-1+4) 
      ans1=0. 
      ans2=0. 
      go to 9000 
!      if (gglo.lt.1.d-24) go to 9000                                   
!      ekt=t*(0.861707)                                                 
!      edelt=rdat1(np1r+2)                                              
!      ans1=(4.1416e-9)*rdat1(np1r)*t**rdat1(np1r+1)*expo(-edelt/ekt)   
!     $        /gglo                                                    
!      write (lun11,*)'ltyp=26',idest1,gglo,edelt,rdat1(np1r),rdat1(np1r
!      idest2=nlev                                                      
!      go to 9000                                                       
!                                                                       
   33 continue 
!     line rates, col                                                   
      idest1=masterdata%idat1(np1i+1) 
      idest2=masterdata%idat1(np1i) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      ggup=leveltemp%rlev(2,masterdata%idat1(np1i+1)) 
      gglo=leveltemp%rlev(2,masterdata%idat1(np1i)) 
      elin=abs(masterdata%rdat1(np1r)) 
      hij=elin*1.d-8 
      if (elin.le.1.d-24) go to 9000 
      nind=4 
      cijpp=masterdata%rdat1(np1r-1+nind) 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      exptmp=expo(-delt) 
      cij=(8.626e-8)*cijpp*exptmp/tsq/gglo 
      cji=(8.626e-8)*cijpp/tsq/ggup 
      ans1=cij*xnx 
      ans2=cji*xnx 
      go to 9000 
!                                                                       
   34 continue 
!     line rates, coll and rad                                          
!           write (lun11,*)'level data'                                 
!           do 1906 ll=1,nlev                                           
!             write (lun11,*)ll,(rlev(mm,ll),mm=1,3),                   
!     $          (ilev(mm,ll),mm=1,3),(klev(mm,ll),mm=1,3)              
! 1906        continue                                                  
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      elin=abs(masterdata%rdat1(np1r)) 
      aij=masterdata%rdat1(np1r+1) 
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         endif 
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
!      ggup=rdat1(np1r-1+4)                                             
!      gglo=rdat1(np1r+2)                                               
      a=masterdata%rdat1(np1r-1+5) 
      hij=elin*1.d-8 
      elammu=elin*1.d-4 
!      flin=aij*hij*hij*ggup/((0.02655)*8.*pi*gglo)                     
      flin=aij*hij*hij*ggup/((0.667274)*gglo) 
      ans1=aij*(ptmp1+ptmp2) 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      sigma=(0.02655)*flin*elin*(1.d-8)/vtherm 
      sigvtherm=sigma 
      ener=12398.41/abs(elin) 
      nb1=nbinc(ener,epi,ncn2) 
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10 
      if (elin.gt.0.99e+9) then 
         ans2=0. 
         sigvtherm=0. 
         endif 
      ans1=ans1+ans2*ggup/(1.d-48+gglo) 
!     notice that opakab does not have abundance in                     
      opakab=sigvtherm 
      delea=0. 
      lfasto=4 
      if (lfasto.ge.4) ans2=0. 
      ans3=ans2*ener*ergsev 
      ans4=ans1*ener*ergsev 
!      if (opakab.gt.1.d-34)                                            
!     $  call linopac(lpri,lun11,opakab,rcem1,rcem2,elin,vturb,t,a,     
!     $               delea,epi,ncn2,opakc,rccemis,
!     $               lfasto)                                           
!      if (lpri.gt.0)                                                   
!     $ write (lun11,*)'ltyp=34',idest1,idest2,elin,flin,ggup,gglo,     
!     $                         a,aij,hij,pi                            
      go to 9000 
!                                                                       
   35 continue 
      lprisv=lpri 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc, ind=15:',lcon                         
      ett=masterdata%rdat1(np1r) 
      if (ett.le.(1.d-24)) go to 9000 
      ntmp=(nrdt-1)/2 
      do ml2=1,ntmp 
        etmpp(ml2)=masterdata%rdat1(np1r-1+1+2*ml2) 
        stmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2) 
        enddo 
      nb1=nbinc(ett,epi,ncn2) 
      gglo=leveltemp%rlev(2,1) 
      ggup=leveltemp%rlev(2,nlevp) 
      if (ggup.le.1.d-24) then 
         write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      numcon2=max(2,ncn2/50) 
      nphint=ncn2-numcon2 
      idest1=masterdata%idat1(np1i-1+6) 
      idest2=masterdata%idat1(np1i-1+5)-masterdata%idat1(np1i-1+7) 
      if (indonly.eq.1) return
      ekt=t*(0.861707) 
      jlo=0 
      ll=nb1 
      lfastl=1 
      do while (ll.le.nphint) 
          epii=epi(ll) 
          efnd=(epii-ett)/13.605692 
          call hunt3(etmpp,ntmp,efnd,jlo,0,lun11) 
          ml2=jlo 
          mlp=ml2+1 
          del1=(efnd-etmpp(ml2))/(etmpp(mlp)-etmpp(ml2)) 
          del2=(efnd-etmpp(mlp))/(etmpp(mlp)-etmpp(ml2)) 
          sg(ll)=-stmpp(ml2)*del2+stmpp(mlp)*del1 
!          if (lpri.gt.1)                                               
!     $    write (lun11,*)ll,epii,sg(ll),ml2,stmpp(ml2),stmpp(mlp),     
!     $              del1,del2                                          
          call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,            &
     &                  ll,nskp,nphint,lrcalc)                          
          ll=ll+nskp 
          enddo 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
      if (lpri.gt.1) then 
        npr=nbinc(ett,epi,ncn2)+2 
        write (lun11,*)'bkh threshold xsection:',                       &
     &         npr,ett,sg(npr)                                          
        endif 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      lpri=lprisv 
      go to 9000 
                                                                        
!                                                                       
   36 continue 
!      photoionization, excited levels:hydrogenic(no l)                 
      lprisv=lpri 
!      lpri=2                                                           
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc, ind=36:',                          &
     &      (masterdata%idat1(np1i-1+mm),mm=1,5)    
      idest1=masterdata%idat1(np1i+nidt-2) 
      ett=leveltemp%rlev(1,nlevp)-leveltemp%rlev(1,idest1) 
      idest2=nlevp 
      if (indonly.eq.1) return
      if (ett.le.1.d-5) go to 9000 
      eth=ett 
      nb1=nbinc(eth,epi,ncn2) 
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,nlevp) 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      ekt=t*(0.861707) 
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in ucalc, ind=36:',                           &
     &   ml,nilin,nelin                                                 
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000 
      mlm=nilin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      nistage=masterdata%idat1(np1i2) 
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      nzel=masterdata%idat1(np1i2) 
      nq=masterdata%idat1(np1i) 
      nq=min(10,nq) 
      zz=float(nzel-nistage+1) 
      sgth=(6.3e-18)*nq*nq/zz/zz 
      if (lpri.gt.1) write (lun11,*)nb1,nq,nzel,nistage,zz,             &
     &                              ett,sgth,idest1,gglo,ggup           
      ll=nb1 
      lfastl=1 
      do while (ll.le.nphint) 
        epii=epi(ll) 
        sg(ll)=sgth*(epii/ett)**(-3) 
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,              &
     &                  ll,nskp,nphint,lrcalc)                          
        ll=ll+nskp 
        enddo 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      lpri=lprisv 
      go to 9000 
!                                                                       
   37 continue 
!     total dr for fe 3pq ions from badnell 2006 Ap. J. Lett 651 L73    
      dirt=0. 
      ekt=0.861707*t 
      t3s2=t**(-1.5) 
      tmr = 1.d-6*t3s2 
      nterm=masterdata%idat1(np1i) 
      do  n = 1,nterm 
        dcfe(n)=masterdata%rdat1(np1r-1+n) 
        defe(n)=masterdata%rdat1(np1r-1+n+4) 
        dirt = dirt + dcfe(n)*expo(-defe(n)/ekt) 
        enddo 
      dirt = dirt*tmr 
      ans1=dirt*xnx 
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
   38 continue 
!     total rr  from badnell http://amdpp.phys.strath.ac.uk/tamoc/DATA/D
      a=masterdata%rdat1(np1r) 
      b=masterdata%rdat1(np1r+1) 
      t0=masterdata%rdat1(np1r+2)/1.e+4 
      t1=masterdata%rdat1(np1r-1+4)/1.e+4 
      if (nrdt.gt.4) then 
        c=masterdata%rdat1(np1r-1+5) 
        t2=masterdata%rdat1(np1r-1+6)/1.e+4 
        b=b+c*exp(-t2/t) 
        endif 
      term1=(T/T0)**(0.5) 
      term2=(1.+(T/T0)**(0.5))**(1.-b) 
      term3=(1.+(T/T1)**(0.5))**(1.+b) 
      rrrt=a/(1.d-48+term1*term2*term3) 
      ans1=rrrt*xnx 
      if (lpri.gt.1) write (lun11,*)a,b,c,t0,t1,t2,                     &
     &         term1,term2,term3,rrrt,ans1                              
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
   39 continue 
!     total dr  from badnell http://amdpp.phys.strath.ac.uk/tamoc/DATA/D
      dirt=0. 
      ekt=0.861707*t 
      t3s2=t**(-1.5) 
      tmr = 1.d-6*t3s2 
      nterm=nrdt/2 
      do  n = 1,nterm 
        dc=masterdata%rdat1(np1r-1+n) 
        dt4=masterdata%rdat1(np1r-1+n+nterm)/1.e+4 
        dirt = dirt + dc*exp(-dt4/t) 
        if (lpri.gt.1) write (lun11,*)n,dc,dirt 
        enddo 
      dirt = dirt*tmr 
      ans1=dirt*xnx 
      if (lpri.gt.1) write (lun11,*)nterm,dirt,ans1 
      idest1=1 
      idest2=0 
      go to 9000 
!                                                                       
   40 continue 
      go to 9000 
!                                                                       
   41 continue 
      go to 9000 
!                                                                       
   42 continue 
      go to 9000 
!                                                                       
   43 continue 
!     total photoionization cross sections tabulated in                 
!     format like 53 (not used)                                         
      go to 9000 
!                                                                       
   44 continue 
      go to 9000 
!                                                                       
   45 continue 
      go to 9000 
!                                                                       
   46 continue 
      go to 9000 
!                                                                       
   47 continue 
      go to 9000 
!                                                                       
   48 continue 
      go to 9000 
!                                                                       
   49 continue 
  499 continue 
!     op pi xsections                                                   
!     old version                                                       
      lprisv=lpri 
!      if (lpri.ge.1) lpri=2                                            
!     these are the initial and final levels and indeces                
!     notice that these are relative to the current ion                 
!     (not relative to the element as a whole)                          
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest4=masterdata%idat1(np1i+nidt-3) 
      idest2=nlevp+max(0,masterdata%idat1(np1i-1+nidt-3))-1 
      if (indonly.eq.1) return
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2,nlevp
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000 
      if (ml.le.0) go to 9000 
      eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
      ett=eth 
      nilin=derivedpointers%npar(ml) 
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml 
      if (nilin.le.0) go to 9000 
      ntmp=nrdt/2 
      do ml2=1,ntmp+1 
        etmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2-1) 
        stmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2)*1.d-18 
        stmpp(ml2)=max(stmpp(ml2),0.d0) 
        enddo 
!      ett=ett+max(0.,13.605692*etmpp(1))                               
      optst=abund1*stmpp(1) 
!      if ((optst.lt.opcrit).and.(lfast.eq.2)) go to 9000               
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1) 
      if (ett.le.0.) go to 9000 
      ntmp2=nptmpdim 
      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11) 
      nb1=nbinc(ett,epi,ncn2) 
      xkt=ett/(0.861707*t) 
      r19=rr/1.e+19 
      mlm=nilin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdti,np1r2,nidti,np1i2,nkdti,np1k2,mlm,                        &
     &  0,lun11)                                                  
      emax=etmpp(ntmp)*13.6+eth 
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,nlevp) 
      idest3=masterdata%idat1(np1i-1+nidti) 
      idest4=idest3+1 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        mllz=derivedpointers%npar(ndtmp) 
        iltmp=0 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(derivedpointers%npar(ndtmp).eq.mllz))                 
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpri.gt.1) then 
             write (lun11,*)nidt2,iltmp,ndtmp 
             write (lun11,*)np1r2,np1i2,np1k2,mlm 
             call dprinto(ltyp2,lrtyp2,lcon2,                        &
     &          nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,                    &
     &          lun11)                 
             endif 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           enddo 
         ggup=masterdata%rdat1(np1r2+1) 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup                       
         endif 
      if (lpri.gt.1) write (lun11,*)'before phint53' 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      if (lpri.gt.1) then 
        write (lun11,*)'type 49 data:',masterdata%idat1(np1i),          &
     &    masterdata%idat1(np1i+nidt-1),t,xnx,                          &
     &    eth,gglo,ggup,swrat                                           
        call dprinto(ndesc,nrdesc,lcon,                              &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)  
        endif 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      tm=t*1.e4 
      bktm=bk*tm/ergsev 
      q2=2.07e-16*xnx*(tm**(-1.5)) 
      emltlv=leveltemp%rlev(2,nlev) 
      rs=q2/emltlv 
      ethion=leveltemp%rlev(1,nlev) 
      emltlv=leveltemp%rlev(2,idest1) 
      rnissel=emltlv*rs
      rnisseu=1.
      rnist=rnissel                                                     &
     &  *exp(-(max(0.d0,13.605692*etmpp(1)))/(0.861707)/t)              &
     &  /(1.e-37+rnisseu)                                                   
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'ett=',ett,etmpp(1),                             &
     &  ett+max(0.d0,13.605692*etmpp(1)),                               &
     &  rnist
      call phint53(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,         &
     &  ans5,ans6,abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,           &
     &  opakc,opakcont,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,       &
     &  lfast,lun11)                                                    
      if (lpri.gt.1) then 
        npr=nb1 
        write (lun11,*)'bautista threshold xsection:',                  &
     &      npr,ett,eth,masterdata%rdat1(np1r),sg(npr),ans2,swrat     
        endif 
      lpri=lprisv 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      go to 9000 
!                                                                       
   50 continue 
!     op line rad. rates                                                
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
!     nb check this out:  no bound-bound decays from continuum          
      if ((idest1.le.0).or.(idest1.ge.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.ge.nlev))                          &
     &      go to 9000                                                  
      elin=abs(masterdata%rdat1(np1r)) 
      if (elin.le.1.d-34) go to 9000 
      aij=masterdata%rdat1(np1r+2) 
      ans1=aij*(ptmp1+ptmp2) 
!      aij=min(aij,1.e+10)                                              
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         endif 
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000 
      flin=(1.d-16)*aij*ggup*elin*elin/((0.667274)*gglo) 
!                                                                       
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      a=masterdata%rdat1(np1r2+1) 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      ener=12398.41/elin 
      dele=ener*vtherm/3.e+10 
      elammu=elin*1.d-4 
      sigma=(0.02655)*flin*elin*(1.d-8)/vtherm 
      sigvtherm=sigma 
      jkkl=derivedpointers%nplini(ml) 
      if (jkkl.le.0) go to 9000 
      ml3=derivedpointers%nplin(jkkl) 
      if (ml3.le.0) go to 9000 
      mlm=ml3
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                              &
     &  0,lun11)                                                  
      elin=abs(masterdata%rdat1(np1r)) 
      ener=12398.41/abs(elin) 
      epiip=ener
      nb1=nbinc(ener,epi,ncn2) 
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10                          &
     &      *flinabs(ptmp1)                                            
!     &      *2.
!     NB this is a fudge to test comparison with prismspect which assumes
!     radiation isotropic in the half space.
!      bremtmpp=bremsa(nb1)/(12.56) 
!      bbnurjp=(min(2.d+4,epiip))**3*(1.571e+22)*2. 
!      ans2=ans1*(bremtmpp/bbnurjp)*gglo/ggup
!                                                                       
!     nb turning off rex                                                   
!      ans2=0. 
      ans2=ans2*max(0.,1.d0-cfrac)
!                                                                       
      if (elin.gt.0.99e+9) then 
         ans2=0. 
         sigvtherm=0. 
         endif 
!      ans1=ans1+ans2*ggup/(1.d-36+gglo)                                
!     note that now opakab does not have abundance in                   
      opakab=sigvtherm 
      lfasto=2 
!      lfasto=4                                                         
      delea=0. 
      lfnd=0 
      lpriu=lpri

!      if (lpri.ge.1) lpriu=3                                           
      call deleafnd(jkion,idest1,                                    &
     &   delea,lfnd,lpriu,lun11)                          
!                                                                       
      if (lfnd.eq.0) delea=masterdata%rdat1(np1r+2)*(4.136e-15) 
      ans4=ans1*ener*ergsev 
      ans3=ans2*ener*ergsev 
      rcem1=abund2*ans4*ptmp1/(1.d-48+ptmp1+ptmp2) 
      rcem2=abund2*ans4*ptmp2/(1.d-48+ptmp1+ptmp2) 
      opakb1=sigvtherm*abund1 
!     this test should prevent calculation when called 
!     from calc_rates_level      
!     since abund1 will be zero                                         
!      lpriu=lpri                                                       
      lpriu=0 
      if ((nrdesc.ne.9).and.(lfasto.le.4).and.(opakb1*delr.gt.1.d-8))   &
     & call linopac(lpriu,lun11,opakb1,                                 &
     &               rcem1,rcem2,elin,vturb,t,a,delea,epi,ncn2,         &
     &               opakc,rccemis,lfasto)              
9801  format (1x,'in ucalc, ind=50:',3i8,10(1pe11.3),5i8,3(1pe11.3), &
     &              i8,3(1pe11.3),i8)
      if (lpri.gt.1)                                                    &
     &  write (lun11,9801)                                              &
     &  ml,nilin,nelin,elin,flin,masterdata%rdat1(np1r+2),gglo,ggup,a,  &
     &  vtherm,vturb,                                                   &
     &  ans1,ans2,idest1,idest2,idest3,idest4,nlev,sigvtherm,           &
     &  bremtmpp,bbnurjp,nb1,abund1,abund2,delea,lfnd                        
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      anstmp=ans1
      ans1=ans2
      ans2=anstmp
!
!     nb testing with suppressed decay!!!
!      ans2=ans2/1.e+8
!      ans1=ans1/1.e+8
!
      if (nrdesc.ne.9) go to 9000 
!                                                                       
!       special for 2 photon                                            
        ansar2=0. 
        em2ph=aij 
        lskp=1 
        emax=ener 
        nbmx=nbinc(emax,epi,ncn2) 
        if (lpri.gt.1)                                                  &
     &  write (lun11,*)'in ucalc, ind=50:',                          &
     &  ml,nilin,nelin,elin,flin,masterdata%rdat1(np1r),gglo,ggup,a,    &
     &  vtherm,ans2,nbmx                                                     
        rcemsum=0. 
        lfastl=0 
        ll=2 
        do while (ll.le.nbmx) 
          ansar2o=ansar2 
          ansar2=epi(ll)*epi(ll)*max(0.d0,(epi(nbmx)-epi(ll))) 
          rcemsum=rcemsum+(ansar2+ansar2o)                              &
     &                   *(epi(ll)-epi(ll-lskp))/2.                     
          call enxt(epi(1),nb1,0,epi,ncn2,t,lfastl,lun11,            &
     &                  ll,lskp,nphint,lrcalc)                          
          ll=ll+lskp 
          enddo 
        rctmp1=0. 
        rctmp2=0. 
        ll=2 
        do while (ll.le.nbmx) 
          ansar2=epi(ll)*epi(ll)*max(0.d0,(epi(nbmx)-epi(ll))) 
          ansar2=ansar2*em2ph*emax/(1.d-24+rcemsum) 
          rctmp1=abund2*ansar2*ptmp1/12.56 
          rctmp2=abund2*ansar2*ptmp2/12.56 
          rccemis(1,ll)=rccemis(1,ll)+rctmp1 
          rccemis(2,ll)=rccemis(2,ll)+rctmp2 
          call enxt(epi(1),nb1,0,epi,ncn2,t,lfastl,lun11,            &
     &                  ll,nskp,nphint,lrcalc)                          
          ll=ll+nskp 
          enddo 
      go to 9000 
!                                                                       
   51 continue 
!     line rates, col, burgess and tully from manuel                    
      idest1=masterdata%idat1(np1i+2) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         eeup=leveltemp%rlev(1,idest1) 
         eelo=leveltemp%rlev(1,idest2) 
         endif 
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      eijry=masterdata%rdat1(np1r) 
      eij=eijry*13.605692 
      elin=12398.41/eij 
      hij=elin*1.d-8 
!      if (lpri.gt.0)                                                   
!     $ write (lun11,*)'type 51 data:',elin                             
      if (elin.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)elin,ekt,delt                                     
!      if (delt.gt.50.) go to 9000                                      
      c=masterdata%rdat1(np1r+1) 
      p1=masterdata%rdat1(np1r+2) 
      p2=masterdata%rdat1(np1r-1+4) 
      p3=masterdata%rdat1(np1r-1+5) 
      p4=masterdata%rdat1(np1r-1+6) 
      p5=masterdata%rdat1(np1r-1+7) 
      tk=t*1.e+4 
!      tk=max(tk,(1.e+4)*12398.54/elin/(0.861707)/50.)                  
      tk=max(tk,2.8777d+6/elin) 
      ik=masterdata%idat1(np1i) 
      cijpp=upsil(ik,eijry,c,p1,p2,p3,p4,p5,tk) 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      cji=(8.626e-8)*cijpp/tsq                                          &
     &      /ggup                                                       
      exptmp=expo(-delt) 
      cij=cji*ggup*exptmp/gglo 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'ltyp=51',c,p1,p2,p3,p4,p5,ik,                    &
     &      eij,idest1,idest2,cij,cji,xnx,cijpp                         
      ans1=cij*xnx 
      ans2=cji*xnx 
      ans6=ans1*eij*ergsev
      ans5=ans2*eij*ergsev
      elin=0. 
      go to 9000 
!                                                                       
   52 continue 
!     same as 59 but rate type 7                                        
      go to 59 
!                                                                       
   53 continue 
  533 continue 
!     op pi xsections                                                   
      lprisv=lpri 
!      if (lpri.ge.1) lpri=2                                            
!     these are the initial and final levels and indeces                
!     notice that these are relative to the current ion                 
!     (not relative to the element as a whole)                          
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=nlevp+masterdata%idat1(np1i-1+nidt-3)-1 
      if (indonly.eq.1) return
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2,nlevp,ml 
      if (lpri.gt.1) write (lun11,*)'bremsa=',bremsa(1),bremsa(10),     &
     &    bremsa(100)
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000 
      if (ml.le.0) go to 9000 
      eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
      eexc=leveltemp%rlev(1,idest1) 
      ett=eth 
      nilin=derivedpointers%npar(ml) 
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml 
      if (nilin.le.0) go to 9000 
      ntmp=nrdt/2 
!      ett=ett+max(0.,13.605692*etmpp(1))                               
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1) 
      if (ett.le.0.) go to 9000 
      nb1=nbinc(ett,epi,ncn2) 
      xkt=ett/(0.861707*t) 
      r19=rr/1.e+19 
      mlm=nilin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      emax=etmpp(ntmp)*13.6+eth 
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,nlevp) 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        if (ndtmp.le.0) go to 9000 
        mllz=derivedpointers%npar(ndtmp) 
        iltmp=0 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(derivedpointers%npar(ndtmp).eq.mllz))                 
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           if (ndtmp.le.0) go to 9000 
           enddo 
!        NB fix to excited level PI and rec                             
         ett=ett+masterdata%rdat1(np1r2) 
         eth=ett 
         ggup=masterdata%rdat1(np1r2+1) 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett                   
         endif 
      sscal=1. 
      do ml2=1,ntmp 
        etmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2-1) 
        stmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2)*1.d-18*sscal 
        stmpp(ml2)=max(stmpp(ml2),0.d0) 
        if (lpri.gt.1) write (lun11,9819)ml2,etmpp(ml2),stmpp(ml2) 
 9819   format (1x,i6,2(1pe11.3)) 
        enddo 
      optst=abund1*stmpp(1) 
!      if ((optst.lt.opcrit).and.(lfast.eq.2)) go to 9000               
      ntmp2=nptmpdim 
!     nb includes extrapolation                                         
!     this is dangerous.  It does the right thing for ground-ground,    
!       but some cross sections should not be extrapolated.             
!      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11)        
      if (lpri.gt.1) write (lun11,*)'before phint53',eexc,eth,lfast 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      if (lpri.gt.1) then 
        write (lun11,*)'type 53 data:',masterdata%idat1(np1i),          &
     &    masterdata%idat1(np1i+nidt-1),t,xnx,                          &
     &    eth,gglo,ggup,swrat                                           
        call dprinto(ndesc,nrdesc,lcon,                              &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)  
        endif 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      tm=t*1.e4 
      bktm=bk*tm/ergsev 
      q2=2.07e-16*xnx*(tm**(-1.5)) 
      emltlv=leveltemp%rlev(2,nlev) 
      rs=q2/emltlv 
      ethion=leveltemp%rlev(1,nlev) 
      emltlv=leveltemp%rlev(2,idest1) 
      rnissel=emltlv*rs
      rnisseu=1.
      rnist=rnissel                                                     &
     &  *exp(-(max(0.d0,13.605692*etmpp(1)))/(0.861707)/t)              &
     &  /(1.e-37+rnisseu)                                                   
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'ett=',ett,etmpp(1),                             &
     &  ett+max(0.d0,13.605692*etmpp(1)),                               &
     &  rnist                                
      call phint53(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,         &
     &  ans5,ans6,abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,           &
     &  opakc,opakcont,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,       &
     &  lfast,lun11)                                                    
!                                                                       
      if (lpri.gt.1) then 
        npr=nb1 
        write (lun11,*)'bautista threshold xsection:',                  &
     &         npr,ett,eth,masterdata%rdat1(np1r),sg(npr),ans2,swrat         
        endif 
      if (lpri.gt.1) then 
        temp=t*1.e+4 
        do ml2=1,ntmp 
          etmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2-1) 
          stmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2) 
          enddo 
        lprim=0 
        call milne(temp,ntmp,etmpp,stmpp,ett/13.6,alphamilne,        &
     &     lun11,lprim)                                                 
        alphamilne=alphamilne*xnx 
        amilnerr=(log10(alphamilne/max(1.d-48,ans2))) 
        if ((abs(amilnerr).gt.0.05)                                     &
     &    .and.((alphamilne.gt.1.d-28).or.(ans2.gt.1.d-28))             &
     &    .and.(lfast.gt.1))                                            &
     &     write (lun11,*)'milne error',alphamilne,ans2,amilnerr        
        endif 
      lpri=lprisv 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      go to 9000 
!                                                                       
   54 continue 
!     h-like cij, bautista (hlike ion)                                  
      idest1=masterdata%idat1(np1i-1+nidt-3) 
      idest2=masterdata%idat1(np1i+nidt-3) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      ans3=0. 
      ans4=0. 
      lprisv=lpri 
!      if (lpri.ge.1) lpri=2                                            
      if (lpri.gt.1) write (lun11,*)'type 54 data:',                    &
     &  masterdata%idat1(np1i-1+nidt-3),masterdata%idat1(np1i+nidt-3)          
      if (leveltemp%rlev(1,idest2).lt.leveltemp%rlev(1,idest1)) then 
        itmp=idest2 
        idest2=idest1 
        idest1=itmp 
        endif 
      if (indonly.eq.1) return
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      dele=abs(eeup-eelo)
      elin=12398.41/abs(eeup-eelo+1.d-24) 
      hij=elin*1.d-8 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
!      if (delt.gt.50.) go to 9000                                      
      ni=leveltemp%ilev(1,idest2) 
      li=leveltemp%ilev(3,idest2) 
      nf=leveltemp%ilev(1,idest1) 
      lf=leveltemp%ilev(3,idest1) 
      if (lpri.gt.1) write (lun11,*)                                    &
     &  eeup,eelo,elin,ni,li,nf,lf                                      
      if (ni.eq.nf) go to 9000 
      if (ni.lt.nf) then 
        ntmp=ni 
        ni=nf 
        nf=ntmp 
        endif 
      iq=masterdata%idat1(np1i+nidt-2) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'before anl1:',ni,nf,li,lf,iq,idest1,idest2,   &
     &  eelo,eeup,masterdata%idat1(np1i-1+nidt-3),                      &
     &  masterdata%idat1(np1i+nidt-3)               
      call anl1(ni,nf,lf,iq,alm,alp,lpri,lun11) 
      ans1=alp 
      if (li.lt.lf) ans1=alm 
!     changing order to agree with type 50
      anstmp=ans2
      ans2=ans1
      ans1=anstmp
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans6
      ans6=-ans5
      ans5=-anstmp
      lpri=lprisv 
      go to 9000 
!                                                                       
   55 continue 
!      hydrogenic pi xsections, bautista format                         
      lprisv=lpri 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc, ind=55:',                          &
     &      (masterdata%idat1(np1i-1+mm),mm=1,5)    
      idest1=masterdata%idat1(np1i+nidt-2) 
      ett=leveltemp%rlev(1,nlevp)-leveltemp%rlev(1,idest1) 
      idest2=nlevp 
      if (indonly.eq.1) return
      if (ett.le.1.d-5) go to 9000 
      eth=ett 
      nb1=nbinc(eth,epi,ncn2) 
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,nlevp) 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      ekt=t*(0.861707) 
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in ucalc, ind=55:',                           &
     &   ml,nilin,nelin                                                 
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000 
      mlm=nilin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      nistage=masterdata%idat1(np1i2) 
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      nzel=masterdata%idat1(np1i2) 
      zz=float(nzel-nistage) 
      sgth=(6.3e-18)/zz/zz 
      ll=nb1 
      lfastl=1 
      do while (ll.le.nphint) 
        epii=epi(ll) 
        sg(ll)=sgth*(epii/ett)**(-3) 
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,              &
     &                  ll,nskp,nphint,lrcalc)                          
        ll=ll+nskp 
        enddo 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)    
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      lpri=lprisv 
      go to 9000 
!                                                                       
   56 continue 
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest2.le.0)                                &
     &   .or.(idest1.gt.nlev).or.(idest2.gt.nlev)) go to 9000      
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &      .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      lprisv=lpri 
!      if (lpri.ge.1) lpri=2                                            
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      dele=abs(eeup-eelo) 
      if (dele.le.1.d-16) go to 9000 
      ntmp=nrdt/2 
      if (ntmp.eq.1) then
        cijpp=masterdata%rdat1(np1r+1) 
        else
        do kl=1,ntmp 
          ttmp(kl)=masterdata%rdat1(np1r-1+kl) 
          enddo 
        tfnd=log10(t*1.e+4) 
        jlo=0 
        call hunt3(ttmp,ntmp,tfnd,jlo,0,lun11) 
        jlo=min(jlo,ntmp-1) 
        nind=ntmp+jlo 
        cijpp=(masterdata%rdat1(np1r-1+nind+1)                          &
     &   -max(1.d-48,masterdata%rdat1(np1r-1+nind)))                    &
     &   *(tfnd-ttmp(jlo))/(ttmp(jlo+1)-ttmp(jlo)+1.d-24)               &
     &     +max(1.d-48,masterdata%rdat1(np1r-1+nind))                          
        endif
      if (lpri.gt.1) write (lun11,*)'type 56:',                         &
     &  idest1,idest2,ggup,gglo,dele,jlo,nind,                          &
     &  masterdata%rdat1(np1r-1+nind),masterdata%rdat1(np1r-1+nind+1),  &
     &  tfnd,ttmp(jlo+1),ttmp(jlo)                                      
!                                                                       
!     NB a fudge for Fe XXIV q line                                     
!      if ((jkion.eq.349).and.(idest1.eq.1).and.(idest2.ge.26)) then    
!         cijpp=cijpp*166.                                              
!         endif                                                         
!                                                                       
      cijpp=max(0.d0,cijpp) 
      ekt=0.861707*t 
      delt=dele/ekt 
      cij=0. 
      exptmp=expo(-delt) 
      if (lpri.gt.1) write (lun11,*)'type 56:',                         &
     &  idest1,idest2,ggup,gglo,dele                                    
      cij=(8.626e-8)*cijpp*exptmp/tsq/gglo 
      cji=(8.626e-8)*cijpp/tsq/ggup 
      ans1=cij*xnx 
      ans2=cji*xnx 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
      if (lpri.gt.1) write (lun11,*)'type 56 data:',                    &
     &  idest1,idest2,dele,cijpp,delt,exptmp,cij,cji,                   &
     &  nind,gglo,ggup,tfnd,ntmp,jlo,ntmp,ttmp(jlo)                     
      lpri=lprisv 
      go to 9000 
!                                                                       
   57 continue 
!     same as  65 (?)                                                   
!     effective charge to be used in coll. ion.                         
      lprisv=lpri 
      lpri=0 
!      if (lprisv.ge.1) lpri=2                                          
      tz=t*1.e+4 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=nlevp 
      if (indonly.eq.1) return
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in ucalc at 57:',idest1,                      &
     &  masterdata%idat1(np1i),masterdata%rdat1(np1r)  
      if ((masterdata%idat1(np1i).le.0)                                 &
     &      .or.(idest1.le.1).or.(idest1.gt.nlevp))                     &
     &        go to 9000                                                
      i57=masterdata%idat1(np1i) 
      eth=max(0.d0,leveltemp%rlev(1,nlevp)-leveltemp%rlev(1,idest1)) 
      ekt=0.861707*t 
!      tz=max(tz,(1.e+4)*eth/(0.861707)/50.)                            
!      tz=max(tz,2.320975e+02*eth)                                      
      e1=leveltemp%rlev(1,idest1) 
      ep=leveltemp%rlev(4,idest1)-e1
      if (ep.le.0.) go to 9000 
      call calt57(tz,xnx,e1,ep,i57,cion,crec,lun11,lpri) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'ltype=57:',cion,crec,gglo,ggup,nlevp,idest1,rinf,&
     &  eth,ekt,ans1,ans2                                               
!                                                                       
!     trying a fudge to test cloudy's ci                                
!      if (lprisv.ge.1) write (lun11,*)'fudging ci for test'            
!      if (i57.eq.2) cion=cion*2.                                       
!      if (i57.eq.2) crec=crec*2.                                       
!      if (i57.eq.3) cion=cion*5.                                       
!      if (i57.eq.3) crec=crec*5.                                       
!      if (i57.eq.4) cion=cion*15.                                      
!      if (i57.eq.4) crec=crec*15.                                      
!                                                                       
      ans1=cion*xnx 
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,idest1) 
!     note that rinf has exponential removed                            
      rinf=gglo/(1.d-48+ggup) 
      ans2=crec*rinf*xnx*xnx 
      ans6=-ans1*eth*ergsev
      ans5=-ans2*eth*ergsev
!     set to zero for ground state because we have more accurate rates  
!     for these levels: types 95 or 25                                  
      if (idest1.eq.1) then 
        ans1=0. 
        ans2=0. 
        endif 
      go to 9000 
!                                                                       
   58 continue 
!      bautista cascade rates. defunct.                                 
      go to 9000 
!                                                                       
   59 continue 
      lprisv=lpri 
      lpril=lpri 
!      if (lpri.ge.1) lpril=2                                           
      if (lpril.gt.1) then
         write (lun11,*)'ltyp=59',ml,derivedpointers%npar(ml) 
         write (lun11,*)(masterdata%rdat1(np1r-1+jj),jj=1,nrdt) 
         write (lun11,*)(masterdata%idat1(np1i-1+jj),jj=1,nidt),nidt 
         write (lun11,*)(masterdata%kdat1(np1k-1+jj),jj=1,nkdt)
         endif
      if (ml.le.0) go to 9000 
!                                                                       
!                                                                       
!     experiment with only vfky                                         
!      if (nrdt.le.6) go to 9000                                        
!                                                                       
      lfastl=1 
      nilin=derivedpointers%npar(ml) 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=masterdata%idat1(np1i+nidt-3) 
!     why was this statement here?                                      
      if (idest4.gt.idest3+1) go to 9000 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=nlevp+masterdata%idat1(np1i-1+nidt-3)-1 
      idest2=max(idest2,1) 
      if (indonly.eq.1) return
!      nb must uncomment these if func2a is called                      
!      if (nrdesc.eq.7) then                                            
!        idest2=nlevp                                                   
!        endif                                                          
      if ((nilin.le.0).or.(nilin.gt.np2)) go to 9000 
      mlm=nilin
      call drd(ltyp2,lrtyp2,lcon2,                                      &
     &  nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                        &
     &  0,lun11)                                                  
      if (lpril.gt.1)                                                   &
     & write (lun11,*)ml,nilin,masterdata%rdat1(np1r),idest1,           &
     & masterdata%rdat1(np1r2),nlevp    
      ett=masterdata%rdat1(np1r2) 
      if ((idest1.gt.nlevp).or.(idest1.le.0)) go to 9000 
      if (ml.le.0) go to 9000 
      if (ett.le.0.) go to 9000 
      nb1=nbinc(ett,epi,ncn2) 
      numcon2=max(2,ncn2/50) 
      nphint=ncn2-numcon2 
      nphint=max(nphint,nb1+1) 
      if (nb1.ge.nphint-1) go to 9000 
      ett=masterdata%rdat1(np1r) 
      nb1=nbinc(ett,epi,ncn2) 
      gglo=leveltemp%rlev(2,1) 
      ggup=leveltemp%rlev(2,nlevp) 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpril.gt.1)                                                 &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpril.gt.1)                                                 &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        if (ndtmp.le.0) go to 9000 
        mllz=derivedpointers%npar(ndtmp) 
        iltmp=0 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(derivedpointers%npar(ndtmp).eq.mllz))                
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           if (ndtmp.le.0) go to 9000 
           enddo 
!        NB fix to excited level PI and rec                             
         ett=ett+masterdata%rdat1(np1r2) 
         eth=ett 
         ggup=masterdata%rdat1(np1r2+1) 
         if (lpril.gt.1)                                                &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett                   
         endif 
      if (lpril.gt.1) write (lun11,*)nlevp,ggup 
      if (ggup.le.1.d-24) then 
        if (lpril.gt.1) write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      ett=masterdata%rdat1(np1r) 
      nb1=nbinc(ett,epi,ncn2) 
      if (lpril.gt.1)                                                   &
     & write (lun11,*)'ett=',ett,nb1,nphint,swrat,gglo,ggup             
      if (nb1.ge.(nphint-1)) go to 9000 
      if ((bremsint(nb1).lt.1.d-20).and.(lpril.gt.1))                   &
     &    write (lun11,*)'skipping 59',                                 &
     &         nb1,bremsint(nb1)                                        
      if (bremsint(nb1).lt.1.d-20) go to 9000 
      if (nrdt.eq.9) then 
          ett=masterdata%rdat1(np1r) 
          emax=masterdata%rdat1(np1r+1) 
          e0=masterdata%rdat1(np1r+2) 
          s0=masterdata%rdat1(np1r-1+4) 
          ya=masterdata%rdat1(np1r-1+5) 
          pp=masterdata%rdat1(np1r-1+6) 
          yw=masterdata%rdat1(np1r-1+7) 
          y0=masterdata%rdat1(np1r-1+8) 
          y1=masterdata%rdat1(np1r-1+9) 
          l2=0 
        else 
          e0=masterdata%rdat1(np1r+1) 
          s0=masterdata%rdat1(np1r+2) 
          ya=masterdata%rdat1(np1r-1+4) 
          pp=masterdata%rdat1(np1r-1+5) 
          yw=masterdata%rdat1(np1r-1+6) 
          y0=0. 
          y1=0. 
          l2=masterdata%idat1(np1i+2) 
        endif 
      ywsq=yw*yw 
      qq=5.5+l2-pp/2. 
      if (lpril.gt.1) write (lun11,*)'qq=',                             &
     &   l2,qq,ya,ywsq,pp,yw,s0                                         
      ll=nb1 
      do while (ll.le.nphint) 
        epii=epi(ll) 
        xx=epii/e0-y0 
        if (nrdt.eq.9) then 
            yy=sqrt(xx*xx+y1*y1) 
          else 
            yy=xx 
          endif 
        yyqq=qq*log(max(1.d-48,yy)) 
        yyqq=exp(-min(60.d0,max(-60.d0,yyqq))) 
        term1=((xx-1.)*(xx-1.)+ywsq) 
        term2=yyqq 
        term3=(1.+sqrt(yy/ya))**(-pp) 
        ff=term1*term2*term3 
        sg(ll)=s0*ff*(1.d-18) 
        if (lpril.gt.1) write (lun11,*)ll,epii,sg(ll),                  &
     &    yy,yyqq,xx,term1,term2,term3,qq,ff                            
        call enxt(ett,nb1,0,epi,ncn2,t,lfastl,lun11,                 &
     &                  ll,nskp,nphint,lrcalc)                          
        ll=ll+nskp 
        enddo 
      ekt=t*(0.861707) 
      lprib=0 
      if (lpril.gt.1) lprib=lpril 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
      if (lpril.gt.1) then 
        npr=nb1 
        write (lun11,*)'verner threshold xsection:',                    &
     &         npr,ett,sg(npr),opakab                                   
        endif 
!     nb this turns off all recombination into excited levels           
!     for type 59...                                                    
      if ((nrdesc.eq.1).or.(idest1.gt.1)) then                         
        ans6=0.                                                        
        ans4=0.                                                        
        ans2=0.                                                        
        endif                                                          
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      lpri=lprisv 
      go to 9000 
!                                                                       
   60 continue 
!      go to 9000                                                       
!     calloway h-like coll. strength                                    
      lpril=0 
!      if (lpri.ge.1) lpril=2                                           
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &    .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      dele=abs(leveltemp%rlev(1,idest2)-leveltemp%rlev(1,idest1)) 
      if (dele.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      delt=dele/ekt 
      temp=t*1.e+4 
      temp=max(temp,0.02*dele*1.e+4/(0.861707)) 
      call calt6062(temp,nrdt,ndesc,np1r,np1i,cijpp) 
!      cijpp=cijpp/2./2.                                                
      cji=(8.626e-8)*cijpp/tsq/(1.d-16+ggup) 
      exptmp=expo(-delt) 
      cij=cji*ggup*exptmp/(1.d-16+gglo) 
      ans1=cij*xnx 
      ans2=cji*xnx 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
      if (lpril.gt.1) then 
        write (lun11,*)'ltyp=60',idest1,idest2,temp,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,                                  &
     &                          (masterdata%rdat1(np1r-1+mm),mm=1,8),jlo 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp,dele,delt 
        endif 
      go to 9000 
!                                                                       
   61 continue 
      go to 9000 
!                                                                       
   62 continue 
      go to 60 
!                                                                       
   63 continue 
!      if (lpri.gt.0) write (lun11,*) 'type 63 data not implemented'    
!      go to 9000                                                       
      lpril=0 
!      if (lpri.ge.1) lpril=2                                           
      idest1=masterdata%idat1(np1i-1+nidt-3) 
      idest2=masterdata%idat1(np1i+nidt-3) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      elin=12398.41/abs(eeup-eelo+1.d-24) 
      dele=abs(eeup-eelo)
      hij=elin*1.d-8 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      if (lpril.ne.0) write (lun11,*)'delt=',delt 
      if (delt.gt.50.) go to 9000 
      ni=leveltemp%ilev(1,idest1) 
      li=leveltemp%ilev(3,idest1) 
      nf=leveltemp%ilev(1,idest2) 
      lf=leveltemp%ilev(3,idest2) 
      sum=0. 
      iq=masterdata%idat1(np1i+nidt-2) 
      if (lpril.ne.0)                                                   &
     & write (lun11,*)'ltyp=63',idest1,idest2,ni,li,nf,lf               
      if (nf.eq.ni) then 
        if (lpril.ne.0)                                                 &
     &     write (lun11,*)'nf=ni',lf,li                                 
        if (abs(lf-li).eq.1) then 
          lff=min(lf,li) 
          lii=max(lf,li) 
                             ! mab                                      
          li1=max(1,lff) 
          do  nn=li1,ni-1 
            if (lpril.ne.0)                                             &
     &        write (lun11,*)'before anl1'                              
           if (lii.ge.1) then 
            if (lpril.ne.0)                                             &
     &        write (lun11,*)'li=1',ni,nn,lii-1,iq                      
            call anl1(ni,nn,lii-1,iq,alm,alp,lpril,lun11) 
!              write (lun11,*)'li=1',ni,nn,lii-1,iq,alp                 
            sum=sum+alp 
           endif 
           if (nn.gt.lii+1) then 
            if (lpril.ne.0)                                             &
     &        write (lun11,*)'nn=li+1',ni,nn,lii+1,iq                   
            call anl1(ni,nn,lii+1,iq,alm,alp,lpril,lun11) 
            sum=sum+alm 
           endif 
          enddo 
          if (lpril.ne.0)                                               &
     &     write (lun11,*)'after anl1',sum                              
          ecm=abs(leveltemp%rlev(1,idest1)                              &
     &             -leveltemp%rlev(1,idest2))*8059.9 
          ecm=0. 
          nnz=masterdata%idat1(np1i-1+4) 
          tbig=t*1.e+4 
          z1=1. 
          rm=1800. 
          il=0 
          psi=0.75/nnz/nnz*lii/(2*lii+1)*ni*ni*(ni*ni-lii*lii) 
          if (lpril.ne.0)                                               &
     &     write (lun11,*)'before amcrs',ecm,ni,lii,sum                 
          call amcrs(ni,lii,tbig,nnz,z1,rm,xnx,sum,ecm,psi,il,cn,       &
     &        lpril,lun11)                                              
          cno=cn 
          iz=masterdata%idat1(np1i-1+4) 
          if (lf.lt.li) then 
            ans1=cn 
            ans2=cn*leveltemp%rlev(2,idest1)/leveltemp%rlev(2,idest2) 
          else 
            ans2=cn 
            ans1=cn*leveltemp%rlev(2,idest2)/leveltemp%rlev(2,idest1) 
          endif 
          if (lpril.ne.0)                                               &
     &     write (lun11,*)'after amcrs',cn,iz,cno,ans1,ans2             
        endif 
      else 
        if (lpril.ne.0)                                                 &
     &     write (lun11,*)'nf.ne.ni'                                    
        aa1=0. 
        if (abs(lf-li).eq.1) then 
          sum=0. 
          nu=max(ni,nf) 
          nll=min(ni,nf) 
          do lff=0,nll-1 
            call anl1(nu,nll,lff,iq,alm,alp,lpri,lun11) 
             sum=sum+alp*(2*lff+3) 
             if (lff.gt.0)  then 
              sum=sum+alm*(2*lff-1) 
             endif 
             if (lff.eq.lf .and. li.gt.lf) aa1=alp 
             if (lff.eq.lf .and. li.lt.lf) aa1=alm 
             if (lpril.ne.0) write (lun11,*)'after anl1',            &
     &           lff,li,lf,sum,alp,alm,aa1                              
          enddo 
          if (lpril.ne.0)                                               &
     &     write (lun11,*)'after anl1',sum,alp,alm,aa1                  
          nnz=masterdata%idat1(np1i-1+4) 
          tbig=t*1.e+4 
          call erc(nll,nu,tbig,nnz,se,sd,sum,lun11,lpril) 
! ***** check if ans1 and ans2 are correct or inverted                  
          ans1=se*(2*lf+1)*aa1/sum 
          ans2=sd*(2*li+1)*aa1/sum 
          if ((nf.gt.ni).or.(lf.gt.li)) then 
           atmp=ans1 
           ans1=ans2 
           ans2=atmp 
          endif 
          if (lpril.ne.0)                                               &
     &     write (lun11,*)'after erc',se,sd,ans1,ans2                   
        endif 
      endif 
!                                                                       
      ans1=ans1*xnx 
      ans2=ans2*xnx 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
      go to 9000 
!                                                                       
   64 continue 
!     hydrogenic pi xsections, bautista format                          
      lprisv=lpri 
      idest1=masterdata%idat1(np1i+nidt-2) 
      if (indonly.eq.1) return
      ett=abs(leveltemp%rlev(1,nlevp)-leveltemp%rlev(1,idest1)) 
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc, ind=64:',                          &
     &      (masterdata%rdat1(np1r-1+mm),mm=1,5)    
      if (ett.le.1.d-5) go to 9000 
      zzz=float(masterdata%idat1(np1i+2)) 
      enn=float(masterdata%idat1(np1i)) 
      eth=ett 
      nb1=nbinc(eth,epi,ncn2) 
      gglo=leveltemp%rlev(2,idest1) 
      swrat=gglo 
      idest2=nlevp 
      ekt=t*(0.861707) 
      ll=nb1 
      lorb=masterdata%idat1(np1i+1) 
      ic=masterdata%idat1(np1i+2) 
      nq=masterdata%idat1(np1i) 
      mm=0 
      lfastl=1 
      do while (ll.le.nphint) 
        mm=mm+1 
        epii=epi(ll) 
        e=epii 
        eth=ett 
        erel=max(0.d0,(e-eth)/13.605692) 
        call hphotx(erel,ic,nq,xsec,lun11,lpri) 
        sg(ll)=xsec(lorb+1)*(1.d-18) 
        stmpp(mm)=xsec(lorb+1) 
        etmpp(mm)=erel 
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,              &
     &                  ll,nskp,nphint,lrcalc)                          
        ll=ll+nskp 
        enddo 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,ans5,ans6,             &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
      lprim=0 
      ntmp=ll-nb1 
      temp=t*1.e+4 
      ntmp=mm 
      call milne(temp,ntmp,etmpp,stmpp,eth/13.6,ans2,lun11,lprim) 
      ans2=ans2*swrat 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      lpri=lprisv 
      go to 9000 
!                                                                       
!                                                                       
   65 continue 
!     effective charge to be used in coll. ion.                         
      tz=t*1.e+4 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=nlevp 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,1) 
      eth=max(0.d0,leveltemp%rlev(1,nlevp)-leveltemp%rlev(1,idest1)) 
      ekt=0.861707*t 
!      if (eth/ekt.gt.50.) go to 9000                                   
      call szirco(masterdata%idat1(np1i),tz,masterdata%rdat1(np1r),  &
     &               cii) 
      ans1=cii*xnx 
!     note that rinf has exponential removed                            
      rinf=(2.08e-22)*gglo/ggup/t/tsq 
      ans2=ans1*rinf*expo(eth/ekt) 
      ans6=ans1*eth*ergsev
      ans5=ans2*eth*ergsev
      go to 9000 
!                                                                       
   66 continue 
!     Like type 69 but, data in fines tructure                          
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &      .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      elin=masterdata%rdat1(np1r) 
      if (elin.le.1.d-24) go to 9000 
      dele=elin
      elin=12398.41/elin 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
!      if (delt.gt.50.) go to 9000                                      
      hij=elin*1.d-8 
      temp=t*1.e+4 
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
      temp=max(temp,2.8777e+6/elin) 
      call calt66(temp,np1r,nrdt,gamma) 
      cijpp=gamma 
      cji=(8.626e-8)*cijpp/tsq/ggup 
        exptmp=expo(-delt) 
        cij=cji*ggup*exptmp/gglo 
      ans1=cij*xnx 
      ans2=cji*xnx 
      if (lpri.gt.1) then 
        write (lun11,*)'ltyp=66',idest1,idest2,elin,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,                                  &
     &   (masterdata%rdat1(np1r-1+mm),mm=1,8),nind,jlo 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp 
        endif 
      elin=0. 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
      go to 9000 
!                                                                       
   67 continue 
!     Effective collision strengths from Keenan et al.                  
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &      .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      elin=abs(masterdata%rdat1(np1r)) 
      hij=elin*1.d-8 
      if (elin.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      dele=elin
      delt=12398.41/elin/ekt 
      temp=t*1.e+4 
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
      temp=max(temp,2.8777e+6/elin) 
      call calt67(temp,np1r,gamma) 
      cijpp=gamma 
      cijpp=max(0.d0,cijpp) 
      cji=(8.626e-8)*cijpp/tsq/ggup 
        exptmp=expo(-delt) 
        cij=cji*ggup*exptmp/gglo 
      ans1=cij*xnx 
      ans2=cji*xnx 
      if (lpri.gt.1) then 
        write (lun11,*)'ltyp=69',idest1,idest2,elin,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,                                  &
     &                  (masterdata%rdat1(np1r-1+mm),mm=1,8),nind,jlo 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp 
        endif 
      elin=0. 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
      go to 9000 
!                                                                       
   68 continue 
!     coll. strength He-like ions by Zhang & Sampason                   
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &      .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      elin=12398.41/abs(eeup-eelo+1.d-24) 
      dele=abs(eeup-eelo+1.d-24) 
      hij=elin*1.d-8 
      if (elin.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      temp=t*1.e+4 
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
      temp=max(temp,2.8777e+6/elin) 
      if (lpri.gt.1) then 
        write (lun11,*)'ltyp=68',idest1,idest2,elin,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,                               &
     &        (masterdata%rdat1(np1r-1+mm),mm=1,8),nind,jlo 
        endif 
      call calt68(temp,np1r,np1i,gamma) 
      cijpp=gamma 
      cijpp=max(cijpp,0.d0) 
      cji=(8.626e-8)*cijpp/tsq/ggup 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
        exptmp=expo(-delt) 
        cij=cji*ggup*exptmp/gglo 
      ans1=cij*xnx 
      ans2=cji*xnx 
      if (lpri.gt.1) then 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp 
        endif 
      elin=0. 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
      go to 9000 
!                                                                       
   69 continue 
!     Kato & Nakazaki (1996) fit to Helike coll. strgt                  
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &      .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      dele=abs(eeup-eelo)
      elin=12398.41/abs(eeup-eelo+1.d-24) 
      hij=elin*1.d-8 
      if (elin.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      m=nrdt 
      temp=t*1.e+4 
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
!      temp=max(temp,2.8777e+6/elin)                                    
      call calt69(temp,m,np1r,gamma,lpri,lun11) 
      cijpp=gamma 
      cijpp=max(cijpp,0.d0) 
      cji=(8.626e-8)*cijpp/tsq/ggup 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
        exptmp=expo(-delt) 
        cij=cji*ggup*exptmp/gglo 
      if (lpri.gt.1) then 
        write (lun11,*)'ltyp=69',idest1,idest2,elin,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,                                 &
     &               (masterdata%rdat1(np1r-1+mm),mm=1,8),nind,jlo 
        endif 
      ans1=cij*xnx 
      ans2=cji*xnx 
      if (lpri.gt.1) then 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp 
        endif 
      elin=0. 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
      go to 9000 
!                                                                       
   99 continue 
!     Coefficients for phot x-section of suplevels                      
!      lfastl=lfast                                                     
      lfastl=3 
      temp=t*1.e+4 
      ans3=0. 
      ans4=0. 
      den=xpx 
      m=1000 
      lpric=0
      mlm=ml
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                              &
     &  0,lun11)                                                  
!      if (lpri.ge.1) lpric=2                                          
      mlion=derivedpointers%npar(ml) 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest1=min(idest1,nlev-1) 
      idest2=nlev+masterdata%idat1(np1i-1+nidt-3)-1 
      idest2=max(idest2,nlev) 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,nlevp) 
      ett=abs(leveltemp%rlev(1,idest1)-leveltemp%rlev(1,nlevp)) 
      if (lpric.gt.1)                                                   &
     & write (lun11,*)'rlev:',idest1,nlevp,                             &
     &    leveltemp%rlev(1,idest1),leveltemp%rlev(1,nlevp) 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpric.gt.1)                                                 &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpric.gt.1)                                                 &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        mllz=derivedpointers%npar(ndtmp) 
        iltmp=0 
        nptmp=mllz 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(nptmp.eq.mllz))                                       
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpric.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           nptmp=0 
           if (ndtmp.ne.0) nptmp=derivedpointers%npar(ndtmp) 
           enddo 
         ggup=masterdata%rdat1(np1r2+1) 
         ett=abs(leveltemp%rlev(1,idest1)+masterdata%rdat1(np1r2)) 
         endif 
       if (lpric.gt.1)                                                  &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett                   
      xkt=ett/(0.861707*t) 
      nb1=nbinc(ett,epi,ncn2) 
      mlm=mlion
      call drd(ltyp2,lrtyp2,lcon2,                                      &
     &  nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                        &
     &  0,lun11)                                                  
      ist=masterdata%idat1(np1i2) 
      ic=ist 
      eth=ett 
      gglo=leveltemp%rlev(2,idest1) 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      if (lpric.gt.1) then 
         write (lun11,*)'type 99 data:',masterdata%idat1(np1i),         &
     &       masterdata%idat1(np1i+nidt-1),t,xnx,eth,gglo,ggup,swrat
         call dprinto(ndesc,nrdesc,lcon,                                &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)  
        endif 
      ettry=ett/13.6 
!     nb eliminating high density for H
!      if (jkion.eq.1) den=min(den,1.d+8)
      m=nrdt
      call calt99(temp,den,ettry,ic,m,np1r,np1i,                        &
     &             ntmp,etmpp,stmpp,rec,al,lun11,lpric,ierr)                 
      if ((ierr.ne.0).and.(lpric.ne.0)) write (lun11,*)'calt99 error'
      if (lpric.gt.1) write (lun11,*)'after  calt70:',rec,stmpp(1) 
      crit53=0.01 
      do mm=1,ntmp 
        stmpp(mm)=stmpp(mm)*1.d-18 
        stmpp(mm)=max(stmpp(mm),0.d0) 
        enddo 
      call phint53hunt(stmpp,etmpp,ntmp,ett,ans1,ans2d,ans3d,ans4s,     &
     & lpric,epi,ncn2,bremsa,t,swrat,xnx,crit53,lfastl,lun11)           
      if (ans2d.le.1.d-48) then 
        ans1=0. 
        ans2=0. 
        go to 9000 
        endif 
      scale=rec*xnx/ans2d 
      ans1=ans1*scale 
!     does the swrat not belong?                                        
!      ans2=rec*xnx*swrat                                               
      ans2=rec*xnx 
!      ans2=ans2d                                                       
      tm=t*1.e4 
      q2=2.07e-16*xnx*(tm**(-1.5)) 
      rs=q2/swrat 
      ans1o=ans1 
!      ans1=min(ans1,ans2/rs)                                           
      if (lpric.ge.2)                                                   &
     & write (lun11,*)'type 99 limit:',ans2,rs,swrat,                   &
     &   xnx,tm,q2,ans1o,ans1,scale,rec                                 
!
!                                                                       
!     nb testing superlevel phot.                                       
!      ans1=0.                                                          
!                                                                       
      go to 9000 
!                                                                       
   71 continue 
!     Transition rates from superlevel to spect. lvls                   
      idest1=masterdata%idat1(np1i-1+nidt-3) 
      idest2=masterdata%idat1(np1i+nidt-3) 
      if (indonly.eq.1) return
      temp=t*1.e+4 
      lpril=0 
      den=xpx 
      m=1000 
      if (lpril.ne.0)                                                   &
     &  write (lun11,*)'before calt71:',masterdata%rdat1(np1r),      &
     &    masterdata%rdat1(np1r+1),masterdata%rdat1(np1r+2)                   
      call calt71(temp,den,ic,m,np1r,np1i,                           &
     &            wav,aij,lun11,lpril)                                  
      if ((idest1.le.0).or.(idest1.gt.nlev).or.                         &
     &   (idest2.le.0).or.(idest2.gt.nlev)) go to 9000                  
      if (lpril.ne.0)                                                   &
     & write (lun11,*)idest1,idest2,aij,wav,ml                          
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      elin=wav 
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      flin=(1.d-16)*aij*ggup*elin*elin/((0.667274)*gglo) 
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      a=masterdata%rdat1(np1r2+1) 
      elammu=elin*1.d-4 
      ans1=aij*(ptmp1+ptmp2) 
!     special fudge for ca i and ca ii                                  
      if ((masterdata%idat1(np1i-1+6).eq.96)                            &
     &   .or.(masterdata%idat1(np1i-1+6).eq.97))                        &
     &   ans1=min(ans1,1.d+10)                                            
!                                                                       
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      sigma=(0.02655)*flin*elin*(1.d-8)/vtherm 
      sigvtherm=sigma 
      ener=12398.41/abs(elin) 
      nb1=nbinc(ener,epi,ncn2) 
      ans2=0. 
!      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10                         
      if (elin.gt.0.99e+9) then 
         ans2=0. 
         sigvtherm=0. 
         endif 
      ans1=ans1+ans2*ggup/(1.d-48+gglo) 
      opakab=sigvtherm 
      ans3=ans2*ener*ergsev 
      ans4=ans1*ener*ergsev 
      if (elin.gt.0.1) then 
        dele=12398.41/(elin+1.d-24) 
        ans4=ans1*dele*(1.602197e-12) 
        endif 
!      ans4=0.                                                          
      if (lpril.ne.0)                                                   &
     & write (lun11,*)' ',vtherm,ans2,ans4,flin                         
!     changing order to allow universal assignment in calc_level_rates_level
      ans2=ans1
      ans1=0.
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      go to 9000 
!                                                                       
   72 continue 
!     Autoinization rates (in s^-1) for satellite lvls                  
      lpril=0 
      idest1=masterdata%idat1(np1i-1+nidt-3) 
      idest2=masterdata%idat1(np1i+nidt-3) 
      if (indonly.eq.1) return
      temp=t*1.e+4 
      call calt72(temp,np1r,nrdt,rate,lun11,lpril) 
      ans2=rate*xnx 
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,1) 
!     note that rinf has exponential removed                            
      rinf=(2.08e-22)*gglo/ggup/t/tsq 
      dele=masterdata%rdat1(np1r+1) 
      ans1=rate*xnx*rinf*xnx*expo(dele/temp) 
      go to 9000 
!                                                                       
   75 continue 
!     Autoinization rates (in s^-1) for satellite lvls                  
!        now including final ion stage                                  
      lpril=0 
!      if (lpri.gt.0) lpril=2                                           
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=masterdata%idat1(np1i+nidt-3) 
      idest2=masterdata%idat1(np1i+nidt-2)+nlev-1 
      idest1=masterdata%idat1(np1i-1+nidt-3) 
      idest1=max(idest1,1) 
      idest2=max(idest2,1) 
      if (indonly.eq.1) return
      temp=t*1.e+4 
      call calt72(temp,np1r,nrdt,rate,lun11,lpril) 
      ans2=rate*xnx 
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,idest1) 
!     note that rinf has exponential removed                            
      tsq=sqrt(t) 
      rinf=(2.08e-22)*gglo/ggup/t/tsq 
      ee=masterdata%rdat1(np1r-1+2) 
      ekt=0.861707*t 
      tt=ekt/ee 
      ans1=0. 
!      ans1=ans2*rinf*xnx/expo(-1./tt)                                  
      if (lpril.ne.0)                                                   &
     & write (lun11,975)                                                &
     &  (masterdata%rdat1(np1r-1+lk),lk=1,3),                           &
     &  (masterdata%idat1(np1i-1+lk),lk=1,4),rate,ans2,                 &
     &  ggup,gglo,rinf,t,ekt,ee,tt,ans1                                 
  975 format (1x,'type 75 calc:',3(1pe11.3),4i4,10(1pe11.3)) 
      go to 9000 
!                                                                       
   73 continue 
!     Fit to coll. strengths satellite lvls Helike ion                  
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &      .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      elin=abs(masterdata%rdat1(np1r)) 
      hij=elin*1.d-8 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      if (elin.le.1.d-24) go to 9000 
      m=1000 
      temp=t*1.e+4 
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
      temp=max(temp,2.8777e+6/elin) 
      crate=0. 
      call calt73(temp,np1r,np1i,crate) 
!      write (lun11,*)'type 73 calc:',                                  
!     $  (rdat1(np1r-1+lk),lk=1,7),(idat1(np1i-1+lk),lk=1,4),crate,     
!     $  gglo,ggup                                                      
      cijpp=crate/gglo 
      cijpp=max(cijpp,0.d0) 
      cji=(8.626e-8)*cijpp/tsq/ggup 
        exptmp=expo(-delt) 
       cij=cji*ggup*exptmp/gglo 
      ans1=cij*xnx 
      ans2=cji*xnx 
      ans6=ans1*elin*ergsev
      ans5=ans2*elin*ergsev
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans6
      ans6=-ans5
      ans5=-anstmp
      if (lpri.gt.1) then 
        write (lun11,*)'ltyp=69',idest1,idest2,elin,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,                                 &
     &                 (masterdata%rdat1(np1r-1+mm),mm=1,8),nind,jlo 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp 
        endif 
      elin=0. 
      go to 9000 
!                                                                       
   74 continue 
!     Delta functions to add to phot. x-sections  DR                    
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=nlevp 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
      if (indonly.eq.1) return
      temp=t*1.e+4 
      den=xpx 
      m=1000 
      rec=0. 
      lprisv=lpri 
!      if (lpri.ge.1) lpri=2                                            
      if (lpri.gt.1) write (lun11,*)'type 74 data:',den,temp,           &
     & (masterdata%rdat1(np1r-1+mm),mm=1,nrdt),                         &
     &  (masterdata%idat1(np1i-1+mm),mm=1,nidt)        
      call calt74(temp,ncn2,epi,bremsa,nrdt,np1r,rate,               &
     &       alpha)                                                     
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,idest2) 
      if (lpri.gt.1) write (lun11,*)'returning from calt74:',        &
     &  rate,alpha,idest1,idest2,gglo,ggup                              
      ans1=rate 
      alpha=alpha*gglo/ggup 
      ans2=alpha 
      lpri=lprisv 
      go to 9000 
!                                                                       
   81 continue 
!     bhatia Fe XIX                                                     
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      if (leveltemp%rlev(1,masterdata%idat1(np1i+1))                    &
     &      .lt.leveltemp%rlev(1,masterdata%idat1(np1i))) then 
        idest2=masterdata%idat1(np1i) 
        idest1=masterdata%idat1(np1i+1) 
        endif 
      if (indonly.eq.1) return
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      elin=12398.41/abs(eeup-eelo+1.d-24) 
      dele=abs(eeup-eelo)
      hij=elin*1.d-8 
      if (elin.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      m=nrdt 
      temp=t*1.e+4 
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
      temp=max(temp,2.8777e+6/elin) 
      if (lpri.gt.1) then 
        write (lun11,*)'ltyp=75',idest1,idest2,elin,flin,ggup,gglo 
        write (lun11,*)'       ',nrdt,                                &
     &        (masterdata%rdat1(np1r-1+mm),mm=1,8),nind,jlo 
        endif 
      cijpp=masterdata%rdat1(np1r) 
      cijpp=max(cijpp,0.d0) 
      cji=(8.626e-8)*cijpp/tsq/ggup 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
        exptmp=expo(-delt) 
        cij=cji*ggup*exptmp/gglo 
      ans1=cij*xnx 
      ans2=cji*xnx 
      ans6=ans1*dele*ergsev
      ans5=ans2*dele*ergsev
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans6
      ans6=-ans5
      ans5=-anstmp
      if (lpri.gt.1) then 
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp 
        endif 
      elin=0. 
      go to 9000 
!                                                                       
   76 continue 
!     2 photon decay (just  like 50)                                    
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000     
      lpril=0                                             
      if (lpri.ge.2) lpril=lpri 
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         endif 
      if (indonly.eq.1) return
      aij=masterdata%rdat1(np1r) 
      elin=12398.41/abs(eeup-eelo) 
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000 
      flin=(1.d-16)*aij*ggup*elin*elin/((0.667274)*gglo) 
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      a=masterdata%rdat1(np1r2+1) 
      elammu=elin*1.d-4 
!      if (flin.le.1.d-10) flin=1.                                      
      ans1=aij 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      ans2=0. 
      ans4=aij*ergsev*12398.41/abs(elin) 
      ansar2=0. 
      em2ph=aij 
      lskp=1 
      emax=12398.41/elin 
      nbmx=nbinc(emax,epi,ncn2) 
      if (lpril.ge.1)                                                   &
     &  write (lun11,*)'in ucalc, ind=76:',                             &
     &  ml,nilin,nelin,elin,flin,masterdata%rdat1(np1r),gglo,ggup,a,    &
     &  vtherm,ans2,nbmx                                                       
        rcemsum=0. 
        lfastl=0 
        ll=1+lskp 
        do while (ll.le.nbmx) 
          ansar2o=ansar2 
          ansar2=epi(ll)*epi(ll)*max(0.d0,(epi(nbmx)-epi(ll))) 
          rcemsum=rcemsum+(ansar2+ansar2o)                              &
     &                   *(epi(ll)-epi(ll-lskp))/2.                     
          call enxt(epi(1),nb1,lpril,epi,ncn2,t,lfastl,lun11,           &
     &                  ll,lskp,nphint,lrcalc)                          
          ll=ll+lskp 
          enddo 
!        rcemsum=(emax**3)/12.                                          
        rctmp1=0. 
        rctmp2=0. 
        ll=2 
        do while (ll.le.nbmx) 
          ansar2=epi(ll)*epi(ll)*max(0.d0,(epi(nbmx)-epi(ll))) 
          ansar2=ansar2*em2ph*emax/(1.d-24+rcemsum) 
          rctmp1=abund2*ansar2*ptmp1/12.56 
          rctmp2=abund2*ansar2*ptmp2/12.56 
          rccemis(1,ll)=rccemis(1,ll)+rctmp1 
          rccemis(2,ll)=rccemis(2,ll)+rctmp2 
          if (lpril.ge.1) write (lun11,*)ll,epi(ll),ansar2,rctmp1,rctmp2
          call enxt(epi(1),nb1,lpril,epi,ncn2,t,lfastl,lun11,           &
     &                  ll,nskp,nphint,lrcalc)                          
          ll=ll+nskp 
          enddo 
        if (lpril.ge.1)                                                  &
     &  write (lun11,*)'in ucalc, ind=76:',                             &
     &  ml,nilin,nelin,elin,flin,masterdata%rdat1(np1r+2),gglo,ggup,a,  &
     &    vtherm,ans2  
        ans4=aij*ergsev*12398.41/abs(elin) 
!       changing order to allow universal assignment in calc_level_rates_level
        anstmp=ans3
        ans3=-ans4
        ans4=-anstmp
        anstmp=ans1
        ans1=ans2
        ans2=anstmp
        go to 9000 
!                                                                       
   77 continue 
!     coll rates from 71                                                
!     Transition rates from superlevel to spect. lvls                   
!      go to 9000                                                       
      den=xpx 
      clu=0. 
      cul=0. 
      idest1=masterdata%idat1(np1i-1+nidt-3) 
      idest2=masterdata%idat1(np1i+nidt-3) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      eup=leveltemp%rlev(1,idest2) 
      elo=leveltemp%rlev(1,idest1) 
      wav=12398.41/(eup-elo+1.d-24) 
      ekt=0.861707*t 
      delt=wav/ekt 
      lprit=0 
!      if (lpri.gt.0) lprit=1                                           
      temp=t*1.e+4 
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
      temp=max(temp,2.8777e+6/wav) 
      call calt77(lprit,lun11,temp,den,np1r,np1i,cul,clu) 
!     this is counterintuitive:  77 should be proportional to density
      ans1=clu
      ans2=cul
      go to 9000 
!                                                                       
   78 continue 
                                                                        
      go to 9000 
!                                                                       
   79 continue 
!     fluorescence lines                                                
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         endif 
      if (indonly.eq.1) return
      elin=abs(masterdata%rdat1(np1r)) 
      flin=masterdata%rdat1(np1r+1) 
!      if (flin.le.1.d-10) flin=1.                                      
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      a=masterdata%rdat1(np1r-1+5) 
      hij=elin*1.d-8 
      elammu=elin*1.d-4 
      aij=(6.67e+7)*gglo*flin/ggup/elammu/elammu 
!     this is a fudge to avoid badnumerics from fine structure.         
      if (flin.le.1.01e-12) aij=1.e+5 
      if (elin.ge.1.e+9) aij=1.e+5 
      ans1=aij*(ptmp1+ptmp2) 
      ans4=aij*(ptmp1+ptmp2)*ergsev*12398.41/abs(elin) 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      sigma=(0.02655)*flin*elin*(1.d-8)/vtherm 
      sigvtherm=sigma 
!     notice that opakab does not have abundance in                     
      opakab=sigvtherm 
!      ans2=(0.02655)*flin*elin*(1.d-8)/vtherm                          
      ans2=0. 
      go to 9000 
!                                                                       
   80 continue 
! Collisional ionization rates gnd of Fe and Ni                         
      go to 9000 
!                                                                       
   82 continue 
!     Fe UTA rad rates                                                  
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
!     nb check this out:  no bound-bound decays from continuum          
      if ((idest1.le.0).or.(idest1.ge.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.ge.nlev))                          &
     &      go to 9000                                                  
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         endif 
      if (indonly.eq.1) return
      gflin=masterdata%rdat1(np1r+2) 
      aij=masterdata%rdat1(np1r-1+4) 
      elin=abs(masterdata%rdat1(np1r)) 
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000 
!      flin=(1.d-16)*aij*ggup*elin*elin/((0.667274)*gglo)               
      flin=gflin 
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      a=masterdata%rdat1(np1r2+1) 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      ener=12398.41/elin 
      dele=ener*vtherm/3.e+10 
      delev=vtherm/(elin*(1.d-8)) 
      delea=masterdata%rdat1(np1r-1+6)*(4.14e-15) 
      elammu=elin*1.d-4 
      ans1=aij*(ptmp1+ptmp2) 
      sigma=(0.02655)*flin/delev 
!      sigvtherm=(0.02655)*flin*elin*(1.d-8)/3.e+10                     
      sigvtherm=sigma 
      jkkl=derivedpointers%nplini(ml) 
      if (jkkl.le.0) go to 9000 
      ml3=derivedpointers%nplin(jkkl) 
      if (ml3.le.0) go to 9000 
      mlm=ml3
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                              &
     &  0,lun11)                                                  
      elin=abs(masterdata%rdat1(np1r)) 
      ener=12398.41/abs(elin) 
      nb1=nbinc(ener,epi,ncn2) 
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10 
!                                                                       
!     turning off rex                                                   
!      ans2=0. 
!                                                                       
!      ans4=ans1*ener*ergsev                                            
!     notice that opakab does not have abundance in                     
      opakab=sigvtherm 
      lfasto=2 
      ans3=ans2*ener*ergsev 
!     this is a cheat.  there is still an error in the 82/83 data that  
!       makes some fluorescence emission                                
!     this test should prevent calculation when called 
!     from calc_rates_level
!     since abund1 will be zero                                         
      opakb1=sigvtherm*abund1 
!      lpriu=lpri                                                       
      lpriu=0 
      rcem1=0. 
      rcem2=0. 
      if (opakb1.gt.1.d-48)                                             &
     & call linopac(lpriu,lun11,opakb1,                                 &
     &               rcem1,rcem2,elin,vturb,t,a,delea,epi,ncn2,         &
     &               opakc,rccemis,lfasto)              
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc, ind=82:',                          &
     &  ml,nilin,nelin,elin,flin,masterdata%rdat1(np1r+2),gglo,ggup,a,  &
     &  vtherm,ans2,                                                    &
     &  idest1,idest2,idest3,idest4,nlev,sigvtherm,bremsa(nb1),nb1      
      anstmp=ans2
      ans2=ans1
      ans1=anstmp
      go to 9000 
!                                                                       
!                                                                       
   83 continue 
!     Fe UTA level data                                                 
      go to 9000 
!                                                                       
   84 continue 
      lprisv=lpri 
      lpril=0 
      go to 9000 
!      if (lpri.ge.1) lpril=2                                           
      if (lpril.gt.1) then
         write (lun11,*)'ltyp=84',ml,derivedpointers%npar(ml) 
         write (lun11,*)(masterdata%rdat1(np1r-1+jj),jj=1,nrdt) 
         write (lun11,*)(masterdata%idat1(np1i-1+jj),jj=1,nidt) 
         write (lun11,*)(masterdata%kdat1(np1k-1+jj),jj=1,nkdt)
         endif
      if (ml.le.0) go to 9000 
      lfastl=lfast 
      nilin=derivedpointers%npar(ml) 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=1 
      if (indonly.eq.1) return
      ntmp=nrdt/2-1 
      ett2=masterdata%rdat1(np1r) 
      ett=masterdata%rdat1(np1r+2)*13.605692 
      ediff=masterdata%rdat1(np1r-1+2*ntmp)*13.605692-ett2 
      scal2=masterdata%rdat1(np1r+1) 
      do ml2=1,ntmp 
        etmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2+1)                     &
     &         -masterdata%rdat1(np1r+2) 
        stmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2+2)*1.d-18*scal2 
        stmpp(ml2)=max(stmpp(ml2),0.d0) 
        if (lpril.gt.1) write (lun11,*)ml2,etmpp(ml2),stmpp(ml2) 
        enddo 
      ett=ett2-(masterdata%rdat1(np1r-1+2*ntmp+1)                       &
     &        -masterdata%rdat1(np1r+2))*13.6 
      ntmp2=nptmpdim 
      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11) 
      nb1=nbinc(ett,epi,ncn2) 
      numcon2=max(2,ncn2/50) 
!         numcon2=200                                                   
      nphint=ncn2-numcon2 
      nphint=max(nphint,nb1+1) 
      if (lpril.gt.1)                                                   &
     & write (lun11,*)'ltyp=84:',ett,ett2,ediff,ntmp,                   &
     &  etmpp(1),stmpp(1),etmpp(ntmp+1),stmpp(ntmp+1),nb1,nphint        
      if (nb1.ge.nphint-1) go to 9000 
       lprib=0 
       lprib=lpril 
!       call phint5384(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,        
!     $   abund1,abund2,ptmp1,ptmp2,xpx,opakab,delr,                    
!     $   opakc,rccemis,lprib,epi,ncn2,bremsa,t,trad,swrat,xnx,crit53,  
!     $    lfast,lun11)                                                 
      ans1=0. 
      ans3=0. 
      ans4=0. 
      ans2=0. 
      lpri=lprisv 
      go to 9000 
                                                                        
!                                                                       
   85 continue 
      lprisv=lpri 
      lpril=0 
      if (lpri.ge.1) lpril=2                                           
      if (lpril.gt.1) then
         write (lun11,*)'ltyp=85',ml,derivedpointers%npar(ml) 
         write (lun11,*)(masterdata%rdat1(np1r-1+jj),jj=1,nrdt) 
         write (lun11,*)(masterdata%idat1(np1i-1+jj),jj=1,nidt) 
         write (lun11,*)(masterdata%kdat1(np1k-1+jj),jj=1,nkdt)
         endif
      if (ml.le.0) go to 9000 
      lfastl=1 
      nilin=derivedpointers%npar(ml) 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=1 
      if (indonly.eq.1) return
      ett2=masterdata%rdat1(np1r+1)*13.605692 
      nmin=masterdata%idat1(np1i) 
      jkk=idest3 
      zc=dfloat(jkk-114) 
      eion=dble(masterdata%rdat1(np1r+1)) 
      kdim=ncn2 
      far=dble(masterdata%rdat1(np1r+2)) 
      gam=dble(masterdata%rdat1(np1r-1+4)) 
      scal=dble(masterdata%rdat1(np1r-1+5)) 
      call pexs(nmin,kdim,zc,eion,far,gam,scal,                      &
     &                etmp8,stmp8,ierr,lpril,lun11)                     
      do mm=1,ncn2 
        stmpp(mm)=(stmp8(mm))*1.d-18 
        enddo 
      call phintfo(stmpp,ett2*0.8,ans1,ans2,ans3,ans4,ans5,ans6,     &
     & abund1,abund2,xpx,opakab,                                        &
     & opakc,opakcont,lpril,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11) 
      opakab=0. 
      ans6=0.
      ans4=0. 
      ans2=0. 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      lpri=lprisv 
      go to 9000 
                                                                        
   86 continue 
!     iron auger data                                                   
!     this statement causes pileup of populations in some superlevels.  
!      if (idat1(np1i-1+nidat-1).ne.idat1(np1i-1+nidat)+1) go to 9000   
      ans1=masterdata%rdat1(np1r+1) 
      ans2=0. 
      idest1=masterdata%idat1(np1i-1+nidt-3) 
!      idest2=nlevp                                                     
      idest2=nlevp+masterdata%idat1(np1i-1+nidt-4)-1 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
!      anstmp=ans1
!      ans1=ans2
!      ans2=anstmp
      go to 9000 
!                                                                       
   87 continue 
      go to 9000 
!                                                                       
   88 continue 
!     op inner shell photoexcitation                                    
      lprisv=lpri 
      idest1=masterdata%idat1(np1i+nidt-2) 
!      idest2=masterdata%idat1(np1i+nidt-3)                            
      idest2=nlevp 
!      if (lpri.ge.1) lpri=2                                            
      if (indonly.eq.1) return
      lunsv=lun11 
      if (lpri.gt.1) write (lun11,*)'ltyp=88,idest1=',idest1,idest2 
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000 
      if (ml.le.0) go to 9000 
      eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
      ett=eth 
      nilin=derivedpointers%npar(ml) 
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml 
      if (nilin.le.0) go to 9000 
      ntmp=nrdt/2 
      do ml2=1,ntmp 
        etmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2-1) 
        stmpp(ml2)=masterdata%rdat1(np1r-1+2*ml2)*1.d-18 
        stmpp(ml2)=max(stmpp(ml2),0.d0) 
        enddo 
      ntmp2=nptmpdim 
      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11) 
!      ett=ett+max(0.,13.605692*etmpp(1))                               
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1) 
      if (ett.le.0.) go to 9000 
      nb1=nbinc(ett,epi,ncn2) 
                                                                        
      r19=rr/1.e+19 
      mlm=nilin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdti,np1r2,nidti,np1i2,nkdti,np1k2,mlm,                        &
     &  0,lun11)                                                  
      emax=etmpp(ntmp)*13.6+eth 
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,idest2) 
      idest3=masterdata%idat1(np1i-1+nidti) 
      idest4=idest3+1 
      if (lpri.gt.1) write (lun11,*)'before phint53',gglo,ggup 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      if (lpri.gt.1) then 
         write (lun11,*)'type 88 data:',masterdata%idat1(np1i),         &
     &      masterdata%idat1(np1i+nidt-1),t,xnx,eth,gglo,ggup,swrat  
        call dprinto(ndesc,nrdesc,lcon,                              &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)  
        endif 
!      if ((lpri.ge.1).and.(idest1.le.4).and.(jkion.eq.29)              
!     $    .and.(abund1.gt.1.d-34))  then                               
!        lun99=99                                                       
!        write (lun99,*)'type 88 data:', idest1, idest2,                
!     $           eth,gglo,ggup,swrat, abund1,abund2                    
!        call dprinto(ndesc,nrdesc,lcon,                                
!     $          nrdt,rdat,nidt,idat,nkdt,kdat,lun99)                   
!        do mm=1,ncn2                                                   
!          opaksv(mm)=opakc(mm)                                         
!          opakc(mm)=0.                                                 
!          enddo                                                        
!        endif                                                          
      lprib=lpri 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      tm=t*1.e4 
      bktm=bk*tm/ergsev 
      q2=2.07e-16*xnx*(tm**(-1.5)) 
      emltlv=leveltemp%rlev(2,nlev) 
      rs=q2/emltlv 
      ethion=leveltemp%rlev(1,nlev) 
      emltlv=leveltemp%rlev(2,idest1) 
      rnissel=emltlv*rs
      rnisseu=1.
      rnist=rnissel                                                     &
     &  *exp(-(max(0.d0,13.605692*etmpp(1)))/(0.861707)/t)              &
     &  /(1.e-37+rnisseu)                                                   
      call phint53(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,         &
     &  ans5,ans6,abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,           &
     &  opakc,opakcont,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,       &
     &  lfast,lun11)                                                    
!      if ((lpri.ge.1).and.(idest1.le.4).and.(jkion.eq.29)              
!     $    .and.(abund1.gt.1.d-34))  then                               
!        nhit=0                                                         
!        do mm=1,ncn2                                                   
!          if ((opakc(mm).gt.1.d-34).and.(epi(mm).gt.500.)              
!     $        .and.(epi(mm).lt.800.)) then                             
!            write (lun99,919)mm,epi(mm),opakc(mm),                     
!     $        opakc(mm)/max(1.d-34,abund1)/xpx,opaksv(mm)+opakc(mm)    
!            nhit=1                                                     
!            endif                                                      
!          opakc(mm)=opaksv(mm)+opakc(mm)                               
!          enddo                                                        
!        if (nhit.eq.0) write (lun99,*)'no cross section'               
!        endif                                                          
      if (lpri.gt.1) then 
        npr=nb1 
        write (lun11,*)'bautista threshold xsection:',                  &
     &         npr,ett,eth,masterdata%rdat1(np1r),sg(npr),ans2,swrat           
        endif 
      ans2=0. 
      ans4=0. 
      ans6=0. 
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans3
      ans3=-ans4
      ans4=-anstmp
      anstmp=ans5
      ans5=-ans6
      ans6=-anstmp
      lpri=lprisv 
      go to 9000 
!                                                                       
   89 continue 
!     saf line rad. rates                                               
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
!     nb check this out:  no bound-bound decays from continuum          
      if ((idest1.le.0).or.(idest1.ge.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.ge.nlev))                          &
     &      go to 9000                                                  
      eeup=leveltemp%rlev(1,idest1) 
      eelo=leveltemp%rlev(1,idest2) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         endif 
      if (indonly.eq.1) return
      aij=masterdata%rdat1(np1r+2) 
      ans1=aij*(ptmp1+ptmp2) 
!      aij=min(aij,1.e+10)                                              
      elin=abs(masterdata%rdat1(np1r)) 
      if (elin.le.1.d-48) go to 9000 
      ggup=leveltemp%rlev(2,idest1) 
      gglo=leveltemp%rlev(2,idest2) 
      if (ml.le.0) go to 9000 
      nilin=derivedpointers%npar(ml) 
      if (nilin.le.0) go to 9000 
      nelin=derivedpointers%npar(nilin) 
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000 
      flin=(1.d-16)*aij*ggup*elin*elin/((0.667274)*gglo) 
!                                                                       
      mlm=nelin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      a=masterdata%rdat1(np1r2+1) 
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5) 
      ener=12398.41/elin 
      dele=ener*vtherm/3.e+10 
      elammu=elin*1.d-4 
      sigma=(0.02655)*flin*elin*(1.d-8)/vtherm 
      sigvtherm=sigma 
      jkkl=derivedpointers%nplini(ml) 
      if (jkkl.le.0) go to 9000 
      ml3=derivedpointers%nplin(jkkl) 
      if (ml3.le.0) go to 9000 
      mlm=ml3
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                              &
     &  0,lun11)                                                  
      elin=abs(masterdata%rdat1(np1r)) 
      ener=12398.41/abs(elin) 
      enerm=ener-eth
      nb1=nbinc(ener,epi,ncn2) 
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10                          &
     &      *flinabs(ptmp1)                                             
!                                                                       
!     turning off rex                                                   
      ans2=0. 
!                                                                       
      if (elin.gt.0.99e+9) then 
         ans2=0. 
         sigvtherm=0. 
         endif 
!      ans1=ans1+ans2*ggup/(1.d-36+gglo)                                
!     note that now opakab does not have abundance in                   
      opakab=sigvtherm 
      lfasto=2 
!      lfasto=4                                                         
      delea=0. 
      lfnd=0 
      lpriu=0 
!      if (lpri.ge.1) lpriu=3                                           
      call deleafnd(jkion,idest1,                                    &
     &   delea,lfnd,lpriu,lun11)                          
!                                                                       
      if (lfnd.eq.0) delea=masterdata%rdat1(np1r+2)*(4.136e-15) 
      ans4=ans1*ener*ergsev 
      ans3=ans2*ener*ergsev 
      ans6=0.
      ans5=0.
      rcem1=abund2*ans4*ptmp1/(1.d-48+ptmp1+ptmp2) 
      rcem2=abund2*ans4*ptmp2/(1.d-48+ptmp1+ptmp2) 
      opakb1=sigvtherm*abund1 
!     this test should prevent calculation when called 
!     from calc_rates_level
!     since abund1 will be zero                                         
!      lpriu=lpri                                                       
      lpriu=0 
      if ((nrdesc.ne.9).and.(lfasto.le.4).and.(opakb1*delr.gt.1.d-8))   &
     & call linopac(lpriu,lun11,opakb1,                                 &
     &               rcem1,rcem2,elin,vturb,t,a,delea,epi,ncn2,         &
     &               opakc,rccemis,lfasto)              
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'in ucalc, ind=89:',                          &
     &  ml,nilin,nelin,elin,flin,masterdata%rdat1(np1r+2),              &
     &  gglo,ggup,a,vtherm,vturb,&
     &  ans1,ans2,idest1,idest2,idest3,idest4,nlev,sigvtherm,           &
     &  bremsa(nb1),nb1,abund1,abund2,delea,lfnd                        
      go to 9000 
!                                                                       
   90 continue 
      go to 9000 
!                                                                       
   91 continue 
!     a values from atomdb.  same as 50.                                
      go to 50 
      go to 9000 
!                                                                       
   92 continue 
!     collision strengths from atomdb                                   
!                                                                       
      lpril=0 
!      if (lpri.ge.1) lpril=2                                           
      if (lpril.gt.1) then
        write (lun11,*)'ltyp=92',ml,derivedpointers%npar(ml) 
        write (lun11,*)(masterdata%rdat1(np1r-1+jj),jj=1,nrdt) 
        write (lun11,*)(masterdata%idat1(np1i-1+jj),jj=1,nidt) 
        write (lun11,*)(masterdata%kdat1(np1k-1+jj),jj=1,nkdt)
        endif
!                                                                       
!                                                                       
!     general stuff                                                     
      lctype=masterdata%idat1(np1i+3-1) 
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      tmin=masterdata%rdat1(np1r) 
      tmax=masterdata%rdat1(np1r+1) 
      do mml=1,20 
        tstr(mml)=masterdata%rdat1(np1r+1+mml) 
        cstr(mml)=masterdata%rdat1(np1r+21+mml) 
        enddo 
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      elin=12398.41/abs(eeup-eelo+1.d-24) 
      hij=elin*1.d-8 
      if (elin.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      temp=t*1.e+4 
      eij=abs(eeup-eelo) 
      eijkev=eij/1.e+3 
      tk=t*1.e+4 
      mlm=derivedpointers%npar(ml)
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      nistage=masterdata%idat1(np1i2) 
      mlm=derivedpointers%npar(mlm)
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      nzel=masterdata%idat1(np1i2) 
      if (lpril.gt.1)                                                   &
     & write (lun11,*)'before calc_maxwell_rates',lctype,tmin,tmax,    &
     &   eijkev,tk,zzz,gglo,ggup                                        
      call calc_maxwell_rates(lun11,lpril,lctype,tmin,tmax,            &
     &  Tstr,cstr, eijkev,  tk, nzel,  gglo,  ggup,  cij, cji, upsilon) 
      ans1=cij*xnx 
      ans2=cji*xnx 
      ans6=ans1*eij*ergsev
      ans5=ans2*eij*ergsev
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans6
      ans6=-ans5
      ans5=-anstmp
      if (lpril.gt.1) then 
        write (lun11,*)'type 92 data',lctype,upsilon,cij,cji,xnx 
        endif 
       go to 9000 
!                                                                       
!      old code for 92 not used                                         
!      delt=12398.41/elin/ekt                                           
!      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)              
!      eijry=eij/13.605692                                              
!      tsq=sqrt(t)                                                      
!      cijpp=0.                                                         
!      temp=max(temp,2.8777e+6/elin)                                    
!      if (lctype.eq.11) then                                           
!!       chianti type 1 (dere et al. 1997)                              
!        cijpp=upsil(1,eijry,cijpp,cstr(1),cstr(2),                     
!     $      cstr(3),cstr(4),cstr(5),tk)                                
!        endif                                                          
!      if (lctype.eq.12) then                                           
!!       chianti type 2 (dere et al. 1997)                              
!        cijpp=upsil(2,eijry,cijpp,cstr(1),cstr(2),                     
!     $      cstr(3),cstr(4),cstr(5),tk)                                
!        endif                                                          
!      if (lctype.eq.13) then                                           
!!       chianti type 3 (dere et al. 1997)                              
!        cijpp=upsil(3,eijry,cijpp,cstr(1),cstr(2),                     
!     $      cstr(3),cstr(4),cstr(5),tk)                                
!        endif                                                          
!      if (lctype.eq.14) then                                           
!!       chianti type 4 (dere et al. 1997)                              
!        cijpp=upsil(4,eijry,cijpp,cstr(1),cstr(2),                     
!     $      cstr(3),cstr(4),cstr(5),tk)                                
!        endif                                                          
!      if (lctype.eq.31) then                                           
!!       sampson goett and clark 1983 type 1                            
!        if (nrdt.lt.7) go to 9000                                      
!        y=eij/ekt                                                      
!        aa=rdat1(np1r+2)                                               
!        co=rdat1(np1r+3)                                               
!        cr=rdat1(np1r+4)                                               
!        crp=rdat1(np1r+5)                                              
!        rr=rdat1(np1r+6)                                               
!        sig=rdat1(np1r+6)                                              
!        z2s=rdat1(np1r+7)                                              
!        zeff=float(idat1(np1i+2))-sig                                  
!        if (y.gt.40.)  go to 9000                                      
!        call expint(y,em1)                                             
!        e1=em1/y*exp(-y)                                               
!        if (y*a+y.le.80) then                                          
!            call eint(y*a+y,ee1,ee2,ee3)                               
!          else                                                         
!            ee1=0.                                                     
!            ee2=0.                                                     
!            ee3=0.                                                     
!          endif                                                        
!        er=0.                                                          
!        er1=0.                                                         
!        if (rr.eq.1.) then                                             
!          er=ee1                                                       
!          er1=ee2                                                      
!          endif                                                        
!        if (rr.eq.2.) then                                             
!          er=ee2                                                       
!          er1=ee3                                                      
!          endif                                                        
!        if (y*a+y.le.40) then                                          
!            qij=co*exp(-y)+1.55*z2s*e1+y*exp(y*a)*(cr*er/(a+1.)**(rr-1.
!     #      +cr1*er1/(a+1.)**rr)                                       
!          else                                                         
!            qij=co*exp(-y)+1.55*z2s*e1                                 
!          endif                                                        
!        cijpp=qij*exp(y)/zeff/zeff                                     
!        endif                                                          
!      if (lctype.eq.32) then                                           
!!       sampson goett and clark 1983 type 2                            
!        endif                                                          
!      if (lctype.eq.33) then                                           
!!       sampson goett and clark 1983 type 3                            
!        endif                                                          
!      if (lctype.eq.41) then                                           
!!       kato and nakazaki 1989 type 1                                  
!        call calt66(temp,np1r+2,rdat1,gamma)                           
!        cijpp=gamma                                                    
!        endif                                                          
!      if (lctype.eq.42) then                                           
!!       kato and nakazaki 1989 type 3                                  
!        go to 9000                                                     
!        endif                                                          
!      if (lctype.gt.100) then                                          
!        ncase=int(lctype/50)                                           
!!       ltype=100 -->ncase=2                                           
!!       ltype=150 -->ncase=3                                           
!!       ltype=200 -->ncase=4                                           
!!       ltype=250 -->ncase=5                                           
!!       ltype=300 -->ncase=6                                           
!!       ltype=350 -->ncase=7                                           
!!       ltype=400 -->ncase=8                                           
!!       ltype=450 -->ncase=9                                           
!!       ltype=500 -->ncase=10                                          
!!       ltype=550 -->ncase=11                                          
!!       ltype=600 -->ncase=12                                          
!!       ltype=650 -->ncase=13                                          
!!       ltype=700 -->ncase=14                                          
!!       ltype=750 -->ncase=15                                          
!!       ltype=800 -->ncase=16                                          
!!       ltype=850 -->ncase=17                                          
!!       ltype=900 or greater -->ncase=18                               
!!       don't do qs                                                    
!        if ((ncase.ge.6).and.(ncase.le.9)) stop 'ncase=6-9'            
!!        if ((ncase.ge.6).and.(ncase.le.9)) go to 9000                 
!        if (ncase.ge.14) stop 'ncase=14'                               
!!        if (ncase.ge.14) go to 9000                                   
!!                                                                      
!        npts=lctype-50*ncase                                           
!        if ((tk.le.tmin).or.(tk.ge.tmax)) go to 9000                   
!        mm=1                                                           
!        do while ((tk.lt.tstr(mm)).and.(mm.lt.npts))                   
!          mm=mm+1                                                      
!          enddo                                                        
!        mm=max(mm-1,1)                                                 
!        cijpp=cstr(mm)+(cstr(mm+1)-cstr(mm))*(tk-tstr(mm))             
!     $                 /(tstr(mm+1)-tstr(mm)+1.d-38)                   
!        endif                                                          
!                                                                       
!      cji=(8.626e-8)*cijpp/tsq/ggup                                    
!      ekt=0.861707*t                                                   
!      delt=12398.41/elin/ekt                                           
!      exptmp=expo(-delt)                                               
!      cij=cji*ggup*exptmp/gglo                                         
!                                                                       
   93 continue 
      go to 9000 
!     op pi xsections                                                   
      lprisv=lpri 
      if (nrdt.gt.3) go to 533 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest2=nlevp+masterdata%idat1(np1i-1+nidt-3)-1 
      if (indonly.eq.1) return
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2,nlevp,ml 
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000 
      if (ml.le.0) go to 9000 
      eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
      eexc=leveltemp%rlev(1,idest1) 
      ett=eth 
      nilin=derivedpointers%npar(ml) 
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml 
      if (nilin.le.0) go to 9000 
      ntmp=nrdt/2 
!      ett=ett+max(0.,13.605692*etmpp(1))                               
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1) 
      if (ett.le.0.) go to 9000 
      nb1=nbinc(ett,epi,ncn2) 
      xkt=ett/(0.861707*t) 
      r19=rr/1.e+19 
      mlm=nilin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,nlevp) 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        if (ndtmp.le.0) go to 9000 
        mllz=derivedpointers%npar(ndtmp) 
        iltmp=0 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(derivedpointers%npar(ndtmp).eq.mllz))                
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           if (ndtmp.le.0) go to 9000 
           enddo 
!        NB fix to excited level PI and rec                             
         ett=ett+masterdata%rdat1(np1r2) 
         eth=ett 
         ggup=masterdata%rdat1(np1r2+1) 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett                   
         endif 
      if (lpri.gt.1) write (lun11,*)'before phint53pl',eexc,eth,lfast 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      if (lpri.gt.1) then 
         write (lun11,*)'type 93 data:',masterdata%idat1(np1i),         &
     &       masterdata%idat1(np1i+nidt-1),t,xnx,eth,gglo,ggup,swrat   
        call dprinto(ndesc,nrdesc,lcon,                              &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)  
        endif 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      sth=1.d-18*masterdata%rdat1(np1r+1) 
      alph=masterdata%rdat1(np1r+2) 
      e1=masterdata%rdat1(np1r) 
      lfastl=1 
      call phint53pl(sth,e1,alph,ett,ans1,ans2,ans3,ans4,            &
     &  abund1,abund2,ptmp1,ptmp2,xpx,opakab,                           &
     &  opakc,opakcont,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,       &
     &  lfastl,lun11)                                                   
      if (lpri.gt.1) then 
        npr=nb1 
        write (lun11,*)'bautista threshold xsection:',                  &
     &         npr,ett,eth,masterdata%rdat1(np1r),sg(npr),ans2,swrat          
        endif 
      lpri=lprisv 
      go to 9000 
!                                                                       
   94 continue 
      go to 9000 
!     op pi xsections                                                   
!     old version                                                       
      lprisv=lpri 
!      if (lpri.ge.1) lpri=2                                            
      if (nrdt.gt.3) go to 499 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest4=masterdata%idat1(np1i+nidt-3) 
      idest2=nlevp+masterdata%idat1(np1i+nidt-4)-1 
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2 
      if (indonly.eq.1) return
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000 
      if (ml.le.0) go to 9000 
      eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
      ett=eth 
      nilin=derivedpointers%npar(ml) 
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml 
      if (nilin.le.0) go to 9000 
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1) 
      if (ett.le.0.) go to 9000 
      nb1=nbinc(ett,epi,ncn2) 
      xkt=ett/(0.861707*t) 
      r19=rr/1.e+19 
      mlm=nilin
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,                           &
     &  0,lun11)                                                  
      ntmp=nrdt/2 
      gglo=leveltemp%rlev(2,idest1) 
      ggup=leveltemp%rlev(2,nlevp) 
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=idest3+1 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        mllz=derivedpointers%npar(ndtmp) 
        nptmp=mllz 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(nptmp.eq.mllz))                                       
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           nptmp=0 
           if (ndtmp.ne.0) nptmp=derivedpointers%npar(ndtmp) 
           enddo 
         ggup=masterdata%rdat1(np1r2+1) 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup                       
         endif 
      if (lpri.gt.1) write (lun11,*)'before phint53pl' 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      if (lpri.gt.1) then 
         write (lun11,*)'type 94 data:',masterdata%idat1(np1i),         &
     &        masterdata%idat1(np1i+nidt-1),t,xnx,eth,gglo,ggup,swrat 
        call dprinto(ndesc,nrdesc,lcon,                              &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)  
        endif 
      lprib=0 
      if (lpri.gt.1) lprib=lpri 
      sth=1.d-18*masterdata%rdat1(np1r+1) 
      alph=masterdata%rdat1(np1r+2) 
      e1=masterdata%rdat1(np1r) 
      lfastl=1 
      call phint53pl(sth,e1,alph,ett,ans1,ans2,ans3,ans4,            &
     &  abund1,abund2,ptmp1,ptmp2,xpx,opakab,                           &
     &  opakc,opakcont,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,       &
     &  lfastl,lun11)                                                   
      if (lpri.gt.1) then 
        npr=nb1 
        write (lun11,*)'bautista threshold xsection:',                  &
     &         npr,ett,eth,masterdata%rdat1(np1r),sg(npr),ans2,swrat          
        endif 
      lpri=lprisv 
      go to 9000 
!                                                                       
   95 continue 
!     bryans ci rates                                                   
      idest4=masterdata%idat1(np1i+nidt-1)+1 
      idest3=masterdata%idat1(np1i+nidt-1) 
      if (nrdesc.eq.5) then 
         idest1=masterdata%idat1(np1i) 
         if (nidt.ge.3) then 
           idest2=nlevp-1+masterdata%idat1(np1i+1) 
           else 
           idest2=nlevp 
           endif 
        else 
         idest2=1 
         idest1=1 
        endif 
      if (indonly.eq.1) return
      ee=masterdata%rdat1(np1r) 
      tmin=masterdata%rdat1(np1r+1) 
      nspline=(nrdt-2)/2 
      ekt=0.861707*t 
      tt=ekt/ee 
!     constant is ln(2)                                                 
      xx=1.-(0.693147)/log(tt+2.) 
      mm=1 
      do while ((mm.lt.nspline).and.(xx.gt.masterdata%rdat1(np1r+1+mm))) 
        if (lpri.gt.1) write (lun11,*)mm,xx,masterdata%rdat1(np1r+1+mm) 
        mm=mm+1 
        enddo 
!      this illustrates the storage scheme for the splines              
!      do mm=1,nspline                                                  
!        tspline(mm)=rdat1(np1r+1+mm)                                   
!        vspline(mm)=rdat1(np1r+1+nspline+mm)                           
!        enddo                                                          
!     linear interpolation                                              
      rho=(masterdata%rdat1(np1r+1+nspline+mm-1)                        &
     &  +(xx-masterdata%rdat1(np1r+1+mm-1))*                            &
     &    (masterdata%rdat1(np1r+1+nspline+mm)                          &
     &   -masterdata%rdat1(np1r+1+nspline+mm-1))                        &
     &    /(masterdata%rdat1(np1r+1+mm)-masterdata%rdat1(np1r+1+mm-1))) 
!      dere equation 7                                                  
      call eint(1./tt,e1,e2,e3) 
      citmp1=1.d-6*e1*rho/sqrt(tt*ee**3) 
      ans1=citmp1*xnx 
      ans2=0. 
      idest1=1 
!      idest2=1                                                         
      ggup=leveltemp%rlev(2,nlevp) 
      gglo=leveltemp%rlev(2,idest1) 
!     note that rinf has exponential removed                            
      tsq=sqrt(t) 
      rinf=(2.08e-22)*gglo/ggup/t/tsq 
      ans2=ans1*rinf*xnx/expo(-1./tt) 
!      idest1=idat1(np1i+nidt-2)                                        
      idest4=masterdata%idat1(np1i+nidt-1)+1 
      idest3=masterdata%idat1(np1i+nidt-1) 
      ans6=-ans1*ee*ergsev
      ans5=-ans2*ee*ergsev
      go to 9000 
!                                                                       
   96 continue 
!     Autoinization rates (in s^-1) for satellite lvls                  
!        from safranova in kato et al. 1997 atndt 67 225                
      lpril=0 
!      if (lpri.gt.0) lpril=2                                           
      idest3=masterdata%idat1(np1i+nidt-1) 
      idest4=masterdata%idat1(np1i+nidt-3) 
      idest2=masterdata%idat1(np1i+nidt-2)+nlev-1 
      idest1=masterdata%idat1(np1i-1+nidt-3) 
      idest1=max(idest1,1) 
      idest2=max(idest2,1) 
      if (indonly.eq.1) return
      temp=t*1.e+4 
      dele=masterdata%rdat1(np1r+2) 
      rs=2.069e-3/(temp**1.5) 
      ekt=0.861707*t 
      rate=rs*expo(-dele/ekt)*masterdata%rdat1(np1r+1) 
      ans2=rate*xnx 
      if (lpril.ne.0)                                                   &
     & write (lun11,979)                                                &
     &  (masterdata%rdat1(np1r-1+lk),lk=1,3),                           &
     &  (masterdata%idat1(np1i-1+lk),lk=1,4),rate,ans2   
  979 format (1x,'type 96 calc:',3(1pe11.3),4i4,10(1pe11.3)) 
      go to 9000 
!                                                                       
   97 continue 
!     ci rate in terms of upsilon                                       
!     to do fexxiv --> fexxv 1s2s(3S) ci from patrick                   
      if (lpri.gt.1)                                                    &
     &  write (lun11,*)'type 97:',ndesc,lcon,nrdt,nidt,nkdt,            &
     &  ml,(masterdata%rdat1(np1r+mm-1),mm=1,nrdt),                     &
     &  (masterdata%idat1(np1i+mm-1),mm=1,nidt),                        &
     &  (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)                              
      idest4=masterdata%idat1(np1i+nidt-1)+1 
      idest3=masterdata%idat1(np1i+nidt-1) 
      if (nrdesc.eq.5) then 
         idest1=masterdata%idat1(np1i) 
         if (nidt.ge.3) then 
           idest2=nlevp-1+masterdata%idat1(np1i+1) 
           else 
           idest2=nlevp 
           endif 
        else 
         idest2=1 
         idest1=1 
        endif 
      if (indonly.eq.1) return
      eth=leveltemp%rlev(4,idest1)-leveltemp%rlev(1,idest1) 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        if (ndtmp.le.0) go to 9000 
        mllz=derivedpointers%npar(ndtmp) 
        iltmp=0 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(derivedpointers%npar(ndtmp).eq.mllz))                 
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           if (ndtmp.le.0) go to 9000 
           enddo 
!        NB fix to excited level PI and rec                             
         ett=ett+masterdata%rdat1(np1r2) 
         eth=ett 
         ggup=masterdata%rdat1(np1r2+1) 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett                   
         endif 
      nspline=nrdt/2 
      ekt=0.861707*t 
      xx=ekt 
      mm=1 
      do while ((mm.lt.nspline)                                         &
     &      .and.(ekt.gt.masterdata%rdat1(np1r+mm-1)))                  
        if (lpri.gt.1) write (lun11,*)mm,xx,masterdata%rdat1(np1r+mm-1) 
        mm=mm+1 
        enddo 
      do mm5=1,nspline 
        rdattmp(mm5)=masterdata%rdat1(np1r+mm5+nspline-1) 
!        NB a fudge for Fe XXIV --> Fe XXV                              
!        if ((jkion.eq.349).and.(idest1.eq.1).and.(idest2.gt.53)) then  
!          rdattmp(mm5)=rdat1(np1r+mm5+nspline-1)*0.5                   
!          endif                                                        
        enddo 
      cijpp=rdattmp(mm)+(xx-masterdata%rdat1(np1r+mm-1))*               &
     &    (rdattmp(mm+1)-rdattmp(mm))                                   &
     &    /(masterdata%rdat1(np1r+mm)-masterdata%rdat1(np1r+mm-1))        
      gglo=leveltemp%rlev(2,idest1) 
      cji=(8.626e-8)*cijpp/tsq                                          &
     &      /ggup                                                       
      delt=eth/ekt 
      exptmp=expo(-delt) 
      cij=cji*ggup*exptmp/gglo 
!                                                                       
!     NB a fudge which boosts this works                                
!      if (idest2.eq.nlevp+1) cij=cij*10.                               
!                                                                       
      ans1=cij*xnx 
!     note that rinf has exponential removed                            
      rinf=(2.08e-22)*gglo/ggup/t/tsq 
      ans2=ans1*rinf*xnx/exptmp 
      ans6=ans1*eth*ergsev
      ans5=ans2*eth*ergsev
!      idest1=idat1(np1i+nidt-2)                                        
      idest4=masterdata%idat1(np1i+nidt-1)+1 
      idest3=masterdata%idat1(np1i+nidt-1) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'ltype=97:',cij,cji,cijpp,mm,rdattmp(mm),         &
     & masterdata%rdat1(np1r+mm-1),xx,idest1,idest2,gglo,ggup,          &
     & delt,exptmp,eth,ekt,ans1,ans2                                         
!                                                                       
      go to 9000 
!                                                                       
   70 continue 
!     old type 70
!     Coefficients for phot x-section of suplevels                      
!      lfastl=lfast                                                     
      lfastl=3 
      temp=t*1.e+4 
      ans3=0. 
      ans4=0. 
      den=xpx 
      m=1000 
      lpric=0
      mlm=ml
      call drd(ltyp,lrtyp,lcon,                                         &
     &  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                              &
     &  0,lun11)                                                  
!      if (lpri.ge.1) lpric=2                                           
      mlion=derivedpointers%npar(ml) 
      idest1=masterdata%idat1(np1i+nidt-2) 
      idest1=min(idest1,nlev-1) 
      idest2=nlev+masterdata%idat1(np1i-1+nidt-3)-1 
      idest2=max(idest2,nlev) 
      ggup=leveltemp%rlev(2,nlevp) 
      ett=abs(leveltemp%rlev(1,idest1)-leveltemp%rlev(1,nlevp)) 
      if (lpric.ge.1)                                                   &
     & write (lun11,*)'rlev:',idest1,nlevp,                             &
     &    leveltemp%rlev(1,idest1),leveltemp%rlev(1,nlevp) 
      if (idest2.gt.nlevp) then 
        jkk3=jkion+1 
        if (lpric.gt.1)                                                 &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        ndtmp=derivedpointers%npfi(13,jkk3) 
        if (lpric.gt.1)                                                 &
     &    write (lun11,*)jkk3,ndtmp,nlevp,idest2                        
        mllz=derivedpointers%npar(ndtmp) 
        iltmp=0 
        nptmp=mllz 
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))          &
     &      .and.(nptmp.eq.mllz))                                       
           mlm=ndtmp
           call drd(ltyp2,lrtyp2,lcon2,                                 &
     &       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                   &
     &       0,lun11)                                             
           iltmp=masterdata%idat1(np1i2+nidt2-2) 
           if (lpric.gt.1) write (lun11,*)nidt2,iltmp,ndtmp 
           ndtmp=derivedpointers%npnxt(ndtmp) 
           nptmp=0 
           if (ndtmp.ne.0) nptmp=derivedpointers%npar(ndtmp) 
           enddo 
         ggup=masterdata%rdat1(np1r2+1) 
         ett=abs(leveltemp%rlev(1,idest1)+masterdata%rdat1(np1r2)) 
         endif 
       if (lpric.ge.1)                                                  &
     &    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett                   
      xkt=ett/(0.861707*t) 
      nb1=nbinc(ett,epi,ncn2) 
      mlm=mlion
      call drd(ltyp2,lrtyp2,lcon2,                                      &
     &  nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,                        &
     &  0,lun11)                                                  
      ist=masterdata%idat1(np1i2) 
      ic=ist 
      eth=ett 
      gglo=leveltemp%rlev(2,idest1) 
      if (ggup.le.1.d-24) then 
        write (lun11,*) 'ggup error' 
        return 
        endif 
      swrat=gglo/ggup 
      if (lpric.ne.0) then 
         write (lun11,*)'type 99 data:',masterdata%idat1(np1i),          &
     &       masterdata%idat1(np1i+nidt-1),t,xnx,eth,gglo,ggup,swrat
         call dprinto(ndesc,nrdesc,lcon,                              &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)  
        endif 
      ettry=ett/13.6 
!     nb eliminating high density for H
      if (jkion.eq.1) den=min(den,1.d+8)
      m=nrdt
      call calt70(temp,den,ettry,ic,m,np1r,np1i,                      &
     &             ntmp,etmpp,stmpp,rec,al,lun11,lpric)                 
      if (lpric.ne.0) write (lun11,*)'after  calt70:',rec,stmpp(1) 
      crit53=0.01 
      do mm=1,ntmp 
        stmpp(mm)=stmpp(mm)*1.d-18 
        stmpp(mm)=max(stmpp(mm),0.d0) 
        enddo 
      call phint53hunt(stmpp,etmpp,ntmp,ett,ans1,ans2d,ans3d,ans4s,  &
     & lpric,epi,ncn2,bremsa,t,swrat,xnx,crit53,lfastl,lun11)           
      if (ans2d.le.1.d-48) then 
        ans1=0. 
        ans2=0. 
        go to 9000 
        endif 
      scale=rec*xnx/ans2d 
      ans1=ans1*scale 
!     does the swrat not belong?                                        
!      ans2=rec*xnx*swrat                                               
      ans2=rec*xnx 
!      ans2=ans2d                                                       
      tm=t*1.e4 
      q2=2.07e-16*xnx*(tm**(-1.5)) 
      rs=q2/swrat 
      ans1o=ans1 
!      ans1=min(ans1,ans2/rs)                                           
      if (lpric.ge.2)                                                   &
     & write (lun11,*)'type 70 limit:',ans2,rs,swrat,                   &
     &   xnx,tm,q2,ans1o,ans1,scale,rec                                 
!                                                                       
!     nb testing superlevel phot.                                       
!      ans1=0.                                                          
!
      go to 9000 
!                                                                       
   98 continue 
!     line rates, col, burgess and tully for chianti                    
      lpril=0 
!      if (lpri.ge.1) lpril=2                                           
      idest1=masterdata%idat1(np1i) 
      idest2=masterdata%idat1(np1i+1) 
      if (indonly.eq.1) return
      if ((idest1.le.0).or.(idest1.gt.nlev)                             &
     &  .or.(idest2.le.0).or.(idest2.gt.nlev))                          &
     &      go to 9000                                                  
      eeup=leveltemp%rlev(1,idest2) 
      eelo=leveltemp%rlev(1,idest1) 
      if (eeup.lt.eelo) then 
         itmp=idest1 
         idest1=idest2 
         idest2=itmp 
         eeup=leveltemp%rlev(1,idest2) 
         eelo=leveltemp%rlev(1,idest1) 
         endif 
      ggup=leveltemp%rlev(2,idest2) 
      gglo=leveltemp%rlev(2,idest1) 
      c=masterdata%rdat1(np1r+2) 
      eijry=masterdata%rdat1(np1r) 
      eij=eijry*13.605692 
      elin=12398.41/eij 
      hij=elin*1.d-8 
      ntem=(nrdt-3)/2 
      if (lpril.gt.1)                                                   &
     & write (lun11,*)'type 98 data:',elin,ntem                         
      if (elin.le.1.d-24) go to 9000 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      if (lpril.gt.1)                                                   &
     & write (lun11,*)elin,ekt,delt                                     
      do mm=1,ntem 
        tstr(mm)=masterdata%rdat1(np1r+2+mm) 
        cstr(mm)=masterdata%rdat1(np1r+2+ntem+mm) 
        if (lpril.gt.1) write (lun11,*)mm,tstr(mm),cstr(mm) 
        enddo 
      tk=t*1.e+4 
      tk=max(tk,2.8777e+6/elin) 
      ik=masterdata%idat1(np1i+nidt-2) 
!     whats going on here                                               
!      cijpp=gglo*upsiln(ik,eijry,c,ntem,cstr,tstr,tk,lpril,lun11)      
      cijpp=upsiln(ik,eijry,c,ntem,cstr,tstr,tk,lpril,lun11) 
      ekt=0.861707*t 
      delt=12398.41/elin/ekt 
      cji=(8.626e-8)*cijpp/tsq                                          &
     &      /ggup                                                       
      exptmp=expo(-delt) 
      cij=cji*ggup*exptmp/gglo 
      if (lpril.gt.1)                                                   &
     & write (lun11,*)'ltyp=98',c,p1,p2,p3,p4,p5,ik,                    &
     &    eij,idest1,idest2,cij,cji,xnx,cijpp,exptmp,delt,gglo,ggup     
      ans1=cij*xnx 
      ans2=cji*xnx 
      ans6=ans1*eij*ergsev
      ans5=ans2*eij*ergsev
!     changing order to allow universal assignment in calc_level_rates_level
      anstmp=ans6
      ans6=-ans5
      ans5=-anstmp
      elin=0. 
      lpril=0 
      go to 9000 
                                                                        
!                                                                       
 9000 continue 
!                                                                       
      call remtms(time2) 
!      write (lun11,*)'ndesc=',ndes!                                    
      tucalc(ndesc)=tucalc(ndesc)+abs(time2-time1)                   
      ncall(ndesc)=ncall(ndesc)+1                                    
!                                                                       
      if (lpri.gt.1)                                                    &
     & write (lun11,9931)krdesc(nrdesc),kdesc(ndesc),ndesc,ans1,ans2,   &
     &     ans3,ans4,idest1,idest2,masterdata%rdat1(np1r)                    
 9931 format (1x,'in ucalc :',a28,a56,i4,4x,4(1pe10.2),2i4,          &
     &     3(1pe10.2)) 
!                                                                       
      return 
      end                                           
