      subroutine dbwk2(linst,abel,np2,nptrt,karg,lpri,lun11,            &
     & nlsvn,ncsvn)                                             
!                                                                       
!     this program manipulates the database                             
!                                                                       
!     functions include:                                                
!        read in and write out                                          
!        add a record (on end)                      (instruction 2)     
!        delete a record                            (instruction 3)     
!        print all records of certain type          (instruction 4)     
!        exit                                       (instruction 5)     
!        create key vectors linking certain records (instruction 7)     
!        sort records by ion number and data type   (instruction 8)     
!                                                                       
!     data structures are:                                              
!      data: the database arrays (integer, real, character)             
!       idat1(nidat1)                                                   
!       rdat1(nrdat1),                                                  
!       kdat1(nkdat1)                                                   
!     descriptions of database entries, and pointers                    
!       nptrs(nptt,ndat2)                                               
!         nptrs(2,nx)=data type                                         
!         nptrs(3,nx)=rate type                                         
!         nptrs(4,nx)=continuation flag                                 
!                       (n=number of continuations to come)             
!         nptrs(5,nx)=number of reals                                   
!         nptrs(6,nx)=number of integers                                
!         nptrs(7,nx)=number of characters                              
!         nptrs(8,nx)=pointer to reals                                  
!         nptrs(9,nx)=pointer to integers                               
!         nptrs(10,nx)=pointer to characters                            
!                                                                       
!       pointers:                                                       
!       next record:                                                    
!         npnxt(ndat2)                                                  
!       parent record (=ion header or element header)                   
!         npar(ndat2)                                                   
!       first record of a given rate type                               
!         npfirst(ntyp)                                                 
!       first record of rate type ntyp for ion nni                      
!         npfi(ntyp,nni)                                                
!       pointer for line data from array containing luminosities        
!         nplin(nnnl)                                                   
!       (inverse) pointer for line data to array containing luminosities
!          from database array                                          
!         nplini(ndat2)                                                 
!       pointer for continuum data (pi xsection) from array containing l
!         npcon(nnml)                                                   
!       pointer to abundance array to first level of element nni        
!         npconi2(ndat2)                                                
!       (inverse) pointer for continuum data (pi xsection) from array co
!           luminosities                                                
!         npconi(ndat2)                                                 
!                                                                       
!       local variables used to construct pointers                      
!       rdat(100)                                                       
!       idat(100)                                                       
!       kdat(100)                                                       
!       nemap(13)                                                       
!       nimap(168)                                                      
!       npnxt2(ndat2)                                                   
!       npfirst2(ntyp)                                                  
!       mlold2(50)                                                      
!                                                                       
      use globaldata
      implicit none 
!
        integer nel(30)                                                         
      integer nion(500)                                                       
      integer, allocatable, dimension(:) :: idat ! integer data
      real(8), allocatable, dimension(:) :: rdat  ! real data
      character(1), allocatable, dimension(:) :: kdat ! character data
      integer, allocatable, dimension(:) :: idato ! integer data
      real(8), allocatable, dimension(:) :: rdato  ! real data
      character(1), allocatable, dimension(:) :: kdato ! character data
      real(8) abel(nl)
      integer nemap(nl),nimap(nni) 
      integer luse(nni),nlines(nni)
      integer nptrt(ndat2),itmp2(ndat2),itmp3(ndat2),                   &
     &          itmp33(ndat2),nptrt3(ndat2),itmp22(ndat2)               
      real(8) xpx,xeltp,t,ergsev,ener,emin,emax,                        &
     &     elin,aij,xee,etst      
      integer linst,nrdt,npassmx,npass,np1r,np1i,np1k,                  &
     &     nidt,nkdt,np1ki,nnz,nnion,nlvtmp2,nlvtmp1,nltot  
      integer                                                           &
     &     nlev,nkdti,ninxt2,ninxt,nilin,netmp,nenxt2,                  &
     &     nenxt,ndel,nchng,mt2,mmtmp,mm1,mm2,mm,mls,mlt,mltyp,mlz      
      integer                                                           &
     &     mlrtr,mlrto,mlrtm,mlrt,mlpar,mlo,mln,mlm,mllz2,mllz,mlleltp, &
     &     mllel,mll2,mll,mlk,mlionr,mliono,mlionm,mlion,mli,mlexp      
      integer                                                           &
     &     mlel1,mlel2,mlel,ml2,ml12,ml1,ml,lup,lun9,lun8,              &
     &     ltypo,lrtyp,lrtypo,lrtyp2,ltyp,lrdat,lrdato,lprisv     
      integer                                                           &
     &     lpril,lprid,lomiti,lomita,llo,ll,lktmp,                      &
     &     lkm,lk,lidat,lidato,lchng,lcon,lcono,kl,klel,klion      
      integer                                                           &
     &     ktt,kkl,kk,jxx,jlk,jkkl,jkk2,jkk,j,jk,jkion,itmp,            &
     &     isum,iso,iqo,ilv,iq,ii,ij,ijk  
      integer                                                           &
     &     ijsv,idest1,idest2,nlsvn,np2,lun11,lpri,nvtot,               &
     &     lenact,mlline,ncsvn,lkdat,lkdato,nlnsv(nni),ntmp,nilino      
      integer miso,miso2,mel,m4,msvtmp(30),mkk,m3 
      integer mlold(ntyp)
!                                                                       
      character(50) kblnk20 
      character(48) kdesc(ntyp) 
      character(30) krdesc(ntyp) 
      character(1) ktst 
      character(1) ku,kd,kt 
      character(12) klab1 
      character(4) ktmp,karg(20),kblnk4 
      character(9) kinam1 
      character(1) kblnk
!                                                                       
      logical ex 
!                                                                       
      data ku/'0'/,kd/'d'/,kt/'t'/ 
      data kblnk4/'    '/,kblnk/' '/ 
      data krdesc(1)/'ground state ionization       '/ 
      data krdesc(2)/'level ionization/recombination'/ 
      data krdesc(3)/'bound-bound collision         '/ 
      data krdesc(4)/'bound-bound radiative         '/ 
      data krdesc(5)/'bound-free collision (level)  '/ 
      data krdesc(6)/'total recombination           '/ 
      data krdesc(8)/'total recombination           '/ 
      data krdesc(7)/'bound-free radiative (level)  '/ 
      data krdesc(9)/'2 photon decay                '/ 
      data krdesc(11)/'element data                  '/ 
      data krdesc(12)/'ion data                      '/ 
      data krdesc(13)/'level data                    '/ 
      data krdesc(23)/'collisional superlevel->spect '/ 
      data krdesc(14)/'radiative superlevel->spect   '/ 
      data krdesc(15)/'CI total rate                 '/ 
      data krdesc(40)/'CI from superlevels           '/ 
      data kdesc(1)/'radiative recombination:  aldrovandi and pequign'/ 
      data kdesc(2)/'charge exch. h0: Kingdon and Ferland            '/ 
      data kdesc(3)/'autoionization: hamilton, sarazin chevalier     '/ 
      data kdesc(4)/'line data radiative: mendosa; raymond and smith '/ 
      data kdesc(5)/'2 photon transition collisional                 '/ 
      data kdesc(6)/'level data                                      '/ 
      data kdesc(7)/'dielectronic recombination: aldrovandi and pequi'/ 
      data kdesc(8)/'dielectronic recombination: arnaud and raymond  '/ 
      data kdesc(9)/'charge exch. H0 Kingdon and Ferland             '/ 
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
      data kdesc(37)/'                                                '/ 
      data kdesc(38)/'                                                '/ 
      data kdesc(39)/'                                                '/ 
      data kdesc(40)/'                                                '/ 
      data kdesc(41)/'                                                '/ 
      data kdesc(42)/'                                                '/ 
      data kdesc(43)/'                                                '/ 
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
      data kdesc(59)/'verner pi xc                                    '/ 
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
!                                                                       
       do mm=1,50 
         write (kblnk20(mm:mm),'(a1)')kblnk 
         enddo 
!                                                                       
       lprid=0 
!                                                                       
!        write (lun11,*)'  ' 
!        write (lun11,*)'in dbwk: istruction index=',linst 
!                                                                       
!       delete record                                                   
        if (linst.eq.3) then 
!
          allocate(idat(2000000))
          allocate(rdat(2000000))
          allocate(kdat(2000000))
!
           read (5,*)ndel 
           do 1011 ml=ndel,np2-1 
             call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,    &
     &                 kdat,ml,lprid,lun11)          
             call dprinto2(ltyp,lrtyp,lcon,                             &
     &        lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                   
!             call dprintn(ltyp,lrtyp,lcon,                             
!     $        lrdat,rdat,lidat,idat,lkdat,kdat,nptr1(ml-1),ml-1)       
 1011        continue 
           np2=np2-1 
!
          deallocate(idat)
          deallocate(rdat)
          deallocate(kdat)
!
           endif 
!                                                                       
!       print out records sorted by isoseq                              
      if (linst.eq.9) then 
!
        allocate(idat(2000000))
        allocate(rdat(2000000))
        allocate(kdat(2000000))
!
!       special loop for element data                                   
        ml=derivedpointers%npfirst(11) 
 2099     continue 
          call dread(ltyp,lrtyp2,lcon,lrdat,rdat,lidat,idat,lkdat,      &
     &                 kdat,ml,lprid,lun11)          
          call dprinto2(ltyp,lrtyp2,lcon,                               &
     &               lrdat,rdat,lidat,idat,lkdat,kdat,lun11)            
          ml=derivedpointers%npnxt(ml) 
          if (ml.ne.0) go to 2099 
        do mls=1,28 
!          write (lun11,*)'isosequence=',mls 
          ml12=derivedpointers%npfirst(12) 
          do mli=1,nni 
            if (nion(mli).ne.0) then 
            mlel=derivedpointers%npar(ml12) 
            call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,     &
     &                 kdat,mlel,lprid,lun11)          
            nnz=idat(1) 
            call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,     &
     &                 kdat,mll2,lprid,lun11)          
            nnion=idat(1) 
            iso=nnz+1-nnion 
            lomiti=0 
!            if ((iso.ge.3).and.(iso.le.8)) lomiti=1                    
!            write (lun11,*)mls,mli,mlel,nnz,nnion,iso                  
            if (iso.eq.mls) then 
!              write (lun11,*)'ion number=',mli                         
              call dprinto2(ltyp,lrtyp2,lcon,                           &
     &               lrdat,rdat,lidat,idat,lkdat,kdat,lun11)            
              do mll=1,ntyp 
                lomita=0 
                if ((lomiti.eq.1).and.((mll.eq.3).or.                   &
     &              (mll.eq.4).or.(mll.eq.5).or.(mll.eq.2).or.          &
     &              (mll.eq.2).or.(mll.eq.13).or.(mll.eq.7))) lomita=1  
                if ((mll.ne.12).and.(mll.ne.11)                         &
     &             .and.(lomita.ne.1)) then                             
                ml1=derivedpointers%npfi(mll,mli) 
                if (ml1.ne.0) then 
!                  write (lun11,*)'rate type=',mll,' ',krdesc(mll),ml1  
                  ml=ml1 
 2097             continue 
                  call dread(ltyp,lrtyp2,lcon,                          &
     &                 lrdat,rdat,lidat,idat,lkdat,                     &
     &                 kdat,ml,lprid,lun11)          
                  call dprinto2(ltyp,lrtyp2,lcon,                       &
     &               lrdat,rdat,lidat,idat,lkdat,kdat,lun11)            
!                  write (lun11,*)ml,mli,npni(ml,mli)                   
                  ml=derivedpointers%npnxt(ml) 
                if (ml.ne.0) go to 2097 
                endif 
              endif 
              enddo 
            endif 
            ml12=derivedpointers%npnxt(ml12) 
          endif 
          enddo 
        enddo 
!
      deallocate(idat)
      deallocate(rdat)
      deallocate(kdat)
!
      endif 
!                                                                       
!      do general sort, new quicksort                                   
      if (linst.eq.20) then 
!
        allocate(idat(2000000))
        allocate(rdat(2000000))
        allocate(kdat(2000000))
!
        write (6,*)'linst=20',np2 
!         step through and unpack into temporary array and key          
          mlt=1 
 7013         continue 
              ml=nptrt(mlt) 
              call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,   &
     &                 kdat,ml,lprid,lun11)          
!             call dprinto(ltyp,lrtyp,lcon,                             
!     $        lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                  
!              write (6,*)'ml=',ml,mlt,ltyp,lrtyp,lidat,idat(lidat),    
!     $           rdat(1)                                               
              if ((ml.le.0).or.(ml.gt.np2).or.(idat(lidat).eq.0)        &
     &          .or.(idat(lidat).gt.nni))                               &
     &          stop 'sort index error'                                 
              if ((ltyp.ne.0).and.(idat(lidat).ne.0)) then 
                mlt=mlt+1 
                itmp2(mlt-1)=ltyp 
                itmp3(mlt-1)=idat(lidat) 
                nptrt(mlt-1)=mlt-1 
                endif 
              if (mlt.lt.np2) go to 7013 
            np2=mlt 
!            write (6,*)'before quicksort',np2                          
!            do ii=1,np2                                                
!              write (6,*)ii,itmp2(ii),itmp3(ii),nptrt(ii)              
!              enddo                                                    
            call sort4(np2,itmp2,itmp3,nptrt) 
!            write (6,*)'after quicksort'                               
!            do ii=1,np2                                                
!              write (6,*)ii,itmp2(ii),itmp3(ii),nptrt(ii)              
!              enddo                                                    
!           now reorder and do according to second index                
            ii=1 
            ij=0 
 8018         continue 
              ii=ii+1 
!              write (6,*)ii,itmp2(ii),itmp2(ii+1),itmp3(ii),nptrt(ii),i
              ij=ij+1 
              itmp33(ij)=itmp3(ii) 
              itmp22(ij)=itmp2(ii) 
              nptrt3(ij)=nptrt(ii) 
              if ((itmp2(ii).ne.itmp2(ii+1)).or.(ii.eq.np2)) then 
!                   write (6,*)'before quicksort 2',ii,itmp2(ii),ij     
!                   do ii2=1,ij                                         
!                    write (6,*)ii2,itmp22(ii2),itmp33(ii2),nptrt3(ii2) 
!                    enddo                                              
                  ijsv=ij 
                  call sort4(ij,itmp33,nptrt3,itmp22) 
!                  write (6,*)'after quicksort 2',ii,itmp2(ii)          
                  do ijk=1,ijsv 
                    nptrt(ii-ijsv+ijk)=nptrt3(ijk) 
                    itmp3(ii-ijsv+ijk)=itmp33(ijk) 
                    itmp2(ii-ijsv+ijk)=itmp22(ijk) 
!                    write (6,*)ijk,itmp22(ijk),itmp33(ijk),nptrt3(ijk) 
                    enddo 
                  ij=0 
                endif 
              if (ii.lt.np2) go to 8018 
!            write (6,*)'after quicksort'                               
             do ii=1,np2-1 
!              write (6,*)ii,itmp2(ii),itmp3(ii),nptrt(ii)              
              if ((itmp2(ii).gt.itmp2(ii+1))                            &
     &             .or.((itmp2(ii).eq.itmp2(ii+1))                      &
     &               .and.(itmp3(ii).gt.itmp3(ii+1)))) then             
                 write (6,*)ii+1,itmp2(ii+1),itmp3(ii+1),nptrt(ii+1) 
                 write (6,*)'sort error!' 
                 stop 
                 endif 
              enddo 
!
        deallocate(idat)
        deallocate(rdat)
        deallocate(kdat)
!
        endif 
!                                                                       
!      do general sort, old bubble                                      
      if (linst.eq.10) then 
!
        allocate(idat(2000000))
        allocate(rdat(2000000))
        allocate(kdat(2000000))
!
        allocate(idato(2000000))
        allocate(rdato(2000000))
        allocate(kdato(2000000))
!
!
        write (6,*)'linst=10',np2 
          ktmp=kblnk4 
          ktmp=karg(1) 
          read (ktmp(1:1),'(a1)')ktst 
          lktmp=lenact(ktmp) 
          isum=0 
          do ml=2,2+lktmp 
            mlm=ml-1 
            read (ktmp(ml:ml),'(i1)')itmp 
            mlexp=lktmp-mlm-1 
            isum=isum+itmp*10**mlexp 
            enddo 
          write (6,*)ktst,isum 
          npass=0 
          npassmx=2*ndat2 
!          npassmx=100                                                  
 4099       npass=npass+1 
            write (6,*)'pass=',npass 
            lchng=0 
            nchng=0 
            mlt=2 
 4013         continue 
              ml=nptrt(mlt) 
              call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,   &
     &                 kdat,ml,lprid,lun11)          
              mlo=nptrt(mlt-1) 
              call dread(ltypo,lrtypo,lcono,                            &
     &                 lrdato,rdato,lidato,idato,lkdato,                &
     &                 kdat,mlo,lprid,lun11)          
              if (ktst.eq.'n') then 
                if (isum.eq.1) then 
                  iq=ltyp 
                  iqo=ltypo 
                elseif(isum.eq.2) then 
                  iq=lrtyp 
                  iqo=lrtypo 
                endif 
              elseif (ktst.eq.'i') then 
                iq=idat(isum) 
                iqo=idato(isum) 
              elseif (ktst.eq.'m') then 
                iq=idat(lidat-isum) 
                iqo=idato(lidato-isum) 
              endif 
              if (iq.gt.iqo) then 
                lchng=1 
                nptrt(mlt)=mlo 
                nptrt(mlt-1)=ml 
!                write (6,*)ml,nptrt(mlt)                               
                nchng=nchng+1 
                endif 
              mlt=mlt+1 
              if (mlt.lt.np2) go to 4013 
!            do lm=1,np2                                                
!              write (6,*)lm,nptrt(lm)                                  
!              enddo                                                    
            write (6,*)'nchng=',nchng 
            if ((npass.lt.npassmx).and.(lchng.ne.0)) go to 4099 
!
        deallocate(idat)
        deallocate(rdat)
        deallocate(kdat)
!
!
        deallocate(idato)
        deallocate(rdato)
        deallocate(kdato)
!
        endif 
!                                                                       
!     print out records after general sort                              
      if (linst.eq.11) then 
!
        allocate(idat(2000000))
        allocate(rdat(2000000))
        allocate(kdat(2000000))
!
        write (6,*)'linst=11',np2 
        mlt=1 
        lun8=6 
        ltyp=0 
 5013     continue 
          ml=nptrt(mlt) 
          ltypo=ltyp 
          ltyp=masterdata%nptrs(2,ml) 
          lrtyp=masterdata%nptrs(3,ml) 
          if ((ltyp.eq.0).or.(lrtyp.eq.0)) then 
!            write (6,*)'skipping',ml,mlt,ltyp,lrtyp                    
            mlt=mlt+1 
            if (mlt.lt.np2) go to 5013 
            endif 
          if (ltypo.ne.ltyp) then 
            write (lun8,*)kdesc(ltyp) 
            write (lun8,*)krdesc(lrtyp) 
            endif 
          call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,       &
     &                 kdat,ml,lprid,lun11)          
          call dprinto2(ltyp,lrtyp,lcon,                                &
     &        lrdat,rdat,lidat,idat,lkdat,kdat,lun8)                    
          mlt=mlt+1 
          if (mlt.lt.np2) go to 5013 
!
        deallocate(idat)
        deallocate(rdat)
        deallocate(kdat)
!
        endif 
!                                                                       
!     print out line list                                               
      if (linst.eq.21) then 
        write (lun11,*)'call pprint(8)'
       endif 
                                                                        
!                                                                       
!     count lines                                                       
      if (linst.eq.24) then 
!     nb now with min and max wavelength                                
      emin=0.1 
      emax=12.39854 
!      emax=0.9e+10                                                     
      ntmp=0 
      do mm=1,nni 
         nlnsv(mm)=0 
         enddo 
      do jlk=1,nlsvn 
         j=jlk 
         ml=derivedpointers%nplin(j) 
         mlline=ml 
         mlm=ml-1 
         call drd(ltyp,lrtyp,lcon,                                      &
     &     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                           &
     &     0,lun11)                                               
         elin=abs(masterdata%rdat1(np1r)) 
         if ((lrtyp.ne.14).and.(abs(elin).gt.0.1)                       &
     &     .and.(elin.ge.emin).and.(elin.le.emax)) then                 
           nilino=nilin 
           nilin=derivedpointers%npar(ml) 
           ergsev=1.602197e-12 
           ener=ergsev*(12398.54)/max(elin,1.d-24) 
           etst=ener/ergsev 
           idest1=masterdata%idat1(np1i) 
           idest2=masterdata%idat1(np1i+1) 
           llo=idest1 
           lup=idest2 
           aij=masterdata%rdat1(np1r+2) 
           if (ltyp.eq.82) aij=masterdata%rdat1(np1r+3) 
           elin=abs(masterdata%rdat1(np1r)) 
!           write (lun11,*)'aij,elin:',aij,elin                         
!                                                                       
!          find ion data                                                
           ml=nilin 
           mlm=ml-1 
           call drd(ltyp,lrtyp,lcon,                                    &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
           jkion=masterdata%idat1(np1i+nidt-1) 
           ntmp=ntmp+1 
           if (nilino.ne.nilin) then 
              nlnsv(jkion)=ntmp 
              ntmp=0 
              endif 
           endif 
         enddo 
         do mm=1,nni 
           write (6,*)mm,nlnsv(mm) 
           enddo 
!        this integer is the isosequence                                
         do miso=1,30 
!          this integer is the ion number                               
           m4=0 
!          this integer is the element                                  
           do mel=1,30 
             msvtmp(mel)=0 
!            count up to the element                                    
             do m3=1,mel 
               m4=m4+1 
               miso2=mel-m3+1 
               if (miso2.eq.miso) then 
                 msvtmp(mel)=nlnsv(m4) 
                 endif 
               enddo 
             enddo 
           write (6,9011)(msvtmp(mkk),mkk=1,30) 
 9011      format (1x,30i7) 
           enddo 
       endif 
!                                                                       
!                                                                       
!     print out levels                                                  
      if (linst.eq.22) then 
                                                                        
!       First look for element data (jk is element index)               
        klel=11 
        lpril=0 
        mlel=derivedpointers%npfirst(klel) 
        jk=0 
!                                                                       
!       step through elements                                           
        do while (mlel.ne.0) 
!                                                                       
!         get element data                                              
          jk=jk+1 
          mt2=mlel-1 
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,                          &
     &      0,lun11)                                              
          mllel=masterdata%idat1(np1i+nidt-1) 
          nnz=masterdata%idat1(np1i) 
          xeltp=masterdata%rdat1(np1r) 
          if (lpril.ne.0)                                               &
     &        write (lun11,*)'element:',jk,mlel,mllel,nnz,              &
     &           (masterdata%kdat1(np1k-1+mm),mm=1,nkdt),xeltp            
!                                                                       
!         now step thru ions (jkk is ion index)                         
          klion=12 
          mlion=derivedpointers%npfirst(klion) 
          jkk=0 
          kl=0 
          do while ((mlion.ne.0).and.(kl.lt.nnz)) 
!                                                                       
              jkk=jkk+1 
!             retrieve ion name from kdati                              
              mlm=mlion-1 
              call drd(ltyp,lrtyp,lcon,                                 &
     &            nrdt,np1r,nidt,np1i,nkdti,np1ki,mlm,                  &
     &            0,lun11)                                        
!                                                                       
!             if not accessing the same element, skip to the next elemen
              mlleltp=masterdata%idat1(np1i+nidt-2) 
!              write (6,*)mlleltp,mllel,mlm,jkk                         
              if (mlleltp.eq.mllel) then 
!                                                                       
                kl=kl+1 
                if (lpril.ne.0)                                         &
     &            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,         &
     &                    (masterdata%kdat1(np1ki+mm-1),mm=1,nkdti)            
                do ktt=1,min(8,nkdti) 
                  write (kinam1(ktt:ktt),'(a1)')                        &
     &               masterdata%kdat1(np1ki-1+ktt) 
                  enddo 
                do ktt=nkdti+1,9 
                  write (kinam1(ktt:ktt),'(a1)')kblnk 
                  enddo 
!                                                                       
!               get level data                                          
                t=10. 
                xee=1. 
                xpx=1. 
                lpril=0 
                call calc_rates_level_lte(jkk,lprid,lun11,t,xee,xpx,    &
     &              nlev)
!                                                                       
!               step thru levels                                        
                do mm2=1,nlev 
!                                                                       
!                 get level pointer                                     
                  mmtmp=derivedpointers%npilev(mm2,jkk) 
                  if (mmtmp.ne.0) then 
!                                                                       
!                   get data                                            
                    write (lun11,9296)mmtmp,kinam1,                     &
     &                   (leveltemp%klev(lk,mm2),lk=1,50),              &
     &                   (leveltemp%rlev(lk,mm2),lk=1,4),               &
     &                   (leveltemp%ilev(lk,mm2),lk=1,5)              
 9296               format (1x,i6,1x,a8,1x,(50a1),4(1pe13.5),5i6) 
!                                                                       
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
!         Go to next element                                            
          if (mlel.ne.0) mlel=derivedpointers%npnxt(mlel) 
          enddo 
!                                                                       
        write (lun11,*)'done with linst=22' 
        endif 
!                                                                       
!     split up records into files                                       
      if (linst.eq.14) then 
        write (6,*)'linst=14',np2 
        mlion=0 
        mlrt=0 
        lun9=9 
        ml=0 
        mlt=1 
        if (nptrt(mlt).ge.np2) mlt=mlt+1 
 1113       continue 
            ml=nptrt(mlt) 
            ltyp=masterdata%nptrs(2,ml) 
            lprid=0 
            mliono=mlion 
            mlrto=mlrt 
!            write (lun11,*)'before dread',ltyp,lrtyp,ml,mlt            
            call drd(ltyp,lrtyp,lcon,                                   &
     &       nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                         &
     &       0,lun11)                                             
!            call dread(ltyp,lrtyp,lcon,                                
!     $        lrdat,rdat,lidat,idat,lkdat,kdat,ml-1,                   
!     $        idat1,rdat1,kdat1,nptrs,lprid,lun11)                     
!            write (lun11,*)'after dread',ltyp,lrtyp,ml,lrdat,          
!     $        lidat,lkdat                                              
!            call dprints(ltyp,lrtyp,lcon,                              
!     $        lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                  
            mlion=masterdata%idat1(np1i+nidt-1) 
            mlrt=ltyp 
            if ((mlion.ne.mliono).or.(mlrt.ne.mlrto)) then 
              close(unit=lun9) 
              write (klab1(1:1),'(a1)')kd 
              mlionr=mlion 
              if (mlion.lt.100.) then 
                  write (klab1(2:2),'(a1)')ku 
                else 
                  mlionm=int(mlion/100) 
                  mlionr=mod(mlion,100) 
                  write (klab1(2:2),'(i1)')mlionm 
                endif 
              if (mlion.lt.10.) then 
                  write (klab1(3:3),'(a1)')ku 
                else 
                  mlionm=int(mlionr/10) 
                  mlionr=mod(mlionr,10) 
                  write (klab1(3:3),'(i1)')mlionm 
                endif 
              write (klab1(4:4),'(i1)')mlionr 
              mlrtr=mlrt 
              write (klab1(5:5),'(a1)')kt 
              if (mlrt.lt.100) then 
                  write (klab1(6:6),'(a1)')ku 
                else 
                  mlrtm=int(mlrt/100) 
                  mlrtr=mod(mlrt,100) 
                  write (klab1(6:6),'(i1)')mlrtm 
                endif 
              if (mlrt.lt.10.) then 
                  write (klab1(7:7),'(a1)')ku 
                else 
                  mlrtm=int(mlrtr/10) 
                  mlrtr=mod(mlrtr,10) 
                  write (klab1(7:7),'(i1)')mlrtm 
                endif 
              write (klab1(8:8),'(i1)')mlrtr 
              do lkm=1,8 
                lk=9-lkm 
                read (klab1(lk:lk),'(a1)')ktmp 
                write (klab1(lk+4:lk+4),'(a1)')ktmp 
                enddo 
              write (klab1(1:4),'(a4)')'tmp/' 
              inquire (file=klab1,exist=ex) 
              if (ex) then 
              open(unit=lun9,file=klab1,access='APPEND') 
              write (6,*)'reopening unit:',klab1 
              else 
              open(unit=lun9,file=klab1) 
              write (6,*)'writing to unit:',klab1 
              endif 
              endif 
            if (ltyp.ne.13)                                             &
     &        call dprinto(ltyp,lrtyp,lcon,                             &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,lun9)   
!     $       call dprinto2(ltyp,lrtyp,lcon,                            
!     $        lrdat,rdat,lidat,idat,lkdat,kdat,lun9)                   
            mlt=mlt+1 
            if (mlt.lt.np2) go to 1113 
        endif 
!                                                                       
!     print out after sort, short print                                 
      if (linst.eq.12) then 
!
        allocate(idat(2000000))
        allocate(rdat(2000000))
        allocate(kdat(2000000))
!
        write (lun11,*)'linst=12',np2 
        mlt=1 
 6013     continue 
          ml=nptrt(mlt) 
!          write (lun11,*)mlt,ml 
          ltyp=masterdata%nptrs(2,ml) 
          call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,       &
     &                 kdat,ml,lprid,lun11)          
          call dprints2(ml,ltyp,lrtyp,lcon,                             &
     &        lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                   
          mlt=mlt+1 
          if (mlt.lt.np2) go to 6013 
!
        deallocate(idat)
        deallocate(rdat)
        deallocate(kdat)
!
        endif 
!                                                                       
!     print pointers
      if (linst.eq.25) then 
!
        allocate(idat(2000000))
        allocate(rdat(2000000))
        allocate(kdat(2000000))
!
        write (lun11,*)'linst=25',np2 
        write (lun11,*)'records pointers:',np2 
        mlt=1 
 6014     continue 
          ml=nptrt(mlt) 
!          write (lun11,*)mlt,ml 
          ltyp=masterdata%nptrs(2,ml) 
          call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,       &
     &                 kdat,ml,lprid,lun11)          
          call dprints2(ml,ltyp,lrtyp,lcon,                             &
     &        lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                   
          write (lun11,*)'ml=',ml,' npar(ml)=',derivedpointers%npar(ml),&
     &       ' npnxt(ml)=', derivedpointers%npnxt(ml),' nplini(ml)=',   &
     &        derivedpointers%nplini(ml)
          mlt=mlt+1 
          if (mlt.lt.np2) go to 6014
!
!       step thru elements
        klel=11 
        mlel=derivedpointers%npfirst(klel) 
        do while (mlel.ne.0) 
          call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,       &
     &                 kdat,mlel,lprid,lun11)          
          call dprints2(ml,ltyp,lrtyp,lcon,                             &
     &        lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                   
          write (lun11,*)'element=',idat(1),' npnxt(mlel)=',            &
     &      derivedpointers%npnxt(mlel)
          mlel=derivedpointers%npnxt(mlel) 
          enddo 
!
!       step thru ions
        klel=11 
        mlel=derivedpointers%npfirst(klel) 
        do while (mlel.ne.0) 
          call drd(ltyp,lrtyp,lcon,                                     &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlel,                       &
     &        0,lun11)                                            
          mllel=masterdata%idat1(np1i+nidt-1) 
          nnz=masterdata%idat1(np1i) 
          klion=12 
          mlion=derivedpointers%npfirst(klion) 
          do while (mlion.ne.0) 
            call drd(ltyp,lrtyp,lcon,                                   &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,mlion,                   &
     &            0,lun11)                                        
!           if not accessing the same element, skip to the next elemen
            mlleltp=masterdata%idat1(np1i+nidt-2) 
            if (mlleltp.eq.mllel) then 
              call dread(ltyp,lrtyp,lcon,nrdt,rdat,nidt,idat,nkdt,      &
     &                 kdat,mlion,lprid,lun11)          
              call dprints2(mlion,ltyp,lrtyp,lcon,                      &
     &        nrdt,rdat,nidt,idat,nkdt,kdat,lun11)                   
              jkk=idat(nidt)
              write (lun11,*)'mlion=',mlion,' npnxt(mlion)=',           &
     &         derivedpointers%npnxt(mlion),lidat
!             step thru rate types
              write (lun11,*)'npfi pointers:'
              do mm=1,ntyp
                if (derivedpointers%npfi(mm,jkk).ne.0) then
                  write (lun11,*)'jkk=',jkk,' type=',mm,' npfi=',       &
     &                derivedpointers%npfi(mm,jkk)
                  endif
                enddo
              write (lun11,*)'number of levels:',                       &
     &          derivedpointers%nlevs(jkk)
              endif
!
!           Go to next ion                                            
            mlion=derivedpointers%npnxt(mlion) 
            enddo 
!
          mlel=derivedpointers%npnxt(mlel) 
          enddo 
!          write (lun11,*)'ion pointers:',nni
          
!
        deallocate(idat)
        deallocate(rdat)
        deallocate(kdat)
!
        endif 
!                                                                       
!                                                                       
!       print out records sorted by ion                                 
        if (linst.eq.8) then 
!
        allocate(idat(2000000))
        allocate(rdat(2000000))
        allocate(kdat(2000000))
!
        do mli=1,nni 
          luse(mli)=0 
          enddo 
        mlel1=derivedpointers%npfirst(11) 
 3099     continue 
          call dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,         &
     &                 kdat,mlel1,lprid,lun11)          
          call dprinto2(ltyp,lrtyp2,lcon,                               &
     &               lrdat,rdat,lidat,idat,lkdat,kdat,lun11)            
          do mli=1,nni 
!            write (lun11,*)'mli,nion(mli)',mli2,nimap(mli2),           
!     $                             nion(mli),mli,npfirst(12)           
            ml12=derivedpointers%npfirst(12) 
            if (ml12.ne.0) then 
 3098         continue 
              mlel2=derivedpointers%npar(ml12) 
!              write (lun11,*)mlel2,ml12                                
              call dread(ltyp,lrtyp2,lcon,lrdat,rdat,lidat,idat,lkdat,  &
     &                 kdat,mlel2,lprid,lun11)          
!              write (lun11,*)mlel2,mlel1,mli,idat(lidat),ml12          
              if ((mlel2.eq.mlel1).and.(mli.eq.idat(lidat))             &
     &              .and.(luse(mli).ne.1)) then                         
                luse(mli)=1 
                call dprinto2(ltyp,lrtyp2,lcon,                         &
     &               lrdat,rdat,lidat,idat,lkdat,kdat,lun11)            
                do mlk=1,13 
                  if ((mlk.ne.12).and.(mlk.ne.11)) then 
                    ml=derivedpointers%npfirst(mlk) 
                    if (ml.ne.0) then 
 3097                 continue 
                      call dread(ltyp,lrtyp2,lcon,lrdat,                &
     &                 rdat,lidat,idat,lkdat,                           &
     &                 kdat,ml,lprid,lun11)          
                      if (mli.eq.idat(lidat))                           &
     &                 call dprinto2(ltyp,lrtyp2,lcon,                  &
     &                 lrdat,rdat,lidat,idat,lkdat,kdat,lun11)          
                      ml=derivedpointers%npnxt(ml) 
                      if (ml.ne.0) go to 3097 
                    endif 
                  endif 
                enddo 
              endif 
              ml12=derivedpointers%npnxt(ml12) 
              if (ml12.ne.0) go to 3098 
            endif 
          enddo 
          mlel1=derivedpointers%npnxt(mlel1) 
          if (mlel1.ne.0) go to 3099 
!
        deallocate(idat)
        deallocate(rdat)
        deallocate(kdat)
!
        endif 
!                                                                       
!       print out records of a given type                               
        if (linst.eq.4) then 
!
          allocate(idat(2000000))
          allocate(rdat(2000000))
          allocate(kdat(2000000))
!
          do mll=1,ntyp 
            ml=derivedpointers%npfirst(mll) 
            ltyp=0 
            if (ml.ne.0) then 
              write (lun11,*)'rate type=',mll,' ',krdesc(mll),ml 
              lrtyp=masterdata%nptrs(3,ml) 
              do mll2=1,ntyp 
                ml2=derivedpointers%npfirst(mll2) 
                if (ml2.ne.0) then 
 1093             continue 
                  ltypo=ltyp 
                  ltyp=masterdata%nptrs(2,ml2) 
                  lrtyp2=masterdata%nptrs(3,ml2) 
                  if (lrtyp2.eq.lrtyp) then 
                    if (ltyp.ne.ltypo)                                  &
     &               write (lun11,*)'data type=',mll2,' ',              &
     &                       kdesc(mll2),ml2                            
                    call dread(ltyp,lrtyp2,lcon,lrdat,                  &
     &                 rdat,lidat,idat,lkdat,                           &
     &                 kdat,ml2,lprid,lun11)          
                    call dprinto2(ltyp,lrtyp2,lcon,                     &
     &               lrdat,rdat,lidat,idat,lkdat,kdat,lun11)            
                  endif 
                  ml2=derivedpointers%npnxt(ml2) 
                  if (ml2.ne.0) go to 1093 
                endif 
              enddo 
            endif 
          enddo 
!
        deallocate(idat)
        deallocate(rdat)
        deallocate(kdat)
!
        endif 
!                                                                       
!                                                                       
!       print out all records                                           
        if (linst.eq.6) then 
!
            allocate(idat(2000000))
            allocate(rdat(2000000))
            allocate(kdat(2000000))

!           write (lun11,*)'linst=6',np2                                
            lun8=8 
            mlt=1 
 1013       continue 
             ml=nptrt(mlt) 
!             write (lun11,981)ml,                                      
!     $         (nptrs(jj,ml),jj=1,10)                                  
!  981        format (1x,'ml,nptrs(3,ml):',11i6) 
             ltyp=masterdata%nptrs(2,ml) 
             lprid=0 
             call dread(ltyp,lrtyp2,lcon,lrdat,                         &
     &                 rdat,lidat,idat,lkdat,                           &
     &                 kdat,ml,lprid,lun11)          
!             write (lun11,*)'ml=',ml                                   
!             write (lun11,*)ltyp,lrtyp,lcon,                           
!     $        lrdat,rdat,lidat,idat,lkdat,kdat,nptrs(3,ml),ml          
             call dprinto2(ltyp,lrtyp,lcon,                             &
     &        lrdat,rdat,lidat,idat,lkdat,kdat,lun8)                    
             mlt=mlt+1 
             if (mlt.lt.np2) go to 1013 
!
           deallocate(idat)
           deallocate(rdat)
           deallocate(kdat)
!
           endif 
!                                                                       
!                                                                       
!       set up pointers
        if (linst.eq.7) then 
!            write (lun11,*)'linst=7 set up pointers',np2              
!
            allocate(idat(2000000))
            allocate(rdat(2000000))
            allocate(kdat(2000000))
!                                                                       
            do 1891 kk=1,ntyp 
               mlold(kk)=0 
 1891          continue 
!                                                                       
!                                                                       
!            pointer structure                                          
!     type    desc         nr  ni  nk      daught  par                  
!     1       rr, a&p      2   1   0               14                   
!     2       hcx          4   1   0               14                   
!     3       ai           2   1   0               14                   
!     4       line dat 1   2   3   0        5      14                   
!     5       line dat 2   4   3   0                4                   
!     6       lev dat  1   4   3   0               14                   
!     7       dr a&p       5   1   0               14                   
!     8       dr a&r       0   0   0               14                   
!     9       hecx         4   1   0               14                   
!     10      lev dat 2    0   2  30                6                   
!     11      2 ph         2   2   0               14                   
!     12      pixc, bpl    5   2   0               14                   
!     13      el           2   2  30       14       0                   
!     14      ion          1   2   8       all     13                   
!     15      pixc bkh 1   5   1   0       20      14                   
!     16      pixc bkh     0   0   0               14                   
!     17      cx: cota     4   3   0               14                   
!     18      rr: cota     3   1   0               14                   
!     19      pixc hullac  0   0   0               14                   
!     20      pixc bkh 2   5   1   0       21      15                   
!     21      pixc bkh 3   4   4  11               20                   
!     22      dr stroey    5   1   0               14                   
!     23      pixc clark   5   2   0       24      14                   
!     24      pixc clark 2 4   4   0               23                   
!     25      ci r&s       0   0   0               14                   
!     26      ci cota      2   2   0               14                   
!                                                                       
           lprisv=lpri 
!                                                                       
!          first zero the pointers                                      
           do ml=1,ntyp 
             derivedpointers%npfirst(ml)=0 
             do ll=1,nni 
               derivedpointers%npfi(ml,ll)=0 
               enddo 
             enddo 
           do ml=1,nl 
             nemap(ml)=0 
             nel(ml)=0 
             enddo 
           do ml=1,nni 
             nion(ml)=0 
             nimap(ml)=0 
             enddo 
           do  ml=1,np2 
              derivedpointers%npar(ml)=0 
              derivedpointers%npnxt(ml)=0 
              enddo 
!          first step through and find all the elements                 
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'the element pointers:'                      
           nenxt=0 
           mlo=0 
           do  ml=1,np2 
             lrtyp=masterdata%nptrs(3,ml) 
             if (lrtyp.eq.11) then 
!               write (lun11,*)'ml,nptrs(3,ml):',ml,                     &
!     &            masterdata%nptrs(3,ml)         
               call dread(ltyp,lrtyp,lcon,                              &
     &          lrdat,rdat,lidat,idat,lkdat,kdat,ml,                    &
     &          0,lun11)                        
!               call dprinto2(ltyp,lrtyp,lcon,                           &
!     &           lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                    
               nenxt2=idat(lidat) 
               nenxt=nenxt+1 
               if (abel(nenxt2).gt.1.e-15) then                        
                 if (derivedpointers%npfirst(lrtyp).eq.0) then 
                    derivedpointers%npfirst(lrtyp)=ml 
                    else 
                    derivedpointers%npnxt(mlo)=ml 
                    endif 
                 mlo=ml 
                 nemap(nenxt2)=nenxt 
                 if (lpri.ne.0)                                         &
     &            write (lun11,*)ml,nenxt,mlo                           
                 nel(nenxt)=ml 
                 endif                                                 
               endif 
             enddo 
!          next step through and put in the ion pointers:               
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'the ion pointers:'                          
           ninxt=0 
           do ml=1,np2 
             lrtyp=masterdata%nptrs(3,ml) 
!             write (lun11,*)ml,masterdata%nptrs(1,ml),                  &
!     &         masterdata%nptrs(2,ml),masterdata%nptrs(3,ml),           &
!     &         masterdata%nptrs(4,ml)
             if (lrtyp.eq.12) then 
               call dread(ltyp,lrtyp,lcon,                              &
     &          lrdat,rdat,lidat,idat,lkdat,kdat,ml,                    &
     &          0,lun11)                        
!               call dprinto2(ltyp,lrtyp,lcon,                           &
!     &           lrdat,rdat,lidat,idat,lkdat,kdat,lun11)                    
               ninxt2=idat(3) 
!               write (lun11,*)idat(2),idat(3)                          
               netmp=nemap(idat(2)) 
!               write (lun11,*)netmp                                    
               if ((netmp.gt.0).and.(netmp.le.nl)) then 
                 ninxt=ninxt+1 
                 nion(ninxt)=ml 
                 nimap(ninxt2)=ninxt 
                 derivedpointers%npar(ml)=nel(netmp) 
                 if (lpri.ne.0)                                         &
     &            write (lun11,*)ml,ninxt,netmp,nel(netmp)              
                 endif 
               endif 
             enddo 
            if (lpri.ne.0) then 
              write (lun11,*)'the ion map:' 
              do mm1=1,168 
                 write (lun11,*)mm1,nimap(mm1) 
                 enddo 
              endif 
!          next step through and put in pointers for others             
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'the other pointers:'                        
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'np2=',np2                                   
           do ml=1,np2 
             lrtyp=masterdata%nptrs(3,ml) 
             if ((lrtyp.ne.11).and.(lrtyp.ne.12).and.(lrtyp.gt.0)) then 
               call dread(ltyp,lrtyp,lcon,                              &
     &          lrdat,rdat,lidat,idat,lkdat,kdat,ml,                    &
     &          0,lun11)                        
               if (lidat.ne.0) then 
                 ninxt=idat(lidat) 
!                 write (lun11,*)ml,lrtyp,ninxt                         
                 if (ninxt.gt.0) then 
!                   write (lun11,*)nimap(ninxt)                         
                   if (nimap(ninxt).gt.0) then 
!                     write (lun11,*)nion(nimap(ninxt))                 
                     derivedpointers%npar(ml)=nion(nimap(ninxt)) 
                     if (lpri.ne.0)                                     &
     &                write (lun11,*)ml,ltyp,lrtyp,lidat,ninxt,         &
     &                nion(nimap(ninxt)),derivedpointers%npar(ml)                       
                     endif 
                   endif 
                 endif 
               endif 
             enddo 
!          now the next pointers                                        
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'filling next pointers'                      
           ml=0 
 1289        ml=ml+1 
             mlo=ml 
             lcon=masterdata%nptrs(4,ml) 
             if (lcon.ne.0) then 
 1292               continue 
                    ml=ml+1 
                    lcon=masterdata%nptrs(4,ml) 
                    if ((lcon.ne.0).and.(ml.lt.np2)) go to 1292 
                  endif 
             lrtyp=masterdata%nptrs(3,mlo) 
             if ((lrtyp.ne.0).and.(derivedpointers%npar(ml).ne.0))      &
     &             then                                                 
               mltyp=mlold(lrtyp) 
               if (lpri.ne.0)                                           &
     &          write (lun11,*)ml,mlo,lrtyp,mltyp,masterdata%nptrs(4,ml)           
               if (mltyp.ne.0) then 
                 derivedpointers%npnxt(mltyp)=mlo 
                 endif 
               mlold(lrtyp)=mlo 
               endif 
             if (ml.lt.np2) go to 1289 
!          now the first pointers                                       
           do kl=1,ntyp 
             if (kl.ne.11)                                              &
     &        derivedpointers%npfirst(kl)=0                                             
             enddo 
           ml=0 
 1023        ml=ml+1 
             lrtyp=masterdata%nptrs(3,ml) 
             if ((lrtyp.ne.0).and.(derivedpointers%npar(ml).ne.0))      &
     &             then                                                 
               if (derivedpointers%npfirst(lrtyp).eq.0)                 &
     &             derivedpointers%npfirst(lrtyp)=ml         
               endif 
             if (ml.lt.np2) go to 1023 
           if (lpri.eq.0) go to 9812 
           write (lun11,*)'the next pointers' 
           do ll=1,ntyp 
             ml=derivedpointers%npfirst(ll) 
             write (lun11,*)'type=',ll,ml 
             if (ml.ne.0) then 
 2022            continue 
                 mln=derivedpointers%npnxt(ml) 
                 write (lun11,*)ml,mln 
                 ml=mln 
                 if (mln.ne.0) go to 2022 
               endif 
             enddo 
 9812        continue 
!          now the next ion pointers                                    
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'the next ion pointers'                      
           jkk=0 
           do jkkl=1,nni 
             jkk2=nimap(jkkl) 
             if (jkk2.ne.0) then 
             jkk=jkk+1 
!             jkk=jkk2                                                  
             do kk=1,ntyp 
                derivedpointers%npfi(kk,jkk)=0 
                enddo 
             mlion=nion(jkk) 
             if (lpri.ne.0)                                             &
     &        write (lun11,*)'jkk=',jkk,nimap(jkkl),mlion               
             do  kk=1,ntyp 
               mlold(kk)=0 
               enddo 
!             write (lun11,*)npar(1032)                                 
             do mlt=1,ntyp 
               ml=derivedpointers%npfirst(mlt) 
               if (ml.ne.0) then 
 2023            continue 
                 lrtyp=masterdata%nptrs(3,ml) 
                 lcon=masterdata%nptrs(4,ml) 
!                 write (lun11,*)ml,lrtyp,mlion,npar(ml)                
!                 write (lun11,*)jkk,mlt                                
                 if ((mlion.eq.derivedpointers%npar(ml))                &
     &             .and.(lrtyp.ne.0)) then 
                   mltyp=mlold(lrtyp) 
                   if (mltyp.gt.0) then 
!                     npni(mltyp,jkk)=ml                                
!                     write (lun11,*)mltyp,ml,jkk                       
                   else 
                     derivedpointers%npfi(lrtyp,jkk)=ml 
                   endif 
                   mlold(lrtyp)=ml 
                   endif 
                 ml=derivedpointers%npnxt(ml) 
                 if (ml.ne.0) go to 2023 
                 endif 
               enddo 
               if (lpri.ne.0) then 
                 do kk=1,ntyp 
                   if (derivedpointers%npfi(kk,jkk).ne.0) then 
                     write (lun11,*)jkk,kk,derivedpointers%npfi(kk,jkk) 
                     ml=derivedpointers%npfi(kk,jkk) 
                     mlz=derivedpointers%npar(ml) 
 1281                  continue 
                       ml=derivedpointers%npnxt(ml) 
                       write (lun11,*)'   ',ml 
                       if ((ml.ne.0).and.                               &
     &                  (mlz.eq.derivedpointers%npar(ml))) go to 1281 
                     endif 
                   enddo 
                 endif 
             endif 
             enddo 
!          now the line pointers                                        
           jkkl=0 
           jkk=0 
          if (lpri.ne.0)                                                &
     &      write (lun11,*)'the line pointers'                          
           do mm=1,nni 
             nlines(mm)=0 
             derivedpointers%nlevs(mm)=0 
             enddo 
           do  ml=1,np2 
             if ((derivedpointers%npar(ml).ne.0).and.                   &
     &           (derivedpointers%npar(ml).le.np2)) then                      
               ltyp=masterdata%nptrs(2,ml) 
               lrtyp=masterdata%nptrs(3,ml) 
               if ((lrtyp.eq.4).or.(lrtyp.eq.14).or.(lrtyp.eq.9)) then 
                 jkkl=jkkl+1 
                 derivedpointers%nplin(jkkl)=ml 
                 derivedpointers%nplini(ml)=jkkl 
                 mlpar=derivedpointers%npar(ml) 
                 call dread(ltyp,lrtyp,lcon,                            &
     &            lrdat,rdat,lidat,idat,lkdat,kdat,mlpar,               &
     &            0,lun11)                      
                 ilv=idat(lidat) 
                 if ((ilv.gt.0).and.(ilv.le.nni))                       &
     &             nlines(ilv)=nlines(ilv)+1                              
                 if (lpri.ne.0)                                         &
     &            write (lun11,*)jkkl,ml                                &
     &            ,derivedpointers%nplin(jkkl),derivedpointers%npar(ml)             
                 endif
               endif 
             enddo 
           nlsvn=jkkl 
!          now the continuum pointers                                   
!          note that problems may occur if pi xsections aren't ordered  
!           the same as levels.                                         
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'the continuum pointers'                     
           jkkl=0 
           ml=1 
!
!          step thru ions
           do jkk=1,nni 
!
             if (lpri.ne.0) write (lun11,*)'jkk=',jkk
!            step thru rate types which have continuum processes
             do jxx=1,2 
!
!              ml points to photoionization data record
               if (jxx.eq.1) then 
                 ml=derivedpointers%npfi(7,jkk) 
                else 
                 ml=derivedpointers%npfi(1,jkk) 
                endif 
!               if (lpri.ne.0) write (lun11,*)'jxx,ml=',jxx,ml
!
!               if (lpri.ne.0)                                           &
!     &           write (lun11,*)jkk,ml,nimap(jkk)                          
!
!              test for if npfi is ok
               if ((ml.ne.0).and.(ml.le.ndat2)) then 
!
!                mllz2 points to parent ion of photoionization data record
                 mllz2=derivedpointers%npar(ml) 
!                 if (lpri.ne.0) write (lun11,*)'ml,mllz2=',ml,mllz2
!
!                loop over pi data
                 do while ((ml.ne.0).and.(ml.le.ndat2))
!
!                  test for if parents agree
!                   if (lpri.ne.0) write (lun11,*)'ml,npar(ml)=',        &
!     &                ml,derivedpointers%npar(ml)
                   if (derivedpointers%npar(ml).eq.mllz2) then
!
!                    get index of associated level                      
                     call dread(ltyp,lrtyp,lcon,                        &
     &                lrdat,rdat,lidat,idat,lkdat,kdat,ml,              &
     &                0,lun11)                      
                     nlvtmp1=idat(lidat-1) 
!                     if (lpri.ne.0) write (lun11,*)'nlvtmp1=',nlvtmp1
                     nlvtmp2=0
!
!                    now search the level list
!                    mll points points to the level data
                     mll=derivedpointers%npfi(13,jkk) 
!
!                    test for if level exists
                     if ((mll.ne.0).and.(mll.le.ndat2)) then
!
                       mllz=derivedpointers%npar(mll) 
!                       if (lpri.ne.0) write (lun11,*)'mll,mllz2)=',     &
!     &                    mll,mllz
                       do while ((derivedpointers%npar(mll).eq.mllz)    &
     &                   .and.(mll.ne.0).and.(nlvtmp1.ne.nlvtmp2))
!
                         call dread(ltyp,lrtyp,lcon,                    &
     &                     lrdat,rdat,lidat,idat,lkdat,kdat,mll,        &
     &                     0,lun11)                    
                         nlvtmp2=idat(lidat-1) 
!
                         if (lpri.ne.0) write (lun11,*)                 &
     &                       'nlvtmp1,nlvtmp2,mll:',                    &
     &                       nlvtmp1,nlvtmp2,mll                             
!
                         if (nlvtmp1.eq.nlvtmp2) then
                           jkkl=jkkl+1 
!                          npcon points from the rrc list to the pi data
                           if (jkkl.gt.nnml) stop 'npcon overflow'
                           derivedpointers%npcon(jkkl)=ml 
!                          npconi points from the level list to the rrc list
                           if (mll.gt.ndat2) stop 'npconi overflow'
                           derivedpointers%npconi(mll)=jkkl 
!                          npconi2 points from the pi data to the rrc list
                           if (ml.gt.ndat2) stop 'npconi2 overflow'
                           derivedpointers%npconi2(ml)=jkkl 
                           mlpar=derivedpointers%npar(ml) 
                           call dread(ltyp,lrtyp,lcon,                  &
     &                       lrdat,rdat,lidat,idat,lkdat,kdat,mlpar,    &
     &                       0,lun11)                      
!                          ilv is the ion index
                           ilv=idat(lidat) 
                           if ((ilv.gt.0).and.(ilv.le.nni))             &
     &                       derivedpointers%nlevs(jkk)                 &
     &                         =max(derivedpointers%nlevs(jkk),nlvtmp2)                  
                           if (lpri.ne.0)                               &
     &                       write (lun11,*)jkkl,ml,                    &
     &                       derivedpointers%npcon(jkkl),               &
     &                       derivedpointers%npar(ml),                  &
     &                       derivedpointers%npconi2(ml),               &
     &                       derivedpointers%npconi(mll),mll,nlvtmp2                        
                           endif
!
!                        end of loop over levels
                         mll=derivedpointers%npnxt(mll) 
                         enddo
!
!                      end of test for level exists
                       endif
!
!                    end of test for if parents agree
                     endif
!
!                  end of loop over pi data
                   ml=derivedpointers%npnxt(ml) 
                   enddo
!
!                end of test if npfi is ok
                 endif
!
!              end of loop over rate types
               enddo
!
             if (lpri.ne.0) write (lun11,*)'jkk=',jkk
             derivedpointers%nlevs(jkk)=0
!            step thru levels
!            mll points points to the level data
             mll=derivedpointers%npfi(13,jkk) 
!
!            test for if level exists
             if ((mll.ne.0).and.(mll.le.ndat2)) then
!
               mllz=derivedpointers%npar(mll) 
               do while ((derivedpointers%npar(mll).eq.mllz)            &
     &                   .and.(mll.ne.0))
!
                 call dread(ltyp,lrtyp,lcon,                            &
     &                     lrdat,rdat,lidat,idat,lkdat,kdat,mll,        &
     &                     0,lun11)                    
                 nlvtmp2=idat(lidat-1) 
                 derivedpointers%nlevs(jkk)                             &
     &              =max(derivedpointers%nlevs(jkk),nlvtmp2)                  
!                 write (lun11,*)'filling nlevs:',jkk,mll,mllz,nlvtmp2,  &
!     &              derivedpointers%npar(mll)
!
!                end of loop over levels
                 mll=derivedpointers%npnxt(mll) 
                 enddo
!
!              end of test for level exists
               endif
!
!            end of loop over ions
             enddo 
!
           ncsvn=jkkl 
           if (lpri.ne.0) then 
             write (lun11,*)'ion, #lines, #levels' 
             nltot=0 
             nvtot=0 
             do mm=1,nni 
               write (lun11,*)mm,nlines(mm),                            &
     &              derivedpointers%nlevs(mm) 
               nltot=nltot+nlines(mm) 
               nvtot=nvtot+derivedpointers%nlevs(mm) 
               enddo 
             write (lun11,*)'totals:',nltot,nvtot 
             endif 
!          now the ion level pointers                                   
           if (lpri.ne.0)                                               &
     &      write (lun11,*)'the ion level pointers'                     
           jkkl=0 
           do jkk=1,nni 
             ml=derivedpointers%npfi(13,jkk) 
             if (ml.ne.0) then 
             mlz=derivedpointers%npar(ml) 
             if (lpri.ne.0)                                             &
     &        write (lun11,*)jkk,ml,nimap(jkk)                          
             kkl=0 
 2235        continue 
             if (ml.le.0) go to 2234 
               if (derivedpointers%npar(ml).ne.mlz) go to 2234 
               kkl=kkl+1 
               jkkl=jkkl+1 
               call dread(ltyp,lrtyp,lcon,                              &
     &          nrdt,rdat,nidt,idat,nkdt,kdat,ml,                       &
     &          0,lun11)                        
               derivedpointers%npilevi(jkkl)=kkl 
               derivedpointers%npilev(kkl,jkk)=jkkl 
               if (lpri.ne.0)                                           &
     &            write (lun11,*)jkkl,ml,idat(nidt-1),kkl,jkk           
               ml=derivedpointers%npnxt(ml) 
               go to 2235 
 2234          continue 
             endif 
             enddo 
          if (lpri.ne.0)                                                &
     &     write (lun11,*)'nlsvn=',nlsvn,', ncsvn=',ncsvn               
!         print out the pointers                                        
          if (lpri.ne.0) then 
            write (lun11,*)'the pointers' 
            do ml=1,np2 
               if (masterdata%nptrs(3,ml).ne.0)                         &
     &         write (lun11,*)ml,masterdata%nptrs(3,ml),                &
     &           derivedpointers%npar(ml),derivedpointers%npnxt(ml) 
               enddo 
            endif 
!
          deallocate(idat)
          deallocate(rdat)
          deallocate(kdat)
!
           endif 
!                                                                       
!                                                                       
      lpri=lprisv 
!                                                                       
!                                                                       
      return 
      END                                           
