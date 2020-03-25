      subroutine setptrs(lun11,lpri,                                    &
     &       np2,ncsvn,nlsvn,                                           &
     &       abcosmic,abel)                         
!                                                                       
!     Name: setptrs.f90  
!     Description:  
!           sets up atomic database pointers
!
!     List of Parameters:
!     Input:
!           lun11: logical unit number for printing
!           lpri: print switch, 1=on, 0=off
!           abel(nl):  element abundances relative to H=1
!      Output
!           np2: atomic data parameter, number of records in atomic database
!           ncsvn: atomic data parameter, number of rrcs in atomic database
!           nlsvn: atomic data parameter, number of lines in atomic database
!           abcosmic(nl):  cosmic element abundances relative to H=1
!                taken from database
!
!     Dependencies:  none
!     Called by:  xstarsetup
!
!     See below for description of database quantities:
!
!     this program set the pointers of the database                     
!       Written by Ke Zhang, Oct.8 2001                                 
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
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      integer, dimension(:), allocatable :: nptrt
      integer, dimension(:), allocatable :: npnxt2
      character(4) karg(20)
      real(8) abcosmic(30)
      real(8) abel(nl) 
      integer melpt(nl) 
      integer mlold(ntyp) 
      integer indx, iion, ilev, icon, iline, i, j 
      integer iel2, iel, lpri, lun11, np2, lrtp 
      integer iilev, mltmpn, mlfnd, nclev, mltst 
      integer mltmp, npartmpn, nlsvn, ncsvn 
      integer lsrt, mml, niter, melptmp, npfirst2, mm
      integer mllo, mlloo, itst, ltyp, lrtyp2, lcon 
      integer nrdt, nidt, nkdt, mll,mlm 
      integer itmp,np1i,np1r,np1k,mm1,lrtyp,mlel,mlion
!                                                                       
      allocate(nptrt(ndat2))
      allocate(npnxt2(ndat2))
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
        indx=1 
!                                                                       
! the main data index                                                   
!      go to 9009 
      if (lpri.ge.1) write (lun11,*)'in setptrs'
!
      if (lpri.gt.2) then 
!       first an experimental print                                     
        write (lun11,*)'np2=',np2 
        do itmp=1,np2 
          CALL DRD(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,      &
     &          itmp,0,Lun11)                                   
          write (lun11,*)'itmp=',itmp 
!          write (lun11,*)'nkdt=',nkdt,(kdat1(np1k-1+mm),mm=1,nkdt)     
          call dprints(ltyp,lrtyp2,lcon,                             &
     &    nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)        
          enddo 
        endif 
! 9009   continue 
                                                                        
! the ion index                                                         
      iion=1 
                                                                        
! the level index                                                       
      ilev=1 
                                                                        
! the continum index                                                    
      icon=1 
                                                                        
! the line index                                                        
      iline=1 
                                                                        
! initialize pointers                                                   
                                                                        
      do i=1,ntyp 
        derivedpointers%npfirst(i)=0 
      enddo 
                                                                        
      do i=1,nni 
        do j=1,ntyp 
          derivedpointers%npfi(j,i)=0 
        enddo 
        derivedpointers%nlevs(i)=0 
      enddo 
                                                                        
      do i=1,ndat2 
        derivedpointers%npnxt(i)=0 
      enddo 
!                                                                       
      do i=1,nl 
        abcosmic(i)=0. 
        enddo 
                                                                        
      mlold(11)=0 
      iel2=0 
      do while ((iel2.le.nl).and.(indx.lt.np2)) 
        iel2=iel2+1 
        iel=iel2 
        if (abel(iel).lt.1.e-15) then                                  
!                                                                       
          if (lpri.gt.1)                                                &
     &     write (lun11,*)'iel=',iel2,iel,abel(iel)                     
!  pass by elements that has neglectable abundance                      
          indx=indx+1                                                  
          do while((masterdata%nptrs(3,indx).ne.11).and.(indx.lt.np2)) 
            indx=indx+1                                                
          enddo                                                        
          iion=iion+iel
!                                                                       
        else                                                           
                                                                        
!  register element record                                              
!     npfirst,npnxt,npar,mlold                                          
                                                                        
          if (lpri.gt.0) write (lun11,*)'npfirst(11):',                 &
     &       derivedpointers%npfirst(11),                               &
     &          indx,mlold(11)                                         
          if (derivedpointers%npfirst(11).eq.0) then 
            derivedpointers%npfirst(11)=indx 
          else 
            derivedpointers%npnxt(mlold(11))=indx 
          endif 
!                                                                       
          CALL DRD(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,      &
     &          indx,0,Lun11)                                   
          iel=masterdata%idat1(np1i) 
          abcosmic(iel)=masterdata%rdat1(np1r) 
                                                                        
          if (lpri.gt.1)                                                &
     &     write (lun11,*)'registering element:',iel,abel(iel),indx,    &
     &                         mlold(11),abcosmic(iel)                  
          mlold(11)=indx 
          derivedpointers%npar(indx)=0 
          indx=indx+1 
                                                                        
!  go through ions                                                      
                                                                        
          ltyp=masterdata%nptrs(2,indx) 
          lrtp=masterdata%nptrs(3,indx) 
          if (lpri.gt.1)                                                &
     &     write (lun11,*)'lrtp=',lrtp,indx                             
          do while(lrtp.eq.12) 
                                                                        
!                                                                       
          if (lpri.gt.1)                                                &
     &     write(lun11,*) iel,iion                                      
                                                                        
!  register ion record                                                  
!  npfirst,npnxt,npar,mlold                                             
                                                                        
            if (derivedpointers%npfirst(12).eq.0) then 
              derivedpointers%npfirst(12)=indx 
            else 
              derivedpointers%npnxt(mlold(12))=indx 
            endif 
            if (lpri.gt.1)                                              &
     &       write (lun11,*)'npfirst(12)=',                             &
     &       derivedpointers%npfirst(12),indx             
            mlold(12)=indx 
            derivedpointers%npar(indx)=mlold(11) 
            indx=indx+1 
                                                                        
!  level records, rate type 13                                          
!  npfirst,npnxt,npar,mlold,npfi,npilev,npilevi                         
                                                                        
            if (masterdata%nptrs(3,indx).eq.13) then 
              derivedpointers%npfi(13,iion)=indx 
              iilev=1 
              if (derivedpointers%npfirst(13).eq.0) then 
                derivedpointers%npfirst(13)=indx 
              else 
                derivedpointers%npnxt(mlold(13))=indx 
              endif 
              if (lpri.gt.1)                                            &
     &         write (lun11,*)'filling npilev:'                         
              do while(masterdata%nptrs(3,indx).eq.13) 
                derivedpointers%npar(indx)=mlold(12) 
                derivedpointers%npnxt(indx)=indx+1 
                derivedpointers%npilev(iilev,iion)=ilev 
                derivedpointers%npilevi(ilev)=iilev 
                if(derivedpointers%npilev(iilev,iion).eq.0)             &
     &                print *, 'AJA **** ', iilev,iion 
                if (lpri.gt.1)                                          &
     &           write (lun11,*)ilev,iilev,indx,iion                    
                iilev=iilev+1 
                ilev=ilev+1 
                indx=indx+1 
              enddo 
              mlold(13)=indx-1 
              derivedpointers%npnxt(indx-1)=0 
            endif 
                                                                        
                                                                        
            do i=1,2 
              if (i.eq.1) then 
                lrtp=7 
              else 
                lrtp=1 
              endif 
              if (masterdata%nptrs(3,indx).eq.lrtp) then 
                derivedpointers%npfi(lrtp,iion)=indx 
                if (derivedpointers%npfirst(lrtp).eq.0) then 
                  derivedpointers%npfirst(lrtp)=indx 
                else 
                  derivedpointers%npnxt(mlold(lrtp))=indx 
                endif 
                if (lpri.gt.1)                                          &
     &           write (lun11,*)'npconi loop',indx                      
                do while(masterdata%nptrs(3,indx).eq.lrtp) 
!                 npcon points from the array of continuum emissivities 
!                    to the photoionization data                        
!                 npconi points from the levels to the arrays of        
!                    array of continuum emissivities                    
!                 npconi2 points from the photoionization data          
!                    to the array of continuum emissivities             
!                    (inverse if npcon)                                 
!                 icon is the index of the continuum emissivity array   
!                    element                                            
!                 indx is the index of the photoionization data         
                  derivedpointers%npar(indx)=mlold(12) 
                  derivedpointers%npnxt(indx)=indx+1 
                  derivedpointers%npcon(icon)=indx 
                  if (lpri.gt.1)                                        &
     &             write (lun11,*)'index into continuum  array:',       &
     &                icon                                              
                  if (lpri.gt.1)                                        &
     &             write (lun11,*)'index of photoionization element:',  &
     &                indx                                              
!                 now search for the level that goes with this          
!                    photoionization data                               
                  mltmpn=derivedpointers%npfi(13,iion) 
                  mlfnd=0 
                  nclev=masterdata%idat1(masterdata%nptrs(6,indx)       &
     &                  +masterdata%nptrs(9,indx)-2) 
                  if (nclev.gt.derivedpointers%nlevs(iion))             &
     &                  derivedpointers%nlevs(iion)=nclev 
                  mltst=nclev 
                  if (lpri.gt.1)                                        &
     &             write (lun11,*)'searching for level:'                
                  mltmp=mltmpn 
                  if (mltmpn.ne.0) then 
                      npartmpn=derivedpointers%npar(mltmpn) 
                    else 
                      npartmpn=0 
                    endif 
                  do while ((mlfnd.ne.mltst).and.(mltmpn.ne.0)          &
     &             .and.(indx.ne.0)                                     &
     &             .and.(npartmpn.eq.derivedpointers%npar(indx)))       
                    mltmp=mltmpn 
                    mlfnd=masterdata%idat1(masterdata%nptrs(6,mltmp)    &
     &                    +masterdata%nptrs(9,mltmp)-2) 
                    mltmpn=derivedpointers%npnxt(mltmp) 
                    if (mltmpn.ne.0) then 
                        npartmpn=derivedpointers%npar(mltmpn) 
                      else 
                        npartmpn=0 
                      endif 
                    if (lpri.gt.1)                                      &
     &              write (lun11,*)mltmp,mlfnd,mltmpn,npartmpn,         &
     &                      derivedpointers%npar(indx),nclev                   
                    enddo 
                  derivedpointers%npconi2(indx)=icon 
!                  npconi(icon)=npfi(13,iion)-1+nclev                   
                  if (mltmp.ne.0) then 
                    derivedpointers%npconi(mltmp)=icon 
                    endif 
                  if (lpri.gt.1)                                        &
     &             write (lun11,*)indx,derivedpointers%npar(indx),      &
     &             icon,nclev,masterdata%nptrs(3,indx),lrtp,mltmp               
                  indx=indx+1 
                  icon=icon+1 
                enddo 
                mlold(lrtp)=indx-1 
                derivedpointers%npnxt(indx-1)=0 
              endif 
           enddo 
!                                                                       
!  lines data and lines pointers, rate type 4, 9 & 14                   
!  npfirst,npnxt,npar,mold,npfi,nplin,nplini                            
                                                                        
            if (lpri.gt.1)                                              &
     &       write (lun11,*)'nplin,nplini,:'                            
            do i=1,3 
              if (i.eq.1) then 
                lrtp=4 
              elseif (i.eq.2) then 
                lrtp=9 
!               I don't think 2 photon should be treated as a line      
!                lrtp=-99                                               
              else 
                lrtp=14 
              endif 
              if (lpri.gt.1) write (lun11,*)' indx=',indx,lrtp,iion 
              if (masterdata%nptrs(3,indx).eq.lrtp) then 
                derivedpointers%npfi(lrtp,iion)=indx 
                if (derivedpointers%npfirst(lrtp).eq.0) then 
                  derivedpointers%npfirst(lrtp)=indx 
                else 
                  derivedpointers%npnxt(mlold(lrtp))=indx 
                endif 
                do while(masterdata%nptrs(3,indx).eq.lrtp) 
                  derivedpointers%npar(indx)=mlold(12) 
                  derivedpointers%npnxt(indx)=indx+1 
                  derivedpointers%nplin(iline)=indx 
                  derivedpointers%nplini(indx)=iline 
                  if (lpri.gt.1)                                        &
     &             write (lun11,*)indx,iline                            
                  indx=indx+1 
                  iline=iline+1 
                enddo 
                mlold(lrtp)=indx-1 
                derivedpointers%npnxt(indx-1)=0 
              endif 
            enddo 
                                                                        
!  pointers for rate types 6,8,3,5,40                                   
!  npfirst,npnxt,npar,mold,npfi                                         
                                                                        
            do i=1,5 
              if (i.eq.1) then 
                lrtp=6 
              elseif (i.eq.2) then 
                lrtp=8 
              elseif (i.eq.3) then 
                lrtp=3 
              elseif (i.eq.4) then 
                lrtp=5 
              else 
                lrtp=40 
              endif 
              if (masterdata%nptrs(3,indx).eq.lrtp) then 
                derivedpointers%npfi(lrtp,iion)=indx 
                if (derivedpointers%npfirst(lrtp).eq.0) then 
                  derivedpointers%npfirst(lrtp)=indx 
                else 
                  derivedpointers%npnxt(mlold(lrtp))=indx 
                endif 
                do while(masterdata%nptrs(3,indx).eq.lrtp) 
                  derivedpointers%npar(indx)=mlold(12) 
                  derivedpointers%npnxt(indx)=indx+1 
                  indx=indx+1 
                enddo 
                mlold(lrtp)=indx-1 
                derivedpointers%npnxt(indx-1)=0 
              endif 
            enddo 
                                                                        
!  pointers for other rate types                                        
!  npfirst,npnxt,npar,mold,npfi                                         
                                                                        
            lrtp=masterdata%nptrs(3,indx) 
            do while((lrtp.ne.12).and.(lrtp.ne.11).and.(lrtp.ne.0)) 
              derivedpointers%npar(indx)=mlold(12) 
              if (derivedpointers%npfirst(lrtp).eq.0) then 
                derivedpointers%npfirst(lrtp)=indx 
              else 
                derivedpointers%npnxt(mlold(lrtp))=indx 
              endif 
              mlold(lrtp)=indx 
              if (derivedpointers%npfi(lrtp,iion).eq.0)                 &
     &              derivedpointers%npfi(lrtp,iion)=indx 
!              write (lun11,*)iion,lrtp,indx,npfi(lrtp,iion)            
              indx=indx+1 
              lrtp=masterdata%nptrs(3,indx) 
            enddo 
                                                                        
!  ionization data and continum pointers, rate type 7 & 1               
!  npfirst,npnxt,npar,mlold,npfi,npcon,npconi,npconi2,nlevs             
                                                                        
                                                                        
            iion=iion+1 
                                                                        
          enddo 
        endif                                                          
      enddo 
                                                                        
      nlsvn=iline-1 
      ncsvn=icon-1 
      write (lun11,*)'number of lines=',nlsvn 
      write (lun11,*)'number of rrcs=',ncsvn 
!
!      if (lpri.eq.0) return
!
!     now do a big print

      go to 9000 
!                                                                       
!     sort the element abundances                                       
      lsrt=0 
      do mml=1,nl 
        melpt(mml)=mml 
      enddo 
      niter=0 
      do while (lsrt.eq.0) 
        lsrt=1 
        niter=niter+1 
        do mml=1,nl-1 
          if (abel(melpt(mml)).lt.abel(melpt(mml+1))) then 
            melptmp=melpt(mml) 
            melpt(mml)=melpt(mml+1) 
            melpt(mml+1)=melptmp 
            lsrt=0 
          endif 
        enddo 
      enddo 
!                                                                       
!                                                                       
!     now redo the element pointers                                     
!     zero the new next pointers                                        
      do mml=1,np2 
        npnxt2(mml)=0 
        enddo 
      npfirst2=0 
      mllo=0 
!     step thru elements                                                
      do mml=1,nl 
        mlloo=mllo 
        mll=derivedpointers%npfirst(11) 
        itst=0 
        do while ((mll.ne.0).and.(itst.ne.melpt(mml))) 
          mlm=mll
          call drd(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,      &
     &      mlm,0,lun11)                                          
          itst=masterdata%idat1(np1i-1+nidt) 
          mllo=mll 
          mll=derivedpointers%npnxt(mll) 
          enddo 
        if (mllo.ne.0) then 
          if (npfirst2.eq.0) then 
              npfirst2=mllo 
            else 
              npnxt2(mlloo)=mllo 
            endif 
          endif 
        enddo 
      npnxt2(mlloo)=0 
      derivedpointers%npfirst(11)=npfirst2 
      do mml=1,np2 
        if ((npnxt2(mml).ne.0).or.(mml.eq.mlloo)) then 
          derivedpointers%npnxt(mml)=npnxt2(mml) 
          endif 
        enddo 
!                                                                       
!                                                                       
 9000  continue 
!                                                                       
!       return 
!                                                                       
!     now print stuff sorted                                            
!      ntptmp=11 
!        mll=derivedpointers%npfirst(ntptmp) 
!        write (lun11,*)'ntptmp=',ntptmp 
!        do while (mll.ne.0) 
!          CALL DRD(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,      &
!     &          mll,0,Lun11)                                    
!          write (lun11,*)'mll=',mll 
!          call dprints(ltyp,lrtyp2,lcon,                             &
!     &    nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)        
!          mll=derivedpointers%npnxt(mll) 
!          enddo 
!     
!        write (lun11,*)'in setptrs filling nptrt'
        mm1=0
        do mm=1,ndat2
          iel=0 
!          write (lun11,*)'mm=',mm 
          CALL DRD(ltyp,lrtyp,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,      &
     &          mm,0,Lun11)                                    
!          call dprints(ltyp,lrtyp,lcon,                                &
!     &      nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)        
          if (lrtyp.ne.0) then
            if (lrtyp.ne.11) then
              mlion=mm
              if (lrtyp.ne.12) then
                mlion=derivedpointers%npar(mm)
                endif
              if (mlion.ne.0) then
                mlel=derivedpointers%npar(mlion)
                if (mlel.ne.0) then
                  CALL DRD(ltyp,lrtyp,lcon,nrdt,np1r,nidt,np1i,nkdt,   &
     &                     np1k,mlel,0,Lun11)                       
                  endif
                endif
              else
                mlel=mm
              endif
            if (mlel.ne.0) then
              iel=masterdata%idat1(np1i)
!              write (lun11,*)'mm=',mm,iel,abel(iel)
              if (iel.eq.0) stop 'iel=0'
!              if (mm.gt.200000) stop
              if (abel(iel).gt.1.e-15) then
                mm1=mm1+1
                nptrt(mm1)=mm
!                write (lun11,*)'mm1=',mm1,mm
!                CALL DRD(ltyp,lrtyp,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k, &
!     &            mm,0,Lun11)                                    
!                call dprints(ltyp,lrtyp,lcon,                           &
!     &            nrdt,np1r,nidt,np1i,nkdt,np1k,lun11)        
                endif
              endif
            endif
          enddo
        np2=mm1
!
      if (lpri.ge.1) write (lun11,*)'done with setptrs'
!        
!       print data
!        call dbwk2(12,abel,np2,nptrt,karg,lpri,lun11,                   &
!     &  nlsvn,ncsvn)                                             
!       print pointers
!        call dbwk2(25,abel,np2,nptrt,karg,lpri,lun11,                   &
!     &  nlsvn,ncsvn)                                             
!       set up pointers (again)
        call dbwk2(7,abel,np2,nptrt,karg,lpri,lun11,                    &
     &  nlsvn,ncsvn)                                             
!       print pointers
!        call dbwk2(25,abel,np2,nptrt,karg,lpri,lun11,                   &
!     &  nlsvn,ncsvn)                                             
!
      deallocate(nptrt)
      deallocate(npnxt2)
!                                                                       
      return 
      END                                           
