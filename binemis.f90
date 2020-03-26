      subroutine binemis(lun11,lpri,xlum,                            &
     &       t,vturbi,epi,ncn2,dpthc,                                   &
     &       nlsvn,                                                     &
     &       eliml,elimh,elum,zrtmp,zremsz)
!                                                                       
!     Name: binemis.f90  
!     Description:  
!           Puts emission lines into continuum bins, using voigt profile.
!           Emission analog of linopac
!
!     List of Parameters:
!     Input:
!           lun11: logical unit number for printing
!           lpri: print switch, 1=on, 0=off
!           xlum:  Input continuum luminosity 1 - 1000 Ry (erg/s/10^38)
!           r: radius in nebula (cm)
!           t:  temperature (10^4 K)
!           vturbi:  ion turbulent speed (km/s)
!           epi(ncn):  continuum energy bins (eV)
!           ncn2:  number of continuum energy bins
!           dpthc(2,ncn): optical depth in continuum bins 
!           nlsvn: atomic data parameter, number of lines in atomic database
!           eliml:  energy lower limit (eV)
!           elimh:  energy upper limit (eV)
!           elum(3,nnnl):  line luminosities (erg/s/10^38)
!           zremsz(ncn):  input spectrum  (erg/s/erg/10^38)
!         Output:
!           zrtmp(4,ncn):  master spectrum array.  (erg/s/erg/10^38)
!
!     Dependencies: voigte
!     Called by:  writespectra
!     
!     revised 3/2/2018 to reduce memory TK
!     note that initially zrtmp contains zrems
!     and that in this routine zrems is a temporary
!     no longer doing dynamic memory
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      integer nbtpp 
      parameter (nbtpp=ncn) 
!                                                                       
!                                                                       
!     passed parameters                                                 
      integer lun11 
!     line luminosities                                                 
      real(8) elum(3,nnnl) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum lum                                                     
      real(8) zrtmp(5,ncn),zremsz(ncn) 
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
      integer kl, nilin, nkdt,nidt,lcon,lrtyp,ltyp,ml 
      integer nlsvn, ln, ll, numcon, lnn, nbtmp, nrdt 
      integer  verbose 
      integer lup,ndtmp,mllz,iion,nitmp,                                &
     &        ltyp2,lrtyp2,lcon2,iltmp,mlpar                
      real(8) eliml, elimh, elmmtpp, dele, etmp, elin, aasmall 
      real(8) egam,profile,deler,delea,vturbi,aatmp 
      real(8) e00,deleepi,etptst,tst,sume,zrsum1,zrsum2,deletpp,         &
     &        deleused,tmpe,zrtp2,zrtp1,bbb,xlum
      integer ml1,mlmin,mlmax,ij,mlm,ldir,ml1m,mlc,ncut,                &
     &        ml2,np1k,np1i,np1r,np1i2,ml1max,ml1min,mm
!     arrays containing printout info                                   
!      integer ilsv(nnnl)
!      real(8) ewsv(nnnl),elsv(nnnl) 
      real(8), allocatable, dimension(:) :: etpp
      real(8), allocatable, dimension(:,:) :: zrtpp2
      real(8), allocatable, dimension(:,:) :: zrtmps
      real(8), allocatable, dimension(:,:) :: zrems
      integer ldon(2) 
!                                                                       
      integer lpri,ncn2,nlsv,nelin 
      real(8) t,delet,deleturb,deleth,e0,vth,vturb,dpcrit 
!                                                                       
!     externally defined functions                                      
      integer nbinc 
      real(8) voigte 
!                                                                       
      data dpcrit/1.e-6/ 
!                                                                       
      allocate(zrtpp2(2,nbtpp))
      allocate(etpp(nbtpp))
      allocate(zrtmps(2,nbtpp))
      allocate(zrems(5,ncn))
!
      verbose=lpri 
!                                                                       
!     open and prepare the fits file for spectral data                  
      if(verbose.gt.0) write (lun11,*)'in binemis:' 
!                                                                       
      nlsv=0 
!                                                                       
!     build spectra data tables                                         
      numcon=ncn2 
      bbb=vturbi 
      do ll=1,ncn2 
        do mm=1,5
          zrems(mm,ll)=zrtmp(mm,ll)
          enddo
        enddo 
      do ll=1,ncn2 
        zrtmps(1,ll)=0. 
        zrtmps(2,ll)=0. 
        zrtpp2(1,ll)=0. 
        zrtpp2(2,ll)=0. 
        enddo
!
      do  lnn=1,nlsvn 
!                                                                       
        ln=lnn 
        ml=derivedpointers%nplin(ln) 
        call drd(ltyp,lrtyp,lcon,                                       &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,ml,                       &
     &          0,lun11)                                          
        elin=abs(masterdata%rdat1(np1r)) 
        egam=masterdata%rdat1(np1r+2) 
        lup=masterdata%idat1(np1i+1) 
        nilin=derivedpointers%npar(ml) 
        call drd(ltyp,lrtyp,lcon,                                       &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,nilin,                    &
     &          0,lun11)                                          
        nelin=derivedpointers%npar(nilin) 
        nilin=masterdata%idat1(np1i+2) 
!       get nuclear mass                                                
        ml=nelin 
        call drd(ltyp,lrtyp,lcon,                                       &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,ml,                       &
     &          0,lun11)                                          
        aatmp=masterdata%rdat1(np1r+1) 
        elmmtpp=(elum(2,ln)+elum(1,ln))/2. 
        if ((lpri.gt.0))                                                &
     &     write (lun11,*)ln,elin,elmmtpp,nilin,nelin,egam,             &
     &      lup,aatmp,ln,nnnl,eliml,elimh,ltyp                                           
        if (((ln.gt.0).and.(ln.le.nnnl))                                &
     &    .and.((elin.gt.eliml).and.(elin.lt.elimh))                    &
     &    .and.(elmmtpp.gt.1.e-31*xlum).and.(aatmp.gt.1.e-24)           &
     &    .and.((nilin.gt.0).and.(nilin.le.nni))                        &
     &    .and.(ltyp.ne.76))                                            &
     &       then                                                       
!                                                                       
!         line parameters                                               
          etmp=12398.42/elin 
          nbtmp=nbinc(etmp,epi,ncn2) 
!                                                                       
!          nlsv=nlsv+1 
!          ilsv(nlsv)=ln 
!          ewsv(nlsv)=-elmmtpp/max(1.e-34,zremsz(nbtmp)) 
!          elsv(nlsv)=elmmtpp 
!          if (lpri.gt.0)                                                &
!     &      write (lun11,*)'nlsv,ilsv(nlsv),elsv(nlsv):',               &
!     &                       nlsv,ilsv(nlsv),elsv(nlsv)                 
!                                                                       
!         find associated type 86 data                                  
          iion=1 
          nitmp=derivedpointers%npfi(13,iion) 
          call drd(ltyp,lrtyp,lcon,                                     &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,nitmp,                    &
     &          0,lun11)                                          
          if (lpri.gt.0)                                                &
     &      write (lun11,*)'searching for ion'                          
          do while ((masterdata%idat1(np1i-1+nidt).ne.nilin)            &
     &          .and.(iion.lt.nni)) 
            iion=iion+1 
            nitmp=derivedpointers%npfi(13,iion) 
            call drd(ltyp,lrtyp,lcon,                                   &
     &          nrdt,np1r,nidt,np1i,nkdt,np1k,nitmp,                    &
     &          0,lun11)                                          
            if (lpri.gt.0)                                              &
    &        write (lun11,*)iion,masterdata%idat1(np1i-1+nidt),         &
    &             nilin,nitmp        
            enddo 
          ndtmp=derivedpointers%npfi(41,iion) 
          delea=0. 
          if (ndtmp.gt.0) then 
            if (lpri.gt.0)                                              &
     &        write (lun11,*)'  found ion',lup,ndtmp                    
            mllz=derivedpointers%npar(ndtmp) 
            call drd(ltyp2,lrtyp2,lcon2,                                &
     &          nrdt,np1r,nidt,np1i2,nkdt,np1k,ndtmp,                   &
     &         0,lun11)                                           
            iltmp=masterdata%idat1(np1i2+1) 
            mlpar=mllz 
            do while ((ndtmp.ne.0).and.(lup.ne.iltmp)                   &
     &         .and.(mlpar.eq.mllz))                                    
               call drd(ltyp2,lrtyp2,lcon2,                             &
     &           nrdt,np1r,nidt,np1i2,nkdt,np1k,ndtmp,                  &
     &          0,lun11)                                          
              iltmp=masterdata%idat1(np1i2+1) 
              if (lpri.gt.0)                                            &
     &           write (lun11,*)'   ',iltmp,ndtmp                 
              ndtmp=derivedpointers%npnxt(ndtmp) 
              mlpar=0 
              if (ndtmp.ne.0) mlpar=derivedpointers%npar(ndtmp) 
              enddo 
            endif 
          if (lup.eq.iltmp) then 
            delea=masterdata%rdat1(np1r+2)*(4.14e-15) 
            egam=masterdata%rdat1(np1r+3) 
            endif 
!                                                                       
!         cheat for narrow line plot                                    
!         delea=0.                                                      
!                                                                       
!         a list of all the deles                                       
!           delea=auger natural width in eV                             
!           deleturb=turbulent width                                    
!           deleth=thermal Doppler width                                
!           dele=thermal+turbulent width                                
!           deler=radiative natural width                               
!           deletpp=goal of resolution of internal grid=dele/8          
!           deleepi=xstar grid spacing                                  
!           deleused=spacing of internal grid=deleepi/int(deleepi/depetp
!           delet=energy offset from line center in units of dele (local
!                                                                       
!         thermal width quantities                                      
          vth=(1.2e+1)*sqrt(t/aatmp) 
          vturb=max(bbb,vth) 
          e0=(12398.42)/max(elin,1.d-49) 
          deleturb=e0*(vturb/3.e+5) 
          deleth=e0*(vth/3.e+5) 
!         old expression                                                
!          dele=deleth+deleturb                                         
!         new expression                                                
          dele=sqrt(deleth*deleth+deleturb*deleturb) 
          deler=egam*(4.14e-15) 
          aasmall=(delea+deler)/(1.e-36+dele)/12.56 
!                                                                       
          ml1=nbtmp 
          if (lpri.ge.1) write (lun11,*)                                &
     &   'e0,elin,elum1,elum2,ml1,deleth,delea:',                       &
     &    e0,elin,elum(1,ln),elum(2,ln),ml1,deleth,delea                
!                                                                       
!         calculate profile on temporary grid                           
          e00=epi(ml1) 
          etmp=e0 
!         deleepi is the grid spacing of the epi grid                   
!         deletpp is the physical energy spacing needed                 
!           for an accurate integration of the voigt profile            
!         ncut is the ratio of these two quantities,                    
!           used for rebinning the calculated voigt profile             
          deleepi=epi(ml1+1)-epi(ml1) 
          deletpp=dele 
          ncut=int(deleepi/deletpp) 
          ncut=max(ncut,1) 
          ncut=min(ncut,nbtpp/10) 
          deleused=deleepi/float(ncut) 
          mlc=0 
          ldir=1 
          ldon(1)=0 
          ldon(2)=0 
          mlmin=nbtpp 
          mlmax=1 
          ml1min=ncn+1 
          ml1max=0 
          ml2=int(nbtpp/2)
          if (lpri.gt.0) write (lun11,*)'ncut=',ncut,deleused,deletpp,  &
     &                                    deleepi                       
!                                                                       
!         calculate profile at continuum bin closest to line center     
          delet=(e00-etmp)/dele 
          if (aasmall.gt.1.e-9) then 
              profile=voigte(abs(delet),aasmall)/1.772 
            else 
              profile=exp(-delet*delet)/1.772 
            endif 
          profile=profile/dele/(1.602197e-12) 
          etpp(ml2)=e00
          zrtpp2(1,ml2)=elum(1,ln)*profile 
          zrtpp2(2,ml2)=elum(2,ln)*profile 
          tst=1. 
!                                                                       
!         now put profile on temporary grid                             
!         work outward in both directions from line center              
          do while ((ldon(1)*ldon(2).eq.0).and.(mlc.lt.int(nbtpp/2)))
!                                                                       
            mlc=mlc+1 
!                                                                       
!           alternate directions                                        
            do ij=1,2 
              ldir=-ldir 
!                                                                       
!             test to see if done in this direction                     
              if (ldon(ij).ne.1) then 
!                                                                       
!               index into temporary grid                               
                mlm=ml2+ldir*mlc 
                mlm=min(nbtpp,max(1,mlm))
                                                                        
                etptst=e00+float(ldir*mlc)*deleused 
!                                                                       
!               test to see if within allowed range                     
                if ((mlm.lt.nbtpp).and.(mlm.gt.1)                       &
     &            .and.(etptst.gt.0.).and.(etptst.lt.epi(ncn2))) then   
!                                                                       
!                 calculate index extremes for later use                
!                 ml1m is index into epi grid                           
!                 ml1min and ml1max are extremes of ml1m                
!                 mlmin and mlmax are extremes of mlm                   
                  mlmin=min(mlm,mlmin) 
                  mlmax=max(mlm,mlmax) 
!                                                                       
!                 store energy bin                                     
                  etpp(mlm)=e00+float(ldir*mlc)*deleused 
!                                                                       
!                 calculate profile                                     
                  delet=(etpp(mlm)-etmp)/dele 
                  if (aasmall.gt.1.e-9) then 
                      profile=voigte(abs(delet),aasmall)/1.772 
                    else 
                      profile=exp(-delet*delet)/1.772 
                    endif 
                  profile=profile/dele/(1.602197e-12) 
!                                                                       
                  zrtpp2(1,mlm)=elum(1,ln)*profile 
                  zrtpp2(2,mlm)=elum(2,ln)*profile 
                  tst=profile 
!                                                                       
!                 print                                                 
                  if (lpri.ge.3) write (lun11,*) 'first write',         &
     &               mlm,etpp(mlm),ij,                                  &
     &               deleused,delet,mlmin,mlmax,ml1,                    &
     &               mlc,profile,zrtpp2(2,mlm)    
!                                                                       
!                 end of test for within range                          
                  endif 
!                                                                       
!               test to see if done in this direction:                  
!                 profile not too small                                 
!                 index within range                                    
!                 energy within range                                   
!                 within specified number of doppler widths (50)        
                if (((tst.lt.dpcrit)                                    &
     &               .or.(mlm.le.1).or.(mlm.ge.nbtpp)                   &
     &               .or.(etptst.le.0.).or.(etptst.ge.epi(ncn2))        &
     &               .or.(mlc.gt.nbtpp)                                 &
     &               .or.(abs(delet).gt.max(50.d0,200.*aasmall)))       &
     &               .and.(ml1min.lt.ml1-2).and.(ml1max.gt.ml1+2)       &
     &               .and.(ml1min.ge.1).and.(ml1max.le.ncn))            &
     &                ldon(ij)=1                                        
!                                                                       
!               end of test for done in this direction                  
                endif 
!                                                                       
!             end of loop over directions                               
              enddo 
!                                                                       
!           end of loop over energies                                   
            enddo 
!                                                                       
!         store into continuum bins                                     
          sume=0. 
          zrsum1=0. 
          zrsum2=0. 
          ml1min=nbinc(etpp(mlmin),epi,ncn2) 
          ml1max=nbinc(etpp(mlmax),epi,ncn2) 
          ml1m=ml1min 
          if (lpri.ge.3) write (lun11,*)'renormalizing profile',        &
     &       ml2,mlmin,mlmax,ml1,ml1min,ml1max                          
          mlmin=max(mlmin,2) 
          mlmax=min(mlmax,nbtpp) 
!                                                                       
!         step through temp grid bins                                   
!         and  sum over intervals                                       
          do mlm=mlmin+1,mlmax 
!                                                                       
            tmpe=abs(etpp(mlm)-etpp(mlm-1)) 
            sume=sume+tmpe 
            zrsum1=zrsum1+(zrtpp2(1,mlm)+zrtpp2(1,mlm-1))*tmpe/2. 
            zrsum2=zrsum2+(zrtpp2(2,mlm)+zrtpp2(2,mlm-1))*tmpe/2. 
            if (lpri.ge.3) write (lun11,*)mlm,etpp(mlm),ml1m,epi(ml1m), &
     &         sume,zrsum1,zrsum2                                      
!                                                                       
!           test to see if you have reached epi grid boundary           
            if (etpp(mlm).gt.epi(ml1m)) then 
!                                                                       
!             store current sum                                         
              if (mlm.eq.mlmax) ml1m=max(1,ml1m-1) 
              if (sume.gt.1.d-24) then 
                zrtp2=zrsum2/sume 
                zrtp1=zrsum1/sume 
                do while ((etpp(mlm).gt.epi(ml1m)).and.(ml1m.lt.ncn2)) 
                  zrtmps(1,ml1m)=zrtp1 
                  zrtmps(2,ml1m)=zrtp2 
                  if (lpri.ge.3) write (lun11,*)mlm,ml1m,               &
     &               epi(ml1m),epi(ml1m+1),etpp(mlm),zrtmps(2,ml1m),    &
     &               zrtmps(1,ml1m)          
                  ml1m=ml1m+1 
                  enddo 
                endif 
!                                                                       
!             reset interval sums                                       
              zrsum2=0. 
              zrsum1=0. 
              sume=0. 
!                                                                       
!             end of test for epi bin boundary                          
              endif 
!                                                                       
!           end of rebinning loop                                       
            enddo 
!                                                                       
          do mlm=mlmin,mlmax 
            zrtpp2(1,mlm)=0.
            zrtpp2(2,mlm)=0.
            enddo

          do ml1m=ml1min,ml1max 
            zrtmp(4,ml1m)=zrtmp(4,ml1m)+zrtmps(2,ml1m) 
            zrtmp(3,ml1m)=zrtmp(3,ml1m)+zrtmps(1,ml1m) 
            if (lpri.ge.3) write (lun11,*)ml1m,                         &
     &           epi(ml1m),zrtmps(1,ml1m),zrtmp(3,ml1m)                 
            enddo 
!
          do ml1m=ml1min,ml1max 
            zrtmps(1,ml1m)=0.
            zrtmps(2,ml1m)=0.
            enddo
!                                                                       
          endif 
        enddo 
!                                                                       
      if (lpri.ge.1) write (lun11,*)'after first binemis loop' 
      do kl=1,numcon 
         if (lpri.ge.1) write (lun11,*)kl,epi(kl),zrems(2,kl),          &
     &          zrems(3,kl),zrtmp(4,kl),zrtmp(5,kl),dpthc(1,kl),        &
     &          zremsz(kl)                                              
         zrtmp(3,kl)=zrtmp(3,kl)+zrems(2,kl) 
         zrtmp(4,kl)=zrtmp(4,kl)+zrems(3,kl) 
         zrtmp(2,kl)=zremsz(kl)*exp(-dpthc(1,kl)) 
         zrtmp(1,kl)=zremsz(kl) 
         zrtmp(5,kl)=zrems(4,kl) 
         enddo 
!
      deallocate(zrtpp2)
      deallocate(etpp)
      deallocate(zrtmps)
      deallocate(zrems)
!         
      return 
      END                                           
