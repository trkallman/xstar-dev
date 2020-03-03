      subroutine linopac(lprie,lun11,optpp,                          &
     &                   rcem1,rcem2,elin,vturbi,t,aatmp,delea,epi,ncn2,&
     &                   opakc,rccemis,lfast)           
!                                                                       
!     Name:  linopac.f90
!     Description:
!         puts line opacity into continuum bins
!     author:  T. Kallman
!
!     List of Parameters:
!         Input:
!         lprie:
!         lun11:
!         optpp:  line center opacity from ucalc (cm^-1)
!         rcem1:  line emissivity outward (erg/s/cm^3)
!         rcem2:  line emissivity inward (erg/s/cm^3)
!         elin: line wavelength (A)
!         vturbi:  ion turbulent speed (km/s)
!         t:  temperature (10^4 K)
!         aatmp: ion mass (amu)
!         delea: natural line width
!         epi(ncn):  continuum energy bins (eV)
!         ncn2:  number of continuum energy bins
!         lfast:  switch to choose between full voigt profile (<=2)
!                 and single bin for line (>2)
!         Output:
!         opakc(ncn): continuum opacity (cm^-1)
!         rccemis(2,ncn): continuum emissivity (erg/s/cm^2)
!                        (in current version this is only modified if
!                        lfast > 2)
!
!     Dependencies:
!         none
!     Called by:
!         ucalc
!
!     this routine puts line opacity into continuum bins                
!     author:  T. Kallman                                               
      use globaldata
!                                                                       
      implicit none 
!                                                                       
      integer nbtpp 
      parameter (nbtpp=20000) 
!                                                                       
      real(8) epi(ncn),opakc(ncn)
      integer ldon(2) 
      real(8)  rccemis(2,ncn)
!                                                                       
      real(8)  prftmp,sum,rcem1,rcem2
!                                                                       
      real(8) vturbi,t,                                                 &
     &  dpcrit,bbb,optpp,delea,aatmp,elin,etmp,vth,                     &
     &  vturb,deleturb,deleth,dele,aasmall,                             &
     &  deleused,deleepi,delet,deletpp,e00,                             &
     &  e0,tst,opsum,optmpo,profile,optp2,                              &
     &  tmpopmx,tmpopo,etptst,opsv4
      integer lpri,lun11,                                               &
     &  ncn2,                                                           &
     &  ij,                                                             &
     &  ldir,mlm,ml1m,ml1,ml2,ml1min,                                   &
     &  ml1max,mlc,mlmin,mlmax,ncut,                                    &
     &  lprie                                       
      integer nbinc,lfast 
      real(8) voigte 
!                                                                       
!     temporary grid for use in calculating profile                     
      real(8) etpp(nbtpp),optpp2(nbtpp)
!                                                                       
!      real(8)  tmpew,tmpewo,tmpop,tmpe,sum,sume                         
      real(8) tmpew,tmpop,tmpe,sume,ergsev 
!                                                                       
      data dpcrit/1.e-6/,ergsev/1.602197e-12/ 
!                                                                       
!                                                                       
      lpri=lprie 
!      lpri=0                                                           
!                                                                       
!     test whether line is in range                                     
      if ((elin.gt.1.e+8).or.(elin.lt.1.)) return 
!                                                                       
!     for scattering model, add in line opacity                         
      bbb=vturbi 
      elin=abs(elin) 
!     thermal width quantities                                          
      vth=(1.29E+1)*sqrt(t/aatmp) 
      vturb=bbb 
!      e0=(12398.41)/max(elin,1.E-24)                                   
      e0=(12398.41)/max(elin,1.d-49) 
      if (e0.le.epi(1)) return 
      deleturb=e0*(vturb/3.E+5) 
      deleth=e0*(vth/3.E+5) 
!     old expression                                                    
!     dele=deleth+deleturb                                              
!     new expression                                                    
      dele=sqrt(deleth*deleth+deleturb*deleturb) 
      aasmall=delea/(1.E-24+dele)/12.56 
!                                                                       
!     continuum bin for line                                            
      ml1=nbinc(e0,epi,ncn2) 
      ml1=max(min(ncn-1,ml1),2) 
!                                                                       
!     here is what we do to get the heating right                       
      prftmp=2./(epi(ml1+1)-epi(ml1-1)) 
      opsv4=optpp*dele 
!                                                                       
!     print line quantities                                             
      if (lpri.ge.1) write (lun11,*)                                    &
     &   'e0,optpp,dpcrit*opakc(ml1),ml1,deleth,delea:',                &
     &    e0,optpp, dpcrit*opakc(ml1),ml1,deleth,delea                  
      if (lpri.ge.1) write (lun11,*)optpp,prftmp,opsv4,dele,            &
     &       opsv4*prftmp,rcem1,rcem2                                   
!                                                                       
!     test for simple calculation                                       
      if (lfast.gt.2) then 
!                                                                       
!         single bin calculation                                        
          opakc(ml1)=opakc(ml1)+opsv4*prftmp 
          rccemis(1,ml1)=rccemis(1,ml1)+rcem1*prftmp/ergsev/12.56 
          rccemis(2,ml1)=rccemis(2,ml1)+rcem2*prftmp/ergsev/12.56 
          return 
!                                                                       
!       full profile calculation                                        
        else 
!                                                                       
!         calculate profile on temporary grid                           
!         set up temporary grid                                         
          e00=epi(ml1) 
          etmp=e0 
!         deleepi is the grid spacing of the epi grid                   
!         deletpp is the physical energy spacing needed                 
!           for an accurate integration of the voigt profile            
!         ncut is the ratio of these two quantities,                    
!           used for rebinning the calculated voigt profile             
          deleepi=epi(ml1+1)-epi(ml1) 
!         expanding step to make broader lines                          
          deletpp=dele 
          ncut=int(deleepi/deletpp) 
          ncut=max(ncut,1) 
          ncut=min(ncut,int(nbtpp/10))
          deleused=deleepi/float(ncut) 
          mlc=0 
          ldir=1 
          ldon(1)=0 
          ldon(2)=0 
          mlmin=nbtpp 
          mlmax=1 
          ml1min=ncn+1 
          ml1max=0 
          ml2=nbtpp/2 
          if (lpri.ge.1) write (lun11,*)'ncut=',ncut,deleused,deletpp,  &
     &                                  deleepi                         
!                                                                       
!         calculate profile at continuum bin closest to line center     
          delet=(e00-etmp)/dele 
          if (aasmall.gt.1.e-6) then 
              profile=voigte(abs(delet),aasmall)/1.772 
            else 
              profile=exp(-delet*delet)/1.772 
            endif 
          etpp(ml2)=e00 
          optpp2(ml2)=optpp*profile 
          tst=1. 
!                                                                       
!         now put profile on temporary grid                             
!         work outward in both directions from line center              
          do while ((ldon(1)*ldon(2).eq.0).and.(mlc.lt.nbtpp/2)) 
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
!                                                                       
!               energy of temporary grid point                          
                etptst=e00+float(ldir*mlc)*deleused 
!                                                                       
!               test to see if within allowed range                     
                if ((mlm.le.nbtpp).and.(mlm.ge.1)                       &
     &           .and.(etptst.gt.0.).and.(etptst.lt.epi(ncn2))) then    
!                                                                       
!                 calculate index extremes for later use                
!                 ml1m is index into epi grid                           
!                 ml1min and ml1max are extremes of ml1m                
!                 mlmin and mlmax are extremes of mlm                   
                  mlmin=min(mlm,mlmin) 
                  mlmax=max(mlm,mlmax) 
!                                                                       
!                 store energy bin!                                     
                  etpp(mlm)=e00+float(ldir*mlc)*deleused 
                                                                        
!                 calculate profile                                     
                  delet=(etpp(mlm)-etmp)/dele 
                  if (aasmall.gt.1.e-9) then 
                      profile=voigte(abs(delet),aasmall)/1.772 
                    else 
                      profile=exp(-delet*delet)/1.772 
                    endif 
!                                                                       
!                 calculate opacity                                     
                  optpp2(mlm)=optpp*profile 
!                  tst=optpp2(mlm)*delr                                 
                  tst=profile 
!                                                                       
!                 print                                                 
                  if (lpri.ge.1) write (lun11,*) 'first write',         &
     &             mlm,etpp(mlm),ij,                                    &
     &             deleused,delet,mlmin,mlmax,ml1,                      &
     &             mlc,profile,optpp2(mlm),                             &
     &             tst                                                  
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
          sum=0. 
          opsum=0. 
          tmpop=0. 
          tmpopmx=0. 
          sume=0. 
          ml1min=nbinc(etpp(mlmin),epi,ncn2) 
          ml1max=nbinc(etpp(mlmax),epi,ncn2) 
          ml1m=ml1min 
          if (lpri.ge.1) write (lun11,*)'renormalizing profile',        &
     &       ml2,mlmin,mlmax,ml1m,ml1min,ml1max                         
          tmpew=0. 
          mlmin=max(mlmin,2) 
          mlmax=min(mlmax,nbtpp) 
!                                                                       
!         step through temp grid bins                                   
!         and  sum over intervals                                       
          do mlm=mlmin+1,mlmax 
!                                                                       
            tmpopo=tmpop 
            tmpop=optpp2(mlm) 
            tmpopmx=max(tmpopmx,tmpop) 
            tmpe=abs(etpp(mlm)-etpp(mlm-1)) 
!                                                                       
!           update interval sum                                         
            sume=sume+tmpe 
            opsum=opsum+(tmpop+tmpopo)*tmpe/2. 
!                                                                       
!           test to see if you have reached epi grid boundary           
            if (etpp(mlm).gt.epi(ml1m)) then 
!                                                                       
!             store current sum                                         
              optmpo=opakc(ml1m) 
              if (sume.gt.1.d-34) then 
                optp2=opsum/sume 
               do while ((etpp(mlm).gt.epi(ml1m)).and.(ml1m.lt.ncn2)) 
                  opakc(ml1m)=opakc(ml1m)+optp2 
!                 print                                                 
                  if (lpri.ge.1) write (lun11,*)mlm,ml1m,               &
     &             epi(ml1m),epi(ml1m+1),etpp(mlm),opakc(ml1m),         &
     &               optmpo,optpp2(mlm),opsum,sume                      
                  ml1m=ml1m+1 
                  enddo 
                endif 
!                                                                       
!             reset interval sums                                       
              tmpopmx=0. 
              opsum=0. 
              sume=0. 
!                                                                       
!             end of test for epi bin boundary                          
              endif 
!                                                                       
!           end of rebinning loop                                       
            enddo 
!                                                                       
!         norm check                                                    
!          rnormchk=ewsv(nlsv)/optpp/dele/(1.e-34+delr)                 
!          if (lpri.ne.0) write (lun11,*)'norm check',nilin,elin,optpp, 
!     $          dele,aasmall,ewsv(nlsv),rnormchk                       
!                                                                       
!       end of test for fast calculation                                
        endif 
!                                                                       
!                                                                       
      return 
      end                                           
