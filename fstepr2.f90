      subroutine fstepr2(unit,hdunum,radin,radout,rdel,temp,pres,abel,  &
     &                xcol,xee,xpx,xi,                                  &
     &                np2,ncsvn,nlsvn,                                  &
     &                rcem,oplin,tau0,                                  &
     &                lun11,lpri,status)                                
!                                                                       
!     Append a FITS extension binary table containing                   
!     nrhs columns and at most nrhdimj rows                             
!     author: T. Bridgman                                               
!                                                                       
!     Parameters:                                                       
!        unit    integer            File unit number                    
!        hdunum  integer            Number of last HDU written          
!        radin   real(8)               inner radius of shell             
!        radout  real(8)               outer radius of shell             
!                                   nb but now it is delr in the call   
!        temp    real(8)               temperature of shell              
!        pres    real(8)               pressure in shell                 
!        nrhdimj  integer            Maximum number of rows             
!        idat1   integer(nidat1)    Needed by the atomic database       
!        rdat1   real(nidat1)       Needed by the atomic database       
!        kdat1   char*nidat1        Needed by the atomic database       
!        nptrs                      Needed by the atomic database       
!        npnxt                      Needed by the atomic database       
!        npfi                       Needed by the atomic database       
!        npfirst                    Needed by the atomic database       
!        npcon                      Needed by the atomic database       
!        npconi                     Needed by the atomic database       
!        npcon2                     Needed by the atomic database       
!        xilev   real(nrhdimj)       Fractional level population array  
!        cemab   real(2,nrhdimj)     Recombination emission             
!        opakab  real(nrhdimj)       Opacity                            
!        tauc    real(2,nrhdimj)     Optical depth                      
!        poptol  real(8)               Tolerance for population level    
!        nzone   integer            Pass number through iteration proces
!        status  integer            Returned status code                
!                                                                       
      use globaldata
      implicit none 
                                                                        
      integer nptmpdim 
      parameter (nptmpdim=500000) 
!                                                                       
!     Allocation for passed parameters                                  
      real(8) tau0(2,nnnl), rcem(2,nnnl) 
      real(4) rtmp 
      real(8) radin, radout,rdel, temp, pres,xcol,xee,xpx,xi 
      integer unit,hdunum, nrows, status
!     line opacities                                                    
      real(8) oplin(nnnl) 
      real(8) abel(nl) 
                                                                        
!     Internal work areas                                               
      real(4) rwrk1(nptmpdim) 
      integer ntptr(nptmpdim) 
      character(10) kion(nptmpdim) 
      character(20) klevl(nptmpdim),klevu(nptmpdim),kblnk20 
      integer tfields,varidat 
      character(16) ttype(10),tform(10),tunit(10) 
      integer colnum,frow,felem,hdutype,ll, ltyp 
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpri 
      integer jkk, nlev
      integer nlplmx,ln,lnn,ml,nlpl,                                    &
     &         nlines,nlsvn,                                            &
     &         nilin,mlm,np2,ncsvn
      integer np1i,np1r,np1k 
      real(8) elin
      real(8) rniss(nd),rnisse(nd) 
      real(8) aij,ergsev,etst,ener,gglo,ggup
      integer idest1,idest2,ilevlo,ilevup,j,ktt,lk 
      character(33) extname 
      character(20) klablo,klabup 
      character(9) kinam1 
                                                                        
!     Database manipulation quantities                                  
      integer nkdt 
      character(1) kblnk
      real(4) elsv(nptmpdim) 
                                                                        
      data kblnk/' '/ 
      data kblnk20/'                    '/ 
!                                                                       
      data tform/'1J','1E','8A','20A','20A','1E','1E','1E',             &
     & '1E','1E'/                                                       
                                                                        
      data ttype/'index','wavelength','ion',                            &
     & 'lower_level','upper_level','emis_inward',                       &
     & 'emis_outward','opacity','tau_in','tau_out'/                     
                                                                        
      data tunit/' ','A',' ',' ',' ','erg/cm^3/s',                      &
     & 'erg/cm^3/s','/cm',' ',' '/                                      
                                                                        
      varidat=0 
!                                                                       
                                                                        
!     Move to the last HDU (hdunum) in the file                         
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr2: Moving to end-of-FITS file'              
      call ftmahd(unit,hdunum,hdutype,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     append a new empty extension after the last HDU                   
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Create the new extension'               
      call ftcrhd(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
!                                                                       
!     Extracting data from the Atomic Database here                     
!                                                                       
      if (lpri.ne.0)                                                    &
     & write (lun11,*)' '                                               
!                                                                       
!     print important lines                                             
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'emission line luminosities (erg/sec/10**38))',   &
     &                 nlsvn                                            
         nlplmx=nptmpdim 
!                                                                       
!     step through lines                                                
      nlpl=0 
      do lnn=1,nlsvn 
!                                                                       
        if ((rcem(1,lnn).gt.1.d-64).or.                                 &
     &      (rcem(2,lnn).gt.1.d-64).or.                                 &
     &      (oplin(lnn).gt.1.d-64)) then                                
!                                                                       
!         get line data                                                 
          ln=lnn 
          ml=derivedpointers%nplin(ln) 
          mlm=ml-1 
          call drd(ltyp,lrtyp,lcon,                                     &
     &      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                          &
     &      0,lun11)                                              
          if (lpri.ne.0) write (lun11,*)ln,ml,masterdata%rdat1(np1r) 
!                                                                       
!         exclude rate type 14                                          
          elin=abs(masterdata%rdat1(np1r)) 
          if ((lrtyp.ne.14).and.(abs(elin).gt.0.1).and.(lrtyp.ne.9)     &
     &       .and.(abs(elin).lt.9.e+9)) then                            
!                                                                       
            ergsev=1.602197e-12 
            ener=ergsev*(12398.41)/max(elin,1.e-24) 
            etst=ener/ergsev 
            idest1=masterdata%idat1(np1i) 
            idest2=masterdata%idat1(np1i+1) 
            aij=masterdata%rdat1(np1r+2) 
            if (lpri.ne.0) write (lun11,*)'line data',elin,ener,etst,   &
     &                       idest1,idest2,aij                          
!                                                                       
!           get ion data                                                
            nilin=derivedpointers%npar(ml) 
            mlm=nilin-1 
            call drd(ltyp,lrtyp,lcon,                                   &
     &        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,                        &
     &        0,lun11)                                            
            do ktt=1,min(8,nkdt) 
               write (kinam1(ktt:ktt),'(a1)')                           &
     &               masterdata%kdat1(np1k-1+ktt) 
              enddo 
            do ktt=nkdt+1,9 
              write (kinam1(ktt:ktt),'(a1)')kblnk 
              enddo 
!                                                                       
!           now find level data                                         
            jkk=masterdata%idat1(np1i+nidt-1) 
            if (lpri.ne.0) write (lun11,*)'ion',kinam1,jkk 
            call func2l(jkk,lpri,lun11,temp,xee,xpx,                    &
     &              rniss,rnisse,nlev)
                                                                        
            ggup=leveltemp%rlev(2,idest1) 
            gglo=leveltemp%rlev(2,idest2) 
            do lk=1,20 
              write (klablo(lk:lk),'(a1)')leveltemp%klev(lk,idest1) 
              write (klabup(lk:lk),'(a1)')leveltemp%klev(lk,idest2) 
              enddo 
            ilevlo=idest1 
            ilevup=idest2 
!                                                                       
            nlpl=nlpl+1 
            nlpl=min(nlpl,nlplmx) 
            ntptr(nlpl)=lnn 
            elsv(nlpl)=elin 
            kion(nlpl)=kinam1 
            klevl(nlpl)=klablo 
            klevu(nlpl)=klabup 
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
            if (lpri.ne.0)                                              &
     &       write (lun11,*)j,elin,oplin(j),rcem(1,j),                  &
     &                      rcem(2,j)                                   
!                                                                       
            endif 
!                                                                       
          endif 
!                                                                       
        enddo 
                                                                        
!      if (nlpl.le.0) return                                            
      nlpl=max(nlpl,1) 
!                                                                       
                                                                        
!     End of atomic database extraction                                 
!----------------------------------------------------------------       
!     define parameters for the binary table (see the above data stateme
      nrows=nlpl 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'before header write'                             
      tfields=10 
!     Build extension name                                              
      extname='XSTAR_RADIAL' 
                                                                        
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Write table headers'                    
!     write the required header parameters for the binary table         
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,         &
     &              varidat,status)                                     
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Add some more keywords'                 
                                                                        
!     Write some model parameters in the extension header               
      call ftpcom(unit,'***********************************',status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      call ftpcom(unit,'Model Keywords',status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     Write values to 3 decimal places                                  
      rtmp=radin 
      call ftpkye(unit,'RINNER',rtmp,3,'[cm] Inner shell radius',       &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=radout 
      call ftpkye(unit,'ROUTER',rtmp,3,'[cm] Outer shell radius',       &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=rdel 
      call ftpkye(unit,'RDEL',rtmp,3,'[cm] distance from face',         &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=temp 
      call ftpkye(unit,'TEMPERAT',rtmp,3,'[10**4K] Shell Temperature',  &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=pres 
      call ftpkye(unit,'PRESSURE',rtmp,3,'[dynes/cm**2] Shell Pressure',&
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      rtmp=xcol 
      call ftpkye(unit,'COLUMN',rtmp,3,'[/cm**2] Column ',              &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      rtmp=xee 
      call ftpkye(unit,'XEE',rtmp,3,'electron fraction',                &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      rtmp=xpx 
      call ftpkye(unit,'DENSITY',rtmp,3,'[/cm**3] Density',             &
     & status)                                                          
      if (status .gt. 0)call printerror(lun11,status) 
!                                                                       
      rtmp=xi 
      call ftpkye(unit,'LOGXI',rtmp,3,                                  &
     & '[erg cm/s] log(ionization parameter)',status)                   
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'after header write'                              
!-------------------------------------------------------------------    
!     Step through the columns and write them to the file               
!                                                                       
!     set 'global' parameters for writing FITS columns                  
      frow=1 
      felem=1 
                                                                        
                                                                        
!     column  1  (Line number)                                          
      colnum=1 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr2: Writing Column ',colnum,nlpl             
      nlines=nlpl 
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  2  (wavelength)                                           
      colnum=2 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr2: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,elsv,status) 
      if (status .gt. 0)call printerror(lun11,status) 
      if (status .gt. 0) return 
                                                                        
                                                                        
!     column  3  (Ion)                                                  
      colnum=3 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr2: Writing Column ',colnum                  
      call ftpcls(unit,colnum,frow,felem,nlines,kion,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  4 (lower Level Designation)                               
      colnum=4 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Writing Column ',colnum                 
      call ftpcls(unit,colnum,frow,felem,nlines,klevl,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  5 (Level Designation)                                     
      colnum=5 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr2: Writing Column ',colnum                  
      call ftpcls(unit,colnum,frow,felem,nlines,klevu,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!----------------------------------------------------------------       
                                                                        
!     column  6                                                         
      colnum=6 
      do ll=1,nlines 
         rwrk1(ll)=0. 
         if (ntptr(ll).ne.0)                                            &
     &      rwrk1(ll)=rcem(1,ntptr(ll))                                 
         if (lpri.ne.0) write (lun11,*)ll,ntptr(ll),rwrk1(ll) 
         enddo 
      if (lpri.ne.0)                                                    &
     & write(lun11,*)'fstepr2: Writing Column ',colnum                  
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  7                                                         
      colnum=7 
      do ll=1,nlines 
         rwrk1(ll)=0. 
         if (ntptr(ll).ne.0)                                            &
     &    rwrk1(ll)=rcem(2,ntptr(ll))                                   
         enddo 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Writing Column ',colnum                 
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  8                                                         
      colnum=8 
      do ll=1,nlines 
         rwrk1(ll)=0. 
         if (ntptr(ll).ne.0)                                            &
     &    rwrk1(ll)=oplin(ntptr(ll))                                    
         enddo 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Writing Column ',colnum                 
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!     column  9                                                         
      colnum=9 
      do ll=1,nlines 
         rwrk1(ll)=0. 
         if (ntptr(ll).ne.0)                                            &
     &    rwrk1(ll)=tau0(1,ntptr(ll))                                   
         enddo 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Writing Column ',colnum                 
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!     column  10                                                        
      colnum=10 
      do ll=1,nlines 
         rwrk1(ll)=0. 
         if (ntptr(ll).ne.0)                                            &
     &    rwrk1(ll)=tau0(2,ntptr(ll))                                   
         enddo 
      if (lpri.ne.0)                                                    &
     & write (lun11,*)'fstepr2: Writing Column ',colnum                 
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
                                                                        
!----------------------------------------------------------------       
!     Compute checksums                                                 
      call ftpcks(unit,status) 
      if (status .gt. 0)call printerror(lun11,status) 
                                                                        
!                                                                       
!                                                                       
      return 
      end                                           
