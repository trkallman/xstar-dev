      subroutine dprinto2(ltyp,lrtyp,lcon,                              &
     &  nrdt,rdat,nidt,idat,nkdt,kdat,lun11)                            
!                                                                       
!     Name: dprinto2.f90  
!     Description:  
!     this  routine printe one record of the database 
!     differs from dprinto because data is passed directrly
!     author:  T. Kallman                                               
!     List of Parameters:
!           Input:
!           ltyp=data type
!           lrtyp=rate type
!           lcon=continuation switch
!           nrdt=number of reals
!           rdat=array of real data
!           nidt=number of integers
!           idat=array of int data
!           nkdt=number of chars
!           kdat=array of char data
!           lun11=logical unit number for printout
!     Dependencies:  none
!     Called by:  ucalc, pprint, setptrs, dprinto2, func1, writespectra2
!                                                                       
!                                                                       
      implicit none 
!                                                                       
      integer nptmpdim 
      parameter (nptmpdim=200000) 
!                                                                       
      real(8) rdat(nptmpdim) 
      integer idat(nptmpdim) 
      character(1) kblnk,kperc,kdat(nptmpdim) 
      integer ltyp,lrtyp,lcon,nrdt,nidt,nkdt,lun11,lsp,ml2,nkd,         &
     &        nkd2,ll2,ml,mm,ll                                         
!                                                                       
      character(400000) kdtt 
      character(1) ktst,kdtt2(400000) 
!                                                                       
      data kblnk/' '/,kperc/'%'/ 
!                                                                       
!      write (lun11,*)ltyp,lrtyp,lcon,nrdt,nidt,nkdt,np1r,np1i,np1k     
!     $  (rdat(mm),mm=1,nrdt),(idat1(np1i-1+mm),mm=1,nidt),             
!     $  kblnk,(kdat1(np1k-1+mm),mm=1,nkdt),kblnk,kper!                 
!      return                                                           
!                                                                       
      write (kdtt(1:18),9823)ltyp,lrtyp,lcon 
      write (kdtt(19:37),9823)nrdt,nidt,nkdt 
!                                                                       
!                                                                       
      lsp=0 
      ml2=0 
      if (lsp.eq.1) go to 9009 
      nkd=38 
      nkd2=39 
      if (nrdt.gt.0) then 
!        nkd2=nkd+nrdt*13                                               
        nkd2=nkd 
        do 303 ml=1,nrdt 
           ml2=nkd+(ml-1)*15 
             write (kdtt(ml2:ml2+14),'(1pe15.7)')rdat(ml) 
             nkd2=nkd2+15 
  303      continue 
        endif 
      nkd=nkd2 
      write (kdtt(nkd:nkd),'(a1)')kblnk 
      nkd=nkd2+1 
      if (nidt.gt.0) then 
        nkd2=nkd+nidt*8 
        do 302 ml=1,nidt 
           ml2=nkd+(ml-1)*8 
           write (kdtt(ml2:ml2+7),'(i8)')idat(ml) 
           ml2=ml2+8 
  302      continue 
        endif 
      nkd=nkd2 
      if (nkdt.gt.0) then 
        write (kdtt(nkd:nkd),'(a1)')kblnk 
        nkd=nkd+1 
        nkd2=nkd+nkdt 
!        write (lun11,*)nkd,nkdt,nkd2,(kdat1(np1k-1+mm),mm=1,nkdt)      
        do 301 ml=1,nkdt 
          ml2=nkd+ml-1 
          write (kdtt(ml2:ml2),'(a1)')kdat(ml) 
           ml2=ml2+1 
  301     continue 
         endif 
 9823    format (3i6) 
!       write (lun11,*)'before write:'                                  
!       write (lun11,*)kdtt                                             
      ml2=ml2-1 
!                                                                       
!                                                                       
      ll2=0 
!     remove spaces                                                     
      ktst=kperc 
      do 3301 ll=1,ml2 
!         ktsto=ktst                                                    
         read(kdtt(ll:ll),'(a1)')ktst 
!         if ((ktst.eq.kblnk).and.(ktsto.eq.kblnk)) go to 3301          
         ll2=ll2+1 
         kdtt2(ll2)=ktst 
 3301    continue 
!                                                                       
      write (lun11,911)(kdtt2(mm),mm=1,ll2),kblnk,kperc 
  911 format (400000a1) 
!                                                                       
      return 
!                                                                       
 9009 continue 
!      call dprints(ltyp,lrtyp,lcon,                                    
!     $  nrdt,rdat,nidt,idat,nkdt,kdat,lun11)                           
!                                                                       
!                                                                       
      return 
      end                                           
