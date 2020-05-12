      subroutine rlbin(jkk1,ebln,emis,opac,n1,emin,emax,epi,ncn2,nlbin, &
     &    lopak,lun11,lpri)
!                                                                       
!     this routine sorts line and continuum emissivities
!     author: T. Kallman       
!     nb1 is a temporary bin number                                         
!                                                                       
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      integer nlbin(nrank,ncn)                                            
!     emissivities                                                 
      real(8) emis(2,n1) 
!     energy bins                                                       
      real(8) epi(ncn) 
!     line opacities                                                    
      real(8) opac(n1) 
!     line bin data
      real(8) ebln(n1) 
      integer nb1,mm,mm2,jkk1,jkk2,nbinc,n1,ncn2,lpri,lun11
      integer lopak
      real(8) ener,emin,emax
      logical done
!
!     in case of garbage skip

      if (lpri.ge.1) write (lun11,*)'testing in rlbin:',jkk1,opac(jkk1),&
     &   emis(1,jkk1)+emis(2,jkk1),ebln(jkk1)
      if (jkk1.le.0) return
!      if ((lopak.eq.0).and.(emis(1,jkk1)+emis(2,jkk1).lt.1.e-37)) return
      if ((opac(jkk1).lt.1.e-37).and.                                   &
     &   (emis(1,jkk1)+emis(2,jkk1).lt.1.e-37)) return
!      if ((lopak.eq.1).and.(opac(jkk1).lt.1.e-37)) return
      if ((ebln(jkk1).gt.emax).or.(ebln(jkk1).lt.emin)) return
!
!     check bin
      ener=12398.41/ebln(jkk1)
      nb1=nbinc(ener,epi,ncn2)
!      do mm=1,ncn2
!         write (lun11,*)mm,epi(mm)
!         enddo      
!
!     now determine rank in bin
      mm=0
      if (lpri.ne.0) then
      write (lun11,*)'in rlbin: jkk1=',jkk1,ebln(jkk1),emin,emax
      write (lun11,*)'emis,opac',emis(1,jkk1)+emis(2,jkk1),opac(jkk1)
      write (lun11,*)'bin:',ebln(jkk1),ener,nb1
      endif
      done=.false.
      do while (.not.done)
        mm=mm+1
        jkk2=nlbin(mm,nb1)
        if (jkk2.eq.0) done=.true.
        if (.not.done) then
          if ((lopak.eq.0).and.                                         &
     &      (emis(1,jkk1)+emis(2,jkk1).gt.emis(1,jkk2)+emis(2,jkk2)))   &
     &     done=.true.
          if ((lopak.eq.1).and.                                         &
     &     (opac(jkk1).gt.opac(jkk2))) done=.true.
          if (mm.ge.nrank) done=.true.
          endif
        if (lpri.ne.0)                                                  &
     &  write (lun11,*)'table:',mm,jkk2,emis(1,jkk2)+emis(2,jkk2),      &
     &   opac(jkk2),done
        enddo
      if (mm.ge.nrank) return
      do mm2=nrank-1,mm,-1
        nlbin(mm2+1,nb1)=nlbin(mm2,nb1)
        enddo
      if (lpri.ne.0) write (lun11,*)'saving:',mm,nb1,jkk1
      nlbin(mm,nb1)=jkk1
!                                                                       
      return
      end
