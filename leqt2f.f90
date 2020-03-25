      subroutine leqt2f(a,m,n,np,b,idgt,wkarea,ier,lun11,lpri) 
!                                                                       
!     Name: leqt2f.f90
!     Description:
!       solve a system of equations using numerical recipes routines
!     Parameters:
!       Input:
!         a(np,np):  coefficient matrix
!         m: number of right hand sides (not used)
!         n: number of unknowns
!         b(np): right hand side vector
!         idgt: (not used)
!         wkarea(1): work area (not used)
!         ier:  error flag
!         lun11: logical unit number for printing
!         lpri: print switch
!     Output:
!          b(np):  vector of answers
!     Dependencies:  ludcmp, lubksb, mprove
!     Called by:  msolvelucy
!
      use globaldata
      implicit none 
!                                                                       
      integer indx(nd),ier,lun11,lpri,n,np,m,npp 
      real(8) a(np,np),b(np),wkarea(1) 
      real(8)  ao(ndss,ndss),bo(ndss),btmp,tmp,sum,errmx,err 
      real(8)  an(ndss,ndss),bn(ndss),d,tmpmx 
      integer mm,ll2,mmmx,mmmxo,jk,idgt,kl 
!                                                                       
!     n had better be less than nd                                      
!                                                                       
!     Not used                                                          
      integer javi 
      real(8) javir 
      javi=m 
!      m=javi                                                           
      javi=ier 
      javir=wkarea(1) 
!      wkarea(1)=javir                                                  
      javi=idgt 
!                                                                       
      do jk=1,n 
        bo(jk)=dble(b(jk)) 
        bn(jk)=dble(b(jk)) 
        do kl=1,n 
           an(jk,kl)=dble(a(jk,kl)) 
           ao(jk,kl)=dble(a(jk,kl)) 
           enddo 
        enddo 
!                                                                       
      npp=ndss 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'before ludcmp',n,npp,np                          
      call ludcmp(an,n,npp,indx,d,lun11,lpri) 
      npp=ndss 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'after ludcmp',n,npp                              
      call lubksb(an,n,npp,indx,bn,lun11,lpri) 
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'after lubksb'                                    
      npp=ndss 
      call mprove(ao,an,n,npp,indx,bo,bn,lun11,lpri) 
      if (lpri.gt.2)                                                    &
     & write (lun11,*)'after mprove',n,npp,np                           
!                                                                       
!        check the solution                                             
         if (lpri.gt.2) write (lun11,*)'checking the solution' 
         errmx=0. 
         do  ll2=1,n 
          sum=0. 
          tmpmx=0. 
          mmmx=0 
          mmmxo=0 
          do  mm=1,n 
            btmp=bn(mm) 
            tmp=dble(a(ll2,mm))*max(0.d0,btmp) 
            if (abs(tmp).ge.tmpmx) then 
              mmmxo=mmmx 
              mmmx=mm 
              tmpmx=max(tmpmx,abs(tmp)) 
              endif 
            sum=sum+tmp 
            enddo 
          sum=sum-dble(b(ll2)) 
          err=sum/max(1.d-24,tmpmx) 
          errmx=max(errmx,abs(err)) 
          if (lpri.gt.2) write (lun11,9246)ll2,bn(ll2),tmpmx,sum,err,   &
     &                                     mmmx,mmmxo                   
 9246     format (1h ,i4,4e12.4,2i4) 
          enddo 
!                                                                       
      do jk=1,n 
         if (lpri.gt.2)                                                 &
     &    write (lun11,*)jk,b(jk)                                       
         if (bn(jk).lt.1.d-36) bn(jk)=0. 
         if (bn(jk).gt.1.d+36) bn(jk)=1.d+36 
         b(jk)=(bn(jk)) 
         enddo 
!                                                                       
         if (lpri.gt.2)                                                 &
     &    write (lun11,*)'leaving leqt'                                 
!      ier=0                                                            
!      wkarea(1)=0.                                                     
!      idgt=0                                                           
!                                                                       
      return 
      end                                           
