      subroutine unsavd(jkstep,ldir,                                    &
     &       lpri,iunit,iunit2,iunit3,iunit4,                           &
     &       t,p,r,rdel,delr,xcol,xee,xpx,zeta,                         &
     &       xilev,rnist,                                               &
     &       rcem,oplin,tau0,                                           &
     &       cemab,cabab,opakab,tauc,                                   &
     &       epi,ncn2,zrems,dpthc,opakc,rccemis,                        &
     &       lunlog,status)                                             
!                                                                       
!     this routine  saves only depths for iterative calculation         
!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      implicit none 
!                                                                       
!                                                                       
!     Allocation for passed parameters                                  
      real(8) r,delr,rdel, t, p, xcol,xee,xpx,zeta 
      integer status
      integer iunit,iunit2,iunit3,iunit4 
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn) 
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
      integer ncn2 
!     line opacities                                                    
      real(8) oplin(nnnl) 
      real(8) tauc(2,nnml) 
      real(8) cemab(2,nnml),opakab(nnml),cabab(nnml) 
      real(8) rccemis(2,ncn)
      real(8) xilev(nnml),rnist(nnml)
      real(8) tau0(2,nnnl), rcem(2,nnnl) 
      real(8) zrems(5,ncn) 
                                                                        
!     continuum optical depths                                          
      integer ldir,lpri,lunlog,jkstep 
      real(8) tau0d(2,nnnl),dpthcd(2,ncn),taucd(2,nnml) 
      integer lind1,lind2,kl,ll,nlyc,nry,nbinc
!                                                                       
      r=0. 
      delr=0. 
      t=0. 
      p=0. 
      if (status .gt. 0)call printerror(lunlog,status) 
      call rstepr(iunit,jkstep,r,delr,rdel,t,p,                         &
     &          xcol,xee,xpx,zeta,                                      &
     &          xilev,rnist,                                            &
     &          lunlog,lpri,status)                                     
      call rstepr2(iunit2,jkstep,r,delr,rdel,t,p,                       &
     &          xcol,xee,xpx,zeta,                                      &
     &          rcem,oplin,tau0d,                                       &
     &          lunlog,lpri,status)                                     
      call rstepr3(iunit3,jkstep,r,delr,rdel,t,p,                       &
     &          xcol,xee,xpx,zeta,                                      &
     &          cemab,cabab,opakab,taucd,                               &
     &          lunlog,lpri,status)                                     
      call rstepr4(iunit4,jkstep,r,delr,rdel,t,p,                       &
     &          xcol,xee,xpx,zeta,                                      &
     &          zrems,dpthcd,opakc,rccemis,                             &
     &          lunlog,lpri,status)                                     
      if (status .gt. 0)call printerror(lunlog,status) 
      lind1=1 
      lind2=2 
      if (ldir.gt.0) lind2=1 
      if (ldir.lt.0) lind1=2 
      do ll=lind1,lind2 
        do kl=1,nnnl 
          tau0(ll,kl)=tau0d(ll,kl) 
          enddo 
        do kl=1,ncn2 
          dpthc(ll,kl)=dpthcd(ll,kl) 
          enddo 
        do kl=1,nnml 
          tauc(ll,kl)=taucd(ll,kl) 
          enddo 
        enddo 
      nlyc=nbinc(13.7d0,epi,ncn2) 
      nry=nlyc+1 
      if (lpri.ne.0) write (lunlog,*)'in unsavd',rdel,t,tauc(1,25),     &
     &                ldir,dpthc(1,nry),dpthc(2,nry)                    
!                                                                       
      return 
      end                                           
