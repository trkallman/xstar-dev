      subroutine savd(jkstep,ldir,                                      &
     &       lpri,iunit,iunit2,iunit3,iunit4,                           &
     &       np2,ncsvn,nlsvn,                                           &
     &       t,p,r,rdel,delr,xcol,xee,xpx,zeta,abel,                    &
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
      real(8) r,delr,rdel, t, p, xcol,xee,xpx,zeta 
      integer iunit,iunit2,iunit3,iunit4 
      integer nlsvn,ncsvn,np2
!     energy bins                                                       
      real(8) epi(ncn) 
!     continuum opacities                                               
      real(8) opakc(ncn) 
!     continuum optical depths                                          
      real(8) dpthc(2,ncn) 
!     continuum emissivities                                            
      real(8) rccemis(2,ncn) 
      integer ncn2 
!     line opacities                                                    
      real(8) oplin(nnnl) 
      real(8) tauc(2,nnml) 
      real(8) cemab(2,nnml),opakab(nnml),cabab(nnml) 
      real(8) xilev(nnml),rnist(nnml)
      real(8) tau0(2,nnnl), rcem(2,nnnl) 
      real(8) zrems(5,ncn) 
      real(8) abel(nl) 
!                                                                        
!     continuum optical depths                                          
      integer ldir,lpri,lunlog,jkstep 
      integer nlyc,nry,nbinc,status
!                                                                       
      if (lpri.ge.1)                                                    &
     & write (lunlog,*)'in savd',jkstep,iunit,iunit2,iunit3,iunit4      
      call fstepr(iunit,jkstep,r,delr,rdel,t,p,abel,                    &
     &          xcol,xee,xpx,zeta,                                      &
     &          np2,ncsvn,nlsvn,                                        &
     &          xilev,rnist,                                            &
     &          lunlog,lpri,status)                                     
      call fstepr2(iunit2,jkstep,r,delr,rdel,t,p,abel,                  &
     &          xcol,xee,xpx,zeta,                                      &
     &          np2,ncsvn,nlsvn,                                        &
     &          rcem,oplin,tau0,                                        &
     &          lunlog,lpri,status)                                     
      call fstepr3(iunit3,jkstep,r,delr,rdel,t,p,abel,                  &
     &          xcol,xee,xpx,zeta,                                      &
     &          np2,ncsvn,nlsvn,                                        &
     &          rnist,cemab,cabab,opakab,tauc,                          &
     &          lunlog,lpri,status)                                     
      call fstepr4(iunit4,jkstep,r,delr,rdel,t,p,abel,                  &
     &          xcol,xee,xpx,zeta,                                      &
     &          np2,ncsvn,nlsvn,                                        &
     &          epi,ncn2,zrems,dpthc,opakc,rccemis,                     &
     &          lunlog,lpri,status)                                     
      if (status .gt. 0)call printerror(lunlog,status) 
      nlyc=nbinc(13.7d0,epi,ncn2) 
      nry=nlyc+1 
      if (lpri.ne.0) write (lunlog,*)'in savd',rdel,t,tauc(1,25),       &
     &                ldir,dpthc(1,nry),dpthc(2,nry)                    
!                                                                       
      return 
      end                                            
