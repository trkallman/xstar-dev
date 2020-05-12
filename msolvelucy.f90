      subroutine msolvelucy(ajisb,cjisb,cjisb2,indb,nindb,nsup,nspmx,&
     &   ipmat,x,ht,cl,ht2,cl2,niter,nit2,nit3,nitmx,nitmx2,lun11,lpri)    
!                                                                       
!     Name: msolvelucy.f90
!     Description:
!       solves a linear system using the iterative technique of 
!       lucy 2001 MNRAS 326 95
!     author:  T. Kallman     
!     Parameters:
!        Input:
!         ajisb(2,ndb)=entries to coefficient matrix
!         cjisb(ndb)=heating-cooling rates
!         indb(2,ndb)=index arracy to coefficient matrix
!         nindb=number of entries to indb
!         nsup=number of superlevels
!         nspmx=max number of superlevels
!         ipmat=number of levels
!         nitmx=max number of iterations
!         nitmx2=max number of superlevel iterations
!         lun11=logical unit number for printing
!         lpri=print switch
!        Output:
!         x(nd)=solution vector
!         ht= heating rates
!         cl= cooling rates
!         niter=number of iterations 
!         nit2=number of iterationson superlevels
!         nit3=number of fixed point iterations
!     Depedencies: leqt2f
!     Called by:  calc_hmc_all
!                                               
!                                                                       
!     solves lucy iteration                                             
      use globaldata
!     author:  T. Kallman                                               
      implicit none 
!                                                                       
      real(8) ajisb(2,ndb),cjisb(ndb),cjisb2(ndb)
      integer indb(2,ndb) 
      real(8) x(nd),xo(nd),xoo(nd) 
      real(8) bmatsup(ndss),ajissup(ndss,ndss),p(ndss),rr(nd) 
      real(8) cjissup(ndss,ndss) 
      real(8) cjissup2(ndss,ndss) 
      real(8) wkarea(1) 
      real(8) tt1, crit, crit2, diff, tt2, diff2, diffs
      real(8) riu(nds),ril(nds),rui(nds),rli(nds),xm,tst,cl 
      real(8) ht, clp, htp,ht2,cl2, eps, eps2
      integer nsup(nd), mm, ipmat, ll 
      integer lpril, lpri, lprisv, idgt, ier 
      integer lun11, niter, nit3, nitmx, nspmx 
      integer nn, nsp, ngood, nspm, nspn, nspcon 
      integer nit2, nitmx2, m2, nindb 
!                                                                       
!                                                                       
!       step thru levels, and form calculate superlevel quantities      
      lprisv=lpri 
      call remtms(tt1) 
!      lpri=0                                                           
      if (lpri.ge.1)                                                    &
     & write (lun11,*)'in msolvelucy',lpri,nindb,lun11                  
!
!     these are parameters which control convergence
!     crit is the value of diff required to converge
      crit=1.d-2                                                    
!     crit2 is the value of diff2 required to converge
      crit2=1.d-2
!     nitmx is passed in and is the maximum number of 
!        lu decomposition iterations allowed
!     nitmx2 is passed in and is the maximum number of 
!        fixed point iterations allowed
!     eps is used in calculating which populations to include in diff
      eps=1.e-30
!     eps2 is used in calculating which populations to include in diff2
      eps2=1.e-30
      diff=1. 
      niter=0 
      nit3=0 
      if (lpri.gt.1) write (lun11,*)'rate matrix'
      do ll=1,nindb 
          mm=min(ipmat,indb(1,ll)) 
          nn=min(ipmat,indb(2,ll))    
          nspm=nsup(mm) 
          nspn=nsup(nn) 
          if (lpri.gt.1) write (lun11,92)ll,mm,nn,nspm,nspn,            &
     &              rr(mm),rr(nn),ajisb(1,ll),ajisb(2,ll)               
   92     format (1x,'     ',5i12,6(1pe13.5)) 
          enddo 
      do while ((diff.gt.crit).and.(niter.lt.nitmx)) 
        niter=niter+1 
        if (lpri.gt.1) write (lun11,*)'iteration=',niter 
        if (lpri.gt.1) write (lun11,*)'initial populations:' 
        do mm=1,ipmat 
          xo(mm)=x(mm) 
          if (lpri.gt.1) write (lun11,*)mm,x(mm),nsup(mm) 
          enddo 
        do mm=1,nspmx 
          p(mm)=0. 
          bmatsup(mm)=0. 
          do nn=1,nspmx 
            ajissup(mm,nn)=0. 
            cjissup(mm,nn)=0. 
            enddo 
          enddo 
        do mm=1,ipmat 
          nsp=nsup(mm) 
          p(nsp)=p(nsp)+x(mm) 
          enddo 
        call remtms(tt2) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)'before constucting matrix',abs(tt2-tt1)       
        tt1=tt2 
        if (lpri.gt.1)                                                  &
     &      write (lun11,*)'constucting the condensed matirx:'          
        ngood=0 
        do mm=1,ipmat 
          nspm=nsup(mm) 
          rr(mm)=x(mm)/(1.d-48+p(nspm)) 
          if (p(nspm).le.1.d-36) rr(mm)=1.
!         nb a test
!          rr(mm)=1.
          enddo 
        do ll=1,nindb 
          mm=min(ipmat,indb(1,ll)) 
          nn=min(ipmat,indb(2,ll))    
          nspm=nsup(mm) 
          nspn=nsup(nn) 
!          if (lpri.ge.1) write (lun11,*)ll,mm,nn,nspm,nspn,ipmat
          if ((nspn.ne.nspm)                                            &
     &         .and.(nspn.ne.0).and.(nspm.ne.0)                         &
     &         .and.((abs(ajisb(1,ll)).gt.1.d-48)                       &
     &           .or.(abs(ajisb(2,ll)).gt.1.d-48))) then                
              ajissup(nspm,nspn)=ajissup(nspm,nspn)                     &
     &            +(ajisb(1,ll))*rr(nn)                                 
              ajissup(nspm,nspm)=ajissup(nspm,nspm)                     &
     &            -(ajisb(2,ll))*rr(mm)                                 
!              ajissup(nspn,nspm)=ajissup(nspn,nspm)                    
!     $            +(ajisb(2,ll))*rr(nn)                                
!              ajissup(nspn,nspn)=ajissup(nspn,nspn)                    
!     $            -(ajisb(1,ll))*rr(mm)                                
              cjissup(nspm,nspn)=cjissup(nspm,nspn)                     &
     &            +(cjisb(ll))*rr(mm)                                   
              cjissup2(nspm,nspn)=cjissup2(nspm,nspn)                   &
     &            +(cjisb2(ll))*rr(mm)                                   
              ngood=ngood+1 
              if (lpri.gt.1) write (lun11,91)ll,mm,nn,nspm,nspn,        &
     &              rr(mm),rr(nn),ajisb(1,ll),ajisb(2,ll),              &
     &              ajissup(nspm,nspn),ajissup(nspm,nspm)
   91         format (1x,'used ',5i12,7(1pe13.5)) 
            endif 
          enddo 
        call remtms(tt2) 
        if (lpri.gt.1)                                                  &
     &     write (lun11,*)'after constucting matrix',abs(tt2-tt1),      &
     &                       ngood                                      
        if (lpri.gt.1) then 
          write (lun11,*)'the condensed populations:' 
          do nsp=1,nspmx 
            write (lun11,*)nsp,p(nsp) 
            enddo 
          write (lun11,*)'the condensed matrix:' 
          do nspm=1,nspmx 
            do nspn=1,nspmx 
              if (abs(ajissup(nspm,nspn)).gt.1.d-48)                    &
     &         write (lun11,*)nspm,nspn,ajissup(nspm,nspn)              
              enddo 
            enddo 
          endif 
!        put in number conservation                                     
!         nspcon=1                                                      
         nspcon=nspmx 
         do mm=1,nspmx 
           ajissup(nspcon,mm)=1. 
           bmatsup(mm)=0. 
           enddo 
        bmatsup(nspcon)=1. 
        lpril=0 
        call remtms(tt1) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)'before leqt',abs(tt2-tt1)                     
        call leqt2f(ajissup,1,nspmx,ndss,bmatsup,idgt,wkarea,ier,    &
     &                      lun11,lpril)                                
         call remtms(tt2) 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*)'after leqt',abs(tt2-tt1)                      
        if (lpri.gt.2) write (lun11,*)'the new condensed populations:' 
        do mm=1,nspmx 
          p(mm)=bmatsup(mm) 
          if (lpri.gt.1) write (lun11,*)mm,p(mm) 
          enddo 
        if (lpri.gt.1) write (lun11,*)'new populations' 
        do mm=1,ipmat 
          nsp=nsup(mm) 
          x(mm)=rr(mm)*p(nsp) 
          if (lpri.gt.1) write (lun11,*)mm,nsp,rr(mm),x(mm) 
          enddo 
        nit2=0 
        diff2=10. 
        do while ((nit2.lt.nitmx2).and.(diff2.ge.crit2)) 
          nit2=nit2+1 
          nit3=nit3+1 
          if (lpri.gt.1) write (lun11,*)'before calculate new x(mm)',   &
     &                                   nit2,nit3                      
          call remtms(tt2) 
          if (lpri.gt.2)                                                &
     &    write (lun11,*)'in diff2 loop',abs(tt2-tt1)                   
          tt1=tt2 
          do mm=1,ipmat 
            riu(mm)=0. 
            rui(mm)=0. 
            ril(mm)=0. 
            rli(mm)=0. 
            enddo 
          if (lpri.gt.3) write (lun11,*)'the riu calculation' 
          do ll=1,nindb 
            mm=indb(1,ll) 
            nn=min(ipmat,indb(2,ll)) 
            if (nn.gt.mm) then 
                riu(mm)=riu(mm)+abs(ajisb(2,ll)) 
                rui(mm)=rui(mm)+abs(ajisb(1,ll))*x(nn) 
                if (lpri.gt.3) write (lun11,*)ll,mm,nn,ajisb(2,ll),     &
     &           ajisb(1,ll),x(nn),riu(mm),rui(mm)                      
              endif 
            enddo 
          if (lpri.gt.3) write (lun11,*)'the ril calculation' 
          do ll=1,nindb 
            mm=indb(1,ll) 
            nn=min(ipmat,indb(2,ll)) 
            if (nn.lt.mm) then 
!               I hope the indeces are in the right order here          
                ril(mm)=ril(mm)+abs(ajisb(2,ll)) 
                rli(mm)=rli(mm)+abs(ajisb(1,ll))*x(nn) 
                if (lpri.gt.3) write (lun11,*)ll,mm,nn,ajisb(2,ll),     &
     &           ajisb(1,ll),x(nn),ril(mm),rli(mm)                      
              endif 
            enddo 
          do mm=1,ipmat 
            xoo(mm)=x(mm) 
            x(mm)=(rli(mm)+rui(mm))/(ril(mm)+riu(mm)+1.d-24) 
            if (lpri.gt.3) write (lun11,*)mm,riu(mm),rui(mm),ril(mm),   &
     &                                     rli(mm),x(mm),xoo(mm)        
            enddo 
          xm=0. 
          do mm=1,ipmat 
            xm=xm+x(mm) 
            enddo 
          if (lpri.gt.3) write (lun11,*)'new and old populations',      &
     &                      xm                                        
          do mm=1,ipmat 
            x(mm)=x(mm)/(1.d-24+xm) 
            enddo 
          m2=1 
          diff2=0. 
          tst=0. 
          do while ((diff2.lt.1.e+3)                                    &
     &           .and.(m2.le.ipmat).and.(tst.lt.1.e+3))                 
            if (lpri.gt.3) write (lun11,*)m2,x(m2),xoo(m2),xo(m2),      &
     &            diff2                                                 
            tst=1. 
            if (x(m2).gt.eps2) tst=xoo(m2)/x(m2) 
            diffs=(tst-1.)*(tst-1.) 
            diff2=diff2+diffs
            if ((lpri.gt.1).and.(diffs.gt.1.))                          &
     &         write (lun11,*)'big diff2:',m2,diffs,diff2,nit2 
            m2=m2+1 
            enddo 
          if (lpri.gt.1) write (lun11,*) 'diff2=',diff2,nit2 
          enddo 
        diff=0. 
        m2=1 
        do while ((m2.le.ipmat).and.(diff.lt.1.e+3)) 
          tst=xo(m2)*(1.d-30) 
          if ((diff.lt.1.e+10).and.(x(m2).gt.eps))                      &
     &      then
            diffs=(min(1.d+10,(xo(m2)-x(m2))/(xo(m2)+x(m2))))**2     
            diff=diff+diffs
            if ((lpri.gt.1).and.(diffs.gt.1.))                          &
     &            write (lun11,*)'big diff:',m2,x(m2),xo(m2),diffs,diff                                                  
            endif
          m2=m2+1 
          enddo 
        if (lpri.gt.2) write (lun11,*) 'diff=',diff 
      enddo 
!                                                                       
      if (lpri.ge.1) write (lun11,*)'heating-cooling in msolvelucy:' 
      cl=0. 
      ht=0. 
      cl2=0. 
      ht2=0. 
      do ll=1,nindb 
        mm=min(indb(1,ll),ipmat) 
        nn=min(indb(2,ll),ipmat) 
        if (cjisb(ll).gt.0.) then 
              cl=cl+x(mm)*cjisb(ll) 
            else 
              ht=ht-x(mm)*cjisb(ll) 
            endif 
        if (cjisb2(ll).gt.0.) then 
              cl2=cl2+x(mm)*cjisb2(ll) 
            else 
              ht2=ht2-x(mm)*cjisb2(ll) 
            endif 
!          if ((lpri.ge.1).and.(abs(cjisb(ll)).gt.1.d-24))               &
          if (lpri.ge.1)                                                &
     &         write (lun11,981)ll,mm,nn,x(mm),cjisb(ll),cjisb2(ll),    &
     &           ht,cl,ht2,cl2          
  981          format (1x,3i6,8(1pe11.3)) 
        enddo 
      if (lpri.ge.1)write (lun11,*)'heating-cooling:',ht-cl, ht2-cl2

      go to 9090 
      if (lpri.gt.2) write (lun11,*)'heating-cooling superlevels:' 
      clp=0. 
      htp=0. 
      do mm=1,nspmx 
        do nn=1,nspmx 
          if (cjissup(mm,nn).gt.0.) then 
              clp=clp+p(mm)*cjissup(mm,nn) 
            else 
              htp=htp-p(mm)*cjissup(mm,nn) 
            endif 
          if ((lpri.gt.2).and.(abs(cjissup(mm,nn)).gt.1.d-24))          &
     &         write (lun11,*)mm,nn,p(mm),cjissup(mm,nn),htp,clp     
          enddo 
        enddo 
      ht=htp 
      cl=clp 
 9090 continue 
!
      lpri=lprisv 
!                                                                       
      return 
      end                                           
