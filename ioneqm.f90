      subroutine ioneqm(z,a,s,n,m,l,lpri,lun11) 
      implicit none 
!                                                                       
!                                                                       
!     Name: ioneqm.f90  
!     Description:  
!          computes ionization equilibrium                      
!          solves a system of ionization equations, attempting               
!          to avoid overflow problems.                                       
!     author:  T. Kallman                                               
!
!     List of Parameters
!         Input:
!         z(m): ionization rates
!         a(m): recombination rates
!         s(m): ion fractions
!         m: length of rate vectors
!         n: number of ions
!         l: minimum ion stage of interest (should be 1)
!         lpri:  print switch
!         lun11:  logical unit number for printing
!
!         Dependencies:
!         none
!
!         Called by:
!         istruc.f90
!                                                                       
      integer m, n 
      real(8)  z(m),a(m),s(n),q(31) 
      real(8)  eps,delt,pl,suml,tst 
      real(8)  sumg, pg 
      integer i,j,jk,jmax,k,l,ll,lpri,mmn 
      integer mmx,lun11 
      integer lprisv 
!                                                                       
      data eps/1.e-6/ 
      data delt/1.e-28/ 
!                                                                       
      lprisv=lpri 
      if (lpri.ge.1) lpri=2 
      if (lpri.ge.2) write (lun11,*)'in ioneqm' 
!                                                                       
!     initialize                                                        
      do jk = 1,n 
         s(jk) = 0. 
         enddo 
!                                                                       
      if ( lpri.ge.2 ) write (lun11,99001) 
!                                                                       
!     form naive ratio                                                  
      do j = 1,m 
         q(j) = a(j)/(z(j)+delt) 
         enddo 
!                                                                       
!     step thru and search for max. q value                             
      jk = l 
  300 jk = jk + 1 
      if ( (jk.lt.n) .and. (q(jk-1).lt.1.) ) goto 300 
      jmax = jk 
!                                                                       
      if ( lpri.ge.2 ) write (lun11,99002) n,m,l,jmax 
!                                                                       
!     step forwards                                                     
      suml = 0. 
      if ( jmax.ne.n ) then 
         pl = 1. 
         mmx = jmax - 1 
  350    mmx = mmx + 1 
         pl = pl/(q(mmx)+delt) 
         suml = suml + pl 
         tst = pl/(suml+delt) 
         if ( lpri.ge.2 ) write (lun11,99003) mmx,q(mmx),suml,pl 
         if ( (tst.gt.eps) .and. (mmx.lt.m) ) goto 350 
      endif 
!                                                                       
!     step backwards                                                    
      sumg = 0. 
      if ( jmax.ne.l ) then 
         pg = 1. 
         mmn = jmax 
  400    mmn = mmn - 1 
         pg = pg*q(mmn) 
         sumg = sumg + pg 
         tst = pg/(sumg+delt) 
         if ( lpri.ge.2 ) write (lun11,99004) mmn,q(mmn),sumg,pg 
         if ( (tst.gt.eps) .and. (mmn.gt.l) ) goto 400 
      endif 
!                                                                       
!                                                                       
      s(jmax) = 1./(1.+suml+sumg) 
      if ( jmax.ne.n ) then 
         do j = jmax,mmx 
            s(j+1) = s(j)/(q(j)+delt) 
            enddo 
        endif 
!                                                                       
      if ( jmax.ne.l ) then 
         k = jmax - mmn 
         do i = 1,k 
            j = jmax - i 
            s(j) = s(j+1)*q(j) 
            enddo 
        endif 
!                                                                       
      if ( lpri.ge.2 ) write (lun11,99005) (ll,q(ll),s(ll),ll=l,n) 
!                                                                       
      lpri=lprisv 
!                                                                       
      return 
99001 format (' ',' in ioneqm ') 
99002 format (' ',' n,m,l,jmax --',4i4) 
99003 format (' ','in greater than loop, j,q,sum,p --',i4,3e12.4) 
99004 format (' ','in less than loop, j,q,sum,p --',i4,3e12.4) 
99005 format (' ',i4,2e12.4) 
      end                                           
