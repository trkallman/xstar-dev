      subroutine anl1(ni,nf,lf,iq,alm,alp,lpri,lun11) 
!                                                                       
!     Name: anl1.f90  
!     Description:  
!       this subroutine is used to calculate the values of the                
!       spontaneous transition rates for electric dipole transitions from     
!       level ni,lf+1 and ni,lf-1 to level nf,lf.                             
!       the transition probabilities (a values) are calculated                
!       using the gordon (1929) formula.                                      
!     author:  M. Bautista                                              
!     Parameters:
!         Input:
!         ni=first principal quantum number
!         nf=final principal quantum number
!         lf=final angular momentum
!         iq=ionic charge           
!         Output:
!         alm: forward a value                                     
!         alp: forward a value                                     
!     Dependencies: hgf, dfact
!     Called by: ucalc
!                                                                       
      implicit none 
!                                                                       
      real(8) y1,y2,x1,x2,x3,x4,x5,t 
      real(8) alm,alp,an,dum 
      integer ni,nf,lf,iq,lpri,lun11 
      integer li,n,l,np,ia1,ia2,ib,ic 
      real(8) x,rev,rn 
!                                                                       
      alm=0. 
!                                                                       
! **** for case a set lower limit of nf=1, for case b nf starts at 2    
!                                                                       
      do 40 li=lf-1,lf+1,2 
      if(li.lt.0) go to 40 
      if(lf.gt.li) go to 100 
      n=ni 
      np=nf 
      l=li 
      go to 101 
  100  n=nf 
      np=ni 
      l=lf 
  101  continue 
!                                                                       
       call dfact(n+l,x1) 
       call dfact(np+l-1,x2) 
       call dfact(2*l-1,x3) 
       call dfact(n-l-1,x4) 
       call dfact(np-l,x5) 
      ia1=-n+l+1 
      ia2=ia1-2 
      ib=-np+l 
      ic=2*l 
      x=-4.*n*np/((n-np)*(n-np)) 
       call hgf(ia1,ib,ic,x,y1) 
       call hgf(ia2,ib,ic,x,y2) 
      rev=abs(n-np) 
      rn=float(n+np) 
      t=(l+1)*log((4.*float(n*np)))+(rn-2*l-2)*log(rev) 
      t=t-log(4.e0)-rn*log(rn) 
      y1=abs((y1-y2*(rev/rn)**2)) 
      y1=log(y1)+t 
      t=2.*y1+x1+x2-2.*x3-x4-x5 
      t=exp(t) 
      an=2.6761e09*iq**4*max(li,lf)*t/(2.*li+1) 
         dum=(1./nf/nf-1./ni/ni)**3. 
      an=dum*an 
      if(li.lt.lf) alm=an 
      if(li.gt.lf) alp=an 
!                                                                       
   40  continue 
!                                                                       
      if (lpri.gt.1) then 
        write (lun11,*)'in anl1:',li,lf,t,ni,nf,iq 
        write (lun11,*) rn,n,np,rev,y1,y2,x 
        write (lun11,*) ia1,ia2,ib,ic,x1,x2,x3,x4,x5,l,an 
        endif 
!                                                                       
      return 
      END                                           
