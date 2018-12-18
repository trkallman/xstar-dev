      subroutine huntf(xx,n,x,jlo,lpri,lun11) 
!                                                                       
!     Name: huntf.f90  
!     Description:  
!           Searches in a list. 
!           assumes logarithmically spaced values
!
!     List of Parameters:
!           Input: 
!           xx(n):  list to be searched
!           n:  search length
!           x: search value
!           lpri:  print switch, 1=on, 0=off
!           lun11: logical unit number to write to
!           Output:
!           jlo:  index of found element
!
!     this version of hunt assumes equally spaced data in log           
!     author:  T. Kallman                                               
!                                                                       
      implicit none 
!                                                                       
      integer n,jlo,lpri,lun11 
      real(8) xx(n),x,xtmp,tst,tst2 
!                                                                       
      xtmp=max(x,xx(2)) 
      jlo=int((n-1)*log(xtmp/xx(1))/log(xx(n)/xx(1)))+1 
      if (jlo.lt.n) then 
        tst=abs(log(x/xx(jlo))) 
        tst2=abs(log(x/xx(jlo+1))) 
        if (tst2.lt.tst) jlo=jlo+1 
        endif 
      jlo=max(1,jlo) 
      jlo=min(n,jlo) 
      if (lpri.ne.0)                                                    &
     &  write (lun11,*)'in huntf',n,xx(1),xx(n),jlo,xx(jlo),x           
!                                                                       
      return 
      end                                           
