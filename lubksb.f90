      subroutine lubksb(a,n,np,indx,b,lun11,lpri) 
!                                                                       
!     Name: lubksb.f90
!     Description:
!       routine to solve system by back substitution
!       from numerical recipes
!     Parameters:
!       Input:
!         a(np,np):  coefficient matrix after lu decomposition
!         np:  first dimension of a
!         n: number of unknowns
!         indx:n
!         b(np): right hand side vector
!         lun11: logical unit number for printing
!         lpri: print switch
!     Output:
!          b(np):  vector of answers
!     Dependencies:  none
!     Called by:  leqt2f
!
      implicit none 
      integer i, ii, j, ll, n, np 
      integer lpri, lun11 
!                                                                       
      integer indx(np) 
      real(8)  a(np,np), b(np), sum 
!                                                                       
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in lubksb',n,np                                  
      ii = 0 
      do 100 i = 1 , n 
         ll = indx(i) 
!         write (lun11,*)'i,ll:',i,ll                                   
         sum = b(ll) 
         b(ll) = b(i) 
         if ( ii.ne.0 ) then 
            do 20 j = ii , i - 1 
               if (lpri.gt.1)                                           &
     &          write (lun11,*)'i,j,ii:',i,j,ii,sum                     
               sum = sum - a(i,j)*b(j) 
   20       continue 
         elseif ( sum.ne.0. ) then 
            ii = i 
         endif 
         if (lpri.gt.1)                                                 &
     &    write (lun11,*)'i,sum:',i,sum                                 
         b(i) = sum 
  100 continue 
      do 200 i = n , 1 , -1 
         sum = b(i) 
         if ( i.lt.n ) then 
            do 120 j = i + 1 , n 
               sum = sum - a(i,j)*b(j) 
  120       continue 
         endif 
         b(i) = sum/a(i,i) 
  200 continue 
      return 
      end                                           
