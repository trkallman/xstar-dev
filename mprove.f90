                                                                        
      subroutine mprove(a,alud,n,np,indx,b,x,lun11,lpri) 
!
!     Name: mprove.f90
!     Description:
!       routine to improve system solution
!       from numerical recipes
!     Parameters:
!       Input:
!         a(np,np):  coefficient matrix
!         a(np,np):  coefficient matrix after lu decomposition
!         np:  first dimension of a
!         n: number of unknowns
!         indx:n
!         b(np): right hand side vector
!         lun11: logical unit number for printing
!         lpri: print switch
!     Output:
!          x(np):  vector of answers
!     Dependencies:  lubksb
!     Called by:  leqt2f
!
      implicit none 
      integer np, n,nl 
      parameter (nl=10000) 
      real(8)  a(np,np), alud(np,np), b(n), x(n), r(nl), sdp 
      integer indx(n), lun11, lpri, i, j 
!
      if (lpri.gt.1)                                                    &
     & write (lun11,*)'in mprove:'                                      
      do 12 i=1,n 
        sdp=-b(i) 
        do 11 j=1,n 
          sdp=sdp+a(i,j)*x(j) 
   11   continue 
        r(i)=sdp 
   12 continue 
      call lubksb(alud,n,np,indx,r,lun11,lpri) 
      do 13 i=1,n 
        x(i)=x(i)-r(i) 
   13 continue 
!                                                                       
      return 
      end                                            
