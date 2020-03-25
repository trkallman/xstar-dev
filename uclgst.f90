      subroutine uclgst(kdum,kresult,ierr)
      implicit none
      integer ierr, ll
      character*(*) kdum,kresult
      character(80) kres2
      character(1) ktmp
!
      ierr=0
!
      read (5,'(a80)')kres2
      do ll=1,80
        read (kres2(ll:ll),'(a1)')ktmp
        write(kresult(ll:ll),'(a1)')ktmp
        enddo
      return
      end
