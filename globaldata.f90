      module globaldata
      implicit none 
!                                                                       
      include './PARAM' 
!                                                                       
!     global xstar data
!     master data
      TYPE :: master_data
        integer :: idat1(nidat1) ! integer data
        real(8) :: rdat1(nrdat1)  ! real data
        integer :: nptrs(nptt,ndat2) ! pointer data
        character(1) :: kdat1(nkdat1) ! character data
      END TYPE master_data
      TYPE(master_data) :: masterdata
      TYPE :: derived_pointers
        integer :: npar(ndat2)   !    pointers to master data
        integer :: npnxt(ndat2)  !    pointers to master data
        integer :: npfirst(ntyp) !    pointers to master data
        integer :: npfi(ntyp,nni)!    pointers to master data
        integer :: nplin(nnnl)   ! pointers to line data
        integer :: nplini(ndat2) ! pointers to line data
        integer :: npcon(nnml)
        integer :: npconi2(ndat2) 
        integer :: npconi(ndat2)
        integer :: npilev(nd,nni)
        integer :: npilevi(nnml)
        integer :: nlevs(nni) 
      END TYPE derived_pointers
      TYPE(derived_pointers) :: derivedpointers
      TYPE :: level_temp
        real(8) :: rlev(10,nd) 
        integer:: ilev(10,nd),nlpt(nd),iltp(nd) 
        character(1) :: klev(100,nd) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
!     compton heating data                                              
      real(8) decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp) 
      end module globaldata

