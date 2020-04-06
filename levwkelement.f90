      subroutine levwkelement(ml_element,lpri,ipmatsv,t,xee,xpx,     &
     &    leveltemp,lun11,rnise,mml,mmu)    
!
!     Name: levwkelement.f90  
!     Description:  
!           Calculates partition function for entire element
!
!     List of Parameters:
!           Input:
!           t: temperature in 10^4K
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           lpri: print switch, 1=on, 0=off
!           lun11: logical unit number for printing
!           nlev: number of levels for the ion
!           output:
!           rniss: lte level population
!           rnisse: lte level population relative to ground 
!                  with exponential removed
!           From Globaldata:
!           rlev(10,ndl):  real data for levels of this ion
!           ilev(10,ndl):  integer data for levels of this ion
!           nlpt(ndl):
!           iltp(ndl):
!

      use globaldata
      implicit none 
!                                                                       
!                                                                       
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,ndl) 
        integer:: ilev(10,ndl),nlpt(ndl),iltp(ndl) 
        character(1) :: klev(100,ndl) 
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
      real(8) rnisi(nd),rnisse(nd),rnise(nd), xileve(nd)
      real(8) ergsev,bk,t,                                              &
     &     bb
      integer lpri,lprisv,nlev,lun11,ml_element, ipmatsv, jk
      integer jkk_ion,lcon,ltyp,lrtyp,np1r,np1i,np1k,nrdt,nidt,nkdt
      integer ml_ion_data_type, ml_element_test, ml_ion, ml
!                                                                       
      real(8) xpx, xee, rnissum, rmidlog, scale,                        &
     &         rnismax, rnismin
      integer mm, klion, nnz, mml, mmu,nnzz,nnnn
                                                                        
      data ergsev/1.602197e-12/ 
      save ergsev
      data bk/1.38062e-16/ 
      save bk
      
      if (lpri.gt.1)                                                    &
     &  write (lun11,901)t,xee,xpx,mml,mmu                              
901   format (1x,'    in levwkelement, inputs:',3(1pe10.3),2i4)
!
      ipmatsv=0
!                                                                       
!     print element information
      call drd(ltyp,lrtyp,lcon,                                         &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_element,             &
     &            0,lun11)                                        
      nnz=masterdata%idat1(np1i+nidt-2)
      if (lpri.gt.1) then 
             write (lun11,902)jk,ml_element,nnz,                        &
     &          (masterdata%kdat1(np1k-1+mm),mm=1,min(8,nkdt))
902           format (1x,'  element:',3(i12,1x),8(1a1))
        endif 

      ml_ion_data_type=12
!     step thru ions
      ml_ion=derivedpointers%npfirst(ml_ion_data_type)
      rnissum=0.
      do while (ml_ion.ne.0)
!
!       test if element belongs to parent of ion
        ml_element_test=derivedpointers%npar(ml_ion)
        if (lpri.gt.1)                                                  &
     &   write (lun11,*)'ml_element_test=',ml_element_test,ml_element
        if (ml_element_test.eq.ml_element) then
!
!         get ion index
          call drd(ltyp,lrtyp,lcon,                                     &
     &            nrdt,np1r,nidt,np1i,nkdt,np1k,ml_ion,                 &
     &            0,lun11)                                        
          jkk_ion=masterdata%idat1(np1i+nidt-1)
          klion=masterdata%idat1(np1i)
          nlev=derivedpointers%nlevs(jkk_ion)
          nnzz=masterdata%idat1(np1i+1)
          nnnn=nnzz-masterdata%idat1(np1i)+1
!
!         test for ion in range
          if ((klion.ge.mml).and.(klion.le.mmu)) then
!
            if (lpri.gt.1)                                              &
     &            write (lun11,903)jkk_ion,ml_ion,klion,                &
     &               (masterdata%kdat1(np1k+mm-1),mm=1,nkdt)
903             format (1x,'      ion:',3(i12,1x),8(1a1))
            lprisv=lpri
            call calc_rates_level_lte(jkk_ion,lpri,lun11,t,xee,xpx,     &
     &              nnzz,nnnn,leveltemp,nlev)
            call levwk(rnisi,rnisse,bb,lpri,nlev,t,xee,xpx,             &
     &              leveltemp,lun11)      
            lpri=lprisv
            if (lpri.gt.1) write (lun11,*)'after calc_rates_level_lte'
            if (lpri.gt.1) write (lun11,*)'filling rnise'
            if (klion.eq.mml) rnise(1+ipmatsv)=rnisi(1)
            do mm=2,nlev
              if (klion.gt.mml) then
                rnise(mm+ipmatsv)=rnise(mm+ipmatsv-1)                   &
     &             *rnisi(mm)/(1.d-97+rnisi(mm-1))
                else
                rnise(mm+ipmatsv)=rnisi(mm)
                endif
              if (lpri.gt.1)                                            &
     &        write (lun11,9022)mm,(leveltemp%klev(ml,mm),ml=1,20)      &
     &          ,leveltemp%rlev(1,mm),leveltemp%rlev(2,mm),             &
     &           xileve(mm+ipmatsv),rnise(mm+ipmatsv),                  &
     &           rnise(mm+ipmatsv-1),rnisi(mm),rnisi(mm-1),ipmatsv
 9022         format (4x,i4,1x,20a1,7(1pe10.3),i4) 
              enddo
!
            else
!
            do mm=1,nlev
              rnise(mm+ipmatsv)=0.
              enddo
!
!           end of test for ion in range
            endif

          ipmatsv=ipmatsv+nlev-1
!
!         end of test if element belongs to parent of ion
          endif
!
!       end of step thru ions
        ml_ion=derivedpointers%npnxt(ml_ion)
        ml_element_test=derivedpointers%npar(ml_ion)
        enddo
!
!     fully stripped
      rnise(ipmatsv+1)=rnise(ipmatsv)                                   &
     &             *rnisi(nlev)/(1.d-97+rnisi(nlev-1))
      mm=1
      if (lpri.gt.1)                                                    &
     &        write (lun11,9022)nlev,(leveltemp%klev(ml,nlev),ml=1,20)  &
     &          ,leveltemp%rlev(1,nlev),leveltemp%rlev(2,nlev),         &
     &           xileve(mm+ipmatsv),rnise(mm+ipmatsv),                  &
     &           rnise(mm+ipmatsv-1),rnisi(nlev),rnisi(nlev-1),ipmatsv
      ipmatsv=ipmatsv+1
!
      rnissum=0.
      do mm=1,ipmatsv
        rnissum=rnissum+rnise(mm)
        enddo
      if (lpri.gt.1) write (lun11,*)'partition function',rnissum
      do mm=1,ipmatsv
        if ((lpri.gt.1).and.(rnise(mm).gt.1.d-97))                      &
     &     write (lun11,*)mm,rnise(mm),rnise(mm)/(1.d-97+rnissum)
        rnise(mm)=rnise(mm)/(1.d-97+rnissum)
        enddo
!
!     not scaling
      return
!
!     now try a fancy trick to avoid too small numbers
      rnismin=1.e+10
      rnismax=0.
      do mm=1,ipmatsv
        rnismin=min(rnismin,rnise(mm))
        rnismax=max(rnismax,rnise(mm))
        enddo
      if ((rnismax.lt.1.d-97).or.(rnismin.lt.1.d-97)) stop 'range error'
      rmidlog=(log10(rnismax)+log10(rnismin))/2.
      scale=10.**rmidlog
      if (lpri.gt.1) write (lun11,*)'scaled partition function',scale
      do mm=1,ipmatsv
        rnise(mm)=rnise(mm)/scale
        if (lpri.gt.1) write (lun11,*)mm,rnise(mm)
        enddo

      return
      end
