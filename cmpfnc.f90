      real(8) function cmpfnc(ee,sxx,lun11,lpri) 
!                                                                       
!     Name: cmpfnc.f90  
!     Description:  
!       this routine is used in the relativistic compton calculation      
!       searches in coheat.dat table and returns compton heating-cooling
!       author:  T. Kallman
!     Parameters:
!           Input:
!           ee=photon energy in m_ec^2
!           sxx=kT/m_ec^2
!           lun11=logical unit number for printing
!           lpri=print switch
!           Output:
!           cmpfnc=Delta(E)/E for compton heating-cooling
!     Dependencies: hunt3
!     Called by:  comp2
!                                                                       
      use globaldata
!                                                                       
      implicit none 
!                                                                       
      real(8) ee,eetp,sxx,ddedsx,ddede,dele,delsx,sxtp 
      integer lun11,lpri,mm,ll,mmm1,llm1,nc2 
!                                                                       
!     Not used                                                          
      integer javi 
                                                                        
      javi=lpri 
      lpri=javi 
!                                                                       
      eetp=ee 
      sxtp=sxx 
      nc2=ncomp 
      cmpfnc=0. 
      if (eetp.gt.1.e-4) then 
        call hunt3(ecomp,nc2,eetp,mm,0,lun11) 
        call hunt3(sxcomp,nc2,sxtp,ll,0,lun11) 
!       if ((mm.gt.1).and.(ll.gt.1)) then                               
               mm=max(2,min(ncomp,mm)) 
               ll=max(2,min(ncomp,ll)) 
               mmm1 = mm - 1 
               llm1 = ll - 1 
               ddedsx = (decomp(ll,mm)-decomp(llm1,mm)                  &
     &                  +decomp(ll,mmm1)-decomp(llm1,mmm1))             &
     &                     /(2.*(sxcomp(ll)-sxcomp(llm1)))              
               ddede = (decomp(ll,mm)-decomp(ll,mmm1)                   &
     &                  +decomp(llm1,mm)-decomp(llm1,mmm1))             &
     &                 /(2.*(ecomp(mm)-ecomp(mmm1)))                    
               dele = ee - ecomp(mmm1) 
               delsx = sxx - sxcomp(llm1) 
               cmpfnc = ddedsx*delsx + ddede*dele + decomp(llm1,mmm1) 
        else 
               cmpfnc = 4.*sxx - ee 
        endif 
<<<<<<< HEAD
!      if (lpri.gt.0)                                                   
=======
!      if (lpri.ne.0)                                                   
>>>>>>> 2d75308c63b9789458ce092c697c7853fcdde44a
!     $ write (lun11,*)'in cmpfnc',ee,sxx,ll,mm,cmpfnc
!                                                                       
      return 
      end                                           
