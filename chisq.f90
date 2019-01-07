      subroutine chisq(ajisb,cjisb,indb,nindb,                          &
     &   ipmat,x,lun11,lpri)                                            
!                                                                       
!     Name: chisq.f90  
!     Description:  
!           evaluate statistical equilibrium solution
!           print out diagonal and largest off diagonal element 
!             in each row
!           no Output returned.  result is printed
!     List of Parameters:
!           Input:
!           ajisb(2,ndb)=elements of rate array (s^-1)
!           cjisb(2,ndb)=elements of heating-cooling rate array (erg s^-1 cm^-3)
!           indb(2,ndb)=indeces for entries of ajisb,cjisb
!           ipmat=number of levels
!           lpri=print switch
!           lun11=logical unit number for printing
!           x=vector of level populations
!           Output:none
!     Dependencies:  none
!     called by:  func
!
      use globaldata
!                                                                       
      real(8) ajisb(2,ndb),cjisb(ndb) 
      integer indb(2,ndb),ll1mx(nd),ll2mx(nd) 
      real(8) x(nd),sum(nd),                                             &
     &       term1mx(nd),term2mx(nd),term1,term2                        
!                                                                       
!        check the solution                                             
         write (lun11,*)'checking the solution' 
        do mm=1,ipmat 
          sum(mm)=0. 
          term1mx(mm)=0. 
          term2mx(mm)=0. 
          ll1mx(mm)=0 
          ll2mx(mm)=0 
          enddo 
        write (lun11,*)'rate equations' 
        do ll=1,nindb 
          mm=min(ipmat,indb(1,ll)) 
          nn=min(ipmat,indb(2,ll)) 
          if (mm.ne.nn) then 
            term1=ajisb(1,ll)*x(nn) 
            term2=ajisb(2,ll)*x(mm) 
            sum(mm)=sum(mm)+term1-term2 
            if (term1.gt.term1mx(mm)) then 
               term1mx(mm)=term1 
               ll1mx(mm)=ll 
               endif 
            if (term2.gt.term2mx(mm)) then 
               term2mx(mm)=term2 
               ll2mx(mm)=ll 
               endif 
            write (lun11,*)ll,mm,nn,ajisb(1,ll),ajisb(2,ll),            &
     &         x(nn),x(mm),sum(mm),term1,term2                          
 9246       format (1h ,i4,4e12.4,2i4) 
            endif 
          enddo 
        write (lun11,*)'sums' 
        do mm=1,ipmat 
          if ((ll1mx(mm).ne.0).and.(ll2mx(mm).ne.0))                    &
     &    write (lun11,*)mm,x(mm),sum(mm),term1mx(mm),                  &
     &       indb(2,ll1mx(mm)),term2mx(mm),indb(2,ll2mx(mm))            
          enddo 
!                                                                       
         if (lpri.gt.2)                                                 &
     &    write (lun11,*)'leaving leqt'                                 
!                                                                       
      return 
      end                                           
