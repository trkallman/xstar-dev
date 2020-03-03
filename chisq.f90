      subroutine chisq(ajisb,indb,nindb,                                &
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
!           indb(2,ndb)=indeces for entries of ajisb
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
      real(8) ajisb(2,ndb)
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
!        write (lun11,*)'rate equations',nindb,ipmat
        do ll=1,nindb 
          mm=min(ipmat,indb(1,ll)) 
          nn=min(ipmat,indb(2,ll)) 
!          write (lun11,*)ll,mm,nn
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
!            write (lun11,92)ll,mm,nn,ajisb(1,ll),ajisb(2,ll),           &
!     &         x(nn),x(mm),sum(mm),term1,term2                          
!   92       format (1x,'     ',3i6,7(1pe13.5)) 
            endif 
          enddo 
        write (lun11,*)'dominant components:'
        write (lun11,*)'level, abundance, error, max rate in, index,',  &
     &     'max rate out, index'
        do mm=1,ipmat 
          if ((ll1mx(mm).ne.0).and.(ll2mx(mm).ne.0))                    &
     &     write (lun11,93)mm,x(mm),sum(mm),term1mx(mm),                &
     &       indb(2,ll1mx(mm)),term2mx(mm),indb(2,ll2mx(mm)),           &
     &       ll1mx(mm),ll2mx(mm)            
   93     format (1x,i6,3(1pe13.5),i6,1pe12.5,3i12) 
          enddo 
!                                                                       
      return 
      end                                           
