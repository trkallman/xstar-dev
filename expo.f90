      real(8) function expo(x) 
!                                                                       
!     Name: expo.f90  
!     Description:  
!        calculates e^x with limits, currently 600
!     Parameters:
!         Input: 
!         x=independent value
!         Output:
!         expo=e^x
!     Dependencies:  none
!     Called by:  ucalc
!
      use globaldata
      implicit none 
!                                                                       
      real(8) x,crit 
!                                                                       
      crit=600.                                                        
!       crit=60. 
!      expo=exp(x)                                                      
!      return                                                           
!      crit=60.                                                         
!      if (x.lt.-crit) then                                             
!        expo=1.e-24                                                    
!      else                                                             
!        xtmp=min(x,crit)                                               
        expo=exp(min(max(x,-crit),crit)) 
!      endif                                                            
!                                                                       
      return 
      END                                           
