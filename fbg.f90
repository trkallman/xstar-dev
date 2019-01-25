      real(8) function fbg(u,gam) 
!                                                                       
!     Name: fbg.f90  
!     Description:  
!     this function computes the free-free gaunt factor.
!     After Kellogg, Baldwin and Koch 1975 Ap J 199, 299.
!     This routine needs to be tested and reevaluated before adoption
!     Parameter:
!         Input:                 
!         u=h nu/kt                                                        
!         gam=z**2 ry/kt, z=charge of scattering ion ry=rydberg constant  
!         Output:
!         fbg=gaunt factor
!     Dependencies:  none
!     Called by: none (was bremem)
!                                                                       
      use globaldata
      implicit none 
!                                                                       
      real(8) a,a1,a2,a3,ai,ak,born,g1,g2,gam,                           &
     &     gam1,gam2,gam3,p,power,t,u,u1,u2,expo                        
      real(8) u4 
      integer m,m1,n 
!      real(8)  t,ai,ak,u4                                               
      dimension a(6,7,3),gam2(6),gam3(6) 
      dimension a1(6,7),a2(6,7),a3(6,7) 
!                                                                       
      equivalence (a1(1,1),a(1,1,1)),(a2(1,1),a(1,1,2)),                &
     &             (a3(1,1),a(1,1,3))                                   
!                                                                       
      data gam2/.7783,1.2217,2.6234,4.3766,20.,70./ 
      data gam3/1.,1.7783,3.,5.6234,10.,30./ 
      data a1/1.001,1.004,1.017,1.036,1.056,1.121,1.001,                &
     &     1.005,1.017,1.046,1.073,1.115,.9991,1.005,                   &
     &     1.030,1.055,1.102,1.176,.9970,1.005,1.035,                   &
     &     1.069,1.134,1.186,.9962,1.004,1.042,1.100,                   &
     &     1.193,1.306,.9874,.9962,1.047,1.156,1.327,                   &
     &     1.485,.9681,.9755,.8363,1.208,1.525,1.955/                   
      data a2/.30290,.16160,.04757,.01300,.00490,-.00320,               &
     &     .49050,.21550,.08357,.02041,.00739,.00029,                   &
     &     .65400,.28330,.08057,.03257,.00759,-.00151,                  &
     &     1.0290,.39100,.12660,.05149,.01274,.00324,                   &
     &     .95690,.48910,.17640,.05914,.01407,-.00024,                  &
     &     1.2690,.75790,.32600,.10770,.02800,.00548,                   &
     &     1.3270,1.0170,1.3980,.20500,.06050,.00187/                   
      data a3/ - 1.3230,-.25400,-.01571,-.001000,-.000184,              &
     &     .00008,-4.7620,-.33860,-.03571,-.001786,-.000300,            &
     &     .00001,-8.3490,-.42060,-.02571,-.003429,-.000234,            &
     &     .00005,-13.231,-.59000,-.04571,-.005714,-.000445,            &
     &     -.00004,-7.6720,-.68520,-.06430,-.005857,-.000420,           &
     &     .00004,-7.1430,-.99470,-.12000,-.010070,-.000851,            &
     &     -.00004,-3.1750,-1.1160,-.84140,-.018210,-.001729,           &
     &     .00023/                                                      
!                                                                       
      gam1 = gam*1000. 
      if ( gam1.gt.100. ) then 
         power = -.134/(gam**.2097) 
         fbg = 1.5*(3.*u)**power 
         return 
      else 
         u2 = u**2 
!                                                                       
!*****compute born approximation gaunt factor                           
!                                                                       
         u1 = u/2. 
         t = u1/3.75 
         u4 = u1/2. 
         if ( u1.gt.2. ) then 
!                                                                       
            ak = 1.2533141 - .07832358/u4 + .02189568/u4**2 -           &
     &           .01062446/u4**3 + .00587872/u4**4 - .00251540/u4**5 +  &
     &           .00053208/u4**6                                        
            ak = ak/(expo(u1)*sqrt(u1)) 
         else 
            ai = 1.0 + 3.5156229*t**2 + 3.0899424*t**4 +                &
     &           1.2067492*t**6 + 0.2659732*t**8 + 0.0360768*t**10 +    &
     &           0.0045813*t**12                                        
            ak = -1.*log(u4)*ai - .57721566 + .42278420*u4**2 +         &
     &           .23069758*u4**4 + .0348859*u4**6 + .00262698*u4**8 +   &
     &           .00010750*u4**10 + .0000074*u4**12                     
         endif 
         born = .5513*expo(u1)*ak 
!                                                                       
!*****compute polymonial factor to multiply born approximation          
!                                                                       
         m=0 
         n=0 
         if ( gam1.ge.1. ) then 
            if ( u.ge..003 ) then 
               if ( u.le..03 ) n = 1 
               if ( (u.le..3) .and. (u.gt..03) ) n = 2 
               if ( (u.le.1.) .and. (u.gt..3) ) n = 3 
               if ( (u.le.5.) .and. (u.gt.1.) ) n = 4 
               if ( (u.le.15.) .and. (u.gt.5.) ) n = 5 
               if ( u.gt.15. ) n = 6 
               if ( gam1.le.1.7783 ) m = 1 
               if ( (gam1.le.3.) .and. (gam1.gt.1.7783) ) m = 2 
               if ( (gam1.le.5.6234) .and. (gam1.gt.3.) ) m = 3 
               if ( (gam1.le.10.) .and. (gam1.gt.5.6234) ) m = 4 
               if ( (gam1.le.30.) .and. (gam1.gt.10.) ) m = 5 
               if ( (gam1.le.100.) .and. (gam1.gt.30.) ) m = 6 
               m1 = m + 1 
               g1 = (a(n,m,1)+a(n,m,2)*u+a(n,m,3)*u2)*born 
               g2 = (a(n,m1,1)+a(n,m1,2)*u+a(n,m1,3)*u2)*born 
               p = (gam1-gam3(m))/gam2(m) 
               fbg = (1.0-p)*g1 + p*g2 
               return 
            endif 
         endif 
      endif 
!      fbg = born 
!                                                                       
      return 
      END                                           
