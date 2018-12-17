      function upsiln(k,eij,cc,ntem,cstr,temps,t,lpri,lun11) 
      use globaldata
      implicit none 
      real(8) cstr(50),temps(50),calc_spline 
      real(8)  eij, cc,      t 
      real(8)   upsiln, xt, y2(50), kte, ups, sups 
      integer k,ntem,lpri,lun11 
!                                                                       
!     this routine calculates upsilons for Burgess and Tully            
!     used by chianti                                                   
!     author:  T. Kallman                                               
!                                                                       
!     t = electron temperature in Kelvin                                
!     p# = spline knot values                                           
!     c = abscissa scale parameter                                      
!     k = transition type                                               
!     eij = transition energy (Ryd)                                     
!                                                                       
!     ;c                                                                
                                                                        
       kte=t/eij/1.57888d5 
                                                                        
                                                                        
      if ((k.EQ.1).OR.(k.EQ.4))                                         &
     &   xt=1 - log(cc)/(log(kte + cc))                                 
      if ((k.EQ.2).OR.(k.EQ.3).OR.(k.EQ.5).OR.(k.EQ.6))                 &
     & xt=kte / (kte +cc)                                               
       if ((ntem.ne.5).and.(lpri.gt.1)) write (lun11,*)'nsplines, ntem=' 
       call prep_spline(temps,cstr,ntem,y2) 
       sups=calc_spline(temps,cstr,y2,ntem,xt,lpri,lun11) 
      if (k.eq.1) ups=sups*log(kte + exp(1.)) 
      if (k.eq.2) ups=sups 
      if (k.eq.3) ups=sups/(kte+1.) 
      if (k.eq.4) ups=sups*log(kte+cc) 
      if (k.eq.5) ups=sups/(kte) 
      if (k.eq.6) ups=10.**sups 
       upsiln=ups 
      if (lpri.ne.0) write (lun11,*)'k=',k,xt 
      if (lpri.ne.0) write (lun11,*)'in upsiln:',t,kte,cc,xt,sups,ups 
!                                                                       
      return 
      END                                           