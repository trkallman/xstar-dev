      real(8) function interpol_huntd(n,x,y,z,jl,ju,lpri,lun11) 

!     Name: interpol_huntd.f90  
!     Description:  
!           does interpolation
!           function needed by calc_maxwell_rates
!
!     List of Parameters:
!           Input: 
!           n:  number of values
!           x(n):  array of x values
!           y(n):  array of y values
!           z:  independent variable value
!           lpri: print switch
!           lun11: logical unit number to print to
!           Output
!           jl: index of lower bound x value
!           jl: index of upper bound x value
!           value:  interpolated value
!
!     Dependencies:
!
!     Called by:
!            calc_maxwell_rates
!
!     author:
!           apec
!
      integer n 
      real(8) x(n) 
      real(8) y(n) 
      real(8) z 
      real(8) grad,d1,df,f,pow 
      integer jl, jm, ju,lpri,lun11 
      logical inc 
!                                                                       
      inc = .false. 
      if (x(n-1) .gt. x(1)) inc=.true. 
      if (( inc .and.((z .lt. x(1)) .or. (z .gt. x(n-1)))) .or.         &
     & (.not.(inc).and.(z .gt. x(1) .or. z .lt. x(n-1)))) then          
       if (lpri.gt.1) then 
         write (lun11,*)"interpol_huntd: Asking for, min is,max is",z,&
     &        x(1),x(n-1)                                               
         write (lun11,*)"interpol_huntd: Cannot extrapolate" 
         endif 
        return 
        endif 
      jl = 0 
      ju = n 
      do while (ju - jl .gt. 1) 
        jm = (ju + jl) / 2 
        if ((z .gt. x(jm)) .eqv. inc) then 
          jl = jm 
        else 
          ju = jm 
        endif 
      enddo 
!     /* ------	Z is sandwiched between JL and JU ------ */             
      if ((x(jl) .gt. 0.).and.(x(ju).gt.0.).and.                        &
     &    (y(jl) .gt. 0.).and.(y(ju) .gt. 0.)) then                     
        grad = (log10(y(ju)) - log10(y(jl))) /                          &
     &      (log10(x(ju)) - log10(x(jl)))                               
        df = grad * (log10(z) - log10(x(jl))) 
        d1 = log10(y(jl)) + df 
        f = pow(10.d0, d1) 
      else 
        f = y(jl)+(y(ju)-y(jl))*((z-x(jl))/(x(ju)-x(jl))) 
      endif 
      interpol_huntd=f 
      if (lpri.gt.1)                                                    &
     &   write (lun11,*)"in interpol_huntd:",z,                       &
     &   x(1),x(n-1),jm,ju,jl,y(ju),y(jl),x(ju),x(jl),grad,d1,f         
      return 
      END                                           
