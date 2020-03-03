!_______________________________________________________________        
                                                                        
!_______________________________________________________________        
      real(8) function calc_sampson_h(om,Z,Te) 
!                                                                        
!     Name: calc_sampson_h.f90  
!     Description:  
!       routines to implement rates from 
!       Sampson, Goett, & Clark, ADNDT 29, 467 
!       author:  apec
!     Parameters:
!       Input:
!       om=array of coefficients
!       z=ion charge
!       Te=temperature in kelvin
!       Output:
!       calc_sampson_h
!     Dependencies:  exintn
!     Called by: calc_maxwell_rates
!
      real(8) om(7) 
      integer Z 
      real(8) Te 
                                                                        
      real(8) result 
      real(8) kT, y,exintn
      real(8) dE, z2s, a, c0,c1,c2,c_sw 
      real(8) a1, a1y, E1y, Ery, Er1y, Zeff2 
      real(8) term 
      real(8) KBOLTZ 
      data KBOLTZ/8.617385e-8  / 
! /*      !< in units of keV/K */                                       
                                                                        
      dE   = om(1) 
      z2s  = om(2) 
      a    = om(3) 
      c0   = om(4) 
      c1   = om(5) 
      c2   = om(6) 
      c_sw = om(7) 
      Zeff2 =  Z*Z 
                                                                        
      kT  = KBOLTZ*Te 
      y   = dE/kT 
      a1  = a+1 
      a1y = a1*y 
      E1y = exintn(y,-1.d0,1) 
      Ery = exintn(a1y,-1.d0,1) 
      Er1y = exintn(a1y,Ery,2) 
!      /* This is E_2(a1y) */                                           
                                                                        
      term = (c1*Ery + c2*Er1y/a1) 
      if (term.gt.0) then 
          result = c0 + 1.333*z2s*E1y*exp(y) + y*exp(a1y)*term 
        else 
          result = c0 + 1.333*z2s*E1y*exp(y) 
        endif 
                                                                        
      result = result*2*c_sw/Zeff2 
                                                                        
      calc_sampson_h=result 
      return 
      END                                           
