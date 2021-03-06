      subroutine calc_maxwell_rates(lun11,lpri,coll_type,min_T,max_T,  &
     &  tarr, om, dE,  T,                                               &
     & z,  degl,  degu,                                                 &
     & exc_rate, dex_rate,upsilon)                                      
!                                                                       
!     Name: calc_maxwell_rates.f90  
!     Description:  
!      Calculates Maxwellian-averaged collision strengths          
!      Calculates Maxwellian-averaged collision strengths for any         
!      transition, given input values taken from the APED (collision type,
!      temperature minimum/maximum, and data vectors Tarr and om.  Also   
!      requires the energy separating the two levels, the electron (or pro
!      temperature, the proton number, and the upper and lower level      
!      degeneracies.  Returns both the excitation and the deexcitation rat
!      author:  apec
!     Parameters:
!        Input:
!        lpri= print switch
!        lun11= logical unit for printing
!        coll_type The collision type, as defined in calc_maxwell_rates
!        min_T=The minimum temperature (in Kelvin) for which the     
!               transition is defined.                              
!        max_T=The maximum temperature (in Kelvin) for which the     
!               transition is defined.                                  
!        Tarr=array of temperatures
!        om=array of collision strengths                                
!        dE=The delta Energy of the transition, in keV               
!        T=Temperature, in Kelvin, of either the electrons or proton
!               involved in the transition.                           
!        Z=The # of protons in the target ion.                      
!        degl=The degeneracy of the lower level in the transition    
!        degu=The degeneracy of the upper level in the transition    
!        Output:
!        exc_rate=The calculated excitation rate coefficient, in cm^3
!        dex_rate=The calculated deexcitation rate coefficient, in cm
!     Called by: ucalc
!     Dependencies:  calc_kato, prepspline, calcspline, errmess, 
!       calc_sampson_s,calc_sampson_p,calc_sampson_h, exint1,expo,
!       interpol_huntd
!_______________________________________________________________        
!                                                                       
      real(8) calcspline 
      real(8) interpol_huntd 
      real(8) calc_kato 
      real(8) calc_sampson_s 
      real(8) calc_sampson_p 
      real(8) calc_sampson_h 
      real(8) exint1,expo 
      integer coll_type 
      real(8)  min_T 
      real(8) max_T 
      real(8) Tarr(100) 
      real(8) om(100) 
      real(8) dE 
      real(8) T 
      integer lun11,lpri 
      integer Z 
      real(8) degl 
      real(8) degu 
      real(8) exc_rate 
      real(8) dex_rate,pow 
!                                                                       
      integer jl, ju 
      integer calc_type 
      integer N_interp, num_spline 
      real(8) a,b,c,d,e 
      real(8) chi, chi_inv, upsilon, rate_coeff,  st, logT, expint_2 
      real(8) xs(5) 
      real(8) xs9(9) 
      data xs/0.00, 0.25, 0.50, 0.75, 1.00/ 
      data xs9/0.00, 0.125, 0.25, 0.375, 0.50, 0.675, 0.80, 0.925, 1.00/ 
      real(8) yp(5), yp9(9), spline_dat(9) 
! /*! \file calc_maxwell_rates.h                                        
!  * \brief Collisional excitation data types                           
!  * Data are stored in a number of different formats;                  
!  * this header lists all the currently-used versions                  
!  * which are found in calc_maxwell_rates.c                            
!  */                                                                   
      integer MAX_UPS 
      data MAX_UPS/20 / 
! /*      !< Length of temp/om arrays in collision structure */         
      integer MAX_CHI 
      data MAX_CHI/200. / 
! /*      !< Do not calculate rate if kT/dE greater than this.*/        
      real(8) KBOLTZ 
      data KBOLTZ/8.617385e-8  / 
! /*      !< in units of keV/K */                                       
      real(8) M_E 
      data M_E/2.7182818284590452354  / 
! /*      !< Euler e */                                                 
      real(8) UPSILON_COLL_COEFF 
      data UPSILON_COLL_COEFF/8.629e-6  / 
! /*      !< sqrt{2 pi / kB} hbar^2/m_e^{1.5} */                        
      integer E_UPSILON 
      data E_UPSILON/1    / 
! /*      !< Electron upsilon values (unitless) */                      
      integer E_RATE_COEFF 
      data E_RATE_COEFF/2 / 
! /*      !< Electron rate coefficient (cm^3/s) */                      
      integer P_UPSILON 
      data P_UPSILON/3    / 
! /*      !< Proton upsilon values (unitless) */                        
      integer P_RATE_COEFF 
      data P_RATE_COEFF/4 / 
! /*      !< Proton rate coefficient (cm^3/s) */                        
                                                                        
      integer BURGESSTULLY 
      data BURGESSTULLY/1  / 
!     /*      !< Burgess-Tully-type data*/                              
      integer CHIANTI_1 
      data CHIANTI_1/11   / 
! /*      !< CHIANTI pre-4.0 type 1 data (5 pt spline) */               
      integer CHIANTI_2 
      data CHIANTI_2/12   / 
! /*      !< CHIANTI pre-4.0 type 2 data (5 pt spline) */               
      integer CHIANTI_3 
      data CHIANTI_3/13   / 
! /*      !< CHIANTI pre-4.0 type 3 data (5 pt spline) */               
      integer CHIANTI_4 
      data CHIANTI_4/14   / 
! /*      !< CHIANTI pre-4.0 type 4 data (5 pt spline) */               
      integer CHIANTI_5 
      data CHIANTI_5/15   / 
! /*      !< CHIANTI pre-4.0 type 5 data (5 pt spline) */               
      integer CHIANTI_6 
      data CHIANTI_6/16   / 
! /*      !< CHIANTI pre-4.0 type 6 data (5 pt spline) */               
                                                                        
      integer CHIANTI4_1 
      data CHIANTI4_1/21   / 
! /*      !< CHIANTI 4.0 type 1 data (9 pt spline) */                   
      integer CHIANTI4_2 
      data CHIANTI4_2/22   / 
! /*      !< CHIANTI 4.0 type 2 data (9 pt spline) */                   
      integer CHIANTI4_3 
      data CHIANTI4_3/23   / 
! /*      !< CHIANTI 4.0 type 3 data (9 pt spline) */                   
      integer CHIANTI4_4 
      data CHIANTI4_4/24   / 
! /*      !< CHIANTI 4.0 type 4 data (9 pt spline) */                   
      integer CHIANTI4_5 
      data CHIANTI4_5/25   / 
! /*      !< CHIANTI 4.0 type 5 data (9 pt spline) */                   
      integer CHIANTI4_6 
      data CHIANTI4_6/26   / 
! /*      !< CHIANTI 4.0 type 6 data (9 pt spline) */                   
                                                                        
      integer SGC_1 
      data SGC_1/31   / 
! /*      !< Sampson, Goett and Clark (1983) S-type He-like data */     
      integer SGC_2 
      data SGC_2/32   / 
! /*      !< Sampson, Goett and Clark (1983) P-type He-like data */     
      integer SGC_3 
      data SGC_3/33   / 
! /*      !< Sampson, Goett and Clark (1983) S-type H-like data */      
      integer KATO_NAKAZAKI_1 
      data KATO_NAKAZAKI_1/41  / 
! /*      !< Kato and Nakazaki (1989), ADNDT 42, 313 */                 
      integer KATO_NAKAZAKI_2 
      data KATO_NAKAZAKI_2/42  / 
! /*      !< Kato and Nakazaki (1989), ADNDT 42, 313 */                 
! /* These must be spaced by at least MAX_UPS */                        
      integer INTERP_E_UPSILON 
      data INTERP_E_UPSILON/100     / 
! /*      !< Include both left & right boundaries */                    
      integer INTERP_P_UPSILON 
      data INTERP_P_UPSILON/200     / 
! /*      !< Include both left & right boundaries */                    
      integer INTERP_E_RATE_COEFF 
      data INTERP_E_RATE_COEFF/300  / 
! /*      !< Include both left & right boundaries */                    
      integer INTERP_P_RATE_COEFF 
      data INTERP_P_RATE_COEFF/400  / 
! /*      !< Include both left & right boundaries */                    
                                                                        
      integer INTERP_E_UPS_OPEN 
      data INTERP_E_UPS_OPEN/150    / 
! /*      !< Include neither boundary */                                
      integer INTERP_P_UPS_OPEN 
      data INTERP_P_UPS_OPEN/250    / 
! /*      !< Include neither boundary */                                
      integer INTERP_E_RATE_OPEN 
      data INTERP_E_RATE_OPEN/350   / 
! /*      !< Include neither boundary */                                
      integer INTERP_P_RATE_OPEN 
      data INTERP_P_RATE_OPEN/450   / 
! /*      !< Include neither boundary */                                
                                                                        
      integer INTERP_E_UPS_INC_MIN 
      data INTERP_E_UPS_INC_MIN/500   / 
! /*      !< Include only minimum; max is out */                        
      integer INTERP_P_UPS_INC_MIN 
      data INTERP_P_UPS_INC_MIN/600   / 
! /*      !< Include only minimum; max is out */                        
      integer INTERP_E_RATE_INC_MIN 
      data INTERP_E_RATE_INC_MIN/700  / 
! /*      !< Include only minimum; max is out */                        
      integer INTERP_P_RATE_INC_MIN 
      data INTERP_P_RATE_INC_MIN/800  / 
! /*      !< Include only minimum; max is out */                        
                                                                        
      integer INTERP_E_UPS_INC_MAX 
      data INTERP_E_UPS_INC_MAX/550   / 
! /*      !< Include only maximum; min is out */                        
      integer INTERP_P_UPS_INC_MAX 
      data INTERP_P_UPS_INC_MAX/650   / 
! /*      !< Include only maximum; min is out */                        
      integer INTERP_E_RATE_INC_MAX 
      data INTERP_E_RATE_INC_MAX/750  / 
! /*      !< Include only maximum; min is out */                        
      integer INTERP_P_RATE_INC_MAX 
      data INTERP_P_RATE_INC_MAX/850  / 
! /*      !< Include only maximum; min is out */                        
      integer PROTON_BT 
      data PROTON_BT/1001 / 
! /*      !< For Burgess-Tully Proton excitation rates */               
                                                                        
      calc_type = -1 
      chi  = dE / (KBOLTZ*T) 
      chi_inv = (KBOLTZ*T) / dE 
                                                                        
      logT = log10(T) 
      st = 0 
      upsilon = 0.0 
      rate_coeff = 0.0 
      if (lpri.ge.2)                                                    &
     &     write (lun11,*)'in calc_maxwell_rates,',t,dE,z,degl,degu,   &
     &                   kboltz,chi,calc_type,coll_type                 
                                                                        
      exc_rate = 0 
      dex_rate = 0 
                                                                        
!  Check to make sure temperature is inside valid range for this */     
!  transition otherwise return defaults set above.  */                  
      if (lpri.ge.2) write (lun11,*)T,min_T,max_T 
      if ((T .lt. min_T).or.(T .gt. max_T)) then 
!         call errmess(lun11,35,"calc_maxwell_rates, Bad transition ")  
        return 
        endif 
                                                                        
!     burgesstully=1                                                    
      if (coll_type .eq. BURGESSTULLY) then 
        expint_2 = exint1(chi,2) 
        a = om(1) 
        b = om(2) 
        c = om(3) 
        d = om(4) 
        e = om(5) 
        upsilon = a + b*chi*expint_2 + c*chi*(1-chi*expint_2)           &
     &     + d*(chi/2)*(1-chi*(1-chi*expint_2)) + e*expint_2            
        calc_type = E_UPSILON 
        if (lpri.ge.2)                                                  &
     &   write (lun11,*)'coll_type=BURGESSTULLY',chi,upsilon,calc_type  
        endif 
                                                                        
!     chianti_1=1                                                       
!     chianti_2=12                                                      
!     chianti_3=13                                                      
!     chianti_4=14                                                      
!     chianti_5=15                                                      
!     chianti_6=16                                                      
      if ((coll_type.eq.CHIANTI_1   ).or.                               &
     &     (coll_type.eq.CHIANTI_2   ).or.                              &
     &     (coll_type.eq.CHIANTI_3   ).or.                              &
     &     (coll_type.eq.CHIANTI_4   ).or.                              &
     &     (coll_type.eq.CHIANTI_5   ).or.                              &
     &     (coll_type.eq.CHIANTI_6   )) then                            
        c = om(1) 
                                                                        
        if (coll_type.eq.CHIANTI_1  ) st = 1 - (log(c)/log(chi_inv+c)) 
        if (coll_type.eq.CHIANTI_2  ) st = chi_inv/(chi_inv+c) 
        if (coll_type.eq.CHIANTI_3  ) st = chi_inv/(chi_inv+c) 
        if (coll_type.eq.CHIANTI_4  ) st = 1 - (log(c)/log(chi_inv+c)) 
        if (coll_type.eq.CHIANTI_5  ) st = chi_inv/(chi_inv+c) 
        if (coll_type.eq.CHIANTI_6  ) st = chi_inv/(chi_inv+c) 
                                                                        
        do mm=1,5 
          spline_dat(mm) = (om(mm)) 
          enddo 
        call prepspline(xs, spline_dat, 5, yp) 
        upsilon = calcspline(xs, spline_dat, yp, 5, st,lpri,lun11) 
                                                                        
        if (coll_type.eq.CHIANTI_1  ) upsilon=upsilon*log(chi_inv+M_E) 
!        if (coll_type .eq. CHIANTI_2)   /* Do nothing */               
        if (coll_type.eq.CHIANTI_3  ) upsilon = upsilon/(chi_inv+1) 
        if (coll_type.eq.CHIANTI_4  ) upsilon = upsilon*log(chi_inv+c) 
        if (coll_type.eq.CHIANTI_5  ) upsilon = upsilon/chi_inv 
        if (coll_type.eq.CHIANTI_6  ) upsilon = pow(10.d0,upsilon) 
        if (coll_type.eq.CHIANTI4_6 ) then 
          calc_type = P_UPSILON 
          else 
          calc_type = E_UPSILON 
          endif 
        if (lpri.ge.2)                                                  &
     &   write (lun11,*)'coll_type=CHIANTI 1-6',chi,upsilon,calc_type   
        endif 
                                                                        
!     chianti4_1=21                                                     
!     chianti4_2=22                                                     
!     chianti4_3=23                                                     
!     chianti4_4=24                                                     
!     chianti4_5=25                                                     
!     chianti4_6=26                                                     
      if ((coll_type.eq.CHIANTI4_1  ).or.                               &
     &     (coll_type.eq.CHIANTI4_2  ).or.                              &
     &     (coll_type.eq.CHIANTI4_3  ).or.                              &
     &     (coll_type.eq.CHIANTI4_4  ).or.                              &
     &     (coll_type.eq.CHIANTI4_5  ).or.                              &
     &     (coll_type.eq.CHIANTI4_6  )) then                            
        num_spline = int(om(1))
        c = om(2) 
        if (lpri.gt.1)                                                  &
     &    write (lun11,*)'in chianti4:',om(1),num_spline,om(2),c        
                                                                        
        if (coll_type.eq.CHIANTI4_1  )st=1 - (log(c)/log(chi_inv+c)) 
        if (coll_type.eq.CHIANTI4_2  )st=chi_inv/(chi_inv+c) 
        if (coll_type.eq.CHIANTI4_3  )st=chi_inv/(chi_inv+c) 
        if (coll_type.eq.CHIANTI4_4  )st=1 - (log(c)/log(chi_inv+c)) 
        if (coll_type.eq.CHIANTI4_5  )st=chi_inv/(chi_inv+c) 
        if (coll_type.eq.CHIANTI4_6  )st=chi_inv/(chi_inv+c) 
                                                                        
        do mm=1,num_spline 
           spline_dat(mm) = (om(2+mm)) 
           if (lpri.gt.1) write (lun11,*)mm,spline_dat(mm) 
           enddo 
        if (num_spline .eq. 5) then 
            call prepspline(xs, spline_dat, num_spline, yp) 
            upsilon = calcspline(xs, spline_dat, yp, num_spline, st, &
     &                            lpri,lun11)                           
          else 
            if (num_spline .eq. 9) then 
                if (lpri.gt.1) write (lun11,*)'before prep_spline:',    &
     &            xs9,spline_dat                                        
                call prepspline(xs9, spline_dat, num_spline, yp9) 
                if (lpri.gt.1) write (lun11,*)'before calc_spline:',    &
     &               st,yp9                                             
                upsilon=calcspline(xs9,spline_dat,yp9,num_spline,st, &
     &                            lpri,lun11)                           
              else 
                call errmess(lun11,71,"calc_maxwell_rates, Chianti  &
     &data with %d values, not 5 or 9 as expected")                          
                stop 
!     $           num_spline)                                           
                return 
              endif 
          endif 
                                                                        
        if (coll_type.eq.CHIANTI4_1) upsilon=upsilon*log(chi_inv+M_E) 
!        if (coll_type .eq. CHIANTI4_2)   /* Do nothing */              
        if (coll_type.eq.CHIANTI4_3) upsilon = upsilon/(chi_inv+1) 
        if (coll_type.eq.CHIANTI4_4) upsilon=upsilon*log(chi_inv+c) 
        if (coll_type.eq.CHIANTI4_5) upsilon = upsilon/chi_inv 
        if (coll_type.eq.CHIANTI4_6) upsilon = pow(10.d0,upsilon) 
        if (coll_type.eq.CHIANTI4_6) then 
            calc_type =P_UPSILON 
          else 
            calc_type=E_UPSILON 
          endif 
        if (lpri.ge.2)                                                  &
     &   write (lun11,*)'coll_type=CHIANTI4 1-6',chi,upsilon,calc_type  
        endif 
!                                                                       
                                                                        
!     SGC_1=31                                                          
      if (coll_type .eq.SGC_1       ) then 
! /* S-type transitions, He-like */                                     
        upsilon = calc_sampson_s(om, Z, T) 
        calc_type =E_UPSILON 
        endif 
                                                                        
!     SGC_2=32                                                          
      if (coll_type .eq. SGC_2) then 
! /* P-type transitions, He-like */                                     
        upsilon = calc_sampson_p(om, Z, T) 
        calc_type =E_UPSILON 
        endif 
                                                                        
!     SGC_3=33                                                          
      if (coll_type .eq. SGC_3) then 
!/* S-type transitions, H-like */                                       
        upsilon = calc_sampson_h(om, Z, T) 
        calc_type =E_UPSILON 
        endif 
                                                                        
!     KATO_NAKAZAKI_1=41                                                
      if (coll_type .eq.KATO_NAKAZAKI_1) then 
        upsilon = calc_kato(1, om, Z, T) 
        calc_type =E_UPSILON 
        endif 
                                                                        
!     KATO_NAKAZAKI_2=42                                                
      if (coll_type .eq.KATO_NAKAZAKI_2) then 
        upsilon = calc_kato(2, om, Z, T) 
        calc_type =E_UPSILON 
        endif 
                                                                        
!     PROTON_BT=1001                                                    
      if (coll_type .eq.PROTON_BT   ) then 
        a = om(1) 
        b = om(2) 
        c = om(3) 
!        /* 0.8 for the np/ne ratio, */                                 
        if ((logT .gt. om(3)).and.( logT .lt. om(4))) then 
          rate_coeff=0.8*pow(10.d0, (a + b*(logT) + c*logT*logT) ) 
          endif 
        calc_type =P_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=PROTON_BT' 
        endif 
                                                                        
!     INTERP_E_UPSILON=100                                              
!     MAX_UPS=20                                                        
      if (lpri.gt.1) write (lun11,*)'coll_type:',coll_type,             &
     &   INTERP_E_UPSILON,INTERP_E_UPSILON+MAX_UPS,MAX_UPS              
      if ((coll_type .ge.INTERP_E_UPSILON) .and.                        &
     &    (coll_type .le.INTERP_E_UPSILON + MAX_UPS)) then              
        N_interp = coll_type - INTERP_E_UPSILON 
        upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,       &
     &       lpri,lun11)                                                
        calc_type = E_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPSILON'      &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_E_UPS_OPEN=150                                             
      if ((coll_type .ge. INTERP_E_UPS_OPEN) .and.                      &
     &     (coll_type .le. INTERP_E_UPS_OPEN + MAX_UPS)) then           
        N_interp = coll_type - INTERP_E_UPS_OPEN 
        if ((T .gt. min_T).and.(T .lt. max_T))                          &
     &      upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,   &
     &       lpri,lun11)                                                
        calc_type = E_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPS_OPEN'     &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_E_UPS_INC_MIN=500                                          
      if ((coll_type .ge. INTERP_E_UPS_INC_MIN) .and.                   &
     &     (coll_type .le. INTERP_E_UPS_INC_MIN + MAX_UPS)) then        
        N_interp = coll_type - INTERP_E_UPS_INC_MIN 
        if (T .lt. max_T)                                               &
     &     upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,    &
     &       lpri,lun11)                                                
        calc_type = E_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPS_INC_MIN'  &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_E_UPS_INC_MAX=550                                          
      if ((coll_type .ge. INTERP_E_UPS_INC_MAX) .and.                   &
     &     (coll_type .le. INTERP_E_UPS_INC_MAX + MAX_UPS)) then        
        N_interp = coll_type - INTERP_E_UPS_INC_MAX 
        if (T .gt. min_T)                                               &
     &       upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,  &
     &       lpri,lun11)                                                
        calc_type = E_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPS_INC_MAX'  &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_UPSILON=200                                              
      if ((coll_type .ge. INTERP_P_UPSILON) .and.                       &
     &     (coll_type .le. INTERP_P_UPSILON + MAX_UPS)) then            
        N_interp = coll_type - INTERP_P_UPSILON 
        upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,       &
     &       lpri,lun11)                                                
        calc_type = P_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPSILON'      &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_UPS_OPEN=250                                             
      if ((coll_type .ge. INTERP_P_UPS_OPEN) .and.                      &
     &     (coll_type .le. INTERP_P_UPS_OPEN + MAX_UPS)) then           
        N_interp = coll_type - INTERP_P_UPS_OPEN 
        if ((T .gt. min_T).and.(T .lt. max_T))                          &
     &     upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,    &
     &       lpri,lun11)                                                
        calc_type = P_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPS_OPEN'     &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_UPS_INC_MIN=600                                          
      if ((coll_type .ge. INTERP_P_UPS_INC_MIN) .and.                   &
     &     (coll_type .le. INTERP_P_UPS_INC_MIN + MAX_UPS)) then        
        N_interp = coll_type - INTERP_P_UPS_INC_MIN 
        if (T .lt. max_T)                                               &
     &     upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,    &
     &       lpri,lun11)                                                
        calc_type = P_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPS_INC_MIN'  &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_UPS_INC_MIN=650                                          
      if ((coll_type .ge. INTERP_P_UPS_INC_MAX) .and.                   &
     &     (coll_type .le. INTERP_P_UPS_INC_MAX + MAX_UPS)) then        
        N_interp = coll_type - INTERP_P_UPS_INC_MAX 
        if (T .gt. min_T)                                               &
     &     upsilon = interpol_huntd(N_interp, Tarr, om, T, jl, ju,    &
     &       lpri,lun11)                                                
        calc_type = P_UPSILON 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPS_INC_MAX'  &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_E_RATE_COEFF=300                                           
      if ((coll_type .ge. INTERP_E_RATE_COEFF) .and.                    &
     &     (coll_type .le. INTERP_E_RATE_COEFF + MAX_UPS)) then         
        N_interp = coll_type - INTERP_E_RATE_COEFF 
        rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl, ju,    &
     &       lpri,lun11)                                                
        calc_type = E_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_COEFF'   &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_E_RATE_OPEN=350                                            
      if ((coll_type .ge. INTERP_E_RATE_OPEN) .and.                     &
     &     (coll_type .le. INTERP_E_RATE_OPEN + MAX_UPS)) then          
        N_interp = coll_type - INTERP_E_RATE_OPEN 
        if ((T .gt. min_T).and.(T .lt. max_T))                          &
     &     rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl, ju, &
     &       lpri,lun11)                                                
        calc_type = E_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_OPEN'    &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_E_RATE_INC_MIN=700                                         
      if ((coll_type .ge. INTERP_E_RATE_INC_MIN) .and.                  &
     &     (coll_type .le. INTERP_E_RATE_INC_MIN + MAX_UPS)) then       
        N_interp = coll_type - INTERP_E_RATE_INC_MIN 
        if (T .lt. max_T)                                               &
     &     rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl, ju, &
     &       lpri,lun11)                                                
        calc_type = E_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_INC_MIN' &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_E_RATE_INC_MAX=750                                         
      if ((coll_type .ge. INTERP_E_RATE_INC_MAX) .and.                  &
     &     (coll_type .le. INTERP_E_RATE_INC_MAX + MAX_UPS)) then       
        N_interp = coll_type - INTERP_E_RATE_INC_MAX 
        if (T .gt. min_T)                                               &
     &     rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl, ju, &
     &       lpri,lun11)                                                
        calc_type = E_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_INC_MAX' &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_RATE_COEFF=400                                           
      if ((coll_type .ge. INTERP_P_RATE_COEFF) .and.                    &
     &     (coll_type .le. INTERP_P_RATE_COEFF + MAX_UPS)) then         
        N_interp = coll_type - INTERP_P_RATE_COEFF 
        rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl, ju,    &
     &       lpri,lun11)                                                
        calc_type = P_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_COEFF'   &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_RATE_OPEN=450                                            
      if ((coll_type .ge. INTERP_P_RATE_OPEN) .and.                     &
     &     (coll_type .le. INTERP_P_RATE_OPEN + MAX_UPS)) then          
        N_interp = coll_type - INTERP_P_RATE_OPEN 
        if ((T .gt. min_T).and.(T .lt. max_T))                          &
     &       rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl,ju,&
     &       lpri,lun11)                                                
        calc_type = P_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_OPEN'    &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_RATE_INC_MIN<800                                         
      if ((coll_type .ge. INTERP_P_RATE_INC_MIN) .and.                  &
     &     (coll_type .le. INTERP_P_RATE_INC_MIN + MAX_UPS)) then       
        N_interp = coll_type - INTERP_P_RATE_INC_MIN 
        if (T .lt. max_T)                                               &
     &     rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl, ju, &
     &       lpri,lun11)                                                
        calc_type = P_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_INC_MIN' &
     &    ,jl,ju                                                        
        endif 
                                                                        
!     INTERP_P_RATE_INC_MAX=850                                         
      if ((coll_type .ge. INTERP_P_RATE_INC_MAX) .and.                  &
     &     (coll_type .le. INTERP_P_RATE_INC_MAX + MAX_UPS)) then       
        N_interp = coll_type - INTERP_P_RATE_INC_MAX 
        if (T .gt. min_T)                                               &
     &       rate_coeff = interpol_huntd(N_interp, Tarr, om, T, jl,ju,&
     &       lpri,lun11)                                                
        calc_type = P_RATE_COEFF 
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_INC_MAX' &
     &    ,jl,ju                                                        
        endif 
                                                                        
      if (calc_type .eq. -1) then 
        call errmess(lun11,43,"calc_trans_rates, Undefined collision &
     &type.")                                                             
        return 
        endif 
                                                                        
      if (calc_type .eq. E_UPSILON) then 
                                                                        
        if (upsilon .le. 0.)  then 
          upsilon = 0. 
! /* Negative upsilon is unphysical. */                                 
          endif 
                                                                        
        if (chi .lt. MAX_CHI) then 
          exc_rate=UPSILON_COLL_COEFF*upsilon*exp(-chi)/(sqrt(T)*degl) 
          dex_rate=UPSILON_COLL_COEFF * upsilon / (sqrt(T)*degu) 
          if (lpri.gt.1) write (lun11,*)'exc_rate,dex_rate:',           &
     &       upsilon,UPSILON_COLL_COEFF,exc_rate,dex_rate               
          endif 
        endif 
                                                                        
      if (calc_type .eq. P_UPSILON) then 
                                                                        
        if (upsilon .le. 0.)  then 
          upsilon = 0. 
! /* Negative upsilon is unphysical. */                                 
          endif 
                                                                        
        exc_rate = 0. 
        dex_rate = 0. 
        call errmess(lun11,59,"calc_rate,  Can't calculate collision &
     &strength for protons.")                                              
        return 
        endif 
                                                                        
      if ((calc_type .eq. P_RATE_COEFF)                                 &
     &  .or.(calc_type .eq. E_RATE_COEFF)) then                         
        if (rate_coeff .lt. 0) then 
          rate_coeff= 0. 
! /* Negative rate coefficient is unphysical. */                        
          endif 
        exc_rate = rate_coeff 
        dex_rate = rate_coeff*expo(chi)*(degl/degu) 
        endif 
!                                                                       
      return 
      END                                           
