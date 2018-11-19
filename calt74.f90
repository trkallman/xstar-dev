    subroutine calt74(temp,np,xse,xss,nddd,np1r,rate,alpha) 
!                                                                       
!   Takes coefficients in data type 74 and any given radiation flux     
!   array given in xse(i) (energiesin eV) and xss(i) (flux values)      
!   and return the the analytic integral over resonances in the         
!   cross sections represented by delta functions.                      
!   The routine also returns the DR recombination coefficient (in       
!   s-1cm-3) for the given value of temp (in Kelvins). alpha MUST       
!   be mutiplied by the stadistical of the recombined state and         
!   divided by that of the recombining state.                           
!                                                                       
!   np is the number of points xse() and nd is the number of real       
!   values in dtype74()                                                 
!      author: M. Bautista                                              
!                                                                       
      use globaldata
       implicit none 
!                                                                       
       integer nddd 
       real(8) temp,rate,alpha,te,ryk,ry,xt,x,hgh,                       &
     &      rm,xsec,factor                                              
       integer np,m,i,ipos,ip,np1r 
       real(8) xse(np),xss(np) 
!                                                                       
       te=temp*1.38066e-16 
       ryk=4.589343e+10 
       ry=13.605692 
       m=(nddd-1)/2 
!                                                                       
       xt=masterdata%rdat1(np1r-1+1) 
       x=masterdata%rdat1(np1r-1+2) 
       hgh=masterdata%rdat1(np1r-1+2+m) 
       alpha=0. 
       if (x/ryk/te.lt.40.) then 
        alpha=exp(-x/ryk/te)*(x+xt)*(x+xt)*hgh 
       endif 
       do i=2,m 
        x=masterdata%rdat1(np1r-1+1+i)*ry 
        hgh=masterdata%rdat1(np1r-1+1+i+m) 
        x=x/ry 
        if (x/ryk/te.lt.40.) then 
         alpha=alpha+exp(-x/ryk/te)*(x+xt)*(x+xt)*hgh 
        endif 
       enddo 
        factor=213.9577e-9 
        alpha=alpha*factor/(te**1.5)/ryk/ryk 
!                                                                       
       if (xse(np).lt.(x+xt)*ry) then 
        rate=0.e0 
        return 
       endif 
!                                                                       
       x=(masterdata%rdat1(np1r-1+2)+masterdata%rdat1(np1r-1+1))*ry 
       i=np/2 
    5  if (xse(i).ge.x) then 
        i=i-1 
        goto 5 
       endif 
       i=i-1 
   10  i=i+1 
       if(xse(i).lt.x.and.xse(i+1).ge.x)then 
        ipos=i 
       else 
        goto 10 
       endif 
!                                                                       
       rm=(xss(ipos+1)-xss(ipos))/(xse(ipos+1)-xse(ipos)) 
       xsec=xss(ipos)+rm*(x-xse(ipos)) 
       rate=xsec*masterdata%rdat1(np1r-1+2+m) 
       do  i=2,m 
        x=(masterdata%rdat1(np1r-1+1+i)+masterdata%rdat1(np1r-1+1))*ry 
        hgh=masterdata%rdat1(np1r-1+1+i+m) 
        ip=ipos 
   20   if (xse(ip).lt.x) then 
         ip=ip+1 
         goto 20 
        endif 
        ip=ip-2 
        ip=ip+1 
        if (ip.le.np) then 
          ipos=ip 
          rm=(xss(ipos+1)-xss(ipos))/(xse(ipos+1)-xse(ipos)) 
          xsec=xss(ipos)+rm*(x-xse(ipos)) 
          rate=rate+xsec*hgh 
          endif 
        enddo 
!                                                                       
        rate=rate*4.752e-22 
!                                                                       
      return 
      END                                           
