      subroutine sort4(n,ira,irb,irc) 
      parameter (ndat2=350000) 
      dimension ira(n),irb(n),irc(n),iwksp(ndat2) 
      real(8) ra(ndat2),rb(ndat2),rc(ndat2),wksp(ndat2) 
      do ll=1,n 
        ra(ll)=dfloat(ira(ll)) 
        rb(ll)=dfloat(irb(ll)) 
        rc(ll)=dfloat(irc(ll)) 
        enddo 
      call sort3(n,ra,rb,rc,wksp,iwksp) 
      do ll=1,n 
        ira(ll)=int(ra(ll)) 
        irb(ll)=int(rb(ll)) 
        irc(ll)=int(rc(ll)) 
        enddo 
      return 
      END                                           
