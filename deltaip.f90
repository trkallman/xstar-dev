      subroutine deltaip(nz,nn,t,xpx,deltae,lpri,lun11)
!
!     this routine calculates the IP shift based on results 
!        from Deprince et al.
!     Parameters are:
!       nz=nuclear charge
!       nn=number of bound electrons
!       t=temperature in 10^4K
!       xpx=nucleon number density
!       deltae=shift in eV
!
      implicit none
      real(8) zeff,xpx,t,emu,deltae
      integer nz,nn,lpri,lun11
!
      zeff=float(nz-nn+1)
      if (zeff.le.0) stop 'in deltaip: zeff error'
      emu=(5.295e-9)/((6.896e+2)*sqrt(t/xpx))
!
      if (emu.lt.1.e-2) return
      emu=min(emu,0.25)
      deltae=-26.98*emu*(zeff)+26.6*emu**1.89
!
      if (lpri.gt.1) write (lun11,*)'in deltaip',nz,nn,t,xpx,emu,deltae
!
      return
      end
