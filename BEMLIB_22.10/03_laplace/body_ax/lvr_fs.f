      subroutine lvr_fs 
     +                  (Iopt
     +                  ,x,s
     +                  ,x0,s0
     +                  ,u,v
     +                  ,psi
     +                  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------------
c Velocity components (u, v) 
c and Stokes streamfunction (psi)
c due to a line vortex ring
c
c If Iopt.eq.1 the subroutine computed the velocity
c If Iopt.ne.1 the subroutine computed the velocity
c              and the Stokes stream function
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi = 3.14159 265358 D0
      pi4 = 4.0*pi

c---
c prepare
c---

      Dx  = x-x0
      Dxs = Dx**2

      rks = 4.0*s*s0/(Dxs+(s+s0)**2)

      call ell_int (rks,F,E)

      RJ30 = E/(1.0-rks)
      RJ31 = (-2.0D0*F + E*(2.0D0-rks)/(1.0D0-rks))/rks

      cf = 4.0D0/dsqrt((Dxs+(s+s0)**2)**3)

      RI30 = cf * RJ30
      RI31 = cf * RJ31

      u = (-s*RI31+s0*RI30)/pi4
      v =  Dx*RI31/pi4

      If(Iopt.eq.1) Go to 99
c---
c compute the Stokes stream function
c---

      rk   = sqrt(rks)
      RJ11 = ((2.0D0-rks)*F-2.0D0*E)/rks

      cf = 4.0D0/sqrt((Dxs+(s+s0)**2))

      RI11 = cf * RJ11
      psi  = rk*s*s0*RI11/pi4

c-----
c Done
c-----

 99   Continue

      Return
      End
