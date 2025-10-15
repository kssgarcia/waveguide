      subroutine lgf_ax_w
     +
     +  (Iopt
     +  ,Ign      ! Green/Neumann flag
     +  ,x,s
     +  ,x0,s0
     +  ,wall
     +  ,G
     +  ,Gx,Gs
     +  )

c==========================================
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the license agreement.
c==========================================

c-----------------------------------------
c Axisymmetric Green's function of Laplace's equation
c in a semi-infinite domain bounded by a plane wall.
c
c The wall is located at x = wall.
c
c    Iopt =  1 compute only the Green's function
c         ne 1 compute Green's function and gradient
c
c    Ign  = 1  Compute the Green's function
c    Ign  = 2  Compute the Neumann  function
c
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi

c-------- 
c prepare
c-------- 

      sign = 1.0D0                   ! Green
      If(Ign.eq.2) sign = -1.0D0     ! Neumann

      Dx   = x-x0
      Dxs  = Dx*Dx
      ss0s = (s+s0)**2

      rks = 4.0D0*s*s0/(Dxs+ss0s)

      call ell_int (rks,F,E)

c-------------------- 
c primary singularity
c------------------- 

      RJ10 = F
      den  = Dsqrt(Dxs+ss0s)
      RI10 = 4.0D0*RJ10/den

      G = RI10/pi4

      If(Iopt.eq.1) Go to 96

c---
c compute I30 and I31
c---

      rksc = 1.0D0-rks
      RJ30 = E/rksc
      RJ31 = (-2.0D0*F+(2.0D0-rks)*E/rksc)/rks
      cf   = 4.0D0/den**3
      RI30 = cf*RJ30
      RI31 = cf*RJ31

c---------
c gradient
c---------

      Gx = - dx * RI30
      Gs = - s*RI30+s0*RI31

      Gx = Gx/pi4
      Gs = Gs/pi4

  96  continue

c-------------
c Image system
c-------------

      Dx  = x-2.0D0*wall+x0
      Dxs = Dx*Dx

      rks = 4.0D0*s*s0/(Dxs+ss0s)

      call ell_int (rks,F,E)

      RJ10 = F
      den  = dsqrt(Dxs+ss0s)
      RI10 = 4.0D0*RJ10/den

      G = G - sign * RI10/pi4

      If(Iopt.eq.1) Go to 99

c--------------------
c compute I30 and I31
c--------------------

      rksc = 1.0D0-rks
      RJ30 = E/rksc
      RJ31 = (-2.0D0*F+(2.0D0-rks)*E/rksc)/rks
      cf   = 4.0D0/den**3
      RI30 = cf*RJ30
      RI31 = cf*RJ31

c---------
c gradient
c---------

      Gxi = - dx * RI30
      Gsi = - s*RI30+s0*RI31

      Gxi = Gxi/pi4
      Gsi = Gsi/pi4

      Gx = Gx - sign * Gxi
      Gs = Gs - sign * Gsi

c-----
c Done
c-----

  99  Continue

      Return
      End
