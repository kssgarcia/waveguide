      subroutine hgf_2dc_fs
     +
     +  (Iopt
     +  ,delta
     +  ,x,y
     +  ,x0,y0
     +  ,Gr,Gi
     +  ,Grx,Gry
     +  ,Gix,Giy
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c-----------------------------------------------------
c Free-space Green's function of the Helmholtz equation
c
c Lapl G = - i delta^2 G
c
c where i is the imaginary unit and delta is real
c
c
c LEGEND
c ------
c
c Iopt =  1 compute only the Green's function
c      ne 1 compute the Green's function and gradient
c-----------------------------------------------------

      Implicit Double Precision (a-h,k,o-z)

c-----------
c  constants
c-----------

      pi  = 3.1415 92653 58979 32384  D0
      pi2 = 2.0D0*pi

c------ 
c begin
c------ 

      dx = x-x0
      dy = y-y0
      r  = Dsqrt(dx*dx+dy*dy)

      rd = r*delta 

      call ker_kei_0 
     +
     +  (Iopt
     +  ,rd
     +  ,ker_0
     +  ,kei_0
     +  ,ker_0_d
     +  ,kei_0_d
     +  )

      Gr = ker_0/pi2
      Gi = kei_0/pi2

      if(Iopt.eq.1) Go to 99

c---------------------
c compute the gradient
c---------------------

      rpi2 = r*pi2

      tx = delta*dx/rpi2
      ty = delta*dy/rpi2

      Grx = tx*ker_0_d
      Gry = ty*ker_0_d
      Gix = tx*kei_0_d
      Giy = ty*kei_0_d

c-----
c Done
c-----

  99  Continue

      Return
      End
