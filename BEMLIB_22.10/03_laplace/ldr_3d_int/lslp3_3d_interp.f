      subroutine lslp3_3d_interp
     +
     +    (Iopt_int
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +    ,x4,y4,z4
     +    ,x5,y5,z5
     +    ,x6,y6,z6
     +
     +    ,ztx1,zty1,ztz1
     +    ,ztx2,zty2,ztz2
     +    ,ztx3,zty3,ztz3
     +    ,ztx4,zty4,ztz4
     +    ,ztx5,zty5,ztz5
     +    ,ztx6,zty6,ztz6
     +
     +    ,al,be,ga
     +    ,xi,eta
     +    ,x,y,z
     +    ,hs
     +    ,zetx,zety,zetz
     +    )

c===============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===============================================

c-------------------------------------------
c Utility of an integrator:
c
c Interpolates over an element for the:
c
c position vector
c the surface metric
c the vectorial density zeta
c  
c Iopt_int = 1 only the position vector
c            2 position vector and rest of variables
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c--------
c prepare
c--------

      al1 = 1.0D0-al
      be1 = 1.0D0-be
      ga1 = 1.0D0-ga

      alal1 = al*al1
      bebe1 = be*be1
      gaga1 = ga*ga1

c---
c compute the basis functions
c---

      ph2 = xi *(xi -al+eta*(al-ga)/ga1)   /al1
      ph3 = eta*(eta-be+xi *(be+ga-1.0D0)/ga)/be1
      ph4 = xi *(1.0D0-xi-eta)/alal1
      ph5 = xi*eta          /gaga1
      ph6 = eta*(1.0D0-xi-eta)/bebe1
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c---
c interpolate the position vector
c---

      x = x1*ph1 + x2*ph2 + x3*ph3 + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3 + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3 + z4*ph4 + z5*ph5 + z6*ph6

c---
c interpolate for the vectorial density
c---

      zetx = ztx1*ph1 + ztx2*ph2 + ztx3*ph3
     +      +ztx4*ph4 + ztx5*ph5 + ztx6*ph6

      zety = zty1*ph1 + zty2*ph2 + zty3*ph3
     +      +zty4*ph4 + zty5*ph5 + zty6*ph6

      zetz = ztz1*ph1 + ztz2*ph2 + ztz3*ph3
     +      +ztz4*ph4 + ztz5*ph5 + ztz6*ph6

c---
c compute xi derivatives of basis functions
c---

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/ga1)/al1
      dph3 =  eta*(be+ga-1.0D0)/(ga*be1)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alal1
      dph5 =  eta/gaga1
      dph6 = -eta/bebe1
      dph1 = -dph2-dph3-dph4-dph5-dph6

c---
c compute xi derivatives of the position vector
c---

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c--
c compute eta derivatives of basis functions
c---

      dph2 =  xi*(al-ga)/(al1*ga1)
      dph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/be1
      dph4 = -xi/alal1
      dph5 =  xi/gaga1
      dph6 =  (1.0D0-xi-2.0D0*eta)/bebe1
      dph1 = -dph2-dph3-dph4-dph5-dph6

c---
c compute eta derivatives of the position vector
c---

      DxDet = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDet = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDet = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c---
c compute the raw normal vector and the
c surface metric hs
c---

      vnxr = DyDxi * DzDet - DyDet * DzDxi
      vnyr = DzDxi * DxDet - DzDet * DxDxi
      vnzr = DxDxi * DyDet - DxDet * DyDxi

      hs  = Dsqrt(vnxr**2+vnyr**2+vnzr**2 )

c-----
c done
c-----

      return
      end
