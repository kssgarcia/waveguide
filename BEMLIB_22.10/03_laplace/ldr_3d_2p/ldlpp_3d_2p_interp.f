      subroutine ldlpp_3d_2p_interp
     +
     +    (Iopt_int
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +    ,x4,y4,z4
     +    ,x5,y5,z5
     +    ,x6,y6,z6
     +
     +    ,vx1,vy1,vz1
     +    ,vx2,vy2,vz2
     +    ,vx3,vy3,vz3
     +    ,vx4,vy4,vz4
     +    ,vx5,vy5,vz5
     +    ,vx6,vy6,vz6
     +
     +    ,q1,q2,q3,q4,q5,q6
     +
     +    ,al,be,ga
     +    ,xi,eta
     +    ,x,y,z
     +    ,vx,vy,vz
     +    ,hs
     +    ,qint
     +    )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------------------
c   Utility of the point-dlp integrator:
c
c   Interpolates over an element for the:
c
c   position vector
c   normal vector
c   the scalar density q
c   the surface metric
c  
c   Iopt_int = 1 only the position vector
c              2 position vector etc
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c prepare
c---

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
c interpolate for the scalar density
c---

      qint = q1*ph1 + q2*ph2 + q3*ph3
     +      +q4*ph4 + q5*ph5 + q6*ph6


      if(Iopt_int.gt.1) then

c---
c interpolate for the normal vector
c---

      vx = vx1*ph1 +vx2*ph2 +vx3*ph3 +vx4*ph4 +vx5*ph5 +vx6*ph6
      vy = vy1*ph1 +vy2*ph2 +vy3*ph3 +vy4*ph4 +vy5*ph5 +vy6*ph6
      vz = vz1*ph1 +vz2*ph2 +vz3*ph3 +vz4*ph4 +vz5*ph5 +vz6*ph6

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
c  compute the raw normal vector and the
c  surface metric hs
c---

      vnxr = DyDxi * DzDet - DyDet * DzDxi
      vnyr = DzDxi * DxDet - DzDet * DxDxi
      vnzr = DxDxi * DyDet - DxDet * DyDxi

      hs  = sqrt(vnxr**2+vnyr**2+vnzr**2 )

      end if

c-----
c done
c-----

      return
      end
