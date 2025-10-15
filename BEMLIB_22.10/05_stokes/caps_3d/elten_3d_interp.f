      subroutine elten_3d_interp
     +
     +     (elst
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     +     ,x4,y4,z4
     +     ,x5,y5,z5
     +     ,x6,y6,z6
     +
     +     ,vnx1,vny1,vnz1
     +     ,vnx2,vny2,vnz2
     +     ,vnx3,vny3,vnz3
     +     ,vnx4,vny4,vnz4
     +     ,vnx5,vny5,vnz5
     +     ,vnx6,vny6,vnz6
     +
     +     ,xr1,yr1,zr1
     +     ,xr2,yr2,zr2
     +     ,xr3,yr3,zr3
     +     ,xr4,yr4,zr4
     +     ,xr5,yr5,zr5
     +     ,xr6,yr6,zr6
     +
     +     ,vnxr1,vnyr1,vnzr1
     +     ,vnxr2,vnyr2,vnzr2
     +     ,vnxr3,vnyr3,vnzr3
     +     ,vnxr4,vnyr4,vnzr4
     +     ,vnxr5,vnyr5,vnzr5
     +     ,vnxr6,vnyr6,vnzr6
     +
     +     ,al,be,ga
     +     ,xi,eta
     +     ,elten
     +     )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------
c Computat the elastic tension (eltenen)
c on a triangle
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  defg(3,3)
      Dimension elten(3,3)

c--------
c prepare
c--------

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c-----------------------------
c  compute the basis functions
c-----------------------------

      ph2 = xi*(xi-al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi*(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0-ph2-ph3-ph4-ph5-ph6

c-----------------------------
c interpolate the
c current normal vector
c----------------------------

      vnx = vnx1*ph1 + vnx2*ph2 + vnx3*ph3 + vnx4*ph4
     +    + vnx5*ph5 + vnx6*ph6

      vny = vny1*ph1 + vny2*ph2 + vny3*ph3 + vny4*ph4
     +    + vny5*ph5 + vny6*ph6

      vnz = vnz1*ph1 + vnz2*ph2 + vnz3*ph3 + vnz4*ph4
     +    + vnz5*ph5 + vnz6*ph6

c-------------------------
c interpolate the 
c reference normal vector
c------------------------

      vnxr = vnxr1*ph1 + vnxr2*ph2 + vnxr3*ph3 + vnxr4*ph4
     +     + vnxr5*ph5 + vnxr6*ph6

      vnyr = vnyr1*ph1 + vnyr2*ph2 + vnyr3*ph3 + vnyr4*ph4
     +     + vnyr5*ph5 + vnyr6*ph6

      vnzr = vnzr1*ph1 + vnzr2*ph2 + vnzr3*ph3 + vnzr4*ph4
     +     + vnzr5*ph5 + vnzr6*ph6

c------------------------------------------
c compute xi derivatives of basis functions
c------------------------------------------

      dph2 =  (2.0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0)/(ga*bec)
      dph4 =  (1.0-2.0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c---------------------------------------
c interpolate xi derivatives of x and xr
c---------------------------------------

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

      DxrDxi = xr1*dph1 + xr2*dph2 + xr3*dph3 + xr4*dph4
     +       + xr5*dph5 + xr6*dph6
      DyrDxi = yr1*dph1 + yr2*dph2 + yr3*dph3 + yr4*dph4
     +       + yr5*dph5 + yr6*dph6
      DzrDxi = zr1*dph1 + zr2*dph2 + zr3*dph3 + zr4*dph4
     +       + zr5*dph5 + zr6*dph6

c-------------------------------------------
c compute eta derivatives of basis functions
c-------------------------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0*eta-be+xi*(be+ga-1.0)/ga)/bec
      pph4 = -xi/alalc
      pph5 =  xi/gagac
      pph6 =  (1.0-xi-2.0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c----------------------------------------
c interpolate eta derivatives of x and xr
c----------------------------------------

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6
      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6
      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

      DxrDet = xr1*pph1 + xr2*pph2 + xr3*pph3 + xr4*pph4
     +       + xr5*pph5 + xr6*pph6
      DyrDet = yr1*pph1 + yr2*pph2 + yr3*pph3 + yr4*pph4
     +       + yr5*pph5 + yr6*pph6
      DzrDet = zr1*pph1 + zr2*pph2 + zr3*pph3 + zr4*pph4
     +       + zr5*pph5 + zr6*pph6

c--------------------------
c solve three systems
c of three linear equations
c for the deformation gradient
c--------------------------

      a11 = DxrDxi
      a12 = DyrDxi
      a13 = DzrDxi

      a21 = DxrDet
      a22 = DyrDet
      a23 = DzrDet

      a31 = vnxr
      a32 = vnyr
      a33 = vnzr

      b1 = DxDxi
      b2 = DxDet
      b3 = 0.0

      call cramer_33
     +
     +   (A11,A12,A13
     +   ,A21,A22,A23
     +   ,A31,A32,A33
     +   ,B1,B2,B3
     +   ,defg(1,1)
     +   ,defg(1,2)
     +   ,defg(1,3)
     +   )

      b1 = DyDxi
      b2 = DyDet
      b3 = 0.0

      call cramer_33
     +
     +   (A11,A12,A13
     +   ,A21,A22,A23
     +   ,A31,A32,A33
     +   ,B1,B2,B3
     +   ,defg(2,1)
     +   ,defg(2,2)
     +   ,defg(2,3)
     +   )

      b1 = DzDxi
      b2 = DzDet
      b3 = 0.0

      call cramer_33
     +
     +   (A11,A12,A13
     +   ,A21,A22,A23
     +   ,A31,A32,A33
     +   ,B1,B2,B3
     +   ,defg(3,1)
     +   ,defg(3,2)
     +   ,defg(3,3)
     +   )

c------------------------------
c compute the projection matrix
c------------------------------

      prjxx = 1.0D0 - vnx*vnx
      prjxy =       - vnx*vny
      prjxz =       - vnx*vnz
      prjyy = 1.0D0 - vny*vny
      prjyz =       - vny*vnz
      prjzz = 1.0D0 - vnz*vnz
      prjyx = prjxy
      prjzx = prjxz
      prjzy = prjyz

c------------------------------
c compute the surface deformation gradient
c------------------------------

      savexx = defg(1,1)
      savexy = defg(1,2)
      savexz = defg(1,3)

      saveyx = defg(2,1)
      saveyy = defg(2,2)
      saveyz = defg(2,3)

      savezx = defg(3,1)
      savezy = defg(3,2)
      savezz = defg(3,3)

      Axx = prjxx*savexx + prjxy*saveyx +  prjxz*savezx
      Axy = prjxx*savexy + prjxy*saveyy +  prjxz*savezy
      Axz = prjxx*savexz + prjxy*saveyz +  prjxz*savezz

      Ayx = prjyx*savexx + prjyy*saveyx +  prjyz*savezx
      Ayy = prjyx*savexy + prjyy*saveyy +  prjyz*savezy
      Ayz = prjyx*savexz + prjyy*saveyz +  prjyz*savezz

      Azx = prjzx*savexx + prjzy*saveyx +  prjzz*savezx
      Azy = prjzx*savexy + prjzy*saveyy +  prjzz*savezy
      Azz = prjzx*savexz + prjzy*saveyz +  prjzz*savezz

c---------------------------------------------
c Compute the left Cauchy-Green stress tensor:
c
c AAT = A . A^T
c----------------------------------------

      AATxx = Axx*Axx + Axy*Axy + Axz*Axz
      AATxy = Axx*Ayx + Axy*Ayy + Axz*Ayz
      AATxz = Axx*Azx + Axy*Azy + Axz*Azz
 
      AATyy = Ayx*Ayx + Ayy*Ayy + Ayz*Ayz
      AATyz = Ayx*Azx + Ayy*Azy + Ayz*Azz

      AATzz = Azx*Azx + Azy*Azy + Azz*Azz

      AATyx = AATxy
      AATzx = AATxz
      AATzy = AATyz

c     FFTxx = defg(1,1)*defg(1,1)
c    +      + defg(1,2)*defg(2,1)
c    +      + defg(1,3)*defg(3,1)

c     FFTxy = defg(1,1)*defg(1,2)
c    +      + defg(1,2)*defg(2,2)
c    +      + defg(1,3)*defg(3,2)

c     FFTxz = defg(1,1)*defg(1,3)
c    +      + defg(1,2)*defg(2,3)
c    +      + defg(1,3)*defg(3,3)

c     FFTyy = defg(2,1)*defg(1,2)
c    +      + defg(2,2)*defg(2,2)
c    +      + defg(2,3)*defg(3,2)

c     FFTyz = defg(2,1)*defg(1,3)
c    +      + defg(2,2)*defg(2,3)
c    +      + defg(2,3)*defg(3,3)

c     FFTzz = defg(3,1)*defg(1,3)
c    +      + defg(3,2)*defg(2,3)
c    +      + defg(3,3)*defg(3,3)

c     FFTyx = FFTxy
c     FFTzx = FFTxz
c     FFTzy = FFTyz

c----------------------------------------
c Compute: 
c
c normal vector:    vn = (DxDxi)x(DxDeta)
c surface metric:   hs = norm(vn)
c----------------------------------------

      vx = DyDxi * DzDet - DyDet * DzDxi
      vy = DzDxi * DxDet - DzDet * DxDxi
      vz = DxDxi * DyDet - DxDet * DyDxi

      hs  = sqrt(vx**2 + vy**2 + vz**2)

      vxr = DyrDxi * DzrDet - DyrDet * DzrDxi
      vyr = DzrDxi * DxrDet - DzrDet * DxrDxi
      vzr = DxrDxi * DyrDet - DxrDet * DyrDxi

      hsr = sqrt(vxr**2 + vyr**2 + vzr**2)

c---
c dilatation
c---

      dilat = hs/hsr

c     write (6,*) "dilatation: ",dilat

c--------------------------
c compute Lamda1 and Lamda2
c--------------------------

      RL1 = log(dilat)
      RL2 = 0.5D0*(AATxx+AATyy+AATzz)-1.0D0

c----
c compute the elastic tensions
c----

      a1  = 0.0
      a2  = elst/1.5
      a3  = elst/3.0

      wl1 = a1-a3+(a1+a2)*RL1       ! dW/dL1
      wl2 = a3                      ! dW/dL2

      fc = 1.0/dilat
 
      elten(1,1) = fc * (wl1*prjxx + wl2*AATxx)
      elten(1,2) = fc * (wl1*prjxy + wl2*AATxy)
      elten(1,3) = fc * (wl1*prjxz + wl2*AATxz)
      elten(2,2) = fc * (wl1*prjyy + wl2*AATyy)
      elten(2,3) = fc * (wl1*prjyz + wl2*AATyz)
      elten(3,3) = fc * (wl1*prjzz + wl2*AATzz)

      elten(2,1) = elten(1,2)
      elten(3,1) = elten(1,3)
      elten(3,2) = elten(2,3)

c-----
c done
c-----

  100 Format (20(1x,f10.5))

      return
      end
