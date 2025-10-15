      subroutine sgf_3d_2p 
     +
     +   (method
     +   ,eps
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,area
     +   ,max1,max2
     +   ,Gxx,Gxy,Gxz
     +   ,Gyx,Gyy,Gyz
     +   ,Gzx,Gzy,Gzz
     +   )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c--------------------------------------
c Three-dimensional doubly periodic
c Green's function of Stokes flow
c
c The point forces are deployed
c parallel to the xy plane
c (perpendicular to the z axis)
c
c One of the point forces is located at: x0, y0, z0
c
c x, y, z: field point coordinates
c
c LEGEND:
c -------
c
c method = 1 use the straight Fourier series method
c method = 2 use the fast summation method of Pozrikidis (1996)
c
c eps: parameter for numerical differentiation
c
c This subroutine computes only the velocity tensor
c not the stress tensor
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double precision k1,k2,k,ks,Lapl

c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      pi2  = 2.0D0*pi
      pi4  = 4.0D0*pi
      srpi = sqrt(pi)

c--------
c prepare
c--------

      dx0 = x-x0
      dy0 = y-y0
      dz0 = z-z0

c---
c initialize
c---

      Gxx = 0.0D0
      Gxy = 0.0D0
      Gxz = 0.0D0
      Gyy = 0.0D0
      Gyz = 0.0D0
      Gzz = 0.0D0

c-------------------------------------------------
c  Sum on real space over the 2D lattice based on
c  equations (3.3) and (4.2) of Pozrikidis (1996)
c
c  This loop is executed only for method = 2
c  Will skip for straight Fourier series (method = 1)
c----------------------------------------------------

c--------
      if(method.eq.2) then
c--------

      fc = 2.0D0/srpi

      Do i1 = -max1,max1
      Do i2 = -max1,max1

        dx = dx0 - i1*a11 - i2*a21 
        dy = dy0 - i1*a12 - i2*a22
        dz = dz0 

        r  = Dsqrt(dx*dx + dy*dy + dz*dz)
        r3 = r*r*r
        w  = ew*r
        w2 = w*w

        expp  = exp(-w)
        ccri  = expp*(1.0D0 -3.0D0*w+w2)/r
        ddri3 = expp*(1.0D0 +      w-w2)/r3

        Gxx = Gxx + ccri + dx*dx * ddri3
        Gxy = Gxy        + dx*dy * ddri3
        Gxz = Gxz        + dx*dz * ddri3
        Gyy = Gyy + ccri + dy*dy * ddri3
        Gyz = Gyz        + dy*dz * ddri3
        Gzz = Gzz + ccri + dz*dz * ddri3

c      if(i1.eq.0) write (6,101) i1,i2,ccri

       End Do
      End Do

c--------
      End If
c--------

c-------------------------------------
c  sum in wavenumber space
c
c  use equations (3.15), (4.5), and (4.5)
c  of Pozrikidis (1996)
c-------------------------------------

      pxx = 0.0D0
      pxy = 0.0D0
      pxz = 0.0D0
      pyy = 0.0D0
      pyz = 0.0D0
      pzz = 0.0D0

      Do i1=-max2,max2
      Do i2=-max2,max2

        k1 = i1*b11 + i2*b21 
        k2 = i1*b12 + i2*b22
 
        ks = k1*k1 + k2*k2
        k  = sqrt(ks)

c--------------------------------------------
c For the pure Fourier method
c use expression (2.13) of Pozrikidis (1996)
c
c Definitions:
c
c      S_3D_2P       :=   A * cos(k1*x+k2*y)
c    d(S_2D_2P)/dz   := - B * sin(k1*x+k2*y)
c  d^2(S_2D_2P)/dz^2 :=   D * cos(k1*x+k2*y)
c--------------------------------------------

c---------
      if(k.gt.0.000001) then  ! skip the zero wave number
c---------

        rho = Dabs(dz0)*k
        fcc = exp(-rho)

c-------------------------
      if(method.eq.1) then      ! straight Fourier series (2.13)
c-------------------------

        A =  (1.0D0+rho)*fcc/k**3
        B = -       rho *fcc/k**2
        D = -(1.0D0-rho)*fcc/k

c------------------------------
      else if(method.eq.2) then     ! accelerated
c------------------------------

        delta = abs(ew*dz0)
        zt    = k/ew

        call sgf_3d_2p_aux 
     +
     +   (delta,zt,ew
     +   ,eps
     +   ,A,B,D
     +   )

c-----------
      end if
c-----------

      if(dz0.lt.0) B = -B

      arg = k1*dx0 + k2*dy0

      fc = Dcos(arg)
      fs = Dsin(arg)

      Lapl = D-ks*A     ! Laplacian of cB

      XX = -(Lapl +k1*k1*A)*fc
      XY =        -k1*k2*A *fc
      XZ =        -k1   *B *fs
      YY = -(Lapl +k2*k2*A)*fc
      YZ =        -k2   *B *fs
      ZZ =         ks   *A *fc

      pxx = pxx + XX
      pxy = pxy + XY
      pxz = pxz + XZ
      pyy = pyy + YY
      pyz = pyz + YZ
      pzz = pzz + ZZ

c     if(i1.eq.0) write (6,101) i1,i2,Lapl,A,B

c------------
      end if       ! skip the zero wavenumber
c------------

c---
      end Do
      end Do
c---

c---------------------------------------
c contribution from the zero wave number
c---------------------------------------

      if(method.eq.2) then

       supp = delta*(delta-2.0D0)*exp(-delta)/ew
       pxx  = pxx - supp
       pyy  = pyy - supp

      end if

c--------------------------------
c add the contributions from real
c and wavenumber space
c--------------------------------

      fc = pi2/area

      Gxx = Gxx + pxx*fc
      Gxy = Gxy + pxy*fc
      Gxz = Gxz + pxz*fc
      Gyy = Gyy + pyy*fc
      Gyz = Gyz + pyz*fc
      Gzz = Gzz + pzz*fc

      Gyx = Gxy
      Gzx = Gxz
      Gzy = Gyz

c------------------------------------
c subtract a simple shear flow
c from the first two diagonal entries
c------------------------------------

      fc = pi4/area

      Gxx = Gxx - fc*abs(dz0)
      Gyy = Gyy - fc*abs(dz0)

c------------------------------
c uncomment the following lines
c to subtract off the Stokeslet
c if desired
c------------------------------
c
c     rr   = sqrt(dx0**2+dy0**2+dz0**2)
c     rri  = 1.0/rr
c     rr3i = 1.0/rr**3
c
c     Gxx = Gxx-rri-dx0*dx0*rr3i
c     Gxy = Gxy    -dx0*dy0*rr3i
c     Gxz = Gxz    -dx0*dz0*rr3i
c     Gyy = Gyy-rri-dy0*dy0*rr3i
c     Gyz = Gyz    -dy0*dz0*rr3i
c     Gzz = Gzz-rri-dz0*dz0*rr3i
c     Gyx = Gxy
c     Gzx = Gxz
c     Gzy = Gyz
c
c------------------------------

c-----
c done
c-----

  100 Format (9(1x,f10.6))
  101 Format (2(1x,i3),9(1x,f15.10))

      return
      end
