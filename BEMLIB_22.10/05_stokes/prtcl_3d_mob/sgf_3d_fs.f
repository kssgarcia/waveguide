      subroutine sgf_3d_fs 
     +
     +  (Iopt
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,Gxx,Gxy,Gxz
     +  ,Gyx,Gyy,Gyz
     +  ,Gzx,Gzy,Gzz
     +  ,px,py,pz
     +  ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +  ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +  ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +  )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c---------------------------------------
c Free-space Green's function: Stokeslet
c
c Pozrikidis (1992, p. 23)
c
c Iopt =  1 generates only the Green's function
c      ne 1 generates the Green's function,
c           pressure, and stress
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      dx = x-x0
      dy = y-y0
      dz = z-z0

      dxx = dx*dx
      dxy = dx*dy
      dxz = dx*dz
      dyy = dy*dy
      dyz = dy*dz
      dzz = dz*dz

      r = Dsqrt(dxx+dyy+dzz)
      r3  = r*r*r
      ri  = 1.0D0/r
      ri3 = 1.0D0/r3

      Gxx = ri + dxx*ri3
      Gxy =      dxy*ri3
      Gxz =      dxz*ri3
      Gyy = ri + dyy*ri3
      Gyz =      dyz*ri3
      Gzz = ri + dzz*ri3

      Gyx = Gxy
      Gzx = Gxz
      Gzy = Gyz


c---
c compute the stress tensor and the pressure
c---

      if(Iopt.gt.1) then

      cf = -6.0D0/(r3*r*r)

      Txxx = dxx*dx * cf
      Txxy = dxy*dx * cf
      Txxz = dxz*dx * cf
      Tyxy = dyy*dx * cf
      Tyxz = dyz*dx * cf
      Tzxz = dzz*dx * cf

      Txyx = Txxy
      Txyy = Tyxy
      Txyz = Tyxz
      Tyyy = dyy*dy * cf
      Tyyz = dyz*dy * cf
      Tzyz = dzz*dy * cf

      Txzx = Txxz
      Txzy = Tyxz
      Txzz = Tzxz
      Tyzy = dyy*dz * cf
      Tyzz = dyz*dz * cf
      Tzzz = dzz*dz * cf

c---------
c pressure
c---------

      cf = 2.0D0*ri3
      px = dx * cf
      py = dy * cf
      pz = dz * cf

      end if

c-----
c done
c-----

      return
      end
