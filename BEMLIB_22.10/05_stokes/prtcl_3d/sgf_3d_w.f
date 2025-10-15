      subroutine sgf_3d_w 
     +
     +   (Iopt
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,wall
     +   ,Gxx,Gxy,Gxz
     +   ,Gyx,Gyy,Gyz
     +   ,Gzx,Gzy,Gzz
     +   ,px,py,pz
     +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +   )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c=========================================
c Green's function for semi-infinite flow
c bounded by a plane wall located at x = wall
c
c Iopt =  1 Compute only the Green's function
c      ne 1 Compute the Green's function,
c           pressure, and stress
c=========================================

      Implicit Double Precision (a-h,o-z)

c-------------
c primary part
c-------------

      dx = x-x0
      dy = y-y0
      dz = z-z0

      dxx = dx*dx
      dxy = dx*dy
      dxz = dx*dz
      dyy = dy*dy
      dyz = dy*dz
      dzz = dz*dz

      r   = Dsqrt(dxx+dyy+dzz)
      r3  = r*r*r
      r5  = r3*r*r
      ri  = 1.0D0/r
      ri3 = 1.0D0/r3
      ri5 = 1.0D0/r5

      Sxx = ri + dxx*ri3
      Sxy =      dxy*ri3
      Sxz =      dxz*ri3
      Syy = ri + dyy*ri3
      Syz =      dyz*ri3
      Szz = ri + dzz*ri3

c------------------------------
      if(Iopt.ne.1) then    ! compute the stress and pressure

      cf = -6.0D0*ri5

      TSxxx = dxx*dx * cf
      TSxxy = dxy*dx * cf
      TSxxz = dxz*dx * cf
      TSyxy = dyy*dx * cf
      TSyxz = dyz*dx * cf
      TSzxz = dzz*dx * cf

      TSyyy = dyy*dy * cf
      TSyyz = dyz*dy * cf
      TSzyz = dzz*dy * cf

      TSyzy = dyy*dz * cf
      TSyzz = dyz*dz * cf
      TSzzz = dzz*dz * cf

      cf = 2.0D0*ri3

      psx = dx * cf
      psy = dy * cf
      psz = dz * cf

      End If
c------------------------------

c-------------
c image system
c-------------

      x0im = 2.0D0*wall-x0
      dx   = x-x0im

      dxx = dx*dx
      dxy = dx*dy
      dxz = dx*dz
      
      r   = Dsqrt(dxx+dyy+dzz)
      r3  = r*r*r
      r5  = r3*r*r
      ri  = 1.0D0/r
      ri3 = 1.0D0/r3
      ri5 = 1.0D0/r5

c----------------
c image stokeslet 
c----------------

      Sxx = Sxx - ri - dxx*ri3
      Sxy = Sxy      - dxy*ri3
      Sxz = Sxz      - dxz*ri3
      Syy = Syy - ri - dyy*ri3
      Syz = Syz      - dyz*ri3
      Szz = Szz - ri - dzz*ri3

      Syx = Sxy
      Szx = Sxz
      Szy = Syz

c-----------
      if(Iopt.ne.1) then   ! compute the stress and pressure

      cf = -6.0D0*ri5

      TSxxx = TSxxx - dxx*dx * cf
      TSxxy = TSxxy - dxy*dx * cf
      TSxxz = TSxxz - dxz*dx * cf
      TSyxy = TSyxy - dyy*dx * cf
      TSyxz = TSyxz - dyz*dx * cf
      TSzxz = TSzxz - dzz*dx * cf

      TSxyx = TSxxy
      TSxyy = TSyxy
      TSxyz = TSyxz
      TSyyy = TSyyy - dyy*dy * cf
      TSyyz = TSyyz - dyz*dy * cf
      TSzyz = TSzyz - dzz*dy * cf

      TSxzx = TSxxz
      TSxzy = TSyxz
      TSxzz = TSzxz
      TSyzy = TSyzy - dyy*dz * cf
      TSyzz = TSyzz - dyz*dz * cf
      TSzzz = TSzzz - dzz*dz * cf

      cf = 2.0D0*ri3

      psx = psx - dx * cf
      psy = psy - dy * cf
      psz = psz - dz * cf

      end if
c-----------

c-----------------------
c image potential dipole
c-----------------------

      PDxx = - ri3 + 3.0D0*dxx*ri5
      PDyx =         3.0D0*dxy*ri5
      PDzx =         3.0D0*dxz*ri5

      PDxy = - PDyx
      PDyy =   ri3 - 3.0D0*dyy*ri5
      PDzy =       - 3.0D0*dyz*ri5

      PDxz = - PDzx
      PDyz =   PDzy
      PDzz =   ri3 - 3.0D0*dzz*ri5

c-----------------------
      if(Iopt.ne.1) then

      r7  = r5*r*r
      ri7 = 1.0D0/r7

      cf  =  6.0D0*ri5
      cf1 = 30.0D0*ri7

      TPDxxx =  (dx+dx+dx) * cf - dxx*dx * cf1
      TPDxxy =   dy        * cf - dxy*dx * cf1
      TPDxxz =   dz        * cf - dxz*dx * cf1
      TPDyxy =   dx        * cf - dyy*dx * cf1
      TPDyxz =                  - dyz*dx * cf1
      TPDzxz =   dx        * cf - dzz*dx * cf1

      TPDxyx = - dy        * cf + dxx*dy * cf1
      TPDxyy = - dx        * cf + dxy*dy * cf1
      TPDxyz =                    dxz*dy * cf1
      TPDyyy = -(dy+dy+dy) * cf + dyy*dy * cf1
      TPDyyz = - dz        * cf + dyz*dy * cf1
      TPDzyz = - dy        * cf + dzz*dy * cf1

      TPDxzx = - dz        * cf + dxx*dz * cf1
      TPDxzy =                    dxy*dz * cf1
      TPDxzz = - dx        * cf + dxz*dz * cf1
      TPDyzy = - dz        * cf + dyy*dz * cf1
      TPDyzz = - dy        * cf + dyz*dz * cf1
      TPDzzz = -(dz+dz+dz) * cf + dzz*dz * cf1

c     no pressure contribution

      end if
c-----------------------

c-----------------------
c image stokeslet dipole
c-----------------------

      SDxx = dx * PDxx
      SDyx = dx * PDyx - dy*ri3
      SDzx = dx * PDzx - dz*ri3

      SDxy = dx * PDxy - dy*ri3
      SDyy = dx * PDyy
      SDzy = dx * PDzy

      SDxz = dx * PDxz - dz*ri3
      SDyz = dx * PDyz
      SDzz = dx * PDzz

c-----------------------
      if(Iopt.ne.1) then

      cf = 6.0D0*ri5

      TSDxxx = dx *TPDxxx 
      TSDxxy = dx *TPDxxy - (   -dxy) * cf
      TSDxxz = dx *TPDxxz - (   -dxz) * cf
      TSDyxy = dx *TPDyxy - (dxx-dyy) * cf
      TSDyxz = dx *TPDyxz - (   -dyz) * cf
      TSDzxz = dx *TPDzxz - (dxx-dzz) * cf

      TSDxyx = dx *TPDxyx + dxy * cf
      TSDxyy = dx *TPDxyy 
      TSDxyz = dx *TPDxyz 
      TSDyyy = dx *TPDyyy + dxy * cf
      TSDyyz = dx *TPDyyz 
      TSDzyz = dx *TPDzyz + dxy * cf

      TSDxzx = dx *TPDxzx + dxz * cf
      TSDxzy = dx *TPDxzy 
      TSDxzz = dx *TPDxzz 
      TSDyzy = dx *TPDyzy + dxz * cf
      TSDyzz = dx *TPDyzz
      TSDzzz = dx *TPDzzz + dxz * cf

      PSDx =  - 2.0D0*ri3 + 6.0D0*dxx*ri5
      PSDy =              - 6.0D0*dxy*ri5
      PSDz =              - 6.0D0*dxz*ri5

      end if
c-----------------------

c---------
c assemble
c---------

      h0   = x0-wall
      h02  = 2.0D0*h0
      h0s2 = 2.0D0*h0*h0

      Gxx = Sxx + h0s2 * PDxx - h02 * SDxx
      Gxy = Sxy + h0s2 * PDxy - h02 * SDxy
      Gxz = Sxz + h0s2 * PDxz - h02 * SDxz

      Gyx = Syx + h0s2 * PDyx - h02 * SDyx
      Gyy = Syy + h0s2 * PDyy - h02 * SDyy
      Gyz = Syz + h0s2 * PDyz - h02 * SDyz

      Gzx = Szx + h0s2 * PDzx - h02 * SDzx
      Gzy = Szy + h0s2 * PDzy - h02 * SDzy
      Gzz = Szz + h0s2 * PDzz - h02 * SDzz

c-----------------------
      if(Iopt.ne.1) then

      Txxx = TSxxx + h0s2 * TPDxxx - h02 * TSDxxx
      Txxy = TSxxy + h0s2 * TPDxxy - h02 * TSDxxy
      Txxz = TSxxz + h0s2 * TPDxxz - h02 * TSDxxz
      Tyxy = TSyxy + h0s2 * TPDyxy - h02 * TSDyxy
      Tyxz = TSyxz + h0s2 * TPDyxz - h02 * TSDyxz
      Tzxz = TSzxz + h0s2 * TPDzxz - h02 * TSDzxz
      
      Txyx = TSxyx + h0s2 * TPDxyx - h02 * TSDxyx
      Txyy = TSxyy + h0s2 * TPDxyy - h02 * TSDxyy
      Txyz = TSxyz + h0s2 * TPDxyz - h02 * TSDxyz
      Tyyy = TSyyy + h0s2 * TPDyyy - h02 * TSDyyy
      Tyyz = TSyyz + h0s2 * TPDyyz - h02 * TSDyyz
      Tzyz = TSzyz + h0s2 * TPDzyz - h02 * TSDzyz

      Txzx = TSxzx + h0s2 * TPDxzx - h02 * TSDxzx
      Txzy = TSxzy + h0s2 * TPDxzy - h02 * TSDxzy
      Txzz = TSxzz + h0s2 * TPDxzz - h02 * TSDxzz
      Tyzy = TSyzy + h0s2 * TPDyzy - h02 * TSDyzy
      Tyzz = TSyzz + h0s2 * TPDyzz - h02 * TSDyzz
      Tzzz = TSzzz + h0s2 * TPDzzz - h02 * TSDzzz

      px = psx - h02 * PSDx
      py = psy - h02 * PSDy
      pz = psz - h02 * PSDz

      end if
c----------------------

c-----
c done
c-----

      return
      end
