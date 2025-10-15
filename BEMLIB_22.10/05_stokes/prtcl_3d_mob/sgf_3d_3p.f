      subroutine sgf_3d_3p 
     +
     +   (Iopt
     +   ,x,y,z
     +   ,x0,y0,z0
     +
     +   ,a11,a12,a13
     +   ,a21,a22,a23
     +   ,a31,a32,a33
     +
     +   ,b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +
     +   ,ew,tau
     +   ,max1,max2
     +
     +   ,Gxx,Gxy,Gxz
     +   ,Gyx,Gyy,Gyz
     +   ,Gzx,Gzy,Gzz
     +
     +   ,Px,Py,Pz
     +
     +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +   )

c==========================================
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------------
c triply periodic Green's function of Stokes flow
c
c FLAGS:
c -----
c
c Iopt =  1 computes G
c Iopt ne 1 computes G, p, T,
c
c SYMBOLS:
c -------
c
c x,y,z:     field point coordinates
c x0,y0,z0:  coordinates of one of the point forces
c
c qxxx:  matrix for summing in wave number space
c vxxx:  matrix for summing in wave number space
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k1,k2,k3,ks

      Dimension qxx(-9:9,-9:9,-9:9)
      Dimension qxy(-9:9,-9:9,-9:9)
      Dimension qxz(-9:9,-9:9,-9:9)
      Dimension qyy(-9:9,-9:9,-9:9)
      Dimension qyz(-9:9,-9:9,-9:9)
      Dimension qzz(-9:9,-9:9,-9:9)

      Dimension vxxx(-9:9,-9:9,-9:9)
      Dimension vxxy(-9:9,-9:9,-9:9)
      Dimension vxxz(-9:9,-9:9,-9:9)
      Dimension vyxy(-9:9,-9:9,-9:9)
      Dimension vyxz(-9:9,-9:9,-9:9)
      Dimension vzxz(-9:9,-9:9,-9:9)

      Dimension vxyx(-9:9,-9:9,-9:9)
      Dimension vxyy(-9:9,-9:9,-9:9)
      Dimension vxyz(-9:9,-9:9,-9:9)
      Dimension vyyy(-9:9,-9:9,-9:9)
      Dimension vyyz(-9:9,-9:9,-9:9)
      Dimension vzyz(-9:9,-9:9,-9:9)

      Dimension vxzx(-9:9,-9:9,-9:9)
      Dimension vxzy(-9:9,-9:9,-9:9)
      Dimension vxzz(-9:9,-9:9,-9:9)
      Dimension vyzy(-9:9,-9:9,-9:9)
      Dimension vyzz(-9:9,-9:9,-9:9)
      Dimension vzzz(-9:9,-9:9,-9:9)

      Dimension ppx(-9:9,-9:9,-9:9)
      Dimension ppy(-9:9,-9:9,-9:9)
      Dimension ppz(-9:9,-9:9,-9:9)

c--------------
c common blocks
c--------------

      common/qqqq_3d/qxx,qxy,qxz,qyy,qyz,qzz

      common/vvvv_3d/vxxx,vxxy,vxxz,vyxy,vyxz,vzxz
     +              ,vxyx,vxyy,vxyz,vyyy,vyyz,vzyz
     +              ,vxzx,vxzy,vxzz,vyzy,vyzz,vzzz
     +              ,ppx,ppy,ppz

c----------
c constants
c----------
 
      pi   = 3.14159 265358 D0
      pi2  = 2.0D0*pi
      pi8  = 8.0D0*pi
      srpi = Dsqrt(pi)

c--------
c prepare
c--------

      dx0 = x-x0
      dy0 = y-y0
      dz0 = z-z0

c-----------
c initialize
c-----------

      Gxx = 0.0D0
      Gxy = 0.0D0
      Gxz = 0.0D0
      Gyy = 0.0D0
      Gyz = 0.0D0
      Gzz = 0.0D0

      if(Iopt.ne.1) then

        Px = 0.0D0
        Py = 0.0D0
        Pz = 0.0D0

        Txxx = 0.0D0
        Txxy = 0.0D0
        Txxz = 0.0D0
        Tyxy = 0.0D0
        Tyxz = 0.0D0
        Tzxz = 0.0D0

        Txyx = 0.0D0
        Txyy = 0.0D0
        Txyz = 0.0D0
        Tyyy = 0.0D0
        Tyyz = 0.0D0
        Tzyz = 0.0D0

        Txzx = 0.0D0
        Txzy = 0.0D0
        Txzz = 0.0D0
        Tyzy = 0.0D0
        Tyzz = 0.0D0
        Tzzz = 0.0D0

      end if

c-------------------
c sum in real space
c-------------------

      fc  = 2.0D0/srpi
      fc2 = 2.0D0*fc
      fc4 = 4.0D0*fc

      Do i1 = -max1,max1
      Do i2 = -max1,max1
      Do i3 = -max1,max1

        dx = dx0 - i1*a11 - i2*a21 - i3*a31
        dy = dy0 - i1*a12 - i2*a22 - i3*a32
        dz = dz0 - i1*a13 - i2*a23 - i3*a33

        r  = Dsqrt(dx**2+dy**2+dz**2)

        r3 = r*r*r
        W  = ew*r
        W2 = W*W

c---
c compute erfc
c---
        T = 1.0D0/(1.0D0+0.5D0*W)

        ERFCC = T*EXP(-W**2-1.26551223
     +         +T*(1.00002368D0+T*( 0.37409196D0
     +         +T*(0.09678418D0+T*(-0.18628806D0
     +         +T*(0.27886807D0+T*(-1.13520398D0
     +         +T*(1.48851587D0+T*(-0.82215223D0
     +         +T*0.17087277D0))))))))
     +        )

        EXPP = Dexp(-W2)
        AUX  = fc*W*EXPP

        ccri  = (ERFCC + (2.0D0*W2 -3.0D0 )*AUX)/r
        ddri3 = (ERFCC + (1.0D0 - 2.0D0*W2 )*AUX)/r3

        Gxx = Gxx + ccri + dx**2 * ddri3
        Gxy = Gxy +        dx*dy * ddri3
        Gxz = Gxz +        dx*dz * ddri3
        Gyy = Gyy + ccri + dy**2 * ddri3
        Gyz = Gyz +        dy*dz * ddri3
        Gzz = Gzz + ccri + dz**2 * ddri3

c-----------
        if(Iopt.ne.1) then
c-----------

        W3 = W2*W
        W4 = W2**2
        r2 = r*r
        r5 = r3*r2

        Ipress = 3;
        Ipress = 2;

c---
        if(Ipress.eq.1) then
         prcf = ERFCC+ fc*W*(1.0D0 + 4.0D0*W2 - 4.0D0/3.0D0*W4)* EXPP
        else if(Ipress.eq.2) then
         prcf = ERFCC + fc*W*(1.0D0 -6.0D0 * W2 + 2.0D0 * W4)* EXPP
        else if(Ipress.eq.3) then
         prcf = ERFCC + fc*W*EXPP
        end if
c---

        Px = Px + 2.0D0*prcf/r3*dx
        Py = Py + 2.0D0*prcf/r3*dy
        Pz = Pz + 2.0D0*prcf/r3*dz

        eeri5 = 2.0D0*(3.0D0*ERFCC + fc*W* 
     +          (3.0D0 + 2.0D0*W2 - 4.0D0*W4) * EXPP )/r5

        ffri3 = fc4 * W3 * EXPP/r3

        Txxx = Txxx - eeri5 * dx*dx*dx - ffri3*(dx+dx+dx)
        Txxy = Txxy - eeri5 * dx*dx*dy - ffri3* dy
        Txxz = Txxz - eeri5 * dx*dx*dz - ffri3* dz
        Tyxy = Tyxy - eeri5 * dy*dx*dy - ffri3* dx
        Tyxz = Tyxz - eeri5 * dy*dx*dz 
        Tzxz = Tzxz - eeri5 * dz*dx*dz - ffri3* dx

        Txyx = Txyx - eeri5 * dx*dy*dx - ffri3* dy
        Txyy = Txyy - eeri5 * dx*dy*dy - ffri3* dx
        Txyz = Txyz - eeri5 * dx*dy*dz
        Tyyy = Tyyy - eeri5 * dy*dy*dy - ffri3*(dy+dy+dy)
        Tyyz = Tyyz - eeri5 * dy*dy*dz - ffri3* dz
        Tzyz = Tzyz - eeri5 * dz*dy*dz - ffri3* dy

        Txzx = Txzx - eeri5 * dx*dz*dx - ffri3* dz
        Txzy = Txzy - eeri5 * dx*dz*dy
        Txzz = Txzz - eeri5 * dx*dz*dz - ffri3* dx
        Tyzy = Tyzy - eeri5 * dy*dz*dy - ffri3* dz
        Tyzz = Tyzz - eeri5 * dy*dz*dz - ffri3* dy
        Tzzz = Tzzz - eeri5 * dz*dz*dz - ffri3*(dz+dz+dz)

c-----------
      end if   ! Iopt
c-----------

c------------
      end Do
      end Do
      end Do
c------------

c-------------------------
c  sum in wavenumber space
c-------------------------

      exx = 0.0D0
      exy = 0.0D0
      exz = 0.0D0
      eyy = 0.0D0
      eyz = 0.0D0
      ezz = 0.0D0

      if(Iopt.ne.1) then

        wwx = 0.0D0   ! for the pressure
        wwy = 0.0D0
        wwz = 0.0D0

        wxxx = 0.0D0
        wxxy = 0.0D0
        wxxz = 0.0D0
        wyxy = 0.0D0
        wyxz = 0.0D0
        wzxz = 0.0D0

        wxyx = 0.0D0
        wxyy = 0.0D0
        wxyz = 0.0D0
        wyyy = 0.0D0
        wyyz = 0.0D0
        wzyz = 0.0D0

        wxzx = 0.0D0
        wxzy = 0.0D0
        wxzz = 0.0D0
        wyzy = 0.0D0
        wyzz = 0.0D0
        wzzz = 0.0D0

      end if

      Do i1=-max2,max2
      Do i2=-max2,max2
      Do i3=-max2,max2

      k1 = i1*b11 + i2*b21 + i3*b31
      k2 = i1*b12 + i2*b22 + i3*b32
      k3 = i1*b13 + i2*b23 + i3*b33

      ks = k1**2 + k2**2 + k3**2

c---
      if(ks.gt.0.0000001) then     ! skip the zero wavenumber
c---

      arg = k1*dx0 + k2*dy0 + k3*dz0

      f = Dcos(arg)

      exx = exx + qxx(i1,i2,i3)*f
      exy = exy + qxy(i1,i2,i3)*f
      exz = exz + qxz(i1,i2,i3)*f
      eyy = eyy + qyy(i1,i2,i3)*f
      eyz = eyz + qyz(i1,i2,i3)*f
      ezz = ezz + qzz(i1,i2,i3)*f

c----------------------
      if(Iopt.ne.1) then
c-----------------------

      f = Dsin(arg)

      wwx = wwx + ppx(i1,i2,i3)*f  ! pressure
      wwy = wwy + ppy(i1,i2,i3)*f  ! pressure
      wwz = wwz + ppz(i1,i2,i3)*f  ! pressure

      wxxx = wxxx + vxxx(i1,i2,i3)*f
      wxxy = wxxy + vxxy(i1,i2,i3)*f
      wxxz = wxxz + vxxz(i1,i2,i3)*f
      wyxy = wyxy + vyxy(i1,i2,i3)*f
      wyxz = wyxz + vyxz(i1,i2,i3)*f
      wzxz = wzxz + vzxz(i1,i2,i3)*f

      wxyx = wxyx + vxyx(i1,i2,i3)*f
      wxyy = wxyy + vxyy(i1,i2,i3)*f
      wxyz = wxyz + vxyz(i1,i2,i3)*f
      wyyy = wyyy + vyyy(i1,i2,i3)*f
      wyyz = wyyz + vyyz(i1,i2,i3)*f
      wzyz = wzyz + vzyz(i1,i2,i3)*f

      wxzx = wxzx + vxzx(i1,i2,i3)*f
      wxzy = wxzy + vxzy(i1,i2,i3)*f
      wxzz = wxzz + vxzz(i1,i2,i3)*f
      wyzy = wyzy + vyzy(i1,i2,i3)*f
      wyzz = wyzz + vyzz(i1,i2,i3)*f
      wzzz = wzzz + vzzz(i1,i2,i3)*f

c-----------
      end if   ! Iopt
c-----------

c----
      end if ! zero wave number
c----

c------------
      end Do
      end Do
      end Do
c------------

      fc = pi8/tau

      Gxx = Gxx + exx*fc
      Gxy = Gxy + exy*fc
      Gxz = Gxz + exz*fc
      Gyx = Gxy
      Gyy = Gyy + eyy*fc
      Gyz = Gyz + eyz*fc
      Gzx = Gxz
      Gzy = Gyz
      Gzz = Gzz + ezz*fc

c--------------
c stress tensor
c-------------

      if(Iopt.ne.1) then

      Txxx = Txxx + wxxx*fc
      Txxy = Txxy + wxxy*fc
      Txxz = Txxz + wxxz*fc
      Tyxy = Tyxy + wyxy*fc
      Tyxz = Tyxz + wyxz*fc
      Tzxz = Tzxz + wzxz*fc

      Txyx = Txyx + wxyx*fc
      Txyy = Txyy + wxyy*fc
      Txyz = Txyz + wxyz*fc
      Tyyy = Tyyy + wyyy*fc
      Tyyz = Tyyz + wyyz*fc
      Tzyz = Tzyz + wzyz*fc

      Txzx = Txzx + wxzx*fc
      Txzy = Txzy + wxzy*fc
      Txzz = Txzz + wxzz*fc
      Tyzy = Tyzy + wyzy*fc
      Tyzz = Tyzz + wyzz*fc
      Tzzz = Tzzz + wzzz*fc

c---
c add the mean pressure gradient
c---

      dx0 = x-x0   ! important: do not remove
      dy0 = y-y0   ! important: do not remove
      dz0 = z-z0   ! important: do not remove

      pgx = dx0*fc
      pgy = dy0*fc
      pgz = dz0*fc

      Txxx = Txxx-pgx  
      Tyxy = Tyxy-pgx
      Tzxz = Tzxz-pgx

      Txyx = Txyx-pgy
      Tyyy = Tyyy-pgy
      Tzyz = Tzyz-pgy

      Txzx = Txzx-pgz
      Tyzy = Tyzy-pgz
      Tzzz = Tzzz-pgz

      Px = Px + wwx*fc + pgx
      Py = Py + wwy*fc + pgy
      Pz = Pz + wwz*fc + pgz

c     write (6,*) Px,Py,Pz

c     Px = -(Txxx+Tyxy+Tzxz)/3.0D0   ! pressure in terms of the trace
c     Py = -(Txyx+Tyyy+Tzyz)/3.0D0
c     Pz = -(Txzx+Tyzy+Tzzz)/3.0D0

c---
      end if ! of Iopt
c---

c-----
c done
c-----

      return
      end
