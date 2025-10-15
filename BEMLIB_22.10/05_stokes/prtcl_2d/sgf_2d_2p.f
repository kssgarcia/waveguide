      subroutine sgf_2d_2p 
     +
     +   (Iopt
     +   ,x,y
     +   ,x0,y0
     +   ,a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,area
     +   ,max1,max2
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,Px,Py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )
 
c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------------------
c Doubly-periodic Green's function  of two-dimensional
c Stokes flow
c
c Iopt = 1 computes G
c      = 2 computes G, p, and T
c
c SYMBOLS
c -------
c
c Gij:     Green's function for the velocity
c Tijk:    Green's function for the stress
c
c a_ij:    lattice base vectors
c b_ij:    wave number base vectors
c
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k1,k2,ks

      Dimension qxx (-15:15,-15:15)
      Dimension qxy (-15:15,-15:15)
      Dimension qyy (-15:15,-15:15)
      Dimension vxxx(-15:15,-15:15)
      Dimension vxxy(-15:15,-15:15)
      Dimension vyxy(-15:15,-15:15)
      Dimension vxyx(-15:15,-15:15)
      Dimension vxyy(-15:15,-15:15)
      Dimension vyyy(-15:15,-15:15)
      Dimension ppx (-15:15,-15:15)
      Dimension ppy (-15:15,-15:15)

c--------------
c common blocks
c--------------

      common/qqqq_2d/qxx,qxy,qyy
      common/vvvv_2d/vxxx,vxxy,vyxy,vxyx,vxyy,vyyy,ppx,ppy

c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      piq  = 0.25D0*pi
      pih  = 0.50D0*pi
      pi2  = 2.00D0*pi
      pi4  = 4.00D0*pi
      pi6  = 6.00D0*pi
      pi8  = 8.00D0*pi
      srpi = sqrt(pi)

c---
c prepare to run
c---

      dxx = x-x0
      dyy = y-y0

      ews = ew**2
 
c-----------------------------------
c move as close to the singular
c point as possible
c by shifting along the base vectors
c-----------------------------------

      Imove = nint(dyy)
      dxx = dxx - Imove*a21   ! move along the second base vector
      dyy = dyy - Imove*a22

      Imove = nint(dxx)
      dxx = dxx - Imove*a11   ! move along the first base vector
      dyy = dyy - Imove*a12

c-----------
c initialize
c-----------

      Gxx = 0.0D0
      Gxy = 0.0D0
      Gyy = 0.0D0

      if(Iopt.ne.1) then

        Px = 0.0D0
        Py = 0.0D0

        Txxx = 0.0D0
        Txxy = 0.0D0
        Tyxx = 0.0D0
        Tyxy = 0.0D0

        Txyx = 0.0D0
        Txyy = 0.0D0
        Tyyx = 0.0D0
        Tyyy = 0.0D0

      end if

c------------------
c sum in real space  
c------------------

      Do i1=-max1,max1
      Do i2=-max1,max1

        dx = dxx-i1*a11-i2*a21
        dy = dyy-i1*a12-i2*a22

        r2 = dx**2+dy**2
        r  = Dsqrt(r2)
        w  = ew*r
        w2 = w**2

        Expp  = exp(-w2)
        cc    = 0.5D0*expint(w2)-Expp
        ddri2 = Expp/r2

        Gxx = Gxx + cc + dx**2 * ddri2
        Gxy = Gxy +      dx*dy * ddri2
        Gyy = Gyy + cc + dy**2 * ddri2

c-------
        if(Iopt.ne.1) then
c-------

        ipress=1;
        ipress=2;

        if(ipress.eq.1) then
         fcpress = 2.0D0*Expp/r2
        elseif(ipress.eq.2) then
         fcpress = 2.0D0*(1.0D0-w2)*Expp/r2
        end if

        px = px + dx*fcpress
        py = py + dy*fcpress

        dx2 = dx**2
        dy2 = dy**2

        fc1 = 2.0D0 * ews * Expp
        fc2 = 4.0D0 * ews/r2 * (1.0D0+1.0D0/w2)*Expp

        Txxx = Txxx + (dx+dx)*fc1 - dx2*dx*fc2
        Txxy = Txxy +  dy    *fc1 - dx2*dy*fc2
        Tyxy = Tyxy               - dy2*dx*fc2

        Txyx = Txyx               - dx2*dy*fc2
        Txyy = Txyy +  dx    *fc1 - dy2*dx*fc2
        Tyyy = Tyyy + (dy+dy)*fc1 - dy2*dy*fc2

c-------
        end if
c-------

      end do
      end do

c---------------------------
c subtract off the Stokeslet
c if desired
c---------------------------
c
c     drs = dxx**2+dyy**2
c     rll = 0.5*log(drs)
c
c     Gxx = Gxx - (-rll+dxx**2 /drs)
c     Gxy = Gxy - (     dxx*dyy/drs)
c     Gyy = Gyy - (-rll+dyy**2 /drs)
c
c---------------------------

c------------------------
c sum in wavenumber space
c------------------------

      exx = 0.0D0
      exy = 0.0D0
      eyy = 0.0D0

      if(Iopt.ne.1) then

        wwx = 0.0D0
        wwy = 0.0D0

        wxxx = 0.0D0
        wxxy = 0.0D0
        wyxy = 0.0D0

        wxyx = 0.0D0
        wxyy = 0.0D0
        wyyy = 0.0D0

      End If

      Do i1=-max2,max2
       Do i2=-max2,max2

        k1 = i1*b11 + i2*b21
        k2 = i1*b12 + i2*b22

        ks = k1**2 + k2**2

        If(ks.gt.0.0000001) then   ! skip the zero wavenumber

        arg = k1*dxx+k2*dyy

        f = Dcos(arg)

        exx = exx + qxx(i1,i2)*f
        exy = exy + qxy(i1,i2)*f
        eyy = eyy + qyy(i1,i2)*f

c---
      if(Iopt.ne.1) then
c---

        f = Dsin(arg)

        wwx = wwx + ppx(i1,i2)*f
        wwy = wwy + ppy(i1,i2)*f

        wxxx = wxxx + vxxx(i1,i2)*f
        wxxy = wxxy + vxxy(i1,i2)*f
        wyxy = wyxy + vyxy(i1,i2)*f

        wxyx = wxyx + vxyx(i1,i2)*f
        wxyy = wxyy + vxyy(i1,i2)*f
        wyyy = wyyy + vyyy(i1,i2)*f

c---
       end if
c---

c---
       end if
c---

c---
      end do
      end do
c---

      fc = pi4/area

      Gxx = Gxx + exx*fc
      Gxy = Gxy + exy*fc
      Gyx = Gxy
      Gyy = Gyy + eyy*fc

c--------------
c stress tensor
c--------------

      if(Iopt.ne.1) then

      dxx = x-x0    ! important: do not remove
      dyy = y-y0    ! important: do not remove

       wwx = (wwx+dxx)*fc   ! add also the mean pressure gradient
       wwy = (wwy+dyy)*fc   ! add also the mean pressure gradient

       Txxx = Txxx + wxxx*fc - wwx
       Txxy = Txxy + wxxy*fc
       Tyxx = Txxy
       Tyxy = Tyxy + wyxy*fc - wwx

       Txyx = Txyx + wxyx*fc - wwy
       Txyy = Txyy + wxyy*fc
       Tyyx = Txyy
       Tyyy = Tyyy + wyyy*fc - wwy

       px = px + wwx
       py = py + wwy

      end if

c-----
c done
c-----

      return
      end
