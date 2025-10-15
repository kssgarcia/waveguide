      subroutine sgf_2d_fs 
     +
     +  (Iopt
     +  ,x,y
     +  ,x0,y0
     +  ,Gxx,Gxy
     +  ,Gyx,Gyy
     +  ,Px,Py
     +  ,Txxx,Txxy,Tyxx,Tyxy
     +  ,Txyx,Txyy,Tyyx,Tyyy
     +  )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c-----------------------------------------
c Free-space Green's function (Stokeslet)
c of two-dimensional Stokes flow
c
c see Pozrikidis (1992, p. 61)
c
c  Iopt =  1 compute only G
c       ne 1 compute G, P, T
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

c--------
c prepare
c--------

      dx = x-x0
      dy = y-y0

      dxx = dx*dx
      dxy = dx*dy
      dyy = dy*dy

      r2  = dxx+dyy
      r   = Dsqrt(r2)
      rl  = log(r)
      ri2 = 1.0D0/r2

      Gxx = -rl + dxx*ri2
      Gxy =       dxy*ri2
      Gyx = Gxy
      Gyy = -rl + dyy*ri2

c--------------------------
c compute the stress tensor
c and the pressure vector
c--------------------------

      if(Iopt.gt.1) then

      r4 = r2*r2
      cf = -4.0D0/r4

      Txxx = dxx*dx * cf
      Txxy = dxy*dx * cf
      Tyxx = Txxy
      Tyxy = dyy*dx * cf

      Txyx = Txxy
      Txyy = Tyxy
      Tyyx = Txyy
      Tyyy = dyy*dy * cf

      cf = 2.0D0*ri2

      Px = dx * cf   ! pressure vector
      Py = dy * cf

      end if

c-----
c done
c-----

      return
      end
