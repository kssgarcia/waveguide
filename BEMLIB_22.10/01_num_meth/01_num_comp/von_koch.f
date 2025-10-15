      subroutine von_koch (RL,m,alpha,mp)

c=================================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=================================================

c----------------------------------------------------
c Compute the vertices of the periodic von Koch line
c
c  SYMBOLS:
c  -------
c xaux, yaux	:	auxiliary points
c----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xfr (4001),yfr (4001),ifr (4001)
      Dimension xaux(4001),yaux(4001),iaux(4001)

      common/Intg/ifr
      common/Real/xfr,yfr

c----------
c constants
c----------

      three  = 3.0D0
      sr3    = Dsqrt(three)
      alphac = 1.0D0 - alpha

c--------------
c m=1 line 
c---------------

      height = Dsqrt( 3.0D0*alpha**2+2.0D0*alpha-1.0D0)

      mp =  4  ! number of vertices

      xfr(4) =  RL 
      yfr(4) =  RL*height
      ifr(4) =  1

      xfr(3) =  RL*alpha
      yfr(3) =  0.0D0
      ifr(3) =  2

      xfr(2) = -xfr(3)
      yfr(2) =  yfr(3)
      ifr(2) =  2

      xfr(1) = -xfr(4)
      yfr(1) =  yfr(4)
      ifr(1) =  1

c------------------------
      if(m.eq.1) Go to 99   ! return
c------------------------

      Do 1 i=2,m

c------------------------
c Save the vertices
c of the previous shape
c into an auxiliary vector
c-------------------------

      kp = mp

      Do j=1,kp
        xaux(j) = xfr(j)
        yaux(j) = yfr(j)
        iaux(j) = ifr(j)
      End Do

      mp = 3*4**(i-1)+1  ! new number of vertices

      xfr(mp) = xaux(kp)  ! last point is the same
      yfr(mp) = yaux(kp)
      ifr(mp) = iaux(kp)

      j = mp

      Do l=1,kp-1

      k  = kp-l
      xa = xaux(k)
      ya = yaux(k)
      xb = xaux(k+1)
      yb = yaux(k+1)

      xc = 0.5D0*(xa+xb)
      yc = 0.5D0*(ya+yb)
      dx = xb-xa
      dy = yb-ya
      ds = Dsqrt(dx**2+dy**2)
      height = Dsqrt(3.0D0*(4.0D0*alpha-1.0))/6.0D0

      j = j-1
      xfr(j) = xc + alphac*dx/3.0D0
      yfr(j) = yc + alphac*dy/3.0D0

      ifr(j) = 2
      j = j-1
      xfr(j) = xc - dy * height 
      yfr(j) = yc + dx * height 

      ifr(j) = 1
      j = j-1
      xfr(j) = xc - alphac*dx/3.0D0
      yfr(j) = yc - alphac*dy/3.0D0

      ifr(j) = 2
      j = j-1
      xfr(j) = xa
      yfr(j) = ya
      ifr(j) = ifr(k)
      
      End Do

  1   Continue

c-----
c done
c-----

  99  Continue

      return
      end
