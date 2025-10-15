      subroutine arc_3d 
     +
     +  (x1,x2,x3
     +  ,y1,y2,y3
     +  ,z1,z2,z3
     +  ,xc,yc,zc
     +  ,a
     +  ,chi1,chi3
     +  )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c----------------------------------------
c  Computes the properties of the circular arc
c  passing through
c  three points x1,x2,x3 in 3D space
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c 1998
c
c see pages: 298-299
c------------------------------------------------

c----------
c SYMBOLS:
c 
c  xc,yc,zc     coordinates of the arc center
c  a            radius of the arc
c  chi1, chi3   angles subtended by the second-first 
c               and second-third point
c               as described in text
c----------

      Implicit Double Precision (a-h,o-z)

c---
c Equations (6.8.8)
c---

      a11 = 2.0*(x1-x2)
      a12 = 2.0*(y1-y2)
      a13 = 2.0*(z1-z2)
      a21 = 2.0*(x3-x2)
      a22 = 2.0*(y3-y2)
      a23 = 2.0*(z3-z2)

      tmp = x2**2+y2**2+z2**2
      b1  = x1**2+y1**2+z1**2 - tmp
      b2  = x3**2+y3**2+z3**2 - tmp

c---
c Equation (6.8.9)
c---

      crx = (y3-y1)*(z2-z1)-(z3-z1)*(y2-y1)
      cry = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1)
      crz = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1)

      a31  = crx
      a32  = cry
      a33  = crz
      b3   = x1*crx + y1*cry + z1*crz 

c     write (6,*)
c     write (6,*) " linear system"
c     write (6,*)
c     write (6,100) a11,a12,a13,b1
c     write (6,100) a21,a22,a23,b2
c     write (6,100) a31,a32,a33,b3
c     write (6,100)

c---
c Solve the 3x3 system
c---

      call cramer_33 
     +
     +   (A11,A12,A13
     +   ,A21,A22,A23
     +   ,A31,A32,A33
     +   ,B1,B2,B3
     +   ,xc,yc,zc
     +   )

c---
c compute arc radius
c---

      as = (x1-xc)**2+(y1-yc)**2+(z1-zc)**2
      a  = sqrt(as)

c---
c Compute the signed angles chi1 and chi3 
c
c ch1 < 0, chi2 = 0, chi3 > 0
c
c---

      prj1 = (x1-xc)*(x2-xc)+(y1-yc)*(y2-yc)+(z1-zc)*(z2-zc)
      prj3 = (x3-xc)*(x2-xc)+(y3-yc)*(y2-yc)+(z3-zc)*(z2-zc)
      chi1 = - acos(prj1/as)
      chi3 =   acos(prj3/as)

c-----
c Done
c-----

 100  Format (10(1x,f15.10))

      Return
      End
