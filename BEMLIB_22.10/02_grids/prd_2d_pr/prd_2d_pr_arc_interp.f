      subroutine prd_2d_pr_arc_interp
     +
     +   (x1,x2,x3,x
     +   ,u1,u2,u3,u
     +   ,v1,v2,v3,v
     +   ,c1,c2,c3,c
     +   ,w1,w2,w3,w
     +   ,z1,z2,z3,z
     +   )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------
c Quadratic Lagrange interpolation
c---------------------------------

      Implicit Double Precision (a-h,o-z)

      p1 = (x-x2)*(x-x3)/((x1-x2)*(x1-x3))
      p2 = (x-x1)*(x-x3)/((x2-x1)*(x2-x3))
      p3 = (x-x1)*(x-x2)/((x3-x1)*(x3-x2))

      u = u1*p1 + u2*p2 + u3*p3
      v = v1*p1 + v2*p2 + v3*p3
      c = c1*p1 + c2*p2 + c3*p3
      w = w1*p1 + w2*p2 + w3*p3
      z = z1*p1 + z2*p2 + z3*p3

c-----
c Done
c-----

      return
      end
