      subroutine prd_2d_pr_splc_interp
     +
     +   (x1,x2,x3,x4
     +   ,x
     +   ,u1,u2,u3,u4,u
     +   ,v1,v2,v3,v4,v
     +   ,c1,c2,c3,c4,c
     +   ,w1,w2,w3,w4,w
     +   ,z1,z2,z3,z4,z
     +   )

c-----------------------------------------
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

      p1 = (x-x2)*(x-x3)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4))
      p2 = (x-x1)*(x-x3)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4))
      p3 = (x-x1)*(x-x2)*(x-x4)/((x3-x1)*(x3-x2)*(x3-x4))
      p4 = (x-x1)*(x-x2)*(x-x3)/((x4-x1)*(x4-x2)*(x4-x3))

      u = u1*p1 + u2*p2 + u3*p3 + u4*p4
      v = v1*p1 + v2*p2 + v3*p3 + v4*p4
      c = c1*p1 + c2*p2 + c3*p3 + c4*p4
      w = w1*p1 + w2*p2 + w3*p3 + w4*p4
      z = z1*p1 + z2*p2 + z3*p3 + z4*p4

c     write (6,*) x1,x2,x3,x4,x
c     write (6,*) u1,u2,u3,u4

c-----
c Done
c-----

      Return
      End
