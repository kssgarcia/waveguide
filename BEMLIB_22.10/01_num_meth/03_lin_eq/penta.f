      subroutine penta
     +
     +  (n
     +  ,a,b,c,d,e,s
     +  ,x
     +  )

c===========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c===========================================

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press
c------------------------------------------------

c----------------------------------------
c Solution of a linear system
c with a pentadiagonal coefficient matrix
c----------------------------------------

      Implicit Double Precision (a-h,o-z)
      Parameter (m=200)
      Dimension e(m),d(m),a(m),b(m),c(m),s(m)
      Dimension d1(m),a1(m),b1(m),c1(m), s1(m)
      Dimension             b2(m),c2(m), s2(m)
      Dimension x(m)

c------------------------------
c Generate the system: Q x = s1 
c------------------------------

      Do i=1,n-2
        c1(i) = c(i)
      End Do

      c1(n-1) = 0.0D0  ! irrelevant
      c1(n)   = 0.0D0  ! irrelevant

      a1(1) = a(1)
      b1(1) = b(1)
      s1(1) = s(1)

      d1(2) = d(2)
      a1(2) = a(2)
      b1(2) = b(2)
      s1(2) = s(2)

      Do i=2,n-1
        i1 = i+1
        w  = e(i1)/d1(i)
        a1(i1) = a(i1) - w*b1(i)
        b1(i1) = b(i1) - w*c1(i)
        d1(i1) = d(i1) - w*a1(i)
        s1(i1) = s(i1) - w*s1(i)
      End Do

c-----------------------------
c  Generate the system T x= s2
c-----------------------------

      b2(1) = b1(1)/a1(1)
      c2(1) = c1(1)/a1(1)
      s2(1) = s1(1)/a1(1)

      Do i=1,n-1
       i1 = i+1
        den =  a1(i1)-d1(i1)*b2(i)
        b2(i1) = (b1(i1)-d1(i1)*c2(i))/den
        c2(i1) =  c1(i1) /den
        s2(i1) = (s1(i1)-d1(i1)*s2(i))/den
      End Do

c-------------------
c  Back substitution
c-------------------

      x(n)  = s2(n)
      x(n-1)= s2(n-1)-b2(n-1)*x(n)

      Do i=n-2,1,-1
        x(i) = s2(i)-b2(i)*x(i+1)-c2(i)*x(i+2)
      End Do

c-----
c Done
c-----

      Return
      End
