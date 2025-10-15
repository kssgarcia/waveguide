      subroutine newton
     +
     +   (xint
     +   ,n,x,f
     +   ,fint
     +   ,a
     +   )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c        C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c----------------------
c Newton interpolation
c----------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(50),f(50)
      Dimension u(50),c(0:50),a(50,50)

c--------
c prepare
c--------

      n1 = n+1

      Do i=1,n+1
        u(i)  = f(i)
        a(i,1)= u(i)
      End Do

c---------------
c newton's table
c---------------

      c(0) = u(1)

      Do m=1,n
        Do i=1,n-m+1
         u(i)=(u(i+1)-u(i))/(x(i+m)-x(i))
         a(i,m+1) = u(i)
        End Do
        c(m) = u(1)
      End Do

c---
c evaluation according to (6.2.31)
c---

      fint = c(n)

      Do i=n,1,-1
        fint = fint*(xint-x(i))+c(i-1)
      End Do

c-----
c done
c-----

      return
      end
