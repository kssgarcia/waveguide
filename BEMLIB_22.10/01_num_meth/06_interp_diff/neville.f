      subroutine Neville (xint,n,x,f,fint,a)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
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
c------------------------------------------------

c--------------------------
c Neville interpolation
c
c Algorithm 6.2.1
c-------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(50),f(50)
      Dimension u(50),a(50,50)

c---
c initialize
c---

      n1 = n+1

      Do i=1,n1
        u(i)   = f(i)
        a(i,1) = u(i)
      End Do

c---
c launch
c---

      Do m=2,n1
       Do i=1,n-m+2

        j = m+i-1

c---- alternative expression
c
c        u(i)=((xint-x(i))*u(i+1)+(x(j)-xint)*u(i))
c    +         /(x(j)-x(i))
c-------------------------------------------------

         u(i)=u(i)+(u(i+1)-u(i))/(1.0+(x(j)-xint)/(xint-x(i)))
        a(i,m)=u(i)
       End Do
      End Do

      fint = u(1)

c----
c Done
c-----

      Return
      End
