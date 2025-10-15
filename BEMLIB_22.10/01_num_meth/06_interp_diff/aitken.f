      subroutine aitken 
     +
     +   (xint
     +   ,n,x
     +   ,f,fint,a
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c----------------------
c  Aitken extrapolation
c----------------------

      Implicit Double Precision (a-h,o-z)
      Dimension x(50),f(50)
      Dimension u(50),a(50,50)

c--------
c prepare
c--------

      n1 = n+1

c-----------
c initialize
c-----------

      Do i=1,n1
        u(i)   = f(i)
        a(i,1) = u(i)
      End Do

c-------
c launch
c-------

      Do m=2,n1
       Do i=m,n1
          u(i)=((xint-x(i))*u(m-1)
     +         -(xint-x(m-1))*u(i))/(x(m-1)-x(i))
          a(i,m)=u(i)
       End Do
      End Do

      fint = u(n1)

c-----
c Done
c-----

      Return
      End
