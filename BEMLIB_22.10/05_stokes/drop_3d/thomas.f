      subroutine thomas 
     +
     +  (N   ! matrix size
     +  ,a   ! diagonal
     +  ,b   ! super-diagonal
     +  ,c   ! sub-diagonal
     +  ,s   ! rhs
     +  ,x   ! solution 
     +  )

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------------
c This program accompanies the book:
c                 C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c          Oxford University Press
c------------------------------------------------

c------------------------------------------
c Thomas algorithm for tridiagonal systems
c
c Coefficient matrix:
c
c  | a1 b1  0   0  ...  0   0    0    |
c  | c2 a2  b2  0  ...  0   0    0    |
c  | 0  c3  a3  b3 ...  0   0    0    |
c  | ..............................   |
c  | 0  0   0   0  ... cn-1 an-1 bn-1 |
c  | 0  0   0   0  ...  0   cn   an   |
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(1026),b(1026),c(1026),s(1026),x(1026)
      Dimension d(1026),y(1026)

      Parameter (tol=0.00000001D0)

c--------
c prepare
c--------

      Na = N-1

c------------------------------
c reduction to upper bidiagonal
c------------------------------

      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)

      Do i=1,Na
       i1 = i+1
       Den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/Den
       y(i1) =(s(i1)-c(i1)*y(i))/Den
      End Do

c------------------
c Back substitution
c------------------

      x(N) = y(N)

      Do i=Na,1,-1
        x(i)= y(i)-d(i)*x(i+1)
      End Do

c-----------------------
c Verification and alarm
c-----------------------

      Res = s(1)-a(1)*x(1)-b(1)*x(2)

      if(abs(Res).gt.tol) write (6,*) " thomas: alarm"

      Do i=2,Na
        Res = s(i)-c(i)*x(i-1)-a(i)*x(i)-b(i)*x(i+1)
        If(abs(Res).gt.tol) write (6,*) " thomas: alarm"
      End Do

      Res = s(N)-c(N)*x(N-1)-a(N)*x(N)

      if(abs(Res).gt.tol) write (6,*) " thomas: alarm"

c-----
c Done
c-----

 100  Format (1x,f15.10)

      Return
      End
