      subroutine thomas_pr
     +
     +  (N   ! matrix size
     +  ,a   ! diagonal 
     +  ,b   ! super-diagonal row
     +  ,c   ! sub-diagonal row
     +  ,s   ! rhs
     +  ,x   ! solution 
     +  )

c-----------------------------------------
c FDLIB, BEMLIB
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
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c      Oxford University Press
c------------------------------------------------

c------------------------------------------
c Thomas algorithm for nearly tridiagonal systems
c with wrapped first and last elements,
c corresponding to a periodic boundary condition
c
c  T . x = s
c 
c Coefficient matrix:
c
c      | a1 b1  0   0  ...  0   0    c1   |
c      | c2 a2  b2  0  ...  0   0    0    |
c      | 0  c3  a3  b3 ...  0   0    0    |
c T =  | ..............................   |
c      | 0  0   0   0  ... cn-1 an-1 bn-1 |
c      | bn 0   0   0  ...  0   cn   an   |
c
c
c The solution is found using a
c modification of Algorithm 3.4.1, p. 143
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(129),b(129),c(129),s(129),x(129)
      Dimension d(129),y(129)

      Dimension x0(129)

      Parameter (tol=0.00000001D0)

c--------
c prepare
c--------

      Na = N-1
      Nb = N-2

      save1  = s(1)
      saveNa = s(Na)

c------------------------------
c First assume that x(N) = 0
c and solve the first N-1 equations
c neglecting the last column
c and the last row
c------------------------------

      x(N) = 0.0D0

c-- REGULAR THOMAS

c-----
c reduction to upper bidiagonal
c-----

      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)

      Do i=1,Nb
       i1 = i+1
       Den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/Den
       y(i1) = (s(i1)-c(i1)*y(i))/Den
      End Do

c----
c back substitution
c----

      x(Na) = y(Na)

      Do i=Nb,1,-1
        x(i)= y(i)-d(i)*x(i+1)
      End Do

c-----
c compute the first residual, defined in terms of the
c last equation
c-----

      R0 = a(N)*x(N) + b(N)*x(1) + c(N)*x(Na) - s(N)

c-----
c save the solution
c-----

      Do i=1,N
        x0(i) = x(i)
      End Do

c------------------------------
c Second, assume that x(N) = 1
c and solve the first N-1 equations
c with a modified RHS
c------------------------------

      x(N) = 1.0D0

      s(1)  = s(1)  - c(1)  * x(N)
      s(Na) = s(Na) - b(Na) * x(N)


c-- REGULAR THOMAS

c-----
c reduction to upper bidiagonal
c-----

      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)

      Do i=1,Nb
       i1 = i+1
       Den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/Den
       y(i1) = (s(i1)-c(i1)*y(i))/Den
      End Do

c-----
c back substitution
c-----

      x(Na) = y(Na)

      Do i=Nb,1,-1
        x(i)= y(i)-d(i)*x(i+1)
      End Do

c------
c compute the second residual:
c-----

      R1 = a(N)*x(N) + b(N)*x(1) + c(N)*x(Na) - s(N)

c-----
c rectify the rhs
c-----

      s(1)  = save1
      s(Na) = saveNa

c---------------------------------
c The residual takes the form:
c
c R = (R1-R0) x(N) + R0
c
c compute the correct value of x(N)
c to make R=0
c---------------------------------

      x(N) = -R0/(R1-R0)

c----------------------------
c compose the solution vector
c using
c
c x = (x1-x0) x(N) + x0
c----------------------------

      Do i=1,Na
       x(i) = (x(i)-x0(i)) * x(N) + x0(i)
      End Do

c-----------------------
c Verification and alarm
c-----------------------

      Res = s(1) - a(1)*x(1) -b(1)*x(2)-c(1)*x(N)

      If(abs(Res).gt.tol) write (6,*) " thomas_pr: alarm, 1",Res

      Do i=2,Na
        Res = s(i)-c(i)*x(i-1)-a(i)*x(i)-b(i)*x(i+1)
        If(abs(Res).gt.tol) write (6,*) " thomas_pr: alarm ",i,Res
      End Do

      Res = s(N)-c(N)*x(Na)-a(N)*x(N)-b(N)*x(1)

      If(abs(Res).gt.tol) write (6,*) " thomas_pr: alarm ",N,Res

c-----
c done
c-----

 100  Format (1x,f15.10)

      return
      end
