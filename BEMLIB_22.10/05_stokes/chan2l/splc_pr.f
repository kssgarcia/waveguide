      subroutine splc_pr
     +
     + (N
     + ,x       ! x-values
     + ,f       ! function values
     + ,a,b,c   ! cubic-spline coefficients
     + )

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c              C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c          Oxford University Press
c------------------------------------------------

c--------------------------------------
c  Cubic spline interpolation of prescribed data 
c  with periodic boundary conditions
c
c  the ith cubic is: f = a x^3 + b x^2 + c x + f_i
c
c  SYMBOLS
c  -------
c
c  N  .... number of intervals
c  x .... x coordinate of prescribed data
c  f .... function value of prescribed data
c	
c  c .... polynomial coefficient related to 1st derivative
c  b .... polynomial coefficient related to 2nd derivative
c  a .... polynomial coefficient related to 3rd derivative
c  h .... interval between prescribed data
c
c  at.... diagonal of tridiagonal matrix
c  bt.... superdiagonal of tridiagonal matrix
c  ct.... subdiagonal of tridiagonal matrix
c
c CAPACITY:
c --------
c
c N = 512 max
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  x(513),f(513)
      Dimension  b(513),a(513),c(513)

c--- local:

      Dimension   h(512)
      Dimension  at(512),bt(512),ct(512),rhs(512)
      Dimension sln(512)

c-----------------------
c compute intervals h(i)
c-----------------------

      Do i=1,N
        h(i) = x(i+1)-x(i)
      End Do 

c--------------------------------
C generate a bordered tridiagonal matrix:
c
c | at1 bt1   0    0  . . . 0        ct1    |
c | ct2 at2   bt2  0  . . . 0        0      |
c | 0   ct3   at3 bt3 . . . 0        0      |
c | 0   0 ...               0        0      |
c | 0   0 ...               bt(N-2)  0      |
c | 0   0 ... ct(N-1)       at(N-1)  bt(N-1)|
c | btN 0 ...  0            ctN      atN    |
c--------------------------------
 
      Do i=1,N
        at(i) = 2.0D0*(h(i)+h(i+1))
        bt(i) = h(i+1)
        ct(i) = h(i)
      End Do 

c-----
c generate the right-hand side
c-----

      Do i=1,N-1
        rhs(i) = 3.0D0*( (f(i+2)-f(i+1))/h(i+1)
     +                  -(f(i+1)-f(i) )/h(i)  )
      End Do 

      rhs(N) = 3.0D0*( (f(2) -f(1))/h(1)
     +                -(f(N+1)-f(N))/h(N) )

c-----
c solve the tridiagonal system for b_i
c-----

        call thomas_pr
     +
     +  (N   ! matrix size
     +  ,at  ! diagonal
     +  ,bt  ! super-diagonal row
     +  ,ct  ! sub-diagonal row
     +  ,rhs ! rhs
     +  ,sln ! solution
     +  )

c-------------------------
c extract the coefficients b
c-------------------------

      Do i=1,N
        b(i+1) = sln(i)
      End Do

      b(1) = b(N+1)

c---------------------------------
c compute the coefficients a and c
c---------------------------------

      Do i=1,N
        a(i) = (b(i+1)-b(i))/(3.0D0*h(i))
        c(i) = (f(i+1)-f(i))/h(i) 
     +        - h(i)*(b(i+1)+2.0D0*b(i))/3.0D0
      End Do

      a(N+1) = a(1)
      c(N+1) = c(1)

c-----
c done
c-----

      return
      end
