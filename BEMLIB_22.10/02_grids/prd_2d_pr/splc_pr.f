      subroutine splc_pr
     +
     + (N
     + ,xp       ! x-values
     + ,fp       ! function values
     + ,a,b,c    ! cubic-spline coefficients
     + )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c------------------------------------------------

c--------------------------------------
c  Cubic spline interpolation of prescribed data 
c  with periodic boundary conditions
c
c  ith cubic: f = a x^3 + b x^2 + c x + f_i
c
c  SYMBOLS
c  -------
c
c  N  .... number of intervals
c  xp .... x coordinate of prescribed data
c  fp .... function value of prescribed data
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
c N=512 max
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  xp(513),fp(513)
      Dimension   b(513), a(513),c(513)

c--- local:

      Dimension   h(512)
      Dimension  at(512),bt(512),ct(512),rhs(512)
      Dimension sln(512)

c--------
c prepare
c--------

      Na = N-1
      N1 = N+1

c-----------------------
c compute intervals h(i)
c-----------------------

      Do i=1,N
        h(i) = xp(i+1)-xp(i)
      End Do 

      h1 = h(1)
      hN = h(N)

c--------------------------------
C Generate the tridiagonal matrix:
c
c | a1 b1 0  .  . . |
c | c1 a2 b2 0  . . |
c | 0  c2 a2 b2 .
c | 0  0 ...
c | ...             |
c--------------------------------
 
      Do i=1,Na
        i1 = i+1
        at(i) = 2.0D0*(h(i)+h(i1))
        bt(i) = h(i1)
        ct(i) = h(i)
      End Do 

      at(N) = 2.0D0*(h(N)+h(1))
      bt(N) = h(1)
      ct(N) = h(N)

c-----
c generate the right-hand side
c-----

      Do i=1,Na
        i1 = i+1
        i2 = i+2
        rhs(i) = 3.0D0*( (fp(i2)-fp(i1))/h(i1)
     +                  -(fp(i1)-fp(i) )/h(i)  )
      End Do 

      rhs(N) = 3.0D0*( (fp(2) -fp(1))/h(1)
     +                -(fp(N1)-fp(N))/h(N) )

c-----
c solve the tridiagonal system for b_i
c-----

        call thomas_pr
     +
     + (N   ! matrix size
     + ,at  ! diagonal
     + ,bt  ! super-diagonal row
     + ,ct  ! sub-diagonal row
     + ,rhs ! rhs
     + ,sln ! solution
     + )

c-------------------------
c extract the coefficients
c-------------------------

      Do i=1,N
        b(i+1) = sln(i)
      End Do

      b(1) = b(N1)

c---------------------------------
c Compute the coefficients a and c
c---------------------------------

      Do i=1,N
        i1 = i+1
        a(i) = (b(i1)-b(i))/(3.0D0*h(i))
        c(i) = (fp(i1)-fp(i))/h(i) 
     +        - h(i)*(b(i1)+2.0D0*b(i))/3.0D0
      End Do

c-----
c Done
c-----

      return
      end
