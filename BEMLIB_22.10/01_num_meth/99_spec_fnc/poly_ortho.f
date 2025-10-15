      program poly_ortho

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
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c------------------------------------------------
c  Generate the coefficients of the first Nmax
c  orthogonal polynomials of a certain family
c  using Gram-Schmidt orthogonalization
c
c  To obtain the coefficients of a desired class 
c  that is not implemented in this code,
c  write a function poly_ortho_xxx
c  reflecting the weighting function w(x).
c  Search for "menu" to identify places where changes need to 
c  be made.
c  In addition, change the normalization according
c  to standard convention.
c
c  Legendre polynomials
c  --------------------
c
c  The orthogonalization integral is defined 
c  in the interval [-1,1]
c  The weighting function is w(x) = 1
c
c  Radau polynomials
c  -----------------
c
c  The orthogonalization integral is defined 
c  over the interval [-1,1]
c  The weighting function is w(x) = 1+x
c
c  Log polynomials
c  ---------------
c
c  The orthogonalization integral is defined 
c  over the interval [0,1]
c  The weighting function is w(x) = -ln(x)
c
c
c  SYMBOLS:
c  -------
c
c  G(n) = <pn|pn>: vector of the norms of each orthogonal polynomial
c  B:	Coeffient matrix defined as follows:
c
c  p0(x) = B(0,0)
c  p1(x) = x   + B(1,0) * p0(x)
c  p2(x) = x^2 + B(2,1) * p1(x) + B(2,0)*p0(x)
c  ...
c
c  pi(x) = x^i + B(i,i-1)* pi-1(x) + ... + B(i,0)*p0(x)
c
c----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension B(0:21,0:21),G(0:21),A(0:21,0:21)

      common/pii/pi,srpi

c----------
c constants
c----------

      pi = 3.14159265358 D0
      srpi = sqrt(pi)

c------------
c preferences
c------------

 97   Continue

      write (6,*)
      write (6,*) " Select the polynomical family"
      write (6,*)
      write (6,*) " Enter"
      write (6,*)
      write (6,*) " 1 for Legendre polynomials"
      write (6,*) " 2 for Radau polynomials"
      write (6,*) " 3 for log polynomials"
      write (6,*) " 4 for polynomials with w(x) = x"
      write (6,*) " 5 for Laguerre polynomials"
      write (6,*) " 6 for Hermite polynomials"
      write (6,*) " 0 to quit"
      write (6,*) "----------" 
      read  (5,*) menu

      if(menu.eq.0) go to 99

 98   Continue

      write (6,*)
      write (6,*) " Enter the maximum polynomial degree desired"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) "----------" 
      read  (5,*) Nmax

      if(Nmax.eq.0) go to 97

c--------------------
c coefficient for p_0
c--------------------

      B(0,0) = 1.0D0

      if(menu.eq.1.or.menu.eq.2) then
        G(0) = 2.0D0
      else if(menu.eq.6) then
        G(0) = srpi
      else
        G(0) = 1.0D0
      end if

c-------------------
c initialize G and B 
c-------------------

      Do i=1,Nmax
        G(i) = 0.0D0
        Do j=0,Nmax
          B(i,j) = 0.0D0
        End Do
      End Do

c----------------------------------
c compute coefficients B(i,j)
c by Gram-Schmidt orthogonalization
c----------------------------------

      Do i=1,Nmax
        B(i,i) = 1.0D0
        Do j=0,i-1
          B(i,j) = - poly_ortho_prj(menu,B,i,j)
          B(i,j) = B(i,j)/G(j)
        End Do
        G(i) = poly_ortho_G(menu,B,G,i)
      End Do
        
c---------------------------
c print the G-S coefficients
c---------------------------

      write (6,*)
      write (6,*) 'The polynomials have the form:'
      write (6,*)
      write (6,*) 'p_i = x^i + B(i,i-1)*p(i-1) + ... + B(i,0)*p(0)'
      write (6,*)
      write (6,*) 'matrix B:'
      write (6,*)

      nstep = -1
      if(menu.eq.1) nstep = -2

      Do i=0,Nmax
        if(i.gt.12) then
          write(6,100) i, (B(i,j),j=11,0,nstep)
          write(6,101)    (B(i,j),j=i,12,nstep)
        else
          write(6,100) i, (B(i,j),j=i,0,nstep)
        end if
      End Do

c----------------------------------
c recover the monomial coefficients
c
c p(i) = A(i,i)*x^i + A(i,i-1)*x^(i-1) +... +A(i,0)'
c----------------------------------

      Do i=0,Nmax
        A(i,i) = B(i,i)
        Do j=0,i-1
          A(i,j) = 0.0D0
          Do l=0,i-1
            A(i,j) = A(i,j)+B(i,l)*A(l,j)
          End Do
        End Do
      End Do

c-------------------------------------
c Legendre and Radau polynomials
c standard normalization is: p_i(x=1) = 1.0
c--------------------------------------

      if(menu.eq.1.or.menu.eq.2) then

        Do i=0,Nmax
          sum = 1.0D0
          Do j=0,i-1
           sum = sum + A(i,j)
          End Do
          Do j=0,i
           A(i,j) = A(i,j)/sum
          End Do
        End Do

      End If

c-------------------------------------
c Laguerre polynomials
c standard normalization is: A(i,i)=(-1)^i
c--------------------------------------

      if(menu.eq.5) then

        fc=1.0D0
        Do i=0,Nmax
          Do j=0,i
           A(i,j) = fc*A(i,j)
          End Do
        fc=-fc
        End Do

      End If

c-------------------------------------
c Hermite polynomials
c standard normalization is: A(i,i)=2^i
c--------------------------------------

      if(menu.eq.6) then

        fc=1.0D0
        Do i=0,Nmax
          Do j=0,i
           A(i,j) = fc*A(i,j)
          End Do
        fc=2.0D0*fc
        End Do

      End If

c---------
c printing
c---------

      write (6,*)
      write (6,*) 'The polynomials have the form:'
      write (6,*)
      write (6,*) 'P(i) = A(i,i)*x^i + A(i,i-1)*x^(i-1) +... +A(i,0)'
      write (6,*)
      write (6,*) "The matrix A is:"
      write (6,*)

      Do i=0,Nmax
        If (i.gt. 12) then
          write(6,100) i,(A(i,j),j=11,0,nstep)
          write(6,101)   (A(i,j),j=i,12,nstep)
        Else
          write(6,100) i,(A(i,j),j=i,0,nstep)
        End If
      End Do

c-------
c repeat
c-------

      go to 98

c-----
c Done
c-----

 99   Continue

 100  Format (i2,12(1x,f10.4))
 101  Format (2x,12(1x,f10.4))

      stop
      end

