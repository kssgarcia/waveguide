      subroutine bcg (n,a,rhs,xsln)

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------------
c This program accompanies the book:
c
c             C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c          Oxford University Press
c------------------------------------------------

c-----------------------------------------------
c  Solution of an arbitrary system
c  by the method of biconjugate gradients
c
c  Algorithm 3.9.2
c
c  SYMBOLS:
c  --------
c
c  a .... matrix
c  n .... size (rows/columns) of matrix a
c  rhs .. right hand side vector (e.g. b, as in Ax=b)
c
c  p  .... search direction
c  pb .... search direction (bar)
c  r ..... residual vector
c  rb..... residual vector (bar)
c  alpha.. scale parameter for solution update
c  beta... scale parameter for search direction
c
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(100,100),rhs(100),x(0:100,100)
      Dimension p(100,100),pb(100,100)
      Dimension r(0:100,100), rb(0:100,100)
      Dimension xsln(100)

c-----------------------------------------------
c initialize vectors x and r (step 0)
c set first-step value of p
c-----------------------------------------------

      Do i=1,n
         x(0,i) = 0.0D0
         p(1,i) = rhs(i)
         r(0,i) = rhs(i)
        pb(1,i) = rhs(i)
        rb(0,i) = rhs(i)
      End Do

c-------------------------------
c compute the sums used in alpha
c-------------------------------

      alpha_n = 0.0D0
      alpha_d = 0.0D0

      Do i=1,n
        alpha_n = alpha_n + r(0,i)*rb(0,i)
        Do j=1,n
          alpha_d = alpha_d + pb(1,i)*a(i,j)*p(1,j)
        End Do
      End Do

c-------------------------------
c set first-step values of alpha
c and vectors: x, r,rb
c-------------------------------

      alpha = alpha_n/alpha_d

      Do i=1,n
          x(1,i) =  x(0,1)+alpha*p(1,i)
          r(1,i) =  r(0,i)
         rb(1,i) = rb(0,i)
         Do j=1,n
           r(1,i) =  r(1,i)-alpha*a(i,j)* p(1,j)
          rb(1,i) = rb(1,i)-alpha*a(j,i)*pb(1,j)
         End Do
      End Do

c------------------------------------------
c loop through the remaining search vectors
c 2 to n, and compute
c alpha, beta, and vectors p, x, and r
c------------------------------------------

      Do 98 k=2,n       ! outer loop over search directions

c---
c sums used in beta
c---
 
      beta_n = 0.0D0
      beta_d = 0.0D0
 
      Do i=1,n
        beta_n = beta_n + r(k-1,i)*rb(k-1,i)
        beta_d = beta_d + r(k-2,i)*rb(k-2,i)
      End Do

      beta = beta_n/beta_d

c---
c compute p and p(bar)
c---

      Do i=1,n
       p (k,i) = r (k-1,i)+beta*p (k-1,i)
       pb(k,i) = rb(k-1,i)+beta*pb(k-1,i)
      End Do

c---
c sums used in alpha
c---
 
      alpha_n = beta_n
      alpha_d = 0.0D0

      Do i=1,n
        Do j=1,n
         alpha_d = alpha_d+pb(k,i)*a(i,j)*p(k,j)
        End Do
      End Do

      alpha=alpha_n/alpha_d

c---
c form the k'th iterations of vectors x and r
c---
 
      Do i=1,n
       x(k,i) = x(k-1,i)+alpha*p(k,i)
      End Do

      Do i=1,n
        r (k,i) = r (k-1,i)
        rb(k,i) = rb(k-1,i)
        Do j=1,n
         r (k,i) = r (k,i)-alpha*a(i,j)*p (k,j)
         rb(k,i) = rb(k,i)-alpha*a(j,i)*pb(k,j)
        End Do
      End Do

  98  Continue              ! end of outer loop

c---------------------
c extract the solution
c---------------------

      Do i=1,n
        xsln(i) = x(n,i)
      End Do

c-----
c Done
c-----

      return
      end
