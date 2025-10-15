      subroutine cg 
     +  
     +   (n
     +   ,a
     +   ,rhs
     +   ,sln
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
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press, 1998
c------------------------------------------------

c-------------------------------------------------------
c  Solution of a linear symmetric positive-definite
c  system by the conjugate gradients method
c
c  Search vectors are chosen by the method of 
c  Hestenes and Steifel (1952)
c
c
c  SYMBOLS:
c  -------
c
c a .... symmetric positive definite matrix
c n .... size (rows/columns) of matrix a
c rhs .. right hand side vector (e.g. b, as in Ax=b)
c x .... evolving solution vector
c sln .. final solution vector
c       
c p .... search directions
c r .... residual vectors
c alpha. scale parameter for solution update
c beta.. scale parameter for search direction
c-------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),rhs(128),x(0:128,128)
      Dimension p(128,128),r(0:128,128)
      Dimension sln(128)

      Dimension prj(0:128,128)

c---------------------------------------------------
c set the initial values of the vectors x and r (step 0)
c set the first-step value of p
c---------------------------------------------------

      Do i=1,n
        x(0,i) = 0.0D0
        r(0,i) = rhs(i)
        p(1,i) = rhs(i)
      End Do

c----------------------------
c form the sums used in alpha
c----------------------------

      alpha_num = 0.0D0
      alpha_den = 0.0D0

      Do i=1,n
        alpha_num = alpha_num + r(0,i)*r(0,i)
        Do j=1,n
         alpha_den = alpha_den + p(1,i)*a(i,j)*p(1,j)
        End Do
      End Do

      alpha = alpha_num/alpha_den

c---------------------------------------------------
c set first step values of alpha and vectors x and r
c---------------------------------------------------

      Do i=1,n
        x(1,i) = x(0,i)+alpha*p(1,i)
        r(1,i) = r(0,i)
        Do j=1,n
          r(1,i) = r(1,i)-alpha*a(i,j)*p(1,j)
        End Do
      End Do

c-------------------------------------------
c loop through the remaining search vectors
c 2 to n, and compute
c alpha, beta, and vectors p, x, and r
c-------------------------------------------

      Do 98 k=2,n        ! outer loop over search directions

c---
c sums used in beta
c---
   
      beta_num = 0.0D0
      beta_den = 0.0D0
 
      Do i=1,n
        beta_num = beta_num + r(k-1,i)**2
        beta_den = beta_den + r(k-2,i)**2
      End Do

      beta = beta_num/beta_den

      Do i=1,n
        p(k,i) = r(k-1,i)+beta*p(k-1,i)
      End Do

c---
c compute the sums used in alpha
c---

      alpha_num = beta_num
      alpha_den = 0.0D0

      Do i=1,n
        Do j=1,n
          alpha_den = alpha_den + p(k,i)*a(i,j)*p(k,j)
        End Do
      End Do

      alpha = alpha_num/alpha_den

c---
c compute the k'th iterations of vectors x and r
c---

      Do i=1,n
        x(k,i) = x(k-1,i)+alpha*p(k,i)
        r(k,i) = r(k-1,i)
        Do j=1,n
          r(k,i)=r(k,i)-alpha*a(i,j)*p(k,j)
        End Do
      End Do

  98  Continue              ! end of outer loop

c---------------------
c Extract the solution
c---------------------

      Do i=1,n
        sln(i) = x(n,i)
      End Do

c---------------------
c OPTIONAL:
c
c Compute and print the matrix of residual projection
c-----------------------------------------------------

      write (6,*)
      write (6,*) " cg: matrix of residual projection:"
      write (6,*) " ----------------------------------"
      write (6,*)

      Do i=0,n
       Do j=0,n
         prj(i,j) = 0.0D0
         Do k=1,n
          prj(i,j) = prj(i,j) + r(i,k)*r(j,k)
         End Do
       End Do
      End Do

      Do i=0,n
        write (6,100) (prj(i,j),j=0,n)
      End Do

c-----
c Done
c-----

 100  Format (10(1x,f10.5))

      Return
      End
