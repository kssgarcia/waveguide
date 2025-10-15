      subroutine krylov (n,a,c,Istop)

c=========================================
c                 FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c       This program accompanies the book:
c              C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c----------------------------------------------------------
c  Coefficients of the characteristic polynomial
c  of an nxn matrix by the method of Krylov sequences. 
c
c  Algorithm (5.4.7)
c
c  It is possible that the coefficient matrix 
c  of the linear system may be ill-conditioned.
c
c  If this occurs, the coefficients will be in error.
c  Exercise caution if the determinant of the
c  linear system is small.
c
c  SYMBOLS:
c  --------
c
c  a .... square matrix
c  n .... size (rows/columns) of matrix a
c  c .... vector of polynomial coefficients
c  Istop. flag for success of linear system solver
c         if Istop = 1, something went wrong
c 
c  x .... vector created by transformation with matrix a
c  b .... iterated matrix for computation of d
c  d .... solution vector from linear system solver
c
c----------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),x(10,0:10),c(0:10)
      Dimension b(10,10),rhs(10),d(10)

c-------------------------------
c establish initial vector x(0,i)
c--------------------------------

      Do j=1,n-2
        x(j,0)=1.0D0
      End Do

      x(n-1,0) = 0.5D0
      x(n,  0) = 0.5D0

c---------------------------------------
c compute the Krylov sequence of vectors
c---------------------------------------

      Do i=1,n
       Do j=1,n
         x(j,i)=0.0D0
           Do k=1,n
            x(j,i) = x(j,i)+a(j,k)*x(k,i-1)
           End Do
        End Do
      End Do

c------------------
c form the matrix b
c------------------

      Do i=1,n
       Do j=1,n
         b(j,i)=x(j,i-1)
       End Do
      End Do

c---
c set vector x(j,n) 
c as the rhs of
c the linear system [b]d=x.
c---

      Do j=1,n
        rhs(j)=x(j,n)
      End Do

c---
c  call the linear solver
c---

      Isym   = 0      ! system is not symmetric
      Iwlpvt = 0      ! pivoting is not enabled

      call gel_krylov
     +
     +  (n,b,rhs,d
     +  ,Isym,Iwlpvt
c    +  ,l,u
     +  ,det
     +  ,Istop
     +  )

      If(Istop.eq.1) then
       write (6,*)
       write (6,*) " krylov: linear solver failed"
       write (6,*)
       Return
      End If

      write (6,*)
      write (6,*) " krylov: linear system and residuals"
      write (6,*)

      Do i=1,n
        res = rhs(i)   ! residual
        Do j=1,n
        res = res-b(i,j)*d(j)
        End Do
        write (6,101) (b(i,j),j=1,n),rhs(i),d(i),res
      End Do

c-----------------------------------------
c generate the standard coefficient vector
c-----------------------------------------

      c(0) = 1.0D0

      Do i=1,n
        c(0) = - c(0)
      End Do

      Do i=1,n
       c(i)=-d(n+1-i)*c(0)
      End Do

c------------
c     Do i=0,n
c      write (6,100) i,c(i),d(i)
c     End Do
c------------

c---
c Done
c---

 100  Format (1x,i3,3(1x,f10.5))
 101  Format (20(1x,f10.5))

      Return
      End
