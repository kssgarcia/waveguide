      subroutine cond_num 
     +
     +   (n,a
     +   ,max,eps
     +   ,cond
     +   )

c=============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=============================================

c-----------------------------------------
c Computes the || ||2 condition bumber of 
c a real matrix using the
c power method with forward or inverse iteration
c
c SYMBOLS:
c -------
c
c a:    Matrix whose condition number is to be found
c n:	Size of the matrix
c max:  maximum number of iterations in the power method
c cond: condition number
c------------------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(512,512),b(512,512),x(512)

c-----------------------------------
c compute the normal matrix b = a aT
c-----------------------------------

      Do i=1,n
       Do j=1,n
         b(i,j) = 0.0D0
         Do l=1,n
           b(i,j) = b(i,j)+a(i,l)*a(j,l)
         End Do
       End Do
      End Do

c------
c print
c------

c     write (6,*)
c     write (6,*) " Original matrix"
c     write (6,*) " ---------------"
c     write (6,*)

c     Do i=1,n
c      write (6,102) (a(i,j),j=1,n)
c     End Do

c     write (6,*)
c     write (6,*) " Normal matrix A * AT"
c     write (6,*) " --------------------"
c     write (6,*)

c     Do i=1,n
c      write (6,102) (b(i,j),j=1,n)
c     End Do

c---------------------------
c invoke the power method
c for the maximum eigenvalue
c---------------------------

c--------------
c     write (6,*)
c     write (6,*) " Size of the matrix :",n
c     write (6,*)
c     write (6,*) " Please enter components of the starting vector"
c     write (6,*) " for the forward iteration"
c     write (6,*) " ----------------------------------------------"
c     read  (5,*) (x(i),i=1,n)
c     write (6,*) 
c     write (6,*) " Thank you"
c
      x(1) = 0.5     ! arbitrary values
      Do i=1,n
       x(i) = 0.0
      End Do
      x(n) = 1.0

c--------------

c------
c forward iteration
c------

      Inv = 1

      call power 
     +
     +   (n,b,x
     +   ,Inv,max,eps
     +   ,ev,icount
     +   ,iflag
     +   )

      evmax = ev

      write (6,*) " -------------------------------"
      write (6,*) " Iterations performed: ",icount
      write (6,*) " Maximum eigenvalue is: ",evmax
      write (6,*)
      write (6,*) " Corresponding eigenvector is"
      write (6,102) (x(i),i=1,n)
      write (6,*) " -------------------------------"

c---------------------------
c invoke the power method 
c with inverse iterations
c for the minimum eigenvalue
c---------------------------

c----------------
c     write (6,*)
c     write (6,*) " Size of the matrix :",n
c     write (6,*)
c     write (6,*) " Please enter components of the starting vector"
c     write (6,*) " for the inverse iteration"
c     write (6,*) " ----------------------------------------------"
c     read  (5,*) (x(i),i=1,n)
c     write (6,*)
c     write (6,*) " Thank you"
c
      x(1) = 0.7     ! arbitrary values
      Do i=1,n
       x(i) = 0.0
      End Do
      x(n) = -1.0
c---------------

c------
c inverse iteration
c------

      Inv = 2

      call power 
     +
     +  (n,b,x
     +  ,Inv,max,eps
     +  ,ev,icount
     +  ,Iflag
     +  )

      evmin = 1.0D0/ev

      write (6,*) " -------------------------------"
      write (6,*)
      write (6,*) " Iterations performed: ",icount
      write (6,*) " Minimum eigenvalue is: ",evmin

      write (6,*)
      write (6,*) " Corresponding eigenvector is"
      write (6,102) (x(i),i=1,n)
      write (6,*)
      write (6,*) " -------------------------------"

c---
c condition number
c---

      cond = Dsqrt( abs(evmax)/abs(evmin) )

c---
c Done
c---

 100  Format (1x,i3,f10.5)
 101  Format (2(1x,f20.10))
 102  Format (20(1x,f12.4))

      return
      end
