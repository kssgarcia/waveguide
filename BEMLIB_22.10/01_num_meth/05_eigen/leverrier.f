      subroutine leverrier (N,A,c)

c======================================
c FDLIB
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c======================================

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c-----------------------------------------
c  Coefficients of the 
c  characteristic polynomial of an nxn matrix 
c  by Leverrier' method
c
c  Algorithm (5.4.3)
c
c  SYMBOLS:
c  --------
c       
c  N .... size (rows/columns) of matrix A
c  A .... square matrix
c  c .... vector of polynomial coefficients
c
c  B .... iterated matrix for computation of c
c  D .... iterated matrix for computation of c
c
c P_N = (-1)^N * (lambda^N + c_1 \lambda^(N-1) + ... + c_N
c
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension A(10,10),B(10,10),D(10,10)
      Dimension c(0:10)

c-----------------------------
c initial elements of matrix B 
c-----------------------------

      Do i=1,N
        Do j=1,N
          B(i,j) = 0.0D0
        End Do
        B(i,i) = 1.0D0
      End Do  

c-----------
c initialize
c-----------

      c(0) = 1.0D0

c---------------------------------------
c outer loop for polynomial coefficients
c---------------------------------------

      Do i=1,N             !  begin outer loop
   
      c(i) = 0.0D0

c---
c update the matrix D
c---

        Do j=1,N
         Do k=1,N
          D(j,k) = 0.0D0
           Do l=1,N
            D(j,k) = D(j,k) + A(j,l)*B(l,k)
           End Do
         End Do
        End Do

c---
c compute c(i)
c---
 
        tr = 0.0D0
 
        Do j=1,N
         tr = tr+D(j,j)
        End Do

        c(i) = -tr/i

c---
c update the matrix B
c---

        Do j=1,N
          Do k=1,N
           B(j,k) = D(j,k)
          End Do
          B(j,j) = B(j,j) + c(i)
        End Do

      End Do !  end of outer loop

c----------------------------------------------
c update coefficients for leading factor (-1)^N
c----------------------------------------------

      fact = (-1.0D0)**N

      Do i=1,N
       c(i) = fact*c(i)
      End Do

c-----
c done
c-----

      Return
      End
