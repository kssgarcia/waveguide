      program house_sym_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c 1998
c------------------------------------------------

c----------------------------------------
c  Driver for transforming a real symmetric matrix a
c  a real symmetric tridiagonal matrix  t
c  by Householder reflections.
c
c  Algorithm (5.7.19)
c
c  Evaluation the characteristic polynomial:
c
c  P(lamda) = Deat(t- lamda*I)
c
c  where I is the identity matrix
c 
c  SYMBOLS:
c  -------
c
c  a: primary matrix
c  t: transformed tridiagonal matrix
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),t(10,10)
      Dimension at(10),bt(10),ct(10),f(0:11)
 
c----------------
c read the matrix
c----------------

      open (8,file="matrix_sym.dat")
        read (8,*) n
        Do i=1,n
         read (8,*) (a(i,j),j=1,n)
        End Do
      close (8)

c---
c printing
c---

      write (6,*)
      write (6,*) " Original matrix"
      write (6,*) " ---------------"
      write (6,*)

      Do i=1,n
         write (6,102) (a(i,j),j=1,n)
      End Do

c---
c call Householder
c---

      call house_sym (n,a,t)

c---
c printing
c---

      write (6,*)
      write (6,*) " Similar tridiagonal matrix"
      write (6,*) " --------------------------"
      write (6,*)
      Do i=1,n
         write (6,102) (t(i,j),j=1,n)
      End Do

c---
c Put the diagonals of t in vectors
c--- 

      Do i= 1,n
        at(i)   = t(i,i)
        bt(i)   = t(i,i+1)
        ct(i+1) = t(i+1,i)
      End Do

c-------------------
c compute the traces
c matrices a and t
c-------------------

      trace_a = a(1,1)
      trace_t = at(1)

      Do i=2,n
        trace_a = trace_a + a (i,i)
        trace_t = trace_t + at(i)
      End Do

      write (6,*)
      write (6,105) trace_a
      write (6,106) trace_t

c---------------------------------------
c evaluate the characteristic polynomial
c of the matrix lambda*I - t
c where I is the identity matrix
c---------------------------------------

      write (6,*)
      write (6,*) " Will evaluate the characteristic polynomial"
      write (6,*) " of the matrix lambda*I - t"
      write (6,*) " where I is the identity matrix"
      write (6,*) " ------------------------------"

  98  Continue

      write (6,*)
      write (6,*) " Please enter the value of lamda"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) rlam

      If(rlam.eq.99) Go to 99

c----
c compute the characteristic polynomial
c---

      f(0) = 1.0D0
      f(1) = at(1)-rlam

      Do i=2,n
        f(i) = (at(i)-rlam)  *f(i-1)
     +         -bt(i-1)*ct(i)*f(i-2)
      End Do

      write (6,*)
      write (6,101) f(n)
      write (6,*)

      Go to 98

c-----
c Done
c-----

  99  Continue

      write (6,*) 
      write (6,*) " Thank you for running me"

 101  Format (" Characteristic polynomial :",f15.10)
 102  Format (20(1x,f8.4))
 105  Format ("   Trace of A = ",f10.5)
 106  Format ("   Trace of T = ",f10.5)

      Stop
      End
