      program lanczos_sym_dr

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------
c Driver for the Lanczos reduction of a real 
c and symmetric matrix into a similar,
c real and symmetric tridiagonal matrix "t"
c-----------------------------------------

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

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),p(10,10),t(10,10)
      Dimension alpha(10),beta(10),f(0:11)
 
c----------------
c read the matrix
c----------------

      open (8,file="matrix_sym.dat")
        read (8,*) n
        Do i=1,n
         read (8,*) (a(i,j),j=1,n)
        End Do
      close (8)

c---------
c printing
c--------

      write (6,*)
      write (6,*) " Matrix:"
      write (6,*) " ------"
      write (6,*)

      Do i=1,n
        write (6,102) (a(i,j),j=1,n)
      End Do

c---------------------
c enter search vector
c--------------------


      write (6,*)
      write (6,*) " Size of the matrix :",n
      write (6,*)

      write (6,*) " Please enter components of the starting vector"
      write (6,*) " ----------------------------------------------"
      read (5,*) (p(1,i),i=1,n)
      write (6,*) 
      write (6,*) " Thank you"

c----------
c Transform
c----------

      call lanczos_s (n,a,p,alpha,beta,istop)

c-------------------------------
c make up the tridiagonal matrix
c-------------------------------

	Do i=1,n
	Do j=1,n
	t(i,j) = 0.0D0
	End Do
	End Do

	Do i=1,n
	t(i,i)   = alpha(i)
	t(i,i+1) =  beta(i)
	t(i+1,i) =  beta(i)
	End Do

c---------
c printing
c---------

      write (6,*)
      write (6,*) " Transformed Tridiagonal Matrix"
      write (6,*) " ------------------------------"
      write (6,*)

      Do i=1,n
        write (6,102) (t(i,j),j=1,n)
      End Do

c-------------------
c compare the traces
c-------------------

      trace_a = a(1,1)
      trace_t = alpha(1)

      Do i=2,n
        trace_a = trace_a + a(i,i)
        trace_t = trace_t + alpha(i)
      End Do

      write (6,105) trace_a
      write (6,106) trace_t

c---
      write (6,*)
  98  write (6,*) " Will evaluate the characteristic polynomial"
      write (6,*)
      write (6,*) " Please enter the value of lamda"
      write (6,*)
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) rlam
      If(rlam.eq.99) Go to 99

      f(0) = 1.0D0
      f(1) = alpha(1)-rlam

      Do i=2,n
        f(i) = (alpha(i)-rlam)*f(i-1)-beta(i-1)**2*f(i-2)
      End Do

      write (6,101) f(n)

      Go to 98

c-----
c Done
c-----

  99  write (6,*) " Thank you for running me"

 101  Format (2(1x,f20.10))
 102  Format (20(1x,f15.10))
 105  Format ("Trace of A = ",f10.5)
 106  Format ("Trace of T = ",f10.5)

      Stop
      End
