      program power_sym_dr

c-----------------------------------------
c FDLIB
c
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
c Oxford University Press, 1998
c------------------------------------------------

c-----------------------------
c  Driver for the power method 
c  for a symmetric matrix
c-----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),b(10,10),x(10)

c------
c input
c------

      open (8,file="matrix_sym.dat")

        read (8,*) n
        Do i=1,n
         read (8,*) (a(i,j),j=1,n)
        End Do

      close (8)

c---
c preferences
c---

  98  Continue

      write (6,*)
      write (6,*) " Enter maximum number of iterations"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) max

      If(max.eq.0) Go to 99

      write (6,*)
      write (6,*) " The diagonals will be shifted"
      write (6,*)
      write (6,*) " Please enter the shift"
      write (6,*) " -----------------------"
      read  (5,*) shift

      Do i=1,n
        Do j=1,n
           b(i,j) = a(i,j)
        End Do
        b(i,i) = b(i,i) - shift
      End Do

      write (6,*)
      write (6,*) " Size of the matrix is :",n
      write (6,*)
      write (6,*) " Please enter components of the starting vector"
      write (6,*) " ----------------------------------------------"
      read  (5,*) (x(i),i=1,n)

      write (6,*) " Thank you"

c---
c call the power method
c---

      call power_sym 
     +
     +   (n,b,x,max
     +   ,ev1,ev2
     +   ,Icount
     +   ,Iflag
     +   )

      ev1 = ev1 + shift
      ev2 = ev2 + shift

      write (6,*)
      write (6,*) " Iterations performed: ",icount
      write (6,*)
      write (6,*) " Dominant eigenvalues are:"
      write (6,*)
      write (6,101) ev1,ev2
      write (6,*)

c-----
c Done
c-----

 100  Format (1x,i3,f10.5)
 101  Format (2(1x,f20.10))

      Go to 98

  99  write (6,*) " Thank you for running me"

      Stop
      End
