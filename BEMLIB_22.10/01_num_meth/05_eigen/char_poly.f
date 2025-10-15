      program char_poly

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c----------------------------------------
c Driver for computing the coefficients
c of the characteristic polynomial
c of a square matrix
c---------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension A(10,10),c(0:10)

c----------------
c Read the matrix
c----------------

      open (8,file="matrix.dat")

        read (8,*) N
        Do i=1,N
         read (8,*) (A(i,j),j=1,N)
        End Do

      close (8)

c------------
c preferences
c------------

 98   Continue

      write (6,*)
      write (6,*) " Choose the method"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*)
      write (6,*) " 1 for the Leverrier method"
      write (6,*) " 2 for the Krylov method"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c----------
c launching
c----------

      if(menu.eq.1) then

         write (6,*)
         write (6,*) " Calling Leverrier"
         write (6,*)

         call leverrier (N,A,c)

      Else

         write (6,*)
         write (6,*) " Calling Krylov"
         write (6,*)

         call krylov (N,A,c,Istop)

      End If

c---------
c printing
c---------

      write (6,*)
      write (6,*) " Coeffients of monomials x**N, x**(N-1), ..."
      write (6,*)

      Do i=0,N
       write (6,100) i,c(i)
      End Do

      Go to 98  ! repeat

c-----
c done
c-----

  99  Continue

      write (6,*) " Thank you for running me"

 100  Format (1x,i3,f10.5)

      Stop
      End
