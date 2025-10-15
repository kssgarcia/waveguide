      program cond_num_dr

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c ``Numerical Computation in Science and Engineering''
c
c Oxford University Press, 1998
c------------------------------------------------

c----------------------------------------
c Driver to compute the || ||2 condition bumber 
c of a real matrix using the
c power method with regular and inverse iteration
c
c SYMBOLS:
c -------
c
c a ... matrix
c n ... matrix size
c max.. maximum number of iterations for the power method
c
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(512,512)

c----------------
c read the matrix
c----------------

  98  Continue

      open (8,file="matrix.dat")

        read (8,*) n
        Do i=1,n
         read (8,*) (a(i,j),j=1,n)
        End Do

      close (8)

c----------
c inquiries
c----------

      write (6,*) 
      write (6,*) " Enter the maximum number of iterations"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) max

      If(max.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter desired accuracy of eigenvalues"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) eps

      If(eps.eq.0) Go to 99

c-----------------------------
c compute the condition number
c-----------------------------

      call cond_num
     +
     +   (n,a
     +   ,max,eps
     +   ,cond
     +   )

c---------
c printing
c---------

      write (6,*)
      write (6,*) " Condition number = ",cond
      write (6,*)

      Go to 98    ! repeat

c-----
c Done
c-----

  99  Continue

      write (6,*) " Thank you for running me"

      Stop
      End
