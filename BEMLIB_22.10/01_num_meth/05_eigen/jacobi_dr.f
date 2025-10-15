      program jacobi_dr

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
c Driver for Jacobi's method
c Diagonalization of a real symmetric matrix
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),b(10,10),u(10,10),uv(10,10)
      Dimension eigen(10)

c----------------
c read the matrix
c----------------

      open (8,file="matrix_sym.dat")

        read (8,*) n
        Do i=1,n
         read (8,*) (a(i,j),j=1,n)
        End Do

      close (8)

      write (6,*) 
      write (6,*) " Please enter maximum number of iterations"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) max

      If(max.eq.0) Go to 99

c-----------------
c print the matrix
c-----------------

      write (6,*)
      write (6,*) " Matrix"
      write (6,*) " ------"
      write (6,*)

      Do i=1,n
         write (6,102) (a(i,j),j=1,n)
      End do

c------------
c call jacobi
c------------

      call jacobi 
     +
     +   (n,a,b,u
     +   ,max,iflag
     +   ,icount
     +   )


      write (6,*) " -------------------------------"
      write (6,*)
      write (6,*) " Iterations performed: ",icount
      write (6,*)
      write (6,*) " ------------------"
      write (6,*) " Transformed matrix"
      write (6,*) " ------------------"
      write (6,*)

      Do i=1,n
        write (6,102) (b(i,j),j=1,n)
      End Do

      write (6,*) " -------------------------------"
      write (6,*)
      write (6,*) " ----------------------"
      write (6,*) " Matrix of eigenvectors"
      write (6,*) " ----------------------"
      write (6,*)

      Do i=1,n
         write (6,102) (u(i,j),j=1,n)
      End Do

      write (6,*) " -------------------------------"

c--------------------
c verify eigenvectors
c--------------------

      Do i=1,n

       Do l=1,n
        eigen(l) = u(l,i)    ! ith eigenvector
       End Do

       Do j=1,n
         uv(j,i) = 0.0D0
         Do l = 1,n
          uv(j,i) = uv(j,i)+a(j,l)*eigen(l)
         End Do
       End Do

       Do j=1,n
        uv(j,i) = uv(j,i)/b(i,i)
       End Do

      End Do

      write (6,*)
      write (6,*) " -------------------------------"
      write (6,*) " Verified matrix of eigenvectors"
      write (6,*) " -------------------------------"
      write (6,*)
      Do i=1,n
        write (6,102) (uv(i,j),j=1,n)
      End Do

c-----
c Done
c-----

  99  Continue

      write (6,*)
      write (6,*) " Thank you for running me"

 100  Format (1x,i3,f10.5)
 101  Format (2(1x,f20.10))
 102  Format (20(1x,f8.4))

      Stop
      End
