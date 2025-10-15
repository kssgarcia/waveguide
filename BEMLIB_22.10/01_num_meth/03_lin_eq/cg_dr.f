      program cg_dr

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
c             C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press, 1998
c------------------------------------------------

c-----------------------------------------
c Solves the linear system Ax = b
c by the method of conjugate gradients
c
c "A" is real symmetric and positive-definite
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),b(128),x(128)

c----------------------------------------
c read the matrix and the right-hand side
c----------------------------------------

      open (3, file="mat_s_vec.dat")

      read (3,*) n

      Do i=1,n
        read  (3,*) (a(i,j),j=1,n)
      End Do

      read (3,*) 

      read (3,*) (b(j),j=1,n)

      close (3)

c--------
c printing
c---------

      write (6,*)
      write (6,*) " matrix A:"
      write (6,*) " ---------"
      write (6,*)

      Do i=1,n
       write (6,101) (a(i,j),j=1,n)
      End Do

      write (6,*)
      write (6,*) " Right-hand side:"
      write (6,*) " ----------------"
      write (6,*)

      write (6,101) (b(j),j=1,n)

c---------------------------
c Do the conjugate gradients
c---------------------------

      call cg (n,a,b,x)

      write (6,*)
      write (6,*) " Solution vector:"
      write (6,*) " ---------------"
      write (6,*) 

      Do i=1,n
        write (6,100) i,x(i)
      End Do
 
c----------
c residuals
c----------

      write (6,*) 
      write (6,*) " Residuals"
      write (6,*) " -----------"
      write (6,*) 

      Do i=1,n
        res = b(i)
        Do j=1,n
         res = res - a(i,j)*x(j)
        End Do
        write (6,100) i,res
      End Do

c-----
c Done
c-----

      write (6,*) "Thank you for running me"

 100  Format (1x,i3,1x,f10.6)
 101  Format (20(1x,f8.4))
 102  Format (f10.5)

      Stop
      End
