      program mapping

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c              C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press, 1998
c------------------------------------------------

c------------------------------------------------
c  Multiply a vector by a square matrix many times
c
c  SYMBOLS:
c  -------
c
c  n...dimension of the matrix
c  a...the matrix
c  b...the vector
c
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(50,50),b(50),c(50)

c------------
c preferences
c------------

      write (6,*)
      write (6,*) " Normalize the vector length "
      write (6,*) "       after each projection ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes  "
      write (6,*) " ---------------------------"

      read  (5,*) norm

c----------------
c read the matrix
c and the vector
c---------------

      open (unit=8,file="matrix_v.dat",status="unknown")

      read (8,*) n
      Do i=1,n
       read  (8, *)  (a(i,j),j=1,n)
       write (6,100) (a(i,j),j=1,n)
      End Do
       read  (8,*) (b(j),j=1,n)

      close (8)
 
c-----------------------
c carry out the mappings
c-----------------------

      Icount = 0
 98   Continue
      Icount = Icount + 1

      Do i=1,n
        c(i) = 0.0D0
        Do j=1,n
         c(i) = c(i) + a(i,j)*b(j)
        End Do
      End Do

c---
c update the vector
c---

      Do i=1,n
        b(i) = c(i)
      End Do

c---
c normalize the vector
c---

      If(norm.eq.1) then

        Rnorm = 0.0D0
        Do i=1,n
          Rnorm = Rnorm + b(i)**2
        End Do
        Rnorm = sqrt(Rnorm)
        Do i=1,n
           b(i) = b(i)/Rnorm
        End Do

      End If
c---

      write (6,*)
      write (6,*)  " Vector at mapping :",Icount
      write (6,*)
      write (6,101) (b(i),i=1,n)

c---

      write (6,*)
      write (6,*) ' One more projection ?'
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes  "
      write (6,*) ' -------------------------'
      read (5,*) more

      If(more.eq.1) Go to 98

      write (6,*) ' Exiting after ',icount,' iterations'

c-----
c Done
c-----

 100  Format (5(2x,f10.5))
 101  Format (20(2x,f10.5))

      Stop
      End
