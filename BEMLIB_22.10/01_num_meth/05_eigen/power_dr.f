      program power_dr

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
c ``Numerical Computation in Science and Engineering''
c Oxford University Press, 1998
c------------------------------------------------
c
c------------------------------
c
c  Driver for the power subroutine
c
c  This program computes the eigenvalues
c  of a real matrix using the power method,
c  with options for deflation, shifting,
c  and inverse iteration
c
c  A cascade of deflations allows the computation
c  of ALL eigenvalues
c
c  SYMBOLS:
c  --------
c
c  a ... targeted matrix (square)
c  n ... size of a
c  h ... Householder matrix for deflation
c  b ... shifted or deflated matrix
c  c ... ancillary matrix
c
c  inv = 1 for regular iterations
c      = 2 for inverse iterations
c
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(512,512),b(512,512),c(512,512),h(512,512)
      Dimension x(512),w(512)

c----------------
c read the matrix
c----------------

      open (8,file="matrix.dat")

        read (8,*) n
        Do i=1,n
         read (8,*) (a(i,j),j=1,n)
        End Do

      close (8)

c------------
c preferences
c------------

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 for regular iterations"
      write (6,*) " 2 for inverse iterations"
      write (6,*) " ------------------------"
      read  (5,*) inv

      If(inv.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter maximum number of iterations"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) max

      If(max.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter desired accuracy"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) eps

      If(eps.eq.0) Go to 99

c------------------
c targer for repeat
c------------------

 96   Continue

      write (6,*)
      write (6,*) " Targeted matrix:"
      write (6,*) " ----------------"
      write (6,*)

      Do i=1,n
       write (6,102) (a(i,j),j=1,n)
      End Do

c---

      write (6,*)
      write (6,*) " The matrix will be shifted diagonally"
      write (6,*)
      write (6,*) " Please enter the shift"
      write (6,*) " ----------------------"
      read  (5,*) shift

      write (6,*)
      write (6,*) " Size of the matrix :",n
      write (6,*)
      write (6,*) " Please enter the components of the starting vector"
      write (6,*) " --------------------------------------------------"
      read  (5,*) (x(i),i=1,n)
      write (6,*)
      write (6,*) " Thank you"

c--------------------
c shift the diagonals
c--------------------

      Do i=1,n
       Do j=1,n
         b(i,j) = a(i,j)
       End do
       b(i,i) = b(i,i)-shift
      End do

c--------------------------
c call the power subroutine
c--------------------------

      call power 
     +
     +   (n,b,x
     +   ,inv,max,eps
     +   ,ev
     +   ,Icount
     +   ,iflag
     +   )

c----------------------
c invert if appropriate
c----------------------

      If(inv.eq.2) ev = 1.0D0/ev

c------------------
c unshift and reset
c------------------

      ev = ev+shift

      Do i=1,n
       b(i,i) = b(i,i)+shift
      End do

      write (6,*) " -------------------------------"
      write (6,*)
      write (6,*) " Iterations performed: ",icount
      write (6,*) " Eigenvalue is:"
      write (6,101) ev 

      write (6,*)
      write (6,*) " Corresponding eigenvector is"
      write (6,102) (x(i),i=1,n)
      write (6,*)
      write (6,*) " -------------------------------"

c-------------------
c proceed to deflate
c-------------------

      write (6,*)
      write (6,*) " Will continue to deflate"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 for Householder deflation"
      write (6,*) "       2 for deflation with a row"
      write (6,*) "       3 to repeat"
      write (6,*) " ---------------------------------"
      read  (5,*) idfl

      If(Idfl.eq.0) Go to 99
      If(Idfl.eq.3) Go to 96

c----------------------
c Householder deflation
c----------------------

      If(Idfl.eq.1) then

c--
c make up the vector w
c---

      sign = 1.0
      If(x(1).lt.0) sign = -1.0

      w(1) = Dsqrt(0.5D0*(1.0+sign*x(1)))

      Do i=2,n
        w(i) = sign*0.5D0*x(i)/w(1)
      End Do

c---
c make up the Householder matrix
c---

      Do i=1,n
        Do j=1,n
         h(i,j) = -2.0*w(i)*w(j)
        End Do
        h(i,i) = h(i,i)+1.0
      End Do

      write (6,*)
      write (6,*) " Householder matrix:"
      write (6,*) " -------------------"
      write (6,*)

      Do i=1,n
       write (6,102) (h(i,j),j=1,n)
      End Do

c--
c carry out the similarity transformation
c---

      Do i=1,n
        Do j=1,n
          c(i,j) = 0.0D0
          Do m=1,n
           c(i,j) = c(i,j)+a(i,m)*h(m,j)
          End Do
        End Do
      End Do 

      Do i=1,n
        Do j=1,n
          b(i,j) = 0.0D0
          Do m=1,n
           b(i,j) = b(i,j)+h(i,m)*c(m,j)
          End Do
        End Do
      End Do 

      write (6,*)
      write (6,*) " Transformed matrix:"
      write (6,*) " -------------------"
      write (6,*)

      Do i=1,n
       write (6,102) (b(i,j),j=1,n)
      End Do

      If(n.eq.2) Go to 99

      write (6,*)
      write (6,*) " Compute one more eigenvalue ?"
      write (6,*)
      write (6,*) " Enter 1 for yes, 0 to quit"
      write (6,*) " ---------------------------"
      read  (5,*) more

      If(more.eq.0) Go to 99

c---
c deflate the matrix
c---

      n = n-1    ! reduce the matrix size

      Do i=1,n
        Do j=1,n
         a(i,j) = b(i+1,j+1)
        End do
      End Do

      Go to 96

      End If

c---------------------
c Deflation with a row
c---------------------

      If(idfl.eq.2) then

       m = 1

c---
c find the maximum element
c of the eigenvector
c and its location
c---

       dfl = Dabs(x(1))

       Do i=2,n
         If(abs(x(i)).gt.dfl) then
           m = i
           dfl = abs(x(i))
         End If
       End Do

c---
c deflate
c---

       dfl = x(m)

       Do i=1,n
        x(i) = x(i)/dfl
        Do j = 1,n
         b(i,j) = a(i,j) - x(i)*a(m,j)
        End Do
       End Do

       write (6,*)
       write (6,*) " Deflated matrix:"
       write (6,*) " ----------------"
       write (6,*)

       Do i=1,n
        write (6,102) (b(i,j),j=1,n)
       End Do

       If(n.eq.2) Go to 99

       write (6,*)
       write (6,*) " Compute one more eigenvalue ?"
       write (6,*)
       write (6,*) " Enter 1 for yes, 0 to quit"
       write (6,*) " ---------------------------"
       read  (5,*) more

       If(more.eq.0) Go to 99

c---
c discard the mth column
c---

        Do i=1,n
         Do j=m,n-1
           b(i,j) = b(i,j+1)
         End Do
        End Do

c---
c discard the mth row
c---

        Do i=m,n-1
         Do j=1,n-1
          b(i,j) = b(i+1,j)
         End Do
        End Do

        Do i=1,n
          Do j=1,n
           a(i,j) = b(i,j)
          End Do
        End Do

        n = n-1             ! reduce the matrix size

        Go to 96

      End If

c-----
c Done
c-----

  99  Continue

      write (6,*) " Thank you for running me"

 100  Format (1x,i3,f10.5)
 101  Format (2(1x,f20.10))
 102  Format (20(1x,f8.4))

      Stop
      End
