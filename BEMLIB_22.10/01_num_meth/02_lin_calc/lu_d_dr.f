      program lu_d_dr

c------------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the 
c stipulations of the licensing agreement
c------------------------------------------

c------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c        Oxford University Press
c------------------------------------------------

c-------------------------------------------
c Driver for the
c Doolittle LU decomposition of a matrix
c Computation of the determinant and inverse 
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double precision l(100,100),linv(100,100)

      Dimension    u(100,100),  a(100,100),v(100,100)
      Dimension uinv(100,100),ainv(100,100)

c------
c input
c------

      write (6,*) 
      write (6,*) " Will perform the Doolittle decomposition"

  97  write (6,*) 
      write (6,*) " Choose the decomposition method"
      write (6,*)
      write (6,*) " Enter"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " 1 to decompose by column"
      write (6,*) " 2 to decompose by row"
      write (6,*) "------------------------"
      read  (5,*) method

      If(method.eq.0) Go to 99

c----------------
c read the matrix
c----------------

      open (unit=8,file='matrix.dat')

      read (8,*) n
      Do i=1,n
        read (8,*) (a(i,j),j=1,n)
      End Do

      close(8)

c------
c print
c------

      write (6,*)
      write (6,*) " Original matrix:"
      write (6,*) " ---------------"
      write (6,*)

      Do i=1,n
        write(6,102) (a(i,j),j=1,n)
      End Do 

c----------
c Decompose
c----------

      call lu_d (method,n,a,l,u)

c------
c print
c------

      write (6,*)
      write (6,*) " Matrix L:"
      write (6,*) " ---------"
      write (6,*)

      Do i=1,n
        write(6,102) (l(i,j),j=1,n)
      End Do 

      write (6,*)
      write (6,*) " matrix U:"
      write (6,*) " --------"
      write (6,*)

      Do i=1,n
        write(6,102) (u(i,j),j=1,n)
      End Do 

c----------------------------
c verify the LU decomposition
c----------------------------

      Do i=1,n
        Do j=1,n
          v(i,j) = 0.0D0
          k = i
          If(j.lt.k) k = j
          Do m=1,k
           v(i,j) = v(i,j) + l(i,m)*u(m,j)
          End Do
        End Do
      End Do

      write (6,*)
      write (6,*) " Matrix LU"
      write (6,*) " --------"
      write (6,*)

      Do i=1,n
        write(6,102) (v(i,j),j=1,n)
      End Do 

c------------
c Determinant
c------------

      Det = 1.0D0

      Do i=1,n
         Det = Det*u(i,i)
      End Do

      write (6,*)
      write (6,101) Det

c-------------
c inverse of A
c-------------

      call inv_l (n,l,linv)
      call inv_u (n,u,uinv)

      Do i=1,n
        Do j=1,n
          ainv(i,j) = 0.0
          k = i
          If(j.lt.k) k = j
            Do m=k,n
             ainv(i,j) = ainv(i,j)+uinv(i,m)*linv(m,j)
            End Do
        End Do
      End Do

      write (6,*)
      write (6,*) " Inverse: "
      write (6,*) " --------"
      write (6,*)

      Do i=1,n
        write(6,102) (ainv(i,j),j=1,n)
      End Do

c---
c Verify the inverse of A
c---

      Do i=1,n
        Do j=1,n
          v(i,j) = 0.D0
            Do m=1,n
             v(i,j) = v(i,j)+a(i,m)*ainv(m,j)
            End Do
        End Do
      End Do

      write (6,*)
      write (6,*) " matrix A*A^(-1)"
      write (6,*) " ---------------"
      write (6,*)

      Do i=1,n
       write(6,102) (v(i,j),j=1,n)
      End Do

      Go to 97

c-----
c Done
c-----

  99  Continue

 102  Format(50(3x,f10.6))
 101  Format(/," Determinant of A = ",f15.9)

      Stop
      End
