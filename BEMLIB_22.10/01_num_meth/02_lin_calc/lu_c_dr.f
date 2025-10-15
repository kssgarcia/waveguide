      program lu_c_dr

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

c-----------------------
c Crout LU decomposition,
c
c Determinant and inverse
c of an arbitrary matrix
c-----------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision l(100,100),linv(100,100)
      Dimension        u(100,100),a(100,100),v(100,100)
      Dimension        uinv(100,100),ainv(100,100)

c----------------
c Read the matrix
c----------------

      open (unit=8,file='matrix.dat')

      read (8,*) n
      Do i=1,n
        read(8,*) (a(i,j),j=1,n)
      End Do

      close(8)

c---
c Print
c---

      write (6,*)
      write (6,*) " Matrix A"
      write (6,*)

      Do i=1,n
        write(6,1000) (a(i,j),j=1,n)
      End Do 


c----------
c Decompose
c----------

      call crout (n,a,l,u)

c---
c Print
c---

      write (6,*)
      write (6,*) " Matrix L"
      write (6,*)

      Do i=1,n
        write(6,1000) (l(i,j),j=1,n)
      End Do 

      write (6,*)
      write (6,*) " Matrix U"
      write (6,*)

      Do i=1,n
         write(6,1000) (u(i,j),j=1,n)
      End Do 

c---
c Verify the LU decomposition
c---

      Do i=1,n
        Do j=1,n
          v(i,j) = 0.0D0
          k = i
          If(j.lt.k) k = j
           Do m=1,k
             v(i,j) = v(i,j)+l(i,m)*u(m,j)
           End Do 
         End Do 
      End Do 

      write (6,*)
      write (6,*) " Matrix LU"
      write (6,*)

      Do i=1,n
        write(6,1000) (v(i,j),j=1,n)
      End Do 

c-----------------
c Determinant of A
c-----------------

      Det = 1.0D0

      Do i=1,n
       Det = Det*l(i,i)
      End Do

      write (6,101) Det

c-------------
c inverse of A
c-------------

      call inv_l (n,l,linv)
      call inv_u (n,u,uinv)

      Do i=1,n
        Do j=1,n
          ainv(i,j) = 0.0D0
           k = i
           If(j.lt.k) k = j
            Do m=k,n
              ainv(i,j) = ainv(i,j)+uinv(i,m)*linv(m,j)
            End Do 
        End Do 
      End Do 
           
      write (6,*)
      write (6,*) " Inverse of A"
      write (6,*)

      Do i=1,n
        write(6,1000) (ainv(i,j),j=1,n)
      End Do 

c------------------------
c Verify the inverse of A 
c------------------------

      Do i=1,n
        Do j=1,n
          v(i,j) = 0.0D0
          Do m=1,n
            v(i,j) = v(i,j)+a(i,m)*ainv(m,j)
          End Do 
        End Do 
      End Do 
           
      write (6,*)
      write (6,*) " matrix A*A(inv)"
      write (6,*)

      Do i=1,n
       write(6,1000) (v(i,j),j=1,n)
      End Do

c-----
c Done
c-----

 101  Format (/," Determinant of A :",f15.10,/)
1000  format (50(3x,f10.6))

      Stop
      End
