      subroutine lu_d
     +  
     +  (Method
     +  ,n,a
     +  ,l,u
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c     Oxford University Press
c------------------------------------------------

c---------------------------------
c Doolittle LU decomposition
c
c If method = 1 use algorithm (2.6.1)
c If method = 2 use algorithm (2.6.2)
c---------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision l(100,100)
      Dimension        u(100,100),a(100,100)

c-----------
c initialize
c-----------

      Do i=1,n
       Do j=1,n
        l(i,j)=0.0D0
        u(i,j)=0.0D0
       End Do
      End Do

c---------------------
c first row and column
c---------------------

      Do i=1,n
       l(i,i)=1.0D0
       u(1,i)=a(1,i)
       l(i,1)=a(i,1)/u(1,1)
      End Do

c----------
c main loop
c----------

      Do k=2,n  ! loop over rows or columns

c-----------
      If(method.eq.1) then
c-----------

       Do j=k,n
        u(k,j)=a(k,j)
        Do m=1,k-1
          u(k,j) = u(k,j) - l(k,m)*u(m,j)
        End Do
       End Do

c-----------
      Else If(method.eq.2) then
c-----------

        Do j=1,k
         u(j,k)=a(j,k)
         Do m=1,j-1
           u(j,k) = u(j,k) - l(j,m)*u(m,k)
         End Do
       End Do

c-----------
      End If
c-----------

      Do i=k+1,n
        l(i,k) = a(i,k)
        Do m=1,k-1
         l(i,k) =l(i,k) - l(i,m)*u(m,k)
        End Do
       l(i,k) = l(i,k)/u(k,k)
      End Do

c-----------
      End Do   ! end of main loop
c-----------

c-----
c Done
c-----

      Return
      End
