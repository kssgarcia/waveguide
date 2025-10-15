      subroutine lu_d_t
     +  
     +    (n,c,a,b
     +    ,l,u
     +    )

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

c------------------------------------------------
c Doolittle decomposition of a tridiagonal matrix
c
c Algorithm (2.6.4)
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double precision l(10,10),u(10,10)
      Dimension a(10),b(10),c(10),d(10),e(10)

c-----------
c initialize
c-----------

      Do i=1,n
       Do j=1,n
        l(i,j) = 0.0D0
        u(i,j) = 0.0D0
       End Do
      End Do

c---
c generate the diagonal vectors
c---

      e(1) = a(1)
      d(2) = c(2)/e(1) 

      Do i=2,n-1
        e(i)   = a(i)-d(i)*b(i-1)
        d(i+1) = c(i+1)/e(i)
      End Do

      e(n) = a(n)-d(n)*b(n-1)

c---
c put the vectors into the matrices L and U
c---

      l(1,1) = 1.0D0
      u(1,1) = e(1)
      u(1,2) = b(1)

      Do i=2,n
        l(i,i)   = 1.0D0
        l(i,i-1) = d(i)
        u(i,i)   = e(i)
        u(i,i+1) = b(i)
      End Do
 
c-----
c Done
c-----

      Return
      End
