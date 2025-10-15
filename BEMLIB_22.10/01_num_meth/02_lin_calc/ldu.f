      subroutine ldu (n,a,l,d,u)

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
c LDU decomposition of a symmetric matrix
c
c Algorithm (2.6.1) with modification discussed
c on page 88
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision l(10,10)

      Dimension a(10,10),u(10,10),d(10)

c-----------
c initialize
c-----------

      Do i=1,n
       Do j=1,n
        l(i,j) = 0.0D0
        u(i,j) = 0.0D0
       End Do
      End Do

c-------
c launch
c-------

      d(1)  = a(1,1)
      Do i=1,n
        u(1,i) = a(1,i)
        l(i,1) = u(1,i)/d(1)
      End Do

      Do k=2,n

       Do j=k,n
        sum=0.0D0
          Do m=1,k-1
            sum = sum + l(k,m)*u(m,j)
	  End Do
	u(k,j) = a(k,j) - sum
       End Do

       d(k) = u(k,k)
       Do j=k,n
         l(j,k) = u(k,j)/d(k)
       End Do
      End Do

c-----
c Done
c-----

      Return
      End
