      subroutine lagrange
     +
     +   (xint
     +   ,n
     +   ,x
     +   ,f
     +   ,fint
     +   )

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

c------------------------
c  Lagrange interpolation
c------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(50),f(50)

      n1 = n+1

c----------------------
c Evaluation of the sum
c----------------------

      sum = 0.0D0

      Do i=1,n1
        prod = 1.0D0
        Do j=1,n1
          If(j.ne.i) prod = prod*(xint-x(j))/(x(i)-x(j))
        End Do
        sum = sum + prod*f(i)
      End Do
	
      fint = sum

c-----
c Done
c-----

      Return
      End
