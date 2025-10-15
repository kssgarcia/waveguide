      subroutine solve_chol_c 
     +
     +  (n,a,rhs
     +  ,x
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------
c
c This routine accomplishes two tasks:
c
c 1) Cholesky decomposition by column
c    Algorithm 2.6.3
c
c 2) Solution of the nxn linear system:
c          a x = rhs
c
c SYMBOLS:
c -------
c
c tol:  tolerance in the residuals for
c       raising a flag
c
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision l(32,32)
      Dimension a(32,32),rhs(32),xx(32),x(32)

      Parameter (tol=0.0000001)

c-----------------------
c Cholesky Decomposition
c-----------------------

      l(1,1) = Dsqrt(a(1,1))

      Do i=2,n
        l(i,1) = a(i,1)/l(1,1)
      End do

      Do j=2,n

        sum = 0.0D0

        Do m=1,j-1
          sum= sum+l(j,m)**2
        End Do

        l(j,j)=sqrt((a(j,j)-sum))

        Do i=j+1,n
          sum = 0.0D0
          Do m=1,j-1
            sum = sum +l(i,m)*l(j,m)
          End Do
          l(i,j)=(a(i,j)-sum)/l(j,j)
        End Do

      End Do

c------------------------
c Display Cholesky matrix
c------------------------
c
c     write (6,*)
c     Do i= 1,n
c      write (6,100) (l(i,j),j=1,i)
c     End Do

c--------------------
c forward substitution
c--------------------

      xx(1) = rhs(1)/l(1,1)
      Do i=2,n
       sum=rhs(i)
       Do j=1,i-1
        sum=sum-l(i,j)*xx(j)
       End Do
       xx(i)=sum/l(i,i)
      End Do

c-----------------
c back substitution
c------------------

      na = n-1
      x(n) = xx(n)/l(n,n)

      Do i=na,1,-1
        sum=xx(i)
        Do j=i+1,n
         sum=sum-l(j,i)*x(j)
        End Do
        x(i)=sum/l(i,i)
      End Do

c-------
c verify
c-------

      Do i=1,n
       res = rhs(i)
       Do j=1,n
        res = res-a(i,j)*x(j)
       End Do
       If(abs(res).gt.tol) then
         write (6,100) res
         write (6,*) " Cholesky alarm"
       End If
      End Do

c-----
c Done
c-----

  100 Format (100(1x,f12.10))

      Return
      End
