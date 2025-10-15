      subroutine fornberg
     +
     +  (x,n,xp
     +  ,m
     +  ,c
     +  )

c===================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===================================

c------------------------------------------------
c This program accompanies the book:
c
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c----------------------------------------
c  Fornberg's algorithm for interpolation 
c  and differentiation
c  of a function defined by N+1 arbitrary 
c  distributed data points.
c
c  Algorithm (6.11.19)
c
c  SYMBOLS:
c  --------
c
c  x .... coordinate for interpolation
c  xp ... x coordinates of prescribed data
c  n .... number of prescribed points is n+1
c  m .... order of highest required derivative
c  c .... fornberg coefficient
c
c  alpha  parameter used in computation of c
c  beta . parameter used in computation of c
c  min .. variable limit
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xp(50),c(0:05,0:10,10)

c-----------
c initialize
c-----------

      c(0,1,1) = 1.0D0

      alpha = 1.0D0

c-----------------
c loop over points
c-----------------

      Do j=2,n+1

       beta = 1.0D0

       min = m 
       if(j-1.lt.m) min = j-1

       Do i=1,j-1
         beta = beta*(xp(j)-xp(i))
         Do k=0,min
          c(k,j,i)=((xp(j)-x)*c(k,j-1,i)
     +           -k*c(k-1,j-1,i))/(xp(j)-xp(i))
         End Do 
       End Do

       Do k=0,min
         c(k,j,i)=(alpha/beta)*(k*c(k-1,j-1,i-1)
     +                 -(xp(j-1)-x)* c(k,j-1,i-1))
       End Do 

       alpha = beta

      End Do

c-----
c done
c-----

      return
      end
