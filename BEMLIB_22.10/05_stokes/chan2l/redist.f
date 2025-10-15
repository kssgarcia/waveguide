      subroutine redist (N,X,Y)

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c---------------------------------
c redistribute nodes evenly with respect to x
c inside one period
c---------------------------------

      Implicit double precision (a-h,o-z)

      Dimension   X(0:513),Y(0:513)

      Dimension  Xint(513),Yint(513)

      Dimension Axint(513),Bxint(513),Cxint(513)
      Dimension Ayint(513),Byint(513),Cyint(513)

      N1 = N+1

c----------------------------------
c define the interpolation variable
c----------------------------------

      Xint(1) = 0.0D0

      Do i=2,N1
       ia = i-1
       Xint(i) = Xint(ia)+Dsqrt((X(i)-X(ia))**2
     +                         +(Y(i)-Y(ia))**2)
      End Do

c----
c spline interpolation for x
c----

      Do i=1,N1
       Yint(i) = X(i)
      End Do

      call splc_pr
     +
     +  (N
     +  ,Xint,Yint
     +  ,Axint,Bxint,Cxint
     +  )

c----
c spline interpolation for y
c----

      Do i=1,N1
       Yint(i) = Y(i)
      End Do

      call splc_pr
     +
     +  (N
     +  ,Xint,Yint
     +  ,Ayint,Byint,Cyint
     +  )

c--------------
c evenly spaced nodex
c nodes 1 and N+1 remain fixed
c--------------

      step = Xint(N1)/N

      Do i=2,N

       XX = step*(i-1.0D0)

       Do j=1,N
        prod = (XX-Xint(j))*(XX-Xint(j+1))
        if(prod.le.0) Go to 3
       End Do

   3   Continue

       D = XX-Xint(j)

       X(i) = X(j) + D*(Cxint(j)+D*(Bxint(j)+D*Axint(j)))
       Y(i) = Y(j) + D*(Cyint(j)+D*(Byint(j)+D*Ayint(j)))
     
      End Do

c-----
c Done
c-----

      Return
      End
