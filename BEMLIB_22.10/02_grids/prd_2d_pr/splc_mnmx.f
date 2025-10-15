      subroutine splc_mnmx
     +
     +   (N,X,Y
     +   ,Z
     +   ,Zmax,Xmax
     +   ,Zmin,Xmin
     +   )

c--------------------------------
c FDLIB, BEMLIB
c
c Computation of the minimum and maximum
c of a function Z defined over a periodic
c line described by the data (X, Y)
c
c The maxima and minima are computed by direct search
c on the cubic-splines interpolated line
c
c SYMBOLS:
c -------
c
c N:  Number of spline intervals
c M: Number of search points within each interval
c
c Axint, Bxint, Cxint: spline coeeficients for X(s)
c Azint, Bzint, Czint: spline coeeficients for Z(s)
c
c s = Xint:  interpolation variable: polygonal arc length
c
c Zmax, Zmin: extreme of Z
c Xmax, Xmin: corresponding x position
c--------------------------------

      Implicit double precision (a-h,o-z)

      Dimension X(0:513),Y(0:513),Z(0:513)

      Dimension Xint(513),Yint(513)

c---
c spline coefficients
c---

      Dimension Axint(513),Bxint(513),Cxint(513)
      Dimension Azint(513),Bzint(513),Czint(513)

      Parameter (M=16)

c     write (6,*) " splc_int: entered",N

c--------
c prepare
c--------

      N1 = N+1
      M1 = M+1

c-----------------------
c interpolation variable
c-----------------------

      Xint(1) = 0.0D0

      Do i=2,N1
       ia = i-1
       Xint(i) = Xint(ia)+Dsqrt((X(i)-X(ia))**2
     +                         +(Y(i)-Y(ia))**2)
      End Do

c------
c interpolated variable: X
c-------

      Do i=1,N1
       Yint(i) = X(i)
      End Do

      call splc_pr
     +
     +  (N
     +  ,Xint,Yint
     +  ,Axint,Bxint,Cxint
     +  )

c------
c interpolated variable: Z
c-------

      Do i=1,N1
       Yint(i) = Z(i)
      End Do

      call splc_pr
     +
     +  (N
     +  ,Xint,Yint
     +  ,Azint,Bzint,Czint
     +  )

c---
c find the maximum
c---

      Zmax = -100.00D0
      Zmin =  100.00D0

      Do i=1,N

       step = (Xint(i+1)-Xint(i))/M

       Do j=1,M1

        D = (j-1.0D0)*step

        ZZ = Z(i) + D*(Czint(i)+D*(Bzint(i)+D*Azint(i)))

        If(ZZ.gt.Zmax) then
          Zmax = ZZ
          Xmax = X(i) + D*(Cxint(i)+D*(Bxint(i)+D*Axint(i)))
        End If
        If(ZZ.lt.Zmin) then
          Zmin = ZZ
          Xmin = X(i) + D*(Cxint(i)+D*(Bxint(i)+D*Axint(i)))
        End If

       End Do
       
      End Do

c-----
c Done
c-----

      return
      end
