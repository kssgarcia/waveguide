      subroutine splc_pr_geo
     +   
     +   (N,X,Y
     +   ,vnx,vny
     +   ,crv
     +   ,s
     +   )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c------------------------------------------

c------------------------------------------
c FDLIB
c
c  compute the normal vector (vnx, vny)
c          curvature (crv)
c          arc length (s)
c
c  along a periodic line
c  using cubic-spline interpolation
c
c  Normal vector points downward:
c
c            crv>0
c            ____
c   \ crv<0 /    \ crv<0 /
c    \_____/      \_____/
c       ||
c       ||
c       \/ 
c
c curvature is positive when the
c surface is downward parabolic
c==================================

      Implicit double precision (a-h,o-z)

      Dimension   X(0:513),  Y(0:513)
      Dimension vnx(0:513),vny(0:513)
      Dimension crv(0:513)
      Dimension   s(0:513)

c---
c spline coefficients
c---

      Dimension Xint(513),Yint(513)

      Dimension Axint(513),Bxint(513),Cxint(513)
      Dimension Ayint(513),Byint(513),Cyint(513)

      Parameter (Nstep=32)

c-------
c common
c-------

      common/splc_xy0/Xint
      common/splc_xy1/Axint,Bxint,Cxint
      common/splc_xy2/Ayint,Byint,Cyint

c---------
c announce
c---------

c     write (6,*) " splc_pr_geo: entered",N

c--------
c prepare
c--------

      N1 = N+1
      N2 = N+2

c----------------------------
c interpolation variable
c is the polygonal arc length
c----------------------------

      Xint(1) = 0.0D0

      Do i=2,N1
       ia = i-1
       Xint(i) = Xint(ia)+Dsqrt((X(i)-X(ia))**2
     +                         +(Y(i)-Y(ia))**2)
      End Do

c--------------
c interpolate X
c--------------

      Do i=1,N1
       Yint(i) = X(i)
      End Do

      call splc_pr
     +
     +  (N
     +  ,Xint,Yint
     +  ,Axint,Bxint,Cxint
     +  )

c--------------
c interpolate Y
c--------------

      Do i=1,N1
       Yint(i) = Y(i)
      End Do

      call splc_pr
     +
     +  (N
     +  ,Xint,Yint
     +  ,Ayint,Byint,Cyint
     +  )

c-----------------------
c compute the arc length
c-----------------------

      s(1) = 0.0D0

      Do i=1,N
 
      XX1 = Xint(i)
      XX2 = Xint(i+1)

      step = (XX2-XX1)/Nstep

c---
c integrate by the trapezoidal rule
c---

      DX  = 0.0D0  ! first point

      Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i)
      Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i)

      sum = 0.50D0*Dsqrt( Dxp**2 + Dyp**2 )

      Do j=2,Nstep
        DX = (j-1.0D0)*step
        Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i)
        Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i)
        sum = sum + Dsqrt( Dxp**2 + Dyp**2 )
      End Do

      DX = XX2-XX1

      Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i)
      Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i)

      sum = sum + 0.50D0*Dsqrt( Dxp**2 + Dyp**2 )

      s(i+1)=s(i)+sum*step

      End Do

c-----------------------------------
c compute the downward normal vector
c and the curvature
c at the nodes
c-----------------------------------

      Do i=1,N

       Den = Dsqrt(Cxint(i)**2+Cyint(i)**2)

       vnx(i) =  Cyint(i)/Den
       vny(i) = -Cxint(i)/Den

       crv(i) = 2.0D0*(Bxint(i)*Cyint(i)
     +                -Byint(i)*Cxint(i))/Den**3
      End Do

c------------
c wrap around
c------------

      s(N2) = s(N1)+s(2)

      vnx(N1) = vnx(1)
      vny(N1) = vny(1)
      crv(N1) = crv(1)

        s(0) = s(N)-s(N1)
      vnx(0) = vnx(N)
      vny(0) = vny(N)
      crv(0) = crv(N)

c-----
c done
c-----

      return
      end
