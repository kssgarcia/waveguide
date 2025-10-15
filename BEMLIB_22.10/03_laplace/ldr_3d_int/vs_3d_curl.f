      subroutine vs_3d_curl
     +
     +   (npts
     +   ,nelm
     +   ,ptl
     +   ,curl
     +   )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------
c Compute the component of the curl of the
c vector potential normal to
c a three-dimensional surface
c at the nodes of a triangular grid
c
c SYMBOLS:
c -------
c
c  ptl(i,j):   vector potential at the ith node, j=1,2,3
c
c  curl(i,1);  normal component of the curl of ptl
c              at the ith node, j=1,2,3
c
c  Dptlpq = d(ptl_q)/dx_p
c
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension ptl(1026,3)

      Dimension Dptlxxa(1026),Dptlyxa(1026),Dptlzxa(1026)
      Dimension Dptlxya(1026),Dptlyya(1026),Dptlzya(1026)
      Dimension Dptlxza(1026),Dptlyza(1026),Dptlzza(1026)

      Dimension curl(1026,3)

      Dimension itally(1026)

      Dimension n    (512,6),nbe (512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6),eet(6)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma
      common/geo2/vna

c-----------
c initialize
c-----------

      Do i=1,npts

        Dptlxxa(i) = 0.0D0
        Dptlyxa(i) = 0.0D0
        Dptlzxa(i) = 0.0D0

        Dptlxya(i) = 0.0D0
        Dptlyya(i) = 0.0D0
        Dptlzya(i) = 0.0D0

        Dptlxza(i) = 0.0D0
        Dptlyza(i) = 0.0D0
        Dptlzza(i) = 0.0D0

        itally(i) = 0

      End Do

c-------------------------------
c compute the Frechet derivative
c of the vector potential at the nodes
c of each triangle
c-------------------------------

      Do k=1,nelm     ! run over elements

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be = beta (k)
      ga = gamma(k)

c---
c triangle coordinates
c of the nodes
c---

      xxi(1) = 0.0D0
      eet(1) = 0.0D0

      xxi(2) = 1.0D0
      eet(2) = 0.0D0

      xxi(3) = 0.0D0
      eet(3) = 1.0D0

      xxi(4) = al
      eet(4) = 0.0D0

      xxi(5) = ga
      eet(5) = 1.0D0-ga

      xxi(6) = 0.0D0
      eet(6) = be

      Do i=1,6    ! run over triangle nodes

        m = n(k,i)   ! global index of local node i on element k

        xi  = xxi(i)
        eta = eet(i)

        call vs_3d_curl_interp
     +
     +     (p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +     ,ptl(i1,1),ptl(i1,2),ptl(i1,3)
     +     ,ptl(i2,1),ptl(i2,2),ptl(i2,3)
     +     ,ptl(i3,1),ptl(i3,2),ptl(i3,3)
     +     ,ptl(i4,1),ptl(i4,2),ptl(i4,3)
     +     ,ptl(i5,1),ptl(i5,2),ptl(i5,3)
     +     ,ptl(i6,1),ptl(i6,2),ptl(i6,3)
     +     ,al,be,ga
     +     ,xi,eta
     +     ,Dptlxx,Dptlyx,Dptlzx
     +     ,Dptlxy,Dptlyy,Dptlzy
     +     ,Dptlxz,Dptlyz,Dptlzz
     +     )

        Dptlxxa(m) = Dptlxxa(m) + Dptlxx
        Dptlyxa(m) = Dptlyxa(m) + Dptlyx
        Dptlzxa(m) = Dptlzxa(m) + Dptlzx

        Dptlxya(m) = Dptlxya(m) + Dptlxy
        Dptlyya(m) = Dptlyya(m) + Dptlyy
        Dptlzya(m) = Dptlzya(m) + Dptlzy

        Dptlxza(m) = Dptlxza(m) + Dptlxz
        Dptlyza(m) = Dptlyza(m) + Dptlyz
        Dptlzza(m) = Dptlzza(m) + Dptlzz

        itally(m) = itally(m)+1

      End Do

      End Do               ! loop over nodes

c--------------------------------------------
c average the Frechet derivative
c over all elements
c
c compute the curl of the potential
c
c project the curl onto the normal direction
c to pick up the normal component
c--------------------------------------------

      Do i=1,npts

        par = float(itally(i))

        Dptlxxa(i) = Dptlxxa(i)/par
        Dptlyxa(i) = Dptlyxa(i)/par
        Dptlzxa(i) = Dptlzxa(i)/par

        Dptlxya(i) = Dptlxya(i)/par
        Dptlyya(i) = Dptlyya(i)/par
        Dptlzya(i) = Dptlzya(i)/par

        Dptlxza(i) = Dptlxza(i)/par
        Dptlyza(i) = Dptlyza(i)/par
        Dptlzza(i) = Dptlzza(i)/par

        curl(i,1) = Dptlyza(i) - Dptlzya(i)
        curl(i,2) = Dptlzxa(i) - Dptlxza(i)
        curl(i,3) = Dptlxya(i) - Dptlyxa(i)

        prj = curl(i,1)*vna(i,1) 
     +      + curl(i,2)*vna(i,2)
     +      + curl(i,3)*vna(i,3)

        curl(i,1) = prj*vna(i,1)
        curl(i,2) = prj*vna(i,2)
        curl(i,3) = prj*vna(i,3)

      End Do

c---
c Done
c---

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
