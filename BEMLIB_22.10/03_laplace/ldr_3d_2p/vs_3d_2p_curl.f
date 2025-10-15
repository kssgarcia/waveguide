      subroutine vs_3d_2p_curl
     +
     +  (npts,nelm
     +  ,ptl
     +  ,curl
     +  )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c=======================================

c----------------------------------
c Computes the normal components of the curl
c of the vector potential
c at the nodes of a triangular grid
c
c SYMBOLS:
c -------
c
c Iedge: interior/edge/corner node index:
c
c     Iedge(i,1) = 0 if global node i is an interior node
c
c     Iedge(i,1) = 1 if global node i is an edge node
c     Iedge(i,2): global index of image node
c
c     Iedge(i,1) = 3 if global node i is a corner node
c     Iedge(i,2): global index of first image node
c     Iedge(i,3): global index of first image node
c     Iedge(i,4): global index of first image node
c
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1090,3)
      Dimension   ne(1090,7)
      Dimension  vna(1090,3)

      Dimension ptl(1090,3)

      Dimension Iedge(1090,4)

      Dimension Dptlxxa(1090),Dptlyxa(1090),Dptlzxa(1090)
      Dimension Dptlxya(1090),Dptlyya(1090),Dptlzya(1090)
      Dimension Dptlxza(1090),Dptlyza(1090),Dptlzza(1090)

      Dimension Dptlxxa_sv(1090),Dptlyxa_sv(1090),Dptlzxa_sv(1090)
      Dimension Dptlxya_sv(1090),Dptlyya_sv(1090),Dptlzya_sv(1090)
      Dimension Dptlxza_sv(1090),Dptlyza_sv(1090),Dptlzza_sv(1090)

      Dimension curl(1090,3)

      Dimension itally(1090)

      Dimension n    (512,6),nbe (512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6),eet(6)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/edgepoints/Iedge

      common/albega/alpha,beta,gamma
      common/geo2/vna

c-----------
c initialize
c-----------

      Do i=1,npts

        Dptlxxa(i) = 0.0
        Dptlyxa(i) = 0.0
        Dptlzxa(i) = 0.0

        Dptlxya(i) = 0.0
        Dptlyya(i) = 0.0
        Dptlzya(i) = 0.0

        Dptlxza(i) = 0.0
        Dptlyza(i) = 0.0
        Dptlzza(i) = 0.0

        itally(i) = 0

      End Do

c-------------------------------
c compute the Frechet derivative
c of the potential at the nodes
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
      be =  beta(k)
      ga = gamma(k)

c---
c triangle coordinates
c of the nodes
c---

      xxi(1) = 0.0
      eet(1) = 0.0

      xxi(2) = 1.0
      eet(2) = 0.0

      xxi(3) = 0.0
      eet(3) = 1.0

      xxi(4) = al
      eet(4) = 0.0

      xxi(5) = ga
      eet(5) = 1.0-ga

      xxi(6) = 0.0
      eet(6) = be

      Do i=1,6

        m = n(k,i)   ! global index of local node i on element k

        xi  = xxi(i)
        eta = eet(i)

        call vs_3d_2p_curl_interp
     +
     +      (p(i1,1),p(i1,2),p(i1,3)
     +      ,p(i2,1),p(i2,2),p(i2,3)
     +      ,p(i3,1),p(i3,2),p(i3,3)
     +      ,p(i4,1),p(i4,2),p(i4,3)
     +      ,p(i5,1),p(i5,2),p(i5,3)
     +      ,p(i6,1),p(i6,2),p(i6,3)
     +
     +      ,ptl(i1,1),ptl(i1,2),ptl(i1,3)
     +      ,ptl(i2,1),ptl(i2,2),ptl(i2,3)
     +      ,ptl(i3,1),ptl(i3,2),ptl(i3,3)
     +      ,ptl(i4,1),ptl(i4,2),ptl(i4,3)
     +      ,ptl(i5,1),ptl(i5,2),ptl(i5,3)
     +      ,ptl(i6,1),ptl(i6,2),ptl(i6,3)
     +
     +      ,al,be,ga
     +      ,xi,eta
     +
     +      ,Dptlxx,Dptlyx,Dptlzx
     +      ,Dptlxy,Dptlyy,Dptlzy
     +      ,Dptlxz,Dptlyz,Dptlzz
     +      )

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

c---
c average the Frechet derivative
c over all elements
c---

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

      End Do

c-----------------------
c Account for the images
c-----------------------

c---
c save the points that have images
c---

      Do 76 i=1,npts

       If(Iedge(i,1).eq.0) Go to 76  ! no images

       Dptlxxa_sv(i) = Dptlxxa(i)
       Dptlyxa_sv(i) = Dptlyxa(i)
       Dptlzxa_sv(i) = Dptlzxa(i)

       Dptlxya_sv(i) = Dptlxya(i)
       Dptlyya_sv(i) = Dptlyya(i)
       Dptlzya_sv(i) = Dptlzya(i)

       Dptlxza_sv(i) = Dptlxza(i)
       Dptlyza_sv(i) = Dptlyza(i)
       Dptlzza_sv(i) = Dptlzza(i)

   76 Continue

c---
c average over the images
c---

      Do 77 i=1,npts

       noim = Iedge(i,1)   ! number of images

       If(noim.eq.0) Go to 77  ! no images

        Do j=2,noim+1     ! Iedge(i,1) images

         image = Iedge(i,j)

         Dptlxxa(i) = Dptlxxa(i)+Dptlxxa_sv(image)
         Dptlyxa(i) = Dptlyxa(i)+Dptlyxa_sv(image)
         Dptlzxa(i) = Dptlzxa(i)+Dptlzxa_sv(image)

         Dptlxya(i) = Dptlxya(i)+Dptlxya_sv(image)
         Dptlyya(i) = Dptlyya(i)+Dptlyya_sv(image)
         Dptlzya(i) = Dptlzya(i)+Dptlzya_sv(image)

         Dptlxza(i) = Dptlxza(i)+Dptlxza_sv(image)
         Dptlyza(i) = Dptlyza(i)+Dptlyza_sv(image)
         Dptlzza(i) = Dptlzza(i)+Dptlzza_sv(image)

        End Do

        noim1 = noim+1

        Dptlxxa(i) = Dptlxxa(i)/noim1
        Dptlyxa(i) = Dptlyxa(i)/noim1
        Dptlzxa(i) = Dptlzxa(i)/noim1

        Dptlxya(i) = Dptlxya(i)/noim1
        Dptlyya(i) = Dptlyya(i)/noim1
        Dptlzya(i) = Dptlzya(i)/noim1

        Dptlxza(i) = Dptlxza(i)/noim1
        Dptlyza(i) = Dptlyza(i)/noim1
        Dptlzza(i) = Dptlzza(i)/noim1

  77  Continue

c---------------------------------
c compute the curl of the potential
c project the curl onto the normal direction
c to pick up the normal component
c---------------------------------

      Do i=1,npts

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
c done
c---

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      Return
      End
