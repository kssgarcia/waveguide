      subroutine sdlp_3d
     +
     +   (npts
     +   ,nelm
     +   ,Intm   ! integration method
     +   ,mint
     +   ,nsp    ! spectral connectivity
     +   ,psp    ! spectral global nodes
     +
     +   ,cx,cy,cz
     +
     +   ,q      ! spectral dipole
     +
     +   ,qsintx ! surface integral of q
     +   ,qsinty 
     +   ,qsintz
     +
     +   ,qsinmx ! surface integral of x ^ q
     +   ,qsinmy 
     +   ,qsinmz
     +
     +   ,qdn     ! surface integral of q . n
     +
     +   ,dlp    ! double layer potential
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c----------------------------------
c Compute the principal value of the
c double-layer potential
c of a scalar function q
c at the nodes of a triangular grid
c on a closed surface
c
c slso compute the surface integral of q
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c spectral:

      Dimension   nsp(512,100)   ! connectivity
      Dimension   psp(20000,3)
      Dimension   q(20000,3),dlp(20000,3)

      Parameter (tol=0.00000001)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/trq/xiq,etq,wq

      common/spectral2/nesp,ngsp

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------------------
c launch the quadrature
c----------------------

c     Do i=1,ngsp
c       write (6,*) q(i,1),q(i,2),q(i,3)
c     End Do

      Do 1 i=1,ngsp   ! loop over the global spectral nodes

c      write (6,*) " Computing the dlp potential at global node",i

       x0 = psp(i,1)
       y0 = psp(i,2)
       z0 = psp(i,3)

       qx0 = q(i,1)
       qy0 = q(i,2)
       qz0 = q(i,3)

c      write (6,100) i,x0,y0,z0,q0

c-----------
c initialize
c-----------

       srf_area = 0.0D0

       ptlx = 0.0D0
       ptly = 0.0D0
       ptlz = 0.0D0

       qdn = 0.0D0

       qsintx = 0.0D0
       qsinty = 0.0D0
       qsintz = 0.0D0

       qsinmx = 0.0D0
       qsinmy = 0.0D0
       qsinmz = 0.0D0

c-----------------------------------------
c Compile the integrals over the triangles
c-----------------------------------------

       Do k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

c---
c apply the quadrature
c---

        call sdlp_3d_integral
     +
     +    (x0,y0,z0
     +    ,nsp
     +    ,psp
     +    ,k
     +    ,Intm
     +    ,mint
     +    ,cx,cy,cz
     +    ,q
     +    ,qx0,qy0,qz0
     +
     +    ,pptlx,pptly,pptlz
     +    ,aarel
     +    ,qqsintx,qqsinty,qqsintz
     +    ,qqsinmx,qqsinmy,qqsinmz
     +    ,qqdn
     +    )

c       write (6,*) q(k,1),q(k,2),q(k,3),qqdn
c       write (6,*) pptlx,pptly,pptlz
c       write (6,*) aarel
c       pause

        ptlx = ptlx + pptlx
        ptly = ptly + pptly
        ptlz = ptlz + pptlz

c       write (6,*) ptlx

        srf_area = srf_area+aarel

        qsintx = qsintx + qqsintx
        qsinty = qsinty + qqsinty
        qsintz = qsintz + qqsintz

        qsinmx = qsinmx + qqsinmx
        qsinmy = qsinmy + qqsinmy
        qsinmz = qsinmz + qqsinmz

        qdn = qdn + qqdn

      End Do

c--------------------
c      If(i.eq.7) then
c       write (6,*)
c       write (6,*) " total surface area computed in ldlp: ",srf_area
c       write (6,*)
c      End If
c--------------------

c---
c account for the principal value
c---

c      write (6,100) i,ptlx,ptly,ptlz

       dlp(i,1) = ptlx - pi4 * qx0
       dlp(i,2) = ptly - pi4 * qy0
       dlp(i,3) = ptlz - pi4 * qz0

c      write (6,100) i,ptlx,ptly,ptlz
c      write (6,100) i,x0,y0,z0,dlp(i,1),dlp(i,2),dlp(i,3),srf_area ! ,qx0,qy0,qz0

  1   Continue               ! loop over nodes
c       pause

c-----
c Done
c-----

  99  Continue

 100  Format (1x,i3,10(f12.4))
 101  Format (f12.8)

      Return
      End
