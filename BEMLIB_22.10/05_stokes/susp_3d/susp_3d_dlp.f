      subroutine susp_3d_dlp
     +
     +   (npts
     +   ,nelm
     +   ,intm   ! integration method
     +   ,mint
     +   ,mpoly
     +   ,necl,ngcl
     +   ,pcl    ! collocation global nodes
     +
     +   ,cx,cy,cz
     +
     +   ,q      ! collocation dipole
     +
     +   ,qsintx ! surface integral of q
     +   ,qsinty 
     +   ,qsintz
     +
     +   ,qsinmx ! surface integral of x ^ q
     +   ,qsinmy 
     +   ,qsinmz
     +
     +   ,qdn        ! surface integral of q . n
     +   ,stresslet  ! stresslet
     +
     +   ,dlp    ! double layer potential
     +   )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------
c Compute the principal value of the
c Stokes double-layer potential
c at the nodes of a triangular grid
c on a closed surface
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension     n(512,6), nbe(512,3)

      Dimension stresslet(3,3),sstresslet(3,3)

c collocation nodes:

      Dimension   pcl(2306,3)
      Dimension     q(2306,3),dlp(2306,3)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------------------
c launch the quadrature
c----------------------

c     Do i=1,ngcl
c       write (6,*) q(i,1),q(i,2),q(i,3)
c     End Do

      Do 1 i=1,ngcl   ! loop over the global collocation nodes

c      write (6,*) " Computing the dlp potential at global node",i

       x0 = pcl(i,1)
       y0 = pcl(i,2)
       z0 = pcl(i,3)

       qx0 = q(i,1)
       qy0 = q(i,2)
       qz0 = q(i,3)

c      write (6,100) i,x0,y0,z0,q0

c-----------
c initialize
c-----------

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

       Do ii=1,3
        Do jj=1,3
         stresslet(ii,jj) = 0.0D0
        End Do
       End Do

       srf_area = 0.0D0

c-----------------------------------------
c Compile the integrals over the triangles
c-----------------------------------------

       Do k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

c---
c apply the quadrature
c---

        call susp_3d_dlp_integral
     +
     +    (x0,y0,z0
     +    ,mpoly
     +    ,pcl
     +    ,k
     +    ,intm
     +    ,mint
     +    ,cx,cy,cz
     +    ,necl,ngcl
     +    ,q
     +    ,qx0,qy0,qz0
     +
     +    ,pptlx,pptly,pptlz
     +    ,area
     +    ,qqsintx,qqsinty,qqsintz
     +    ,qqsinmx,qqsinmy,qqsinmz
     +    ,qqdn
     +    ,sstresslet
     +    )

c       write (6,*) pptlz,aarel

        ptlx = ptlx + pptlx
        ptly = ptly + pptly
        ptlz = ptlz + pptlz

c       write (6,*) ptlx

        qdn = qdn+qqdn

        qsintx = qsintx + qqsintx
        qsinty = qsinty + qqsinty
        qsintz = qsintz + qqsintz

        qsinmx = qsinmx + qqsinmx
        qsinmy = qsinmy + qqsinmy
        qsinmz = qsinmz + qqsinmz

        Do ii=1,3
         Do jj=1,3
          stresslet(ii,jj) = stresslet(ii,jj) + sstresslet(ii,jj) 
         End Do
        End Do

        srf_area = srf_area+area

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

       dlp(i,1) = ptlx - pi4*qx0
       dlp(i,2) = ptly - pi4*qy0
       dlp(i,3) = ptlz - pi4*qz0

c      write (6,100) i,dlp(i,1),dlp(i,2),dlp(i,3),srf_area ! ,qx0,qy0,qz0

  1   Continue               ! loop over nodes

c-----
c Done
c-----

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
