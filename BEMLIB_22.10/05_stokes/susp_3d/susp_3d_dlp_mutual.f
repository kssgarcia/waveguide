      subroutine susp_3d_dlp_mutual
     +
     +   (npts
     +   ,nelm
     +   ,intm      ! integration method
     +   ,mint
     +   ,mpoly
     +   ,necl,ngcl 
     +   ,pcl
     +   ,q        ! dipole
     +
     +   ,nev  ! number of evaluation points
     +   ,pev  ! evaluation points
     +   ,dlp  ! double layer potential
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
c compute the double-layer potential
c of a scalar function q
c over a particle
c at the evaluatoin points pev
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     q(2306,3),dlp(2306,3)
      Dimension   ncl(512,100)   ! connectivity
      Dimension   pcl(2306,3)

      Dimension   pev(2306,3)    ! evaluation points

c--------------
c common blocks
c--------------

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8
c----------------------
c launch the quadrature
c----------------------

c     Do i=1,ngcl
c       write (6,*) q(i,1),q(i,2),q(i,3)
c     End Do

c-----------------
      Do i=1,nev   ! loop over evaluation points
c-----------------

c      write (6,*) " Computing the dlp at global node",i

       x0 = pev(i,1)
       y0 = pev(i,2)
       z0 = pev(i,3)

c-----------
c initialize
c-----------

       dlp(i,1) = 0.0D0
       dlp(i,2) = 0.0D0
       dlp(i,3) = 0.0D0

       srf_area = 0.0D0

       testxx = 0.0D0
       testxy = 0.0D0
       testxz = 0.0D0
       testyx = 0.0D0
       testyy = 0.0D0
       testyz = 0.0D0
       testzx = 0.0D0
       testzy = 0.0D0
       testzz = 0.0D0

c-----------------------------------------
c compile the integrals over the triangles
c-----------------------------------------

       Do k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

c---
c apply the quadrature
c---

c     write (6,*) x0,y0,z0
c     write (6,*) npts,nelm,intm,mint,mpoly,necl,ngcl

        call susp_3d_dlp_integral_mutual
     +
     +    (x0,y0,z0
     +    ,mpoly
     +    ,pcl
     +    ,k
     +    ,intm
     +    ,mint
     +    ,necl,ngcl
     +    ,q
     +    ,ptlx,ptly,ptlz
     +    ,area
     +    ,tstxx,tstxy,tstxz
     +    ,tstyx,tstyy,tstyz
     +    ,tstzx,tstzy,tstzz
     +    )

c       write (6,*) pptlz,aarel

        dlp(i,1) = dlp(i,1) + ptlx
        dlp(i,2) = dlp(i,2) + ptly
        dlp(i,3) = dlp(i,3) + ptlz

        srf_area = srf_area+area

        testxx = testxx+tstxx
        testxy = testxy+tstxy
        testxz = testxz+tstxz
        testyx = testyx+tstyx
        testyy = testyy+tstyy
        testyz = testyz+tstyz
        testzx = testzx+tstzx
        testzy = testzy+tstzy
        testzz = testzz+tstzz

       End Do ! loop over elements

c      write (6,100) i,dlp(i,1),dlp(i,2),dlp(i,3),srf_area/pi4

c      write (6,100) i,testxx,testxy,testxz
c      write (6,100) i,testyx,testyy,testyz
c      write (6,100) i,testzx,testzy,testzz

c-------------
      End Do         ! loop over evaluation points
c-------------

c-----
c done
c-----

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
