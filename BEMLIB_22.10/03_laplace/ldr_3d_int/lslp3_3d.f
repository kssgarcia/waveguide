      subroutine lslp3_3d
     +
     +   (npts,nelm
     +   ,mint,NGL
     +   ,zeta
     +   ,slp3
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
c Compute the single-layer potential
c of a vector function zeta
c at the nodes of a triangular grid
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension zeta(1026,3)

      Dimension slp3(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

      Parameter (tol=0.00000001)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/var1/Iflow,Ign,wall

      common/geo2/vna

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

c----------------------
c launch the quadrature
c----------------------

      Do 1 i=1,npts

c      write (6,*) " Computing the single-layer potential at point",i

       x0 = p(i,1)
       y0 = p(i,2)
       z0 = p(i,3)

c-----
c initialize
c-----

       srf_area = 0.0D0
       ptlx     = 0.0D0
       ptly     = 0.0D0
       ptlz     = 0.0D0

c----------------------
c Compile the integrals
c over the triangles
c---------------------

       Do 2 k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

       i1 = n(k,1)
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

c---
c integration will be performed only
c if zeta is nonzero
c---

      test1 = abs(zeta(i1,1))+abs(zeta(i2,1))+abs(zeta(i3,1))
     +     +  abs(zeta(i4,1))+abs(zeta(i5,1))+abs(zeta(i6,1))

      test2 = abs(zeta(i1,2))+abs(zeta(i2,2))+abs(zeta(i3,2))
     +     +  abs(zeta(i4,2))+abs(zeta(i5,2))+abs(zeta(i6,2))

      test3 = abs(zeta(i1,3))+abs(zeta(i2,3))+abs(zeta(i3,3))
     +     +  abs(zeta(i4,3))+abs(zeta(i5,3))+abs(zeta(i6,3))

      if(test1.le.tol.and.test2.le.tol
     +               .and.test3.le.tol) Go to 2

c----------------------------------
c non-singular integration:
c use a regular triangle quadrature
c-----------------------------------

        if(    i.ne.n(k,1).and.i.ne.n(k,2)
     +    .and.i.ne.n(k,3).and.i.ne.n(k,4)
     +    .and.i.ne.n(k,5).and.i.ne.n(k,6)
     +    ) then

         call lslp3_3d_integrate
     +
     +     (x0,y0,z0
     +     ,i,k
     +     ,mint
     +     ,zeta
     +     ,pptlx,pptly,pptlz
     +     ,arelm
     +     )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

         Go to 2

        end if

c-------------------------------------------
c Singular integration:
c
c If the point i is a vertex node,
c integrate over the flat triangle
c defined by the node
c using the polar integration rule
c
c If the point i is a mid-node,
c breakup the curved triangle into four
c flat triangles and integrate over the flat triangles
c using the polar integration rule
c  
c   Iopt_int = 1 only the position vector
c              2 position vector and rest of variables
c
c--------------------------------------------

c--------------------------------------------
c singular element with singularity at node 1
c this is a vertex node
c--------------------------------------------

       if(i.eq.n(k,1)) then

          x1 =    p(i1,1)
          y1 =    p(i1,2)
          z1 =    p(i1,3)
        ztx1 = zeta(i1,1)
        zty1 = zeta(i1,2)
        ztz1 = zeta(i1,3)

          x2 =    p(i2,1)
          y2 =    p(i2,2)
          z2 =    p(i2,3)
        ztx2 = zeta(i2,1)
        zty2 = zeta(i2,2)
        ztz2 = zeta(i2,3)

          x3 =    p(i3,1)
          y3 =    p(i3,2)
          z3 =    p(i3,3)
        ztx3 = zeta(i3,1)
        zty3 = zeta(i3,2)
        ztz3 = zeta(i3,3)

         call lslp3_3d_integrate_sing
     +
     +       (NGL
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,ztx1,zty1,ztz1
     +       ,ztx2,zty2,ztz2
     +       ,ztx3,zty3,ztz3
     +       ,pptlx,pptly,pptlz
     +       ,arelm
     +       )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 2
c vertex node
c--------------------------------------------

        else if(i.eq.n(k,2)) then

          x1 =    p(i2,1)
          y1 =    p(i2,2)
          z1 =    p(i2,3)
        ztx1 = zeta(i2,1)
        zty1 = zeta(i2,2)
        ztz1 = zeta(i2,3)

          x2 =    p(i3,1)
          y2 =    p(i3,2)
          z2 =    p(i3,3)
        ztx2 = zeta(i3,1)
        zty2 = zeta(i3,2)
        ztz2 = zeta(i3,3)

          x3 =    p(i1,1)
          y3 =    p(i1,2)
          z3 =    p(i1,3)
        ztx3 = zeta(i1,1)
        zty3 = zeta(i1,2)
        ztz3 = zeta(i1,3)

         call lslp3_3d_integrate_sing
     +
     +       (NGL
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,ztx1,zty1,ztz1
     +       ,ztx2,zty2,ztz2
     +       ,ztx3,zty3,ztz3
     +       ,pptlx,pptly,pptlz
     +       ,arelm
     +       )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 3
c vertex node
c--------------------------------------------

        else if(i.eq.n(k,3)) then

          x1 =    p(i3,1)
          y1 =    p(i3,2)
          z1 =    p(i3,3)
        ztx1 = zeta(i3,1)
        zty1 = zeta(i3,2)
        ztz1 = zeta(i3,3)

          x2 =    p(i1,1)
          y2 =    p(i1,2)
          z2 =    p(i1,3)
        ztx2 = zeta(i1,1)
        zty2 = zeta(i1,2)
        ztz2 = zeta(i1,3)

          x3 =    p(i2,1)
          y3 =    p(i2,2)
          z3 =    p(i2,3)
        ztx3 = zeta(i2,1)
        zty3 = zeta(i2,2)
        ztz3 = zeta(i2,3)

         call lslp3_3d_integrate_sing
     +
     +      (NGL
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     +      ,ztx1,zty1,ztz1
     +      ,ztx2,zty2,ztz2
     +      ,ztx3,zty3,ztz3
     +      ,pptlx,pptly,pptlz
     +      ,arelm
     +      )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 4
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.n(k,4)) then

          x1 =    p(i4,1)
          y1 =    p(i4,2)
          z1 =    p(i4,3)
        ztx1 = zeta(i4,1)
        zty1 = zeta(i4,2)
        ztz1 = zeta(i4,3)

          x2 =    p(i6,1)
          y2 =    p(i6,2)
          z2 =    p(i6,3)
        ztx2 = zeta(i6,1)
        zty2 = zeta(i6,2)
        ztz2 = zeta(i6,3)

          x3 =    p(i1,1)
          y3 =    p(i1,2)
          z3 =    p(i1,3)
        ztx3 = zeta(i1,1)
        zty3 = zeta(i1,2)
        ztz3 = zeta(i1,3)

         call lslp3_3d_integrate_sing
     +
     +      (NGL
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     +      ,ztx1,zty1,ztz1
     +      ,ztx2,zty2,ztz2
     +      ,ztx3,zty3,ztz3
     +      ,pptlx,pptly,pptlz
     +      ,arelm
     +      )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i4,1)
          y1 =    p(i4,2)
          z1 =    p(i4,3)
        ztx1 = zeta(i4,1)
        zty1 = zeta(i4,2)
        ztz1 = zeta(i4,3)

          x2 =    p(i3,1)
          y2 =    p(i3,2)
          z2 =    p(i3,3)
        ztx2 = zeta(i3,1)
        zty2 = zeta(i3,2)
        ztz2 = zeta(i3,3)

          x3 =    p(i6,1)
          y3 =    p(i6,2)
          z3 =    p(i6,3)
        ztx3 = zeta(i6,1)
        zty3 = zeta(i6,2)
        ztz3 = zeta(i6,3)

         call lslp3_3d_integrate_sing
     +
     +       (NGL
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,ztx1,zty1,ztz1
     +       ,ztx2,zty2,ztz2
     +       ,ztx3,zty3,ztz3
     +       ,pptlx,pptly,pptlz
     +       ,arelm
     +       )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i4,1)
          y1 =    p(i4,2)
          z1 =    p(i4,3)
        ztx1 = zeta(i4,1)
        zty1 = zeta(i4,2)
        ztz1 = zeta(i4,3)

          x2 =    p(i5,1)
          y2 =    p(i5,2)
          z2 =    p(i5,3)
        ztx2 = zeta(i5,1)
        zty2 = zeta(i5,2)
        ztz2 = zeta(i5,3)

          x3 =    p(i3,1)
          y3 =    p(i3,2)
          z3 =    p(i3,3)
        ztx3 = zeta(i3,1)
        zty3 = zeta(i3,2)
        ztz3 = zeta(i3,3)

         call lslp3_3d_integrate_sing
     +
     +       (NGL
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,ztx1,zty1,ztz1
     +       ,ztx2,zty2,ztz2
     +       ,ztx3,zty3,ztz3
     +       ,pptlx,pptly,pptlz
     +       ,arelm
     +       )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i4,1)
          y1 =    p(i4,2)
          z1 =    p(i4,3)
        ztx1 = zeta(i4,1)
        zty1 = zeta(i4,2)
        ztz1 = zeta(i4,3)

          x2 =    p(i2,1)
          y2 =    p(i2,2)
          z2 =    p(i2,3)
        ztx2 = zeta(i2,1)
        zty2 = zeta(i2,2)
        ztz2 = zeta(i2,3)

          x3 =    p(i5,1)
          y3 =    p(i5,2)
          z3 =    p(i5,3)
        ztx3 = zeta(i5,1)
        zty3 = zeta(i5,2)
        ztz3 = zeta(i5,3)

         call lslp3_3d_integrate_sing
     +
     +       (NGL
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,ztx1,zty1,ztz1
     +       ,ztx2,zty2,ztz2
     +       ,ztx3,zty3,ztz3
     +       ,pptlx,pptly,pptlz
     +       ,arelm
     +       )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 5
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.n(k,5)) then

          x1 =    p(i5,1)
          y1 =    p(i5,2)
          z1 =    p(i5,3)
        ztx1 = zeta(i5,1)
        zty1 = zeta(i5,2)
        ztz1 = zeta(i5,3)

          x2 =    p(i4,1)
          y2 =    p(i4,2)
          z2 =    p(i4,3)
        ztx2 = zeta(i4,1)
        zty2 = zeta(i4,2)
        ztz2 = zeta(i4,3)

          x3 =    p(i2,1)
          y3 =    p(i2,2)
          z3 =    p(i2,3)
        ztx3 = zeta(i2,1)
        zty3 = zeta(i2,2)
        ztz3 = zeta(i2,3)

         call lslp3_3d_integrate_sing
     +
     +      (NGL
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     +      ,ztx1,zty1,ztz1
     +      ,ztx2,zty2,ztz2
     +      ,ztx3,zty3,ztz3
     +      ,pptlx,pptly,pptlz
     +      ,arelm
     +      )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i5,1)
          y1 =    p(i5,2)
          z1 =    p(i5,3)
        ztx1 = zeta(i5,1)
        zty1 = zeta(i5,2)
        ztz1 = zeta(i5,3)

          x2 =    p(i1,1)
          y2 =    p(i1,2)
          z2 =    p(i1,3)
        ztx2 = zeta(i1,1)
        zty2 = zeta(i1,2)
        ztz2 = zeta(i1,3)

          x3 =    p(i4,1)
          y3 =    p(i4,2)
          z3 =    p(i4,3)
        ztx3 = zeta(i4,1)
        zty3 = zeta(i4,2)
        ztz3 = zeta(i4,3)

         call lslp3_3d_integrate_sing
     +
     +      (NGL
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     +      ,ztx1,zty1,ztz1
     +      ,ztx2,zty2,ztz2
     +      ,ztx3,zty3,ztz3
     +      ,pptlx,pptly,pptlz
     +      ,arelm
     +      )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i5,1)
          y1 =    p(i5,2)
          z1 =    p(i5,3)
        ztx1 = zeta(i5,1)
        zty1 = zeta(i5,2)
        ztz1 = zeta(i5,3)

          x2 =    p(i6,1)
          y2 =    p(i6,2)
          z2 =    p(i6,3)
        ztx2 = zeta(i6,1)
        zty2 = zeta(i6,2)
        ztz2 = zeta(i6,3)

          x3 =    p(i1,1)
          y3 =    p(i1,2)
          z3 =    p(i1,3)
        ztx3 = zeta(i1,1)
        zty3 = zeta(i1,2)
        ztz3 = zeta(i1,3)

         call lslp3_3d_integrate_sing
     +
     +        (NGL
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,ztx1,zty1,ztz1
     +        ,ztx2,zty2,ztz2
     +        ,ztx3,zty3,ztz3
     +        ,pptlx,pptly,pptlz
     +        ,arelm
     +        )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i5,1)
          y1 =    p(i5,2)
          z1 =    p(i5,3)
        ztx1 = zeta(i5,1)
        zty1 = zeta(i5,2)
        ztz1 = zeta(i5,3)

          x2 =    p(i3,1)
          y2 =    p(i3,2)
          z2 =    p(i3,3)
        ztx2 = zeta(i3,1)
        zty2 = zeta(i3,2)
        ztz2 = zeta(i3,3)

          x3 =    p(i6,1)
          y3 =    p(i6,2)
          z3 =    p(i6,3)
        ztx3 = zeta(i6,1)
        zty3 = zeta(i6,2)
        ztz3 = zeta(i6,3)

         call lslp3_3d_integrate_sing
     +
     +      (NGL
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     +      ,ztx1,zty1,ztz1
     +      ,ztx2,zty2,ztz2
     +      ,ztx3,zty3,ztz3
     +      ,pptlx,pptly,pptlz
     +      ,arelm
     +      )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 6
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.n(k,6)) then

          x1 =    p(i6,1)
          y1 =    p(i6,2)
          z1 =    p(i6,3)
        ztx1 = zeta(i6,1)
        zty1 = zeta(i6,2)
        ztz1 = zeta(i6,3)

          x2 =    p(i1,1)
          y2 =    p(i1,2)
          z2 =    p(i1,3)
        ztx2 = zeta(i1,1)
        zty2 = zeta(i1,2)
        ztz2 = zeta(i1,3)

          x3 =    p(i4,1)
          y3 =    p(i4,2)
          z3 =    p(i4,3)
        ztx3 = zeta(i4,1)
        zty3 = zeta(i4,2)
        ztz3 = zeta(i4,3)

         call lslp3_3d_integrate_sing
     +
     +     (NGL
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     +     ,ztx1,zty1,ztz1
     +     ,ztx2,zty2,ztz2
     +     ,ztx3,zty3,ztz3
     +     ,pptlx,pptly,pptlz
     +     ,arelm
     +     )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

        srf_area = srf_area+arelm

          x1 =    p(i6,1)
          y1 =    p(i6,2)
          z1 =    p(i6,3)
        ztx1 = zeta(i6,1)
        zty1 = zeta(i6,2)
        ztz1 = zeta(i6,3)

          x2 =    p(i4,1)
          y2 =    p(i4,2)
          z2 =    p(i4,3)
        ztx2 = zeta(i4,1)
        zty2 = zeta(i4,2)
        ztz2 = zeta(i4,3)

          x3 =    p(i2,1)
          y3 =    p(i2,2)
          z3 =    p(i2,3)
        ztx3 = zeta(i2,1)
        zty3 = zeta(i2,2)
        ztz3 = zeta(i2,3)

         call lslp3_3d_integrate_sing
     +
     +        (NGL
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,ztx1,zty1,ztz1
     +        ,ztx2,zty2,ztz2
     +        ,ztx3,zty3,ztz3
     +        ,pptlx,pptly,pptlz
     +        ,arelm
     +        )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i6,1)
          y1 =    p(i6,2)
          z1 =    p(i6,3)
        ztx1 = zeta(i6,1)
        zty1 = zeta(i6,2)
        ztz1 = zeta(i6,3)

          x2 =    p(i2,1)
          y2 =    p(i2,2)
          z2 =    p(i2,3)
        ztx2 = zeta(i2,1)
        zty2 = zeta(i2,2)
        ztz2 = zeta(i2,3)

          x3 =    p(i5,1)
          y3 =    p(i5,2)
          z3 =    p(i5,3)
        ztx3 = zeta(i5,1)
        zty3 = zeta(i5,2)
        ztz3 = zeta(i5,3)

         call lslp3_3d_integrate_sing
     +
     +        (NGL
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,ztx1,zty1,ztz1
     +        ,ztx2,zty2,ztz2
     +        ,ztx3,zty3,ztz3
     +        ,pptlx,pptly,pptlz
     +        ,arelm
     +        )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

          x1 =    p(i6,1)
          y1 =    p(i6,2)
          z1 =    p(i6,3)
        ztx1 = zeta(i6,1)
        zty1 = zeta(i6,2)
        ztz1 = zeta(i6,3)

          x2 =    p(i5,1)
          y2 =    p(i5,2)
          z2 =    p(i5,3)
        ztx2 = zeta(i5,1)
        zty2 = zeta(i5,2)
        ztz2 = zeta(i5,3)

          x3 =    p(i3,1)
          y3 =    p(i3,2)
          z3 =    p(i3,3)
        ztx3 = zeta(i3,1)
        zty3 = zeta(i3,2)
        ztz3 = zeta(i3,3)

         call lslp3_3d_integrate_sing
     +
     +        (NGL
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,ztx1,zty1,ztz1
     +        ,ztx2,zty2,ztz2
     +        ,ztx3,zty3,ztz3
     +        ,pptlx,pptly,pptlz
     +        ,arelm
     +        )

         ptlx = ptlx + pptlx
         ptly = ptly + pptly
         ptlz = ptlz + pptlz

         srf_area = srf_area+arelm

c------------
       End If    ! done integrating over a singular element
c------------

  2   Continue

       slp3(i,1) = ptlx
       slp3(i,2) = ptly
       slp3(i,3) = ptlz

c---------------------
c      If(i.eq.7) then
c       write (6,*)
c       write (6,*) " total surface area computed in lslp3: ",srf_area
c       write (6,*)
c      End If
c---------------------

  1   Continue               ! loop over nodes

c---
c Done
c---

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      Return
      End
