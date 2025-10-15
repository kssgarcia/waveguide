      subroutine lslp3_3d_2p
     +
     +                    (npts,nelm
     +                    ,mint,NGL
     +                    ,zeta
     +                    ,Iflow
     +                    ,slp3
     +                    )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c----------------------------------
c Computes the Laplace single-layer potential
c of a vector function zeta
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

      Dimension     p(1090,3)
      Dimension    ne(1090,7)
      Dimension   vna(1090,3)
      Dimension  zeta(1090,3)
      Dimension Iedge(1090,4)

      Dimension slp3(1090,3)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),  beta(512),gamma(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

      Parameter (tol=0.000001)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe
      common/edgepoints/Iedge
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

c----------------------
c launch the quadrature
c----------------------

      Do 1 i = 1,npts

c      write (6,*) " Computing the vector potential at point",i

       x0 = p(i,1)
       y0 = p(i,2)
       z0 = p(i,3)

c-----
c initialize
c-----

       srf_area = 0.0
       ptlx     = 0.0
       ptly     = 0.0
       ptlz     = 0.0

c---------------------
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

      If(test1.le.tol.and.test2.le.tol
     +               .and.test3.le.tol) Go to 2

c---
c test for singular integration
c---

        ii = i   ! test node

        If(   ii.eq.n(k,1).or.ii.eq.n(k,2)
     +    .or.ii.eq.n(k,3).or.ii.eq.n(k,4)
     +    .or.ii.eq.n(k,5).or.ii.eq.n(k,6)
     +    ) then

         Go to 77     ! singular triangle

        End If

c---
c interrogate the images
c
c noim: number of images
c---

      noim = Iedge(i,1)

      If(noim.ne.0) then

       Do j=2,noim+1
        ii = Iedge(i,j)   ! jth image
        If(   ii.eq.n(k,1).or.ii.eq.n(k,2)
     +    .or.ii.eq.n(k,3).or.ii.eq.n(k,4)
     +    .or.ii.eq.n(k,5).or.ii.eq.n(k,6)
     +    ) then
         Go to 77     ! singular triangle
        End If

       End Do
      End If

c----------------------------------
c non-singular integration:
c use a regular triangle quadrature
c-----------------------------------

        call Intr_quad_lslp3
     +
     +                    (x0,y0,z0
     +                    ,i,k
     +                    ,mint,Iflow
     +                    ,pptlx,pptly,pptlz
     +                    ,nelm
     +                    ,arelm
     +                    )

        ptlx = ptlx+pptlx
        ptly = ptly+pptly
        ptlz = ptlz+pptlz

        srf_area = srf_area+arelm

        Go to 2

c------------------

  77  Continue

c-----------------------
c
c Singular integration:
c
c If the point i is a vertex node,
c integrate over the flat triangle
c defined by the node
c using the polar integration rule
c
c If the point i is a mid node,
c breakup the curved triangle into four
c flat triangles
c and integrate over the flat triangles
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

       If(ii.eq.n(k,1)) then

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )


        ptlx = ptlx + pptlx
        ptly = ptly + pptly
        ptlz = ptlz + pptlz

        srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 2
c vertex node
c--------------------------------------------

        Else If(ii.eq.n(k,2)) then

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

        ptlx = ptlx + pptlx
        ptly = ptly + pptly
        ptlz = ptlz + pptlz

        srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 3
c vertex node
c--------------------------------------------

        Else If(ii.eq.n(k,3)) then

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

        ptlx = ptlx + pptlx
        ptly = ptly + pptly
        ptlz = ptlz + pptlz

        srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 4
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        Else If(ii.eq.n(k,4)) then

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )


        ptlx = ptlx + pptlx
        ptly = ptly + pptly
        ptlz = ptlz + pptlz

        srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 5
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        Else If(ii.eq.n(k,5)) then

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

        ptlx = ptlx + pptlx
        ptly = ptly + pptly
        ptlz = ptlz + pptlz

        srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 6
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        Else If(ii.eq.n(k,6)) then

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

         call Intr_lin_sing_lslp3
     +
     +                      (NGL,Iflow
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,pptlx,pptly,pptlz
     +                      ,arelm
     +                      )

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

c--------------------
c      If(i.eq.7) then
c       write (6,*)
c       write (6,*) " total surface area: ",srf_area
c       write (6,*)
c      End If
c--------------------

  1   Continue               ! loop over nodes

c---
c Done
c---

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      Return
      End

c==============================================

      subroutine intr_quad_lslp3
     +
     +                         (x0,y0,z0
     +                         ,ipoint,k
     +                         ,mint,Iflow
     +                         ,ptlx,ptly,ptlz
     +                         ,nelm
     +                         ,area
     +                         )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c--------------------------------------
c Compute the Biot-Savart integral
c for the vector potential
c over a non-singular triangle
c numbered k
c using a triangle quadrature
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1090,3)
      Dimension   ne(1090,7)
      Dimension  vna(1090,3)
      Dimension zeta(1090,3)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),beta(512),gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma
      common/geo2/vna
      common/zzetaa/zeta

      common/trq/xiq,etq,wq

      common/lgfi/method_lgf,Max1,max2,Max3
      common/lgfr/a11,a12,a21,a22,b11,b12,b21,b22,ew,area_lgf

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c---
c initialize
c---

      Iopt_int = 2   ! need the surface metric
      Iopt_lgf = 1   ! need only G 

      area = 0.0
      flow = 0.0

c---
c triangle nodes
c---

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be = beta (k)
      ga = gamma(k)

      ptlx = 0.0
      ptly = 0.0
      ptlz = 0.0

c---
c loop over integration points
c---

      Do 1 i = 1,mint

        xi  = xiq(i)
        eta = etq(i)

        call interp_lslp3
     +
     +                  (Iopt_int
     +                  ,p(i1,1),p(i1,2),p(i1,3)
     +                  ,p(i2,1),p(i2,2),p(i2,3)
     +                  ,p(i3,1),p(i3,2),p(i3,3)
     +                  ,p(i4,1),p(i4,2),p(i4,3)
     +                  ,p(i5,1),p(i5,2),p(i5,3)
     +                  ,p(i6,1),p(i6,2),p(i6,3)
     +                  ,zeta(i1,1),zeta(i1,2),zeta(i1,3)
     +                  ,zeta(i2,1),zeta(i2,2),zeta(i2,3)
     +                  ,zeta(i3,1),zeta(i3,2),zeta(i3,3)
     +                  ,zeta(i4,1),zeta(i4,2),zeta(i4,3)
     +                  ,zeta(i5,1),zeta(i5,2),zeta(i5,3)
     +                  ,zeta(i6,1),zeta(i6,2),zeta(i6,3)
     +                  ,al,be,ga
     +                  ,xi,eta
     +                  ,x,y,z
     +                  ,hs
     +                  ,zet1,zet2,zet3
     +                  )

c---
c compute the Green's function
c---

       call lgf_3d_2p
     +                (Iopt_lgf
     +                ,method_lgf
     +                ,x,y,z
     +                ,x0,y0,z0
     +                ,a11,a12,a21,a22
     +                ,b11,b12,b21,b22
     +                ,ew,area_lgf
     +                ,Max1,Max2,Max3
     +                ,G
     +                ,Gx,Gy,Gz
     +                )

c---
c apply the triangle quadrature
c---

      area = area + 0.5*hs*wq(i)

      cf = G*0.5*hs*wq(i)

      ptlx = ptlx + zet1*cf 
      ptly = ptly + zet2*cf 
      ptlz = ptlz + zet3*cf 

  1   Continue

c---
c Done
c---

 100  Format (4(1x,f10.5))

      Return
      End

c==============================================

      subroutine Intr_lin_sing_lslp3
     +
     +                        (NGL,Iflow
     +                        ,x1,y1,z1
     +                        ,x2,y2,z2
     +                        ,x3,y3,z3
     +                        ,ztx1,zty1,ztz1
     +                        ,ztx2,zty2,ztz2
     +                        ,ztx3,zty3,ztz3
     +                        ,uxel,uyel,uzel
     +                        ,area
     +                        )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c--------------------------------------------------------
c Computes the laplace single-layer potential
c for a vector density function over
c a linear (flat) triangle defined by three points 1-2-3
c
c SYMBOLS:
c -------
c
c asm: triangle area computed by numerical integration
c-------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension zz(20),ww(20)

      common/zwl/zz,ww

      common/lgfi/method_lgf,Max1,max2,Max3
      common/lgfr/a11,a12,a21,a22,b11,b12,b21,b22,ew,area_lgf

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c---
c flags
c---

      Iopt_lgf = 1    ! need G only

c---
c compute triangle area and surface metric
c---

      dx = sqrt( (x2-x1)**2+(y2-y1)**2+(z2-z1)**2 )
      dy = sqrt( (x3-x1)**2+(y3-y1)**2+(z3-z1)**2 )

      vnx = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1)
      vny = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1)
      vnz = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)

      area = 0.5*sqrt( vnx**2+vny**2+vnz**2 )
      hs   = 2.0*area

      vnx = vnx/hs
      vny = vny/hs
      vnz = vnz/hs

c---
c initialize
c---

      asm = 0.0

      uux = 0.0
      uuy = 0.0
      uuz = 0.0

c---
c Double Gaussian quadrature
c---

      Do 1 i=1,NGL

      ph    = piq*(1.0+zz(i))
      cph   = cos(ph)
      sph   = sin(ph)
      rmax  = 1.0/(cph+sph)
      rmaxh = 0.5*rmax

      bsm = 0.0

      bix = 0.0
      biy = 0.0
      biz = 0.0

      Do 2 j=1,NGL

       r = rmaxh*(1.0+zz(j))

       xi = r*cph
       et = r*sph
       zt = 1.0-xi-et

       x = x1*zt + x2*xi + x3*et
       y = y1*zt + y2*xi + y3*et
       z = z1*zt + z2*xi + z3*et

       zetx = ztx1*zt + ztx2*xi + ztx3*et
       zety = zty1*zt + zty2*xi + zty3*et
       zetz = ztz1*zt + ztz2*xi + ztz3*et

       call lgf_3d_2p
     +                (Iopt_lgf
     +                ,method_lgf
     +                ,x,y,z
     +                ,x1,y1,z1
     +                ,a11,a12,a21,a22
     +                ,b11,b12,b21,b22
     +                ,ew,area_lgf
     +                ,Max1,Max2,Max3
     +                ,G
     +                ,Gx,Gy,Gz
     +                )

         cf = r*ww(j)

         bsm = bsm+cf

         cf = cf*G
         bix = bix + zetx*cf
         biy = biy + zety*cf
         biz = biz + zetz*cf

  2     Continue

        cf = ww(i)*rmaxh

        asm = asm + bsm*cf

        uux = uux + bix*cf
        uuy = uuy + biy*cf
        uuz = uuz + biz*cf

  1   Continue

c---
c finish up the quadrature
c---

      cf = piq*hs

      asm = asm*cf

      uxel = uux*cf
      uyel = uuy*cf
      uzel = uuz*cf

c-----------------------------
c  if all went well,
c  asm should be equal to area
c
c     write (6,100) i,area,asm
c-----------------------------

 100  Format (1x,i3,2(f10.5))

      Return
      End

c=============================================

      subroutine interp_lslp3
     +
     +                      (Iopt_int
     +                      ,x1,y1,z1
     +                      ,x2,y2,z2
     +                      ,x3,y3,z3
     +                      ,x4,y4,z4
     +                      ,x5,y5,z5
     +                      ,x6,y6,z6
     +                      ,ztx1,zty1,ztz1
     +                      ,ztx2,zty2,ztz2
     +                      ,ztx3,zty3,ztz3
     +                      ,ztx4,zty4,ztz4
     +                      ,ztx5,zty5,ztz5
     +                      ,ztx6,zty6,ztz6
     +                      ,al,be,ga
     +                      ,xi,eta
     +                      ,x,y,z
     +                      ,hs
     +                      ,zetx,zety,zetz
     +                      )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c-------------------------------------------
c   Utility of the Laplace single-layer itntegrator
c
c   Interpolates over an element for the:
c
c   position vector
c   the surface metric
c   the vectorial density zeta
c  
c   Iopt_int = 1 only the position vector
c              2 position vector and rest of variables
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c prepare
c---

      alc = 1.0-al
      bec = 1.0-be
      gac = 1.0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c----------------------------
c compute the basis functions
c----------------------------

      ph2 = xi *(xi -al+eta*(al-ga)/gac)   /alc
      ph3 = eta*(eta-be+xi *(be+ga-1.0)/ga)/bec
      ph4 = xi *(1.0-xi-eta)/alalc
      ph5 = xi*eta          /gagac
      ph6 = eta*(1.0-xi-eta)/bebec
      ph1 = 1.0-ph2-ph3-ph4-ph5-ph6

c--------------------------------
c interpolate the position vector
c--------------------------------

      x = x1*ph1 + x2*ph2 + x3*ph3
     +  + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3
     +  + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3
     +  + z4*ph4 + z5*ph5 + z6*ph6

c--------------------------------------
c interpolate for the vectorial density
c--------------------------------------

      zetx = ztx1*ph1 + ztx2*ph2 +ztx3*ph3
     +      +ztx4*ph4 + ztx5*ph5 +ztx6*ph6

      zety = zty1*ph1 + zty2*ph2 +zty3*ph3
     +      +zty4*ph4 + zty5*ph5 +zty6*ph6

      zetz = ztz1*ph1 + ztz2*ph2 +ztz3*ph3
     +      +ztz4*ph4 + ztz5*ph5 +ztz6*ph6

c------------------------------------------
c compute xi derivatives of basis functions
c------------------------------------------

      dph2 =  (2.0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0)/(ga*bec)
      dph4 =  (1.0-2.0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c----------------------------------------------
c compute xi derivatives of the position vector
c----------------------------------------------

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c-------------------------------------------
c compute eta derivatives of basis functions
c-------------------------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0*eta-be+xi*(be+ga-1.0)/ga)/bec
      pph4 = -xi/alalc
      pph5 =  xi/gagac
      pph6 =  (1.0-xi-2.0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c-----------------------------------------------
c compute eta derivatives of the position vector
c-----------------------------------------------

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6
      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6
      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

c---------------------------------------
c  compute the raw normal vector and the
c  surface metric hs
c---------------------------------------

      vnxr = DyDxi * DzDet - DyDet * DzDxi
      vnyr = DzDxi * DxDet - DzDet * DxDxi
      vnzr = DxDxi * DyDet - DxDet * DyDxi

      hs = sqrt(vnxr**2+vnyr**2+vnzr**2 )

c---
c Done
c---

      Return
      End
