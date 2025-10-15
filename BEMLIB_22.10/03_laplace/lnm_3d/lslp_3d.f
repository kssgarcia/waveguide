      subroutine lslp_3d
     +
     +   (npts,nelm
     +   ,mint,NGL
     +   ,f
     +   ,Iflow
     +   ,slp
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
c Compute the single-layer potential
c of the scalar function f
c at the nodes of a triangular grid
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension    f(1026)
      Dimension  slp(1026)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

      Parameter (tol = 0.0000001)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

c----------------------
c launch the quadrature
c----------------------

c---
      Do i=1,npts  ! loop over nodes
c---

c      write (6,*) " Computing the single-layer potential at point",i

       x0 = p(i,1)
       y0 = p(i,2)
       z0 = p(i,3)

c-----------
c initialize
c-----------

       srf_area = 0.0D0
       ptl      = 0.0D0

c----------------------
c Compile the integrals
c over the triangles
c---------------------

       Do k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

       i1 = n(k,1)
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

c---
c integration will be performed only
c only if f is nonzero at least at one node
c---

      test = dabs(f(i1))+dabs(f(i2))+dabs(f(i3))
     +     + dabs(f(i4))+dabs(f(i5))+dabs(f(i6))

      if(test.gt.tol) then   ! if test= 0 skip the integration
                             !            over this triangle

c-------------------------------------------
c singular integration:
c --------------------
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
c non-singular integration:
c ------------------------
c
c use a regular triangle quadrature
c--------------------------------------------

c--------------------------------------------
c singular element with singularity at node 1
c this is a vertex node
c--------------------------------------------

       if(i.eq.i1) then

          x1 = p(i1,1)
          y1 = p(i1,2)
          z1 = p(i1,3)
          f1 = f(i1)

          x2 = p(i2,1)
          y2 = p(i2,2)
          z2 = p(i2,3)
          f2 = f(i2)

          x3 = p(i3,1)
          y3 = p(i3,2)
          z3 = p(i3,3)
          f3 = f(i3)

         call lslp_3d_integral_sing
     +
     +      (NGL
     +      ,Iflow
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     +      ,f1,f2,f3
     +      ,pptl
     +      ,arelm
     +      )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 2
c vertex node
c--------------------------------------------

        else if(i.eq.i2) then

          x1 = p(i2,1)
          y1 = p(i2,2)
          z1 = p(i2,3)
          f1 = f(i2)

          x2 = p(i3,1)
          y2 = p(i3,2)
          z2 = p(i3,3)
          f2 = f(i3)

          x3 = p(i1,1)
          y3 = p(i1,2)
          z3 = p(i1,3)
          f3 = f(i1)

          call lslp_3d_integral_sing
     +
     +        (NGL
     +        ,Iflow
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,f1,f2,f3
     +        ,pptl
     +        ,arelm
     +        )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 3
c vertex node
c--------------------------------------------

        else if(i.eq.i3) then

          x1 = p(i3,1)
          y1 = p(i3,2)
          z1 = p(i3,3)
          f1 = f(i3)

          x2 = p(i1,1)
          y2 = p(i1,2)
          z2 = p(i1,3)
          f2 = f(i1)

          x3 = p(i2,1)
          y3 = p(i2,2)
          z3 = p(i2,3)
          f3 = f(i2)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 4
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.i4) then

          x1 = p(i4,1)
          y1 = p(i4,2)
          z1 = p(i4,3)
          f1 = f(i4)

          x2 = p(i6,1)
          y2 = p(i6,2)
          z2 = p(i6,3)
          f2 = f(i6)

          x3 = p(i1,1)
          y3 = p(i1,2)
          z3 = p(i1,3)
          f3 = f(i1)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i4,1)
          y1 = p(i4,2)
          z1 = p(i4,3)
          f1 = f(i4)

          x2 = p(i3,1)
          y2 = p(i3,2)
          z2 = p(i3,3)
          f2 = f(i3)

          x3 = p(i6,1)
          y3 = p(i6,2)
          z3 = p(i6,3)
          f3 = f(i6)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i4,1)
          y1 = p(i4,2)
          z1 = p(i4,3)
          f1 = f(i4)

          x2 = p(i5,1)
          y2 = p(i5,2)
          z2 = p(i5,3)
          f2 = f(i5)

          x3 = p(i3,1)
          y3 = p(i3,2)
          z3 = p(i3,3)
          f3 = f(i3)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i4,1)
          y1 = p(i4,2)
          z1 = p(i4,3)
          f1 = f(i4)

          x2 = p(i2,1)
          y2 = p(i2,2)
          z2 = p(i2,3)
          f2 = f(i2)

          x3 = p(i5,1)
          y3 = p(i5,2)
          z3 = p(i5,3)
          f3 = f(i5)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 5
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.i5) then

          x1 = p(i5,1)
          y1 = p(i5,2)
          z1 = p(i5,3)
          f1 = f(i5)

          x2 = p(i4,1)
          y2 = p(i4,2)
          z2 = p(i4,3)
          f2 = f(i4)

          x3 = p(i2,1)
          y3 = p(i2,2)
          z3 = p(i2,3)
          f3 = f(i2)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i5,1)
          y1 = p(i5,2)
          z1 = p(i5,3)
          f1 = f(i5)

          x2 = p(i1,1)
          y2 = p(i1,2)
          z2 = p(i1,3)
          f2 = f(i1)

          x3 = p(i4,1)
          y3 = p(i4,2)
          z3 = p(i4,3)
          f3 = f(i4)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i5,1)
          y1 = p(i5,2)
          z1 = p(i5,3)
          f1 = f(i5)

          x2 = p(i6,1)
          y2 = p(i6,2)
          z2 = p(i6,3)
          f2 = f(i6)

          x3 = p(i1,1)
          y3 = p(i1,2)
          z3 = p(i1,3)
          f3 = f(i1)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i5,1)
          y1 = p(i5,2)
          z1 = p(i5,3)
          f1 = f(i5)

          x2 = p(i3,1)
          y2 = p(i3,2)
          z2 = p(i3,3)
          f2 = f(i3)

          x3 = p(i6,1)
          y3 = p(i6,2)
          z3 = p(i6,3)
          f3 = f(i6)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

c--------------------------------------------
c singular element with singularity at node 6
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.i6) then

          x1 = p(i6,1)
          y1 = p(i6,2)
          z1 = p(i6,3)
          f1 = f(i6)

          x2 = p(i1,1)
          y2 = p(i1,2)
          z2 = p(i1,3)
          f2 = f(i1)

          x3 = p(i4,1)
          y3 = p(i4,2)
          z3 = p(i4,3)
          f3 = f(i4)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

        srf_area = srf_area+arelm

          x1 = p(i6,1)
          y1 = p(i6,2)
          z1 = p(i6,3)
          f1 = f(i6)

          x2 = p(i4,1)
          y2 = p(i4,2)
          z2 = p(i4,3)
          f2 = f(i4)

          x3 = p(i2,1)
          y3 = p(i2,2)
          z3 = p(i2,3)
          f3 = f(i2)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i6,1)
          y1 = p(i6,2)
          z1 = p(i6,3)
          f1 = f(i6)

          x2 = p(i2,1)
          y2 = p(i2,2)
          z2 = p(i2,3)
          f2 = f(i2)

          x3 = p(i5,1)
          y3 = p(i5,2)
          z3 = p(i5,3)
          f3 = f(i5)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

          x1 = p(i6,1)
          y1 = p(i6,2)
          z1 = p(i6,3)
          f1 = f(i6)

          x2 = p(i5,1)
          y2 = p(i5,2)
          z2 = p(i5,3)
          f2 = f(i5)

          x3 = p(i3,1)
          y3 = p(i3,2)
          z3 = p(i3,3)
          f3 = f(i3)

          call lslp_3d_integral_sing
     +
     +       (NGL
     +       ,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,f1,f2,f3
     +       ,pptl
     +       ,arelm
     +       )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

c--------------------
c regular integration
c---------------------

         else

         call lslp_3d_integral
     +
     +    (x0,y0,z0
     +    ,k
     +    ,mint
     +    ,Iflow
     +    ,f
     +    ,pptl
     +    ,arelm
     +    )

         ptl = ptl + pptl

         srf_area = srf_area+arelm

c------------
       end if    ! done integrating over an element
c------------

c----
       end if ! of Iskip
      End Do  ! loop over triangles
c----

       slp(i) = ptl

c      write (6,*) "slp integrator activated",ptl

c---------------------
c       write (6,*) " total surface area computed in lslp: ",srf_area
c---------------------

c----
      End Do        ! loop over nodes
c----

c---
c done
c---

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
