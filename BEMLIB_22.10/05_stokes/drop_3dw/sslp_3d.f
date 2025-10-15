      subroutine sslp_3d
     +
     +  (nelm
     +  ,npts
     +  ,mint
     +  ,NGL
     +  ,slp
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------
c Compute the single-layer potential
c at the nodes of a triangular grid
c
c SYMBOLS:
c -------
c
c slp(k,i): ith component of slp at
c           the kth node
c
c srtn:  surface tension
c Dsrtn: surface gradient of the surface tension
c------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension crvm(1026)

      Dimension  nvel(1026)
      Dimension  srtn(1026)
      Dimension Dsrtn(1026,3)
      Dimension    Df(1026,3)
      Dimension   slp(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)
      Dimension  Dfel(512,4)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo2/vna
      common/geo3/crvm

      common/tension/srtn

      common/var/shrt,wall

      common/veloc1/nvelt,nvel

      common/dff/Df
      common/dffel/Dfel

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------
c  Go to 99    ! to bypass the computation of the slp
c----------

c-----------------------------------------------
c  Compute Df averaged over each element
c
c  Dfel(k,i) is the ith component of the
c            average value of df over element k
c            for i=1,2,3
c
c  Dfel(k,4) contains the magnitude of Dfel
c
c----
c
c     call dfel_3d (npts,nelm)
c
c     write (6,*)
c     write (6,*) " Average of df over the elements"
c     write (6,*) " -------------------------------"
c     write (6,*)
c
c     Do k=1,nelm
c       write (6,100) k,(dfel(k,i),i=1,4)
c     End Do
c
c-----------------------------------------------

c-----------------------
c compute Df at the nodes
c------------------------

c-----
c mean curvature of the interface
c at the nodes
c-----

      call crvm_3d (nelm,npts)

c---
c surface gradient
c of the surface tension
c at the nodes: Dsrtn
c---

      call sgrad_3d
     +
     +  (npts
     +  ,nelm
     +  ,srtn
     +  ,Dsrtn
     +  )

c-----------
c compute Df
c-----------

c     write (6,*)
c     write (6,*)  " drop_3d_slp: Df, crvm:"
c     write (6,*)

      Do i=1,npts

       cf = 2.0D0*crvm(i)*srtn(i)

       Df(i,1) = cf*vna(i,1) - Dsrtn(i,1)
       Df(i,2) = cf*vna(i,2) - Dsrtn(i,2)
       Df(i,3) = cf*vna(i,3) - Dsrtn(i,3)

c      prj = vna(i,1)*Df(i,1)+vna(i,2)*Df(i,2)
c    +                      + vna(i,3)*Df(i,3)
c
c      write (6,100) i,srtn(i),Df(i,1),Df(i,2),Df(i,3),crvm(i),prj
c    +                ,Dsrtn(i,1),Dsrtn(i,2),Dsrtn(i,3)

      End Do

c-------------------------------------
c loop over velocity evaluation points
c-------------------------------------

      Do 1 node=1,nvelt

       i = nvel(node)  ! i is the global node label

c      write (6,*) " drop_3d_slp: computing the slp at point ",i

       x0 = p(i,1)
       y0 = p(i,2)
       z0 = p(i,3)

       us = 0.0D0
       vs = 0.0D0
       ws = 0.0D0

       Do 2 k=1,nelm     ! run over elements

        uxel = 0.0D0
        uyel = 0.0D0
        uzel = 0.0D0

c----------------------------------
c non-singular element:
c use the Gauss triangle quadrature
c----------------------------------

        if(    i.ne.n(k,1).and.i.ne.n(k,2)
     +    .and.i.ne.n(k,3).and.i.ne.n(k,4)
     +    .and.i.ne.n(k,5).and.i.ne.n(k,6)
     +    ) then

        call sslp_3d_integral
     +
     +   (x0,y0,z0
     +   ,k
     +   ,mint
     +   ,uxel,uyel,uzel
     +   )

        Go to 3                ! Do another element

        end if

c-------------------------------------------
c Singular element:
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
c--------------------------------------------

      i1 = n(k,1)  ! global node label
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

c--------------------------------------------
c vertex node
c singular element with singularity at node 1
c--------------------------------------------

       if(i.eq.n(k,1)) then

          x1 =  p(i1,1)
          y1 =  p(i1,2)
          z1 =  p(i1,3)
         fx1 = Df(i1,1)
         fy1 = Df(i1,2)
         fz1 = Df(i1,3)

          x2 =  p(i2,1)
          y2 =  p(i2,2)
          z2 =  p(i2,3)
         fx2 = Df(i2,1)
         fy2 = Df(i2,2)
         fz2 = Df(i2,3)

          x3 =  p(i3,1)
          y3 =  p(i3,2)
          z3 =  p(i3,3)
         fx3 = Df(i3,1)
         fy3 = Df(i3,2)
         fz3 = Df(i3,3)

         call sslp_3d_integral_sing
     +
     +    (NGL
     +
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +
     +    ,fx1,fy1,fz1
     +    ,fx2,fy2,fz2
     +    ,fx3,fy3,fz3
     +
     +    ,uxel,uyel,uzel
     +    )

c--------------------------------------------
c vertex node
c singular element with singularity at node 2
c--------------------------------------------

        else if(i.eq.n(k,2)) then

          x1 =  p(i2,1)
          y1 =  p(i2,2)
          z1 =  p(i2,3)
         fx1 = Df(i2,1)
         fy1 = Df(i2,2)
         fz1 = Df(i2,3)

          x2 =  p(i3,1)
          y2 =  p(i3,2)
          z2 =  p(i3,3)
         fx2 = Df(i3,1)
         fy2 = Df(i3,2)
         fz2 = Df(i3,3)

          x3 =  p(i1,1)
          y3 =  p(i1,2)
          z3 =  p(i1,3)
         fx3 = Df(i1,1)
         fy3 = Df(i1,2)
         fz3 = Df(i1,3)

         call sslp_3d_integral_sing
     +
     +    (NGL
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +
     +    ,fx1,fy1,fz1
     +    ,fx2,fy2,fz2
     +    ,fx3,fy3,fz3
     +
     +    ,uxel,uyel,uzel
     +    )

c--------------------------------------------
c vertex node
c singular element with singularity at node 3
c--------------------------------------------

        else if(i.eq.n(k,3)) then

          x1 =  p(i3,1)
          y1 =  p(i3,2)
          z1 =  p(i3,3)
         fx1 = Df(i3,1)
         fy1 = Df(i3,2)
         fz1 = Df(i3,3)

          x2 =  p(i1,1)
          y2 =  p(i1,2)
          z2 =  p(i1,3)
         fx2 = Df(i1,1)
         fy2 = Df(i1,2)
         fz2 = Df(i1,3)

          x3 =  p(i2,1)
          y3 =  p(i2,2)
          z3 =  p(i2,3)
         fx3 = Df(i2,1)
         fy3 = Df(i2,2)
         fz3 = Df(i2,3)

         call sslp_3d_integral_sing
     +
     +    (NGL
     +
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +
     +    ,fx1,fy1,fz1
     +    ,fx2,fy2,fz2
     +    ,fx3,fy3,fz3
     +
     +    ,uxel,uyel,uzel
     +    )

c--------------------------------------------
c edge node
c singular element with singularity at node 4
c integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.n(k,4)) then

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)
         fy1 = Df(i4,2)
         fz1 = Df(i4,3)

          x2 =  p(i6,1)
          y2 =  p(i6,2)
          z2 =  p(i6,3)
         fx2 = Df(i6,1)
         fy2 = Df(i6,2)
         fz2 = Df(i6,3)

          x3 =  p(i1,1)
          y3 =  p(i1,2)
          z3 =  p(i1,3)
         fx3 = Df(i1,1)
         fy3 = Df(i1,2)
         fz3 = Df(i1,3)

         call sslp_3d_integral_sing
     +
     +    (NGL
     + 
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     + 
     +    ,fx1,fy1,fz1
     +    ,fx2,fy2,fz2
     +    ,fx3,fy3,fz3
     + 
     +    ,uxel,uyel,uzel
     +    )

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)
         fy1 = Df(i4,2)
         fz1 = Df(i4,3)

          x2 =  p(i3,1)
          y2 =  p(i3,2)
          z2 =  p(i3,3)
         fx2 = Df(i3,1)
         fy2 = Df(i3,2)
         fz2 = Df(i3,3)

          x3 =  p(i6,1)
          y3 =  p(i6,2)
          z3 =  p(i6,3)
         fx3 = Df(i6,1)
         fy3 = Df(i6,2)
         fz3 = Df(i6,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)
         fy1 = Df(i4,2)
         fz1 = Df(i4,3)

          x2 =  p(i5,1)
          y2 =  p(i5,2)
          z2 =  p(i5,3)
         fx2 = Df(i5,1)
         fy2 = Df(i5,2)
         fz2 = Df(i5,3)

          x3 =  p(i3,1)
          y3 =  p(i3,2)
          z3 =  p(i3,3)
         fx3 = Df(i3,1)
         fy3 = Df(i3,2)
         fz3 = Df(i3,3)

         call sslp_3d_integral_sing
     +     
     +      (NGL
     + 
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     + 
     +      ,fx1,fy1,fz1
     +      ,fx2,fy2,fz2
     +      ,fx3,fy3,fz3
     + 
     +      ,uxel,uyel,uzel
     +      )

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)
         fy1 = Df(i4,2)
         fz1 = Df(i4,3)

          x2 =  p(i2,1)
          y2 =  p(i2,2)
          z2 =  p(i2,3)
         fx2 = Df(i2,1)
         fy2 = Df(i2,2)
         fz2 = Df(i2,3)

          x3 =  p(i5,1)
          y3 =  p(i5,2)
          z3 =  p(i5,3)
         fx3 = Df(i5,1)
         fy3 = Df(i5,2)
         fz3 = Df(i5,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

c--------------------------------------------
c edge node
c singular element with singularity at node 5
c integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.n(k,5)) then

          x1 =  p(i5,1)
          y1 =  p(i5,2)
          z1 =  p(i5,3)
         fx1 = Df(i5,1)
         fy1 = Df(i5,2)
         fz1 = Df(i5,3)

          x2 =  p(i4,1)
          y2 =  p(i4,2)
          z2 =  p(i4,3)
         fx2 = Df(i4,1)
         fy2 = Df(i4,2)
         fz2 = Df(i4,3)

          x3 =  p(i2,1)
          y3 =  p(i2,2)
          z3 =  p(i2,3)
         fx3 = Df(i2,1)
         fy3 = Df(i2,2)
         fz3 = Df(i2,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

          x1 = p(i5,1)
          y1 = p(i5,2)
          z1 = p(i5,3)
         fx1 = Df(i5,1)
         fy1 = Df(i5,2)
         fz1 = Df(i5,3)

          x2 = p(i1,1)
          y2 = p(i1,2)
          z2 = p(i1,3)
         fx2 = Df(i1,1)
         fy2 = Df(i1,2)
         fz2 = Df(i1,3)

          x3 =  p(i4,1)
          y3 =  p(i4,2)
          z3 =  p(i4,3)
         fx3 = Df(i4,1)
         fy3 = Df(i4,2)
         fz3 = Df(i4,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

          x1 =  p(i5,1)
          y1 =  p(i5,2)
          z1 =  p(i5,3)
         fx1 = Df(i5,1)
         fy1 = Df(i5,2)
         fz1 = Df(i5,3)

          x2 =  p(i6,1)
          y2 =  p(i6,2)
          z2 =  p(i6,3)
         fx2 = Df(i6,1)
         fy2 = Df(i6,2)
         fz2 = Df(i6,3)

          x3 =  p(i1,1)
          y3 =  p(i1,2)
          z3 =  p(i1,3)
         fx3 = Df(i1,1)
         fy3 = Df(i1,2)
         fz3 = Df(i1,3)

         call sslp_3d_integral_sing
     +
     +      (NGL
     + 
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     + 
     +      ,fx1,fy1,fz1
     +      ,fx2,fy2,fz2
     +      ,fx3,fy3,fz3
     + 
     +      ,uxel,uyel,uzel
     +      )

          x1 =  p(i5,1)
          y1 =  p(i5,2)
          z1 =  p(i5,3)
         fx1 = Df(i5,1)
         fy1 = Df(i5,2)
         fz1 = Df(i5,3)

          x2 =  p(i3,1)
          y2 =  p(i3,2)
          z2 =  p(i3,3)
         fx2 = Df(i3,1)
         fy2 = Df(i3,2)
         fz2 = Df(i3,3)

          x3 =  p(i6,1)
          y3 =  p(i6,2)
          z3 =  p(i6,3)
         fx3 = Df(i6,1)
         fy3 = Df(i6,2)
         fz3 = Df(i6,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

c--------------------------------------------
c edge node
c singular element with singularity at node 6
c integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.n(k,6)) then

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)
         fy1 = Df(i6,2)
         fz1 = Df(i6,3)

          x2 =  p(i1,1)
          y2 =  p(i1,2)
          z2 =  p(i1,3)
         fx2 = Df(i1,1)
         fy2 = Df(i1,2)
         fz2 = Df(i1,3)

          x3 =  p(i4,1)
          y3 =  p(i4,2)
          z3 =  p(i4,3)
         fx3 = Df(i4,1)
         fy3 = Df(i4,2)
         fz3 = Df(i4,3)

         call sslp_3d_integral_sing
     +
     +      (NGL
     + 
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     + 
     +      ,fx1,fy1,fz1
     +      ,fx2,fy2,fz2
     +      ,fx3,fy3,fz3
     + 
     +      ,uxel,uyel,uzel
     +      )

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)
         fy1 = Df(i6,2)
         fz1 = Df(i6,3)

          x2 =  p(i4,1)
          y2 =  p(i4,2)
          z2 =  p(i4,3)
         fx2 = Df(i4,1)
         fy2 = Df(i4,2)
         fz2 = Df(i4,3)

          x3 =  p(i2,1)
          y3 =  p(i2,2)
          z3 =  p(i2,3)
         fx3 = Df(i2,1)
         fy3 = Df(i2,2)
         fz3 = Df(i2,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)
         fy1 = Df(i6,2)
         fz1 = Df(i6,3)

          x2 =  p(i2,1)
          y2 =  p(i2,2)
          z2 =  p(i2,3)
         fx2 = Df(i2,1)
         fy2 = Df(i2,2)
         fz2 = Df(i2,3)

          x3 =  p(i5,1)
          y3 =  p(i5,2)
          z3 =  p(i5,3)
         fx3 = Df(i5,1)
         fy3 = Df(i5,2)
         fz3 = Df(i5,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)
         fy1 = Df(i6,2)
         fz1 = Df(i6,3)

          x2 =  p(i5,1)
          y2 =  p(i5,2)
          z2 =  p(i5,3)
         fx2 = Df(i5,1)
         fy2 = Df(i5,2)
         fz2 = Df(i5,3)

          x3 =  p(i3,1)
          y3 =  p(i3,2)
          z3 =  p(i3,3)
         fx3 = Df(i3,1)
         fy3 = Df(i3,2)
         fz3 = Df(i3,3)

         call sslp_3d_integral_sing
     +
     +     (NGL
     + 
     +     ,x1,y1,z1
     +     ,x2,y2,z2
     +     ,x3,y3,z3
     + 
     +     ,fx1,fy1,fz1
     +     ,fx2,fy2,fz2
     +     ,fx3,fy3,fz3
     + 
     +     ,uxel,uyel,uzel
     +     )

c------------
       end if    ! done integrating over a singular element
c------------

  3    Continue

       us = us + uxel
       vs = vs + uyel
       ws = ws + uzel

  2   Continue

      slp(i,1) = - us/pi8
      slp(i,2) = - vs/pi8
      slp(i,3) = - ws/pi8

  1   Continue               ! loop over nodes

c-----
c done
c-----

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
