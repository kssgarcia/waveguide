      subroutine caps_3d_dfel
     +
     +  (npts
     +  ,nelm
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------------------------
c Compute:
c
c    1) the surface average of Df over each element
c
c    2) the value of Df at nodes
c       by averaging over hosting element values
c
c SYMBOLS:
c --------
c
c Dfel(k,i)  ith component of the surface average
c            value of df over at the kth element
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3),srtn(1026)
      Dimension   ne(1026,7)
      Dimension   vna(1026,3)
c     Dimension   df(1026,3)

      Dimension      n(512,6),nbe (512,3)
      Dimension  alpha(512),  beta(512), gamma(512)
      Dimension   arel(512)
      Dimension   Dfel(512,4)

      Dimension  elten(1026,3,3)     ! in plane elastic tensions
      Dimension     tst(1026,3)      ! transverse shear tension

      Dimension  xxi(6), eet(6)
      Dimension DxDx(6),DyDx(6),DzDx(6)
      Dimension DxDe(6),DyDe(6),DzDe(6)
      Dimension   vx(6),  vy(6),  vz(6)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna

      common/tension/srtn

      common/tenten/elten,tst
      common/dffel/Dfel

c-----------
c initialize
c-----------

      Ichoose = 2    ! for interp_p

c----------------------------------
c perform contour integration along
c the element sides
c----------------------------------

c     write (6,*)
c     write (6,*) "Dfel:"
c     write (6,*)

      Do 1 k=1,nelm

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be =  beta(k)
      ga = gamma(k)

      alc = 1.0-al
      bec = 1.0-be
      gac = 1.0-ga

c---
c six nodes
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
      eet(5) = gac

      xxi(6) = 0.0
      eet(6) = be

      Do i=1,6

        xi  = xxi(i)
        eta = eet(i)

c------
c interpolate the position vector,
c dx/dxi, dx/deta, normal vector
c line and surface metrics
c------

        call interp_p 
     +
     +      (Ichoose
     +      ,p(i1,1),p(i1,2),p(i1,3)
     +      ,p(i2,1),p(i2,2),p(i2,3)
     +      ,p(i3,1),p(i3,2),p(i3,3)
     +      ,p(i4,1),p(i4,2),p(i4,3)
     +      ,p(i5,1),p(i5,2),p(i5,3)
     +      ,p(i6,1),p(i6,2),p(i6,3)
     +      ,al,be,ga
     +      ,xi,eta
     +      ,x,y,z
     +      ,DxDx(i),DyDx(i),DzDx(i)
     +      ,DxDe(i),DyDe(i),DzDe(i)
     +      ,vx(i),vy(i),vz(i)
     +      ,hxi,het,hs
     +      )

c------
c interpolate the averaged
c unit normal vector
c-----

        call caps_3d_dfel_interp
     +
     +      (vna(i1,1),vna(i1,2),vna(i1,3)
     +      ,vna(i2,1),vna(i2,2),vna(i2,3)
     +      ,vna(i3,1),vna(i3,2),vna(i3,3)
     +      ,vna(i4,1),vna(i4,2),vna(i4,3)
     +      ,vna(i5,1),vna(i5,2),vna(i5,3)
     +      ,vna(i6,1),vna(i6,2),vna(i6,3)
     +      ,al,be,ga
     +      ,xi,eta
     +      ,vx(i),vy(i),vz(i)
     +      )

      End Do

c-----------------------------------
c Reset the tensions at the mid-nodes
c by linear interpolation,
c per Saroja's suggestion
c-----------------------------------

      Do j=1,3

       tst(i4,j) = alc*tst(i1,j) + al*tst(i2,j)
       tst(i5,j) = gac*tst(i3,j) + ga*tst(i2,j)
       tst(i6,j) = bec*tst(i1,j) + be*tst(i3,j)

       Do l=1,3
        elten(i4,j,l) = alc*elten(i1,j,l) + al*elten(i2,j,l)
        elten(i5,j,l) = gac*elten(i3,j,l) + ga*elten(i2,j,l)
        elten(i6,j,l) = bec*elten(i1,j,l) + be*elten(i3,j,l)
       End Do

      End Do

c--------------------------------
c initialize the contour integral
c--------------------------------

      cvx = 0.0
      cvy = 0.0
      cvz = 0.0

c----------------------------------
c line integral along segment 1-4-2
c
c side of eta = 0
c----------------------------------

      bvx1 = vy(1)*DzDx(1)-vz(1)*DyDx(1)
      bvy1 = vz(1)*DxDx(1)-vx(1)*DzDx(1)
      bvz1 = vx(1)*DyDx(1)-vy(1)*DxDx(1)

      bvx2 = vy(4)*DzDx(4)-vz(4)*DyDx(4)
      bvy2 = vz(4)*DxDx(4)-vx(4)*DzDx(4)
      bvz2 = vx(4)*DyDx(4)-vy(4)*DxDx(4)

      bvx3 = vy(2)*DzDx(2)-vz(2)*DyDx(2)
      bvy3 = vz(2)*DxDx(2)-vx(2)*DzDx(2)
      bvz3 = vx(2)*DyDx(2)-vy(2)*DxDx(2)

c---
c project binormal on tension
c---

      prj1 = bvx1*tst(i1,1)+bvy1*tst(i1,2)+bvz1*tst(i1,3)
      prj2 = bvx2*tst(i4,1)+bvy2*tst(i4,2)+bvz2*tst(i4,3)
      prj3 = bvx3*tst(i2,1)+bvy3*tst(i2,2)+bvz3*tst(i2,3)

      tx1 =
     +      bvx1 * elten(i1,1,1)
     +    + bvy1 * elten(i1,2,1)
     +    + bvz1 * elten(i1,3,1)
     +    + prj1 *   vna(i1,1)
     +    + bvx1 *  srtn(i1)

      ty1 =
     +      bvx1 * elten(i1,1,2)
     +    + bvy1 * elten(i1,2,2)
     +    + bvz1 * elten(i1,3,2)
     +    + prj1 *   vna(i1,2)
     +    + bvy1 *   srtn(i1)

      tz1 =
     +      bvx1 * elten(i1,1,3)
     +    + bvy1 * elten(i1,2,3)
     +    + bvz1 * elten(i1,3,3)
     +    + prj1 *   vna(i1,3)
     +    + bvz1 *  srtn(i1)


      tx2 =
     +      bvx2 * elten(i4,1,1)
     +    + bvy2 * elten(i4,2,1)
     +    + bvz2 * elten(i4,3,1)
     +    + prj2 *   vna(i4,1)
     +    + bvx2 *  srtn(i4)

      ty2 =
     +      bvx2 * elten(i4,1,2)
     +    + bvy2 * elten(i4,2,2)
     +    + bvz2 * elten(i4,3,2)
     +    + prj2 *   vna(i4,2)
     +    + bvy2 *  srtn(i4)

      tz2 =
     +      bvx2 * elten(i4,1,3)
     +    + bvy2 * elten(i4,2,3)
     +    + bvz2 * elten(i4,3,3) 
     +    + prj2 *   vna(i4,3)
     +    + bvz2 *  srtn(i4)


      tx3 =
     +      bvx3 * elten(i2,1,1)
     +    + bvy3 * elten(i2,2,1)
     +    + bvz3 * elten(i2,3,1)
     +    + prj3 *   vna(i2,1)
     +    + bvx3 *  srtn(i2)

      ty3 = 
     +      bvx3 * elten(i2,1,2)
     +    + bvy3 * elten(i2,2,2)
     +    + bvz3 * elten(i2,3,2)
     +    + prj3 *   vna(i2,2)
     +    + bvy3 *  srtn(i2)

      tz3 =
     +      bvx3 * elten(i2,1,3)
     +    + bvy3 * elten(i2,2,3)
     +    + bvz3 * elten(i2,3,3)
     +    + prj3 *   vna(i2,3)
     +    + bvz3 *  srtn(i2)

      cvx = cvx + al*tx1 + tx2 + alc*tx3
      cvy = cvy + al*ty1 + ty2 + alc*ty3
      cvz = cvz + al*tz1 + tz2 + alc*tz3

c----------------------------------
c line integral along segment 2-5-3
c
c mixed side
c----------------------------------

      bvx1 = vy(2)*DzDx(2)-vz(2)*DyDx(2)
      bvy1 = vz(2)*DxDx(2)-vx(2)*DzDx(2)
      bvz1 = vx(2)*DyDx(2)-vy(2)*DxDx(2)

      bvx2 = vy(5)*DzDx(5)-vz(5)*DyDx(5)
      bvy2 = vz(5)*DxDx(5)-vx(5)*DzDx(5)
      bvz2 = vx(5)*DyDx(5)-vy(5)*DxDx(5)

      bvx3 = vy(3)*DzDx(3)-vz(3)*DyDx(3)
      bvy3 = vz(3)*DxDx(3)-vx(3)*DzDx(3)
      bvz3 = vx(3)*DyDx(3)-vy(3)*DxDx(3)

      prj1 = bvx1*tst(i2,1)+bvy1*tst(i2,2)+bvz1*tst(i2,3)
      prj2 = bvx2*tst(i5,1)+bvy2*tst(i5,2)+bvz2*tst(i5,3)
      prj3 = bvx3*tst(i3,1)+bvy3*tst(i3,2)+bvz3*tst(i3,3)

c---
c project binormal on tension
c---

      tx1 =
     +      bvx1 * elten(i2,1,1)
     +    + bvy1 * elten(i2,2,1)
     +    + bvz1 * elten(i2,3,1)
     +    + prj1 *   vna(i2,1)
     +    + bvx1 *  srtn(i2)

      ty1 =
     +      bvx1 * elten(i2,1,2)
     +    + bvy1 * elten(i2,2,2)
     +    + bvz1 * elten(i2,3,2)
     +    + prj1 *   vna(i2,2)
     +    + bvy1 *  srtn(i2)

      tz1 =
     +      bvx1 * elten(i2,1,3)
     +    + bvy1 * elten(i2,2,3)
     +    + bvz1 * elten(i2,3,3)
     +    + prj1 *   vna(i2,3)
     +    + bvz1 *   srtn(i2)


      tx2 =
     +      bvx2 * elten(i5,1,1)
     +    + bvy2 * elten(i5,2,1)
     +    + bvz2 * elten(i5,3,1)
     +    + prj2 *   vna(i5,1)
     +    + bvx2 *  srtn(i5)

      ty2 =
     +      bvx2 * elten(i5,1,2)
     +    + bvy2 * elten(i5,2,2)
     +    + bvz2 * elten(i5,3,2)
     +    + prj2 *   vna(i5,2)
     +    + bvy2 *  srtn(i5)

      tz2 =
     +      bvx2 * elten(i5,1,3)
     +    + bvy2 * elten(i5,2,3)
     +    + bvz2 * elten(i5,3,3)
     +    + prj2 *   vna(i5,3)
     +    + bvz2 *   srtn(i5)


      tx3 =
     +      bvx3 * elten(i3,1,1)
     +    + bvy3 * elten(i3,2,1)
     +    + bvz3 * elten(i3,3,1)
     +    + prj3 *   vna(i3,1)
     +    + bvx3 *  srtn(i3)

      ty3 =
     +      bvx3 * elten(i3,1,2)
     +    + bvy3 * elten(i3,2,2)
     +    + bvz3 * elten(i3,3,2)
     +    + prj3 *   vna(i3,2)
     +    + bvy3 *   srtn(i3)

      tz3 =
     +      bvx3 * elten(i3,1,3)
     +    + bvy3 * elten(i3,2,3)
     +    + bvz3 * elten(i3,3,3)
     +    + prj3 *   vna(i3,3)
     +    + bvz3 *  srtn(i3)

      cvx = cvx - (gac*tx1 + tx2 + ga*tx3)
      cvy = cvy - (gac*ty1 + ty2 + ga*ty3)
      cvz = cvz - (gac*tz1 + tz2 + ga*tz3)

c---

      bvx1 = vy(2)*DzDe(2)-vz(2)*DyDe(2)
      bvy1 = vz(2)*DxDe(2)-vx(2)*DzDe(2)
      bvz1 = vx(2)*DyDe(2)-vy(2)*DxDe(2)

      bvx2 = vy(5)*DzDe(5)-vz(5)*DyDe(5)
      bvy2 = vz(5)*DxDe(5)-vx(5)*DzDe(5)
      bvz2 = vx(5)*DyDe(5)-vy(5)*DxDe(5)

      bvx3 = vy(3)*DzDe(3)-vz(3)*DyDe(3)
      bvy3 = vz(3)*DxDe(3)-vx(3)*DzDe(3)
      bvz3 = vx(3)*DyDe(3)-vy(3)*DxDe(3)

      prj1 = bvx1*tst(i2,1)+bvy1*tst(i2,2)+bvz1*tst(i2,3)
      prj2 = bvx2*tst(i5,1)+bvy2*tst(i5,2)+bvz2*tst(i5,3)
      prj3 = bvx3*tst(i3,1)+bvy3*tst(i3,2)+bvz3*tst(i3,3)

c---
c project binormal on tension
c---

      tx1 = 
     +    + bvx1 * elten(i2,1,1)
     +    + bvy1 * elten(i2,2,1)
     +    + bvz1 * elten(i2,3,1)
     +    + prj1 *   vna(i2,1)
     +    + bvx1 *  srtn(i2)

      ty1 =
     +    + bvx1 * elten(i2,1,2)
     +    + bvy1 * elten(i2,2,2)
     +    + bvz1 * elten(i2,3,2)
     +    + prj1 *   vna(i2,2)
     +    + bvy1 *  srtn(i2)

      tz1 =
     +    + bvx1 * elten(i2,1,3)
     +    + bvy1 * elten(i2,2,3)
     +    + bvz1 * elten(i2,3,3)
     +    + prj1 *   vna(i2,3)
     +    + bvz1 *   srtn(i2)


      tx2 =
     +    + bvx2 * elten(i5,1,1)
     +    + bvy2 * elten(i5,2,1)
     +    + bvz2 * elten(i5,3,1)
     +    + prj2 *   vna(i5,1)
     +    + bvx2 *  srtn(i5)

      ty2 =
     +    + bvx2 * elten(i5,1,2)
     +    + bvy2 * elten(i5,2,2)
     +    + bvz2 * elten(i5,3,2)
     +    + prj2 *   vna(i5,2)
     +    + bvy2 *   srtn(i5)

      tz2 =
     +    + bvx2 * elten(i5,1,3)
     +    + bvy2 * elten(i5,2,3)
     +    + bvz2 * elten(i5,3,3)
     +    + prj2 *   vna(i5,3)
     +    + bvz2 *  srtn(i5)


      tx3 =
     +    + bvx3 * elten(i3,1,1)
     +    + bvy3 * elten(i3,2,1)
     +    + bvz3 * elten(i3,3,1)
     +    + prj3 *   vna(i3,1)
     +    + bvx3 *  srtn(i3)

      ty3 =
     +    + bvx3 * elten(i3,1,2)
     +    + bvy3 * elten(i3,2,2)
     +    + bvz3 * elten(i3,3,2)
     +    + prj3 *   vna(i3,2)
     +    + bvy3 *  srtn(i3)

      tz3 =
     +    + bvx3 * elten(i3,1,3)
     +    + bvy3 * elten(i3,2,3)
     +    + bvz3 * elten(i3,3,3)
     +    + prj3 *   vna(i3,3)
     +    + bvz3 *   srtn(i3)

      cvx = cvx + gac*tx1 + tx2 + ga*tx3
      cvy = cvy + gac*ty1 + ty2 + ga*ty3
      cvz = cvz + gac*tz1 + tz2 + ga*tz3

c--------------------------------------
c line integral along segment 3-6-1
c
c side of xi = 0
c--------------------------------------

      bvx1 = vy(1)*DzDe(1)-vz(1)*DyDe(1)
      bvy1 = vz(1)*DxDe(1)-vx(1)*DzDe(1)
      bvz1 = vx(1)*DyDe(1)-vy(1)*DxDe(1)

      bvx2 = vy(6)*DzDe(6)-vz(6)*DyDe(6)
      bvy2 = vz(6)*DxDe(6)-vx(6)*DzDe(6)
      bvz2 = vx(6)*DyDe(6)-vy(6)*DxDe(6)

      bvx3 = vy(3)*DzDe(3)-vz(3)*DyDe(3)
      bvy3 = vz(3)*DxDe(3)-vx(3)*DzDe(3)
      bvz3 = vx(3)*DyDe(3)-vy(3)*DxDe(3)

      prj1 = bvx1*tst(i1,1)+bvy1*tst(i1,2)+bvz1*tst(i1,3)
      prj2 = bvx2*tst(i6,1)+bvy2*tst(i6,2)+bvz2*tst(i6,3)
      prj3 = bvx3*tst(i3,1)+bvy3*tst(i3,2)+bvz3*tst(i3,3)

c---
c project binormal on tension
c---

      tx1 =
     +      bvx1 * elten(i1,1,1)
     +    + bvy1 * elten(i1,2,1)
     +    + bvz1 * elten(i1,3,1)
     +    + prj1 *   vna(i1,1)
     +    + bvx1 *  srtn(i1)

      ty1 =
     +      bvx1 * elten(i1,1,2)
     +    + bvy1 * elten(i1,2,2)
     +    + bvz1 * elten(i1,3,2)
     +    + prj1 *   vna(i1,2)
     +    + bvy1 *  srtn(i1)

      tz1 = 
     +      bvx1 * elten(i1,1,3)
     +    + bvy1 * elten(i1,2,3)
     +    + bvz1 * elten(i1,3,3)
     +    + prj1 *   vna(i1,3)
     +    + bvz1 *  srtn(i1)


      tx2 = 
     +      bvx2 * elten(i6,1,1)
     +    + bvy2 * elten(i6,2,1)
     +    + bvz2 * elten(i6,3,1)
     +    + prj2 *   vna(i6,1)
     +    + bvx2 *  srtn(i6)

      ty2 = 
     +      bvx2 * elten(i6,1,2)
     +    + bvy2 * elten(i6,2,2)
     +    + bvz2 * elten(i6,3,2)
     +    + prj2 *   vna(i6,2)
     +    + bvy2 *  srtn(i6)

      tz2 = 
     +      bvx2 * elten(i6,1,3)
     +    + bvy2 * elten(i6,2,3)
     +    + bvz2 * elten(i6,3,3)
     +    + prj2 *   vna(i6,3)
     +    + bvz2 *  srtn(i6)


      tx3 = 
     +      bvx3 * elten(i3,1,1)
     +    + bvy3 * elten(i3,2,1)
     +    + bvz3 * elten(i3,3,1)
     +    + prj3 *   vna(i3,1)
     +    + bvx3 *   srtn(i3)

      ty3 =
     +      bvx3 * elten(i3,1,2)
     +    + bvy3 * elten(i3,2,2)
     +    + bvz3 * elten(i3,3,2)
     +    + prj3 *   vna(i3,2)
     +    + bvy3 *  srtn(i3)

      tz3 = 
     +      bvx3 * elten(i3,1,3)
     +    + bvy3 * elten(i3,2,3)
     +    + bvz3 * elten(i3,3,3)
     +    + prj3 *   vna(i3,3)
     +    + bvz3 *  srtn(i3)

      cvx = cvx - (be*tx1 + tx2 + bec*tx3)
      cvy = cvy - (be*ty1 + ty2 + bec*ty3)
      cvz = cvz - (be*tz1 + tz2 + bec*tz3)

c-----------------------------------------------
c Done computing the line integral
c
c the factor 0.5 comes from the trapezoidal rule
c-----------------------------------------------

      cf = 0.50/arel(k)

      Dfel(k,1) = cf*cvx
      Dfel(k,2) = cf*cvy
      Dfel(k,3) = cf*cvz

      Dfel(k,4) = sqrt(Dfel(k,1)**2+Dfel(k,2)**2
     +                +Dfel(k,3)**2)

c     write (6,100) k,(Dfel(k,l),l=1,4)

  1   Continue

c--------------------------------------------------
c  compute:
c          Df at the nodes (points) 
c          by averaging values from adjacent elements
c--------------------------------------------------
c
c     Do i=1,npts
c
c      Df(i,1) = 0.0
c      Df(i,2) = 0.0
c      Df(i,3) = 0.0
c
c      Do j=2,ne(i,1)+1
c        m = ne(i,j)
c        df(i,1) = df(i,1) + Dfel(m,1)
c        df(i,2) = df(i,2) + Dfel(m,2)
c        df(i,3) = df(i,3) + Dfel(m,3)
c      End Do
c
c      cf = ne(i,1)+1.0-1.0
c
c      df(i,1) = df(i,1)/cf
c      df(i,2) = df(i,2)/cf
c      df(i,3) = df(i,3)/cf
c
c     End Do
c
c-----
c done
c-----

 100  Format (1x,i3,10(1x,f10.5))

      return
      end
