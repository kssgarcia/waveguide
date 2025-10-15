      subroutine elm_geom 
     +
     +   (nelm,npts
     +   ,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c----------------------------------------
c Compute:
c
c  the surface area of the individual elements
c  x, y, and z moments over each element
c
c  the total particle surface area and volume
c
c  the mean curvature of each element
c
c  Also compute the average value of the
c  normal vector at each node.
c  This is done by computing the normal vector 
c  at the six nodes of each triangle, 
c  and then averaging the contributions
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p     (1026,3)
      Dimension ne    (1026,7)
      Dimension vna   (1026,3)
      Dimension itally(1026)

      Dimension  n    (512,6),nbe (512,3)
      Dimension  alpha(512),  beta(512), gamma(512)
      Dimension  arel (512),  xmom(512),  ymom(512),zmom(512)
      Dimension crvmel(512)

      Dimension xxi (6),eet (6)
      Dimension DxDx(6),DyDx(6),DzDx(6)
      Dimension DxDe(6),DyDe(6),DzDe(6)
      Dimension vx  (6),vy  (6),vz  (6)

      Dimension zz (20),ww (20)
      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo3/vna
      common/geo6/crvmel

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

c-----------
c initialize
c-----------

      Ichoose = 2     ! for the point interpolation subroutine

      area = 0.0D0
      vlm  = 0.0D0

      cx = 0.0D0
      cy = 0.0D0
      cz = 0.0D0

      Do i=1,npts
         vna(i,1) = 0.0D0
         vna(i,2) = 0.0D0
         vna(i,3) = 0.0D0
       itally(i)   = 0
      End Do

c-----

      Do 1 k=1,nelm

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al  = alpha(k)
      be  = beta (k)
      ga  = gamma(k)

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      arel(k) = 0.0D0

      xmom(k) = 0.0D0
      ymom(k) = 0.0D0
      zmom(k) = 0.0D0

c------------------------------------------------
c compute surface area of the individual elements
c         x, y, and z moments over each element
c         total particle surface area and volume
c------------------------------------------------

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call interp_p 
     +
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +    ,al,be,ga
     +    ,xi,eta
     +    ,x,y,z
     +    ,DxDxi,DyDxi,DzDxi
     +    ,DxDet,DyDet,DzDet
     +    ,vnx,vny,vnz
     +    ,hxi,het,hs
     +    ,Ichoose
     +    )

        cf = hs*wq(i)

        arel(k) = arel(k) + cf
        xmom(k) = xmom(k) + cf*x
        ymom(k) = ymom(k) + cf*y
        zmom(k) = zmom(k) + cf*z
        vlm = vlm + (x*vnx+y*vny+z*vnz)*cf

      End Do

      arel(k) = 0.5D0*arel(k)

      xmom(k) = 0.5D0*xmom(k)
      ymom(k) = 0.5D0*ymom(k)
      zmom(k) = 0.5D0*zmom(k)

      area = area +arel(k)
      cx   = cx   +xmom(k)
      cy   = cy   +ymom(k)
      cz   = cz   +zmom(k)

c------------------------------------------------------
c compute the average value of the normal vector
c
c         the mean curvature as a contour integral
c         using the nifty formula (4.2.10) 
c         of Pozrikidis (1997)
c------------------------------------------------------

      xxi(1) = 0.0D0
      eet(1) = 0.0D0

      xxi(2) = 1.0D0
      eet(2) = 0.0D0

      xxi(3) = 0.0D0
      eet(3) = 1.0D0

      xxi(4) = al
      eet(4) = 0.0D0

      xxi(5) = ga
      eet(5) = gac

      xxi(6) = 0.0D0
      eet(6) = be

      Do i=1,6

        xi  = xxi(i)
        eta = eet(i)

        call interp_p 
     +
     +   (p(i1,1),p(i1,2),p(i1,3)
     +   ,p(i2,1),p(i2,2),p(i2,3)
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,p(i4,1),p(i4,2),p(i4,3)
     +   ,p(i5,1),p(i5,2),p(i5,3)
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +   ,al,be,ga
     +   ,xi,eta
     +   ,x,y,z
     +   ,DxDx(i),DyDx(i),DzDx(i)
     +   ,DxDe(i),DyDe(i),DzDe(i)
     +   ,vx(i),vy(i),vz(i)
     +   ,hxi,het,hs
     +   ,Ichoose
     +   )

        m = n(k,i)

        vna(m,1) = vna(m,1) + vx(i)
        vna(m,2) = vna(m,2) + vy(i)
        vna(m,3) = vna(m,3) + vz(i)

        itally(m) = itally(m)+1

      End Do

c----------------------------------
c computation of curvature
c line integral along segment 1-4-2
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

      crvx = al*bvx1 + bvx2 + alc*bvx3
      crvy = al*bvy1 + bvy2 + alc*bvy3
      crvz = al*bvz1 + bvz2 + alc*bvz3

c---
c computation of curvature
c line integral along segment 2-5-3
c---

      bvx1 = vy(2)*DzDx(2)-vz(2)*DyDx(2)
      bvy1 = vz(2)*DxDx(2)-vx(2)*DzDx(2)
      bvz1 = vx(2)*DyDx(2)-vy(2)*DxDx(2)

      bvx2 = vy(5)*DzDx(5)-vz(5)*DyDx(5)
      bvy2 = vz(5)*DxDx(5)-vx(5)*DzDx(5)
      bvz2 = vx(5)*DyDx(5)-vy(5)*DxDx(5)

      bvx3 = vy(3)*DzDx(3)-vz(3)*DyDx(3)
      bvy3 = vz(3)*DxDx(3)-vx(3)*DzDx(3)
      bvz3 = vx(3)*DyDx(3)-vy(3)*DxDx(3)

c---
      crvx = crvx - gac*bvx1 - bvx2 - ga*bvx3
      crvy = crvy - gac*bvy1 - bvy2 - ga*bvy3
      crvz = crvz - gac*bvz1 - bvz2 - ga*bvz3
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

      crvx = crvx + gac*bvx1 + bvx2 + ga*bvx3
      crvy = crvy + gac*bvy1 + bvy2 + ga*bvy3
      crvz = crvz + gac*bvz1 + bvz2 + ga*bvz3

c--
c computation of curvature
c line integral along segment 3-6-1
c---

      bvx1 = vy(1)*DzDe(1)-vz(1)*DyDe(1)
      bvy1 = vz(1)*DxDe(1)-vx(1)*DzDe(1)
      bvz1 = vx(1)*DyDe(1)-vy(1)*DxDe(1)

      bvx2 = vy(6)*DzDe(6)-vz(6)*DyDe(6)
      bvy2 = vz(6)*DxDe(6)-vx(6)*DzDe(6)
      bvz2 = vx(6)*DyDe(6)-vy(6)*DxDe(6)

      bvx3 = vy(3)*DzDe(3)-vz(3)*DyDe(3)
      bvy3 = vz(3)*DxDe(3)-vx(3)*DzDe(3)
      bvz3 = vx(3)*DyDe(3)-vy(3)*DxDe(3)

      crvx = crvx - be*bvx1 - bvx2 - bec*bvx3
      crvy = crvy - be*bvy1 - bvy2 - bec*bvy3
      crvz = crvz - be*bvz1 - bvz2 - bec*bvz3

c---
c will project curvature vector onto the normal
c vector at the centroid
c---

      xi  = 1.0D0/3.0D0
      eta = 1.0D0/3.0D0

      call interp_p 
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +  ,alpha(k),beta(k),gamma(k)
     +  ,xi,eta
     +  ,x,y,z
     +  ,DxDxi,DyDxi,DzDxi
     +  ,DxDet,DyDet,DzDet
     +  ,vnx,vny,vnz
     +  ,hxi,het,hs
     +  ,Ichoose
     +  )

      crvmel(k) = 0.25D0*(crvx*vnx+crvy*vny+crvz*vnz)/arel(k)

1     Continue

c-----------------------
c averaged normal vector
c-----------------------

      Do i=1,npts

        par = float(itally(i))

        vna(i,1) = vna(i,1)/par
        vna(i,2) = vna(i,2)/par
        vna(i,3) = vna(i,3)/par

        par = sqrt( vna(i,1)**2+vna(i,2)**2+vna(i,3)**2 )

        vna(i,1) = vna(i,1)/par
        vna(i,2) = vna(i,2)/par
        vna(i,3) = vna(i,3)/par

      End Do

c---
c final computation of the surface-centroid
c and volume
c---

      cx = cx/area
      cy = cy/area
      cz = cz/area

      vlm = 0.5*vlm/3.0

c-----
c Done
c-----

      return
      end
