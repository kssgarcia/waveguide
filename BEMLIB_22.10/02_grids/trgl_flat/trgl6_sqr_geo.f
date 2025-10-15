      subroutine printel (k,Index,c)

c==========================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c==========================================

c-------------------------------------
c print successive nodes of an element
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  p(1090,3)
      Dimension ne(1090,9)
      Dimension c(1090)

      Dimension n(512,6),nbe(512,3)

      common/points/p,ne
      common/elmnts/n,nbe

c----------
c constants
c----------

      nfour  = 4
      nseven = 7

c---------------------------------
c Draw one 6-point curved triangle
c---------------------------------

      if(Index.eq.1) then

c     write (6,100) nseven
c     write (1,100) nseven
        i = n(k,1)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,4)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,2)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,5)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,3)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,6)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,1)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c---------------------------------
c Draw four 3-point flat triangles
c---------------------------------

      elseif(Index.eq.2) then

c--- first

c     write (1,100) nfour
      i = n(k,1)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,1)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c--- second

c     write (1,100) nfour
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,2)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,5)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c--- third

c     write (1,100) nfour
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,5)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c--- fourth

c     write (1,100) nfour
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,5)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,3)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c----
      end if
c----

c---
c Done
c---

  100 Format(1x,i4,10(1x,f12.5))
  101 Format(10(1x,f12.5))

      Return
      End

c==============================================

      subroutine printel_stl (k)

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c-------------------------------------------
c Print element k in file unit 8 (trgl6.stl)
c The element divided into four subelements
c------------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension  p(1090,3)
      Dimension ne(1090,9)

      Dimension  n(512,6),nbe(512,3)

      common/points/p,ne
      common/elmnts/n,nbe

c------------------------------------
c draw four three-node flat triangles
c------------------------------------

c--- first

      apos1x = p(n(k,4),1)-p(n(k,1),1)
      apos1y = p(n(k,4),2)-p(n(k,1),2)
      apos1z = p(n(k,4),3)-p(n(k,1),3)
      apos2x = p(n(k,6),1)-p(n(k,1),1)
      apos2y = p(n(k,6),2)-p(n(k,1),2)
      apos2z = p(n(k,6),3)-p(n(k,1),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm
      vnormz = vnormz/vnormm

      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,1)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,4)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,6)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c--- second

      apos1x = p(n(k,2),1)-p(n(k,4),1)
      apos1y = p(n(k,2),2)-p(n(k,4),2)
      apos1z = p(n(k,2),3)-p(n(k,4),3)
      apos2x = p(n(k,5),1)-p(n(k,4),1)
      apos2y = p(n(k,5),2)-p(n(k,4),2)
      apos2z = p(n(k,5),3)-p(n(k,4),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm
      vnormz = vnormz/vnormm
      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,4)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,2)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,5)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c--- third

      apos1x = p(n(k,5),1)-p(n(k,4),1)
      apos1y = p(n(k,5),2)-p(n(k,4),2)
      apos1z = p(n(k,5),3)-p(n(k,4),3)
      apos2x = p(n(k,6),1)-p(n(k,4),1)
      apos2y = p(n(k,6),2)-p(n(k,4),2)
      apos2z = p(n(k,6),3)-p(n(k,4),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm
      vnormz = vnormz/vnormm

      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,4)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,5)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,6)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c--- fourth

      apos1x = p(n(k,5),1)-p(n(k,6),1)
      apos1y = p(n(k,5),2)-p(n(k,6),2)
      apos1z = p(n(k,5),3)-p(n(k,6),3)
      apos2x = p(n(k,3),1)-p(n(k,6),1)
      apos2y = p(n(k,3),2)-p(n(k,6),2)
      apos2z = p(n(k,3),3)-p(n(k,6),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm

      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,6)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,5)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,3)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c-----
c done
c-----

   91 Format(" facet normal",10(1x,f12.5))
   92 Format("   outer loop")
   93 Format("    vertex",10(1x,f12.5))
   94 Format("    vertex",10(1x,f12.5))
   95 Format("    vertex",10(1x,f12.5))
   96 Format("   endloop")
   97 Format(" endfacet")

  100 Format(1x,i4,10(1x,f12.5))

      Return
      End


c=========================================

      subroutine abc 
     +
     +   (x1,y1,z1
     +   ,x2,y2,z2
     +   ,x3,y3,z3
     +   ,x4,y4,z4
     +   ,x5,y5,z5
     +   ,x6,y6,z6
     +   ,al,be,ga
     +   )

c==========================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement
c==========================================

c--------------------------------
c Compute the xi and eta mapping
c coefficients alpha, beta, gamma
c--------------------------------

      Implicit Double Precision (a-h,o-z)

      d42 = sqrt( (x4-x2)**2 + (y4-y2)**2 + (z4-z2)**2 )
      d41 = sqrt( (x4-x1)**2 + (y4-y1)**2 + (z4-z1)**2 )
      d63 = sqrt( (x6-x3)**2 + (y6-y3)**2 + (z6-z3)**2 )
      d61 = sqrt( (x6-x1)**2 + (y6-y1)**2 + (z6-z1)**2 )
      d52 = sqrt( (x5-x2)**2 + (y5-y2)**2 + (z5-z2)**2 )
      d53 = sqrt( (x5-x3)**2 + (y5-y3)**2 + (z5-z3)**2 )

      al = 1.0D0/(1.0D0+d42/d41)
      be = 1.0D0/(1.0D0+d63/d61)
      ga = 1.0D0/(1.0D0+d52/d53)

c---
c Done
c---

      Return
      End

c=========================================

      subroutine interp_vn_hs 
     +
     +    (x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +    ,x4,y4,z4
     +    ,x5,y5,z5
     +    ,x6,y6,z6
     +    ,alpha,beta,gamma
     +    ,xi,eta
     +    ,x,y,z
     +    ,vnx,vny,vnz
     +    ,hs
     +    )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c------------------------------------
c  compute the normal vector
c  and the surface metric coefficient
c------------------------------------

      Implicit Double Precision (a-h,o-z)
 
c---
c preparations
c---

      alc = 1.0D0-alpha
      bec = 1.0D0-beta
      gmc = 1.0D0-gamma

c---
c  compute position vector (x, y, z)
c---

      ph2 = xi*(xi-alpha+eta*(alpha-gamma)/gmc)/alc
      ph3 = eta*(eta-beta+xi*(beta+gamma-1.0D0)/gamma)/bec
      ph4 = xi*(1.0D0-xi-eta)/(alpha*alc)
      ph5 = xi*eta/(gamma*gmc)
      ph6 = eta*(1.0D0-xi-eta)/(beta*bec)
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

      x = x1*ph1 + x2*ph2 + x3*ph3 + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3 + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3 + z4*ph4 + z5*ph5 + z6*ph6

c---
c compute dx/dxi from xi derivatives of phi
c---

      dph2 = (2.0D0*xi-alpha+eta*(alpha-gamma)/gmc)/alc
      dph3 =  eta*(beta+gamma-1.0D0)/(gamma*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/(alpha*alc)
      dph5 =  eta/(gamma*gmc)
      dph6 = -eta/(beta*bec)
      dph1 = -dph2-dph3-dph4-dph5-dph6

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c---
c  compute dx/deta from eta derivatives of phi
c---

      dph2 = xi*(alpha-gamma)/(alc*gmc)
      dph3 = (2.0*eta-beta+xi*(beta+gamma-1.0)/gamma)/bec
      dph4 =-xi/(alpha*alc)
      dph5 = xi/(gamma*gmc)
      dph6 = (1.0-xi-2.0*eta)/(beta*bec)
      dph1 = -dph2-dph3-dph4-dph5-dph6

      DxDet = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDet = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDet = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c---
c  normal vector:       vn = (DxDxi)x(DxDeta)
c  and surface metric:  hs=norm(vn)
c---

      vnx = DyDxi*DzDet - DyDet*DzDxi
      vny = DzDxi*DxDet - DzDet*DxDxi
      vnz = DxDxi*DyDet - DxDet*DyDxi

      hs = sqrt(vnx*vnx + vny*vny + vnz*vnz)

      vnx = vnx/hs
      vny = vny/hs
      vnz = vnz/hs

c---
c Done
c---

      return
      end

c=========================================================

      subroutine elm_geom
     +
     +   (nelm,npts,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c----------------------------------------
c Compute:
c
c (a) the surface area of the individual elements.
c (b) the x, y, and z surface moments over each element.
c (c) the total surface area and volume.
c (d) the mean curvature of each element.
c (e) the average value of the normal vector at each node;
c     this is done by computing the normal vector at the 6 nodes
c     of each triangle, and then averaging the contributions.
c
c  SYMBOLS:
c  --------
c
c  area:     total surface area
c  vlm:      total volume enclosed by the surface area
c  cx,cy,cz: surface centroid
c
c  itally(i): number of elements sharing a node
c             for averaging
c
c  Iedge: interior/edge/corner node index:
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
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension      p(1090,3)
      Dimension     ne(1090,9)
      Dimension    vna(1090,3)
      Dimension vna_sv(1090,3)
      Dimension itally(1090)     ! internal

      Dimension  Iedge(1090,4)

      Dimension       n(512,6), nbe(512,3)
      Dimension   alpha(512),  beta(512), gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)

      Dimension  crvmel(512)

      Dimension  xxi(6), eet(6)

      Dimension DxDx(6),DyDx(6),DzDx(6)
      Dimension DxDe(6),DyDe(6),DzDe(6)
      Dimension   vx(6),  vy(6),  vz(6)

      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/edgepoints/Iedge

      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna
      common/geo6/crvmel

      common/trq/xiq,etq,wq

c---
c initialize
c---

      Ichoose = 2   ! for subroutine interp_p

      area = 0.0
      vlm  = 0.0
      cx   = 0.0
      cy   = 0.0
      cz   = 0.0

      Do i=1,npts

        vna(i,1) = 0.0
        vna(i,2) = 0.0
        vna(i,3) = 0.0

        itally(i) = 0    ! for averaging at a node

      End Do

c------------------
      Do 1 k=1,nelm    ! loop over elements
c------------------

       i1 = n(k,1)     ! global node labels
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be = beta (k)
       ga = gamma(k)

       al1 = 1.0-al
       be1 = 1.0-be
       ga1 = 1.0-ga

c---
c initialize
c element variables
c---

       arel(k) = 0.0   ! element area

       xmom(k) = 0.0   ! element moments
       ymom(k) = 0.0
       zmom(k) = 0.0

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
       eet(5) = ga1

       xxi(6) = 0.0
       eet(6) = be

c------------------------------------------------
c compute:
c         surface area of the individual elements
c         x, y, and z surface moments
c         total surface area and enclosed volume
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

      arel(k) = 0.5*arel(k)
      xmom(k) = 0.5*xmom(k)
      ymom(k) = 0.5*ymom(k)
      zmom(k) = 0.5*zmom(k)

      area = area +arel(k)

      cx = cx + xmom(k)
      cy = cy + ymom(k)
      cz = cz + zmom(k)

c------------------------------------------------------
c compute:
c
c        (a) the average value of the normal vector
c        (b) the mean curvature as a contour integral
c            using the nifty formula (4.2.10)
c            of Pozrikidis (1997)
c------------------------------------------------------

c---
c triangle coordinates
c of the nodes
c---

      Do i=1,6

        xi  = xxi(i)
        eta = eet(i)

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
     +    ,DxDx(i),DyDx(i),DzDx(i)
     +    ,DxDe(i),DyDe(i),DzDe(i)
     +    ,vx(i),vy(i),vz(i)
     +    ,hxi,het,hs
     +    ,Ichoose
     +    )

        m = n(k,i)   ! global index of local node i on element k

        vna(m,1) = vna(m,1) + vx(i)
        vna(m,2) = vna(m,2) + vy(i)
        vna(m,3) = vna(m,3) + vz(i)

        itally(m) = itally(m)+1

      End Do

c-----------------------------------
c computation of mean curvature:
c line integral along segment 1-4-2
c-----------------------------------

      bvx1 = vy(1)*DzDx(1)-vz(1)*DyDx(1)
      bvy1 = vz(1)*DxDx(1)-vx(1)*DzDx(1)
      bvz1 = vx(1)*DyDx(1)-vy(1)*DxDx(1)

      bvx2 = vy(4)*DzDx(4)-vz(4)*DyDx(4)
      bvy2 = vz(4)*DxDx(4)-vx(4)*DzDx(4)
      bvz2 = vx(4)*DyDx(4)-vy(4)*DxDx(4)

      bvx3 = vy(2)*DzDx(2)-vz(2)*DyDx(2)
      bvy3 = vz(2)*DxDx(2)-vx(2)*DzDx(2)
      bvz3 = vx(2)*DyDx(2)-vy(2)*DxDx(2)

      crvx = al*bvx1 + bvx2 + al1*bvx3
      crvy = al*bvy1 + bvy2 + al1*bvy3
      crvz = al*bvz1 + bvz2 + al1*bvz3

c---
c computation of mean curvature:
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
      crvx = crvx - ga1*bvx1 - bvx2 - ga*bvx3
      crvy = crvy - ga1*bvy1 - bvy2 - ga*bvy3
      crvz = crvz - ga1*bvz1 - bvz2 - ga*bvz3
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

      crvx = crvx + ga1*bvx1 + bvx2 + ga*bvx3
      crvy = crvy + ga1*bvy1 + bvy2 + ga*bvy3
      crvz = crvz + ga1*bvz1 + bvz2 + ga*bvz3

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

      crvx = crvx - be*bvx1 - bvx2 - be1*bvx3
      crvy = crvy - be*bvy1 - bvy2 - be1*bvy3
      crvz = crvz - be*bvz1 - bvz2 - be1*bvz3

      cf   = 0.25/arel(k)

c--------------------------------------
c one way to compute the mean curvature
c is to consider the norm of the contour
c integral around the edges:
c
c     crvx = cf*crvx
c     crvy = cf*crvy
c     crvz  = cf*crvz
c     crvmel(k) = sqrt(crvx**2+crvy**2+crvz**2)
c--------------------------------------

c---------------------------------------------
c another way to compute the mean curvature is to
c project the curvature vector onto the normal
c vector at the element centroid
c---------------------------------------------

      xi  = 1.0/3.0
      eta = 1.0/3.0

      call interp_p
     +
     +   (p(i1,1),p(i1,2),p(i1,3)
     +   ,p(i2,1),p(i2,2),p(i2,3)
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,p(i4,1),p(i4,2),p(i4,3)
     +   ,p(i5,1),p(i5,2),p(i5,3)
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +   ,alpha(k),beta(k),gamma(k)
     +   ,xi,eta
     +   ,x,y,z
     +   ,DxDxi,DyDxi,DzDxi
     +   ,DxDet,DyDet,DzDet
     +   ,vnx,vny,vnz
     +   ,hxi,het,hs
     +   ,Ichoose
     +   )

      crvmel(k) = cf*(crvx*vnx + crvy*vny + crvz*vnz)


c-------------

  1   Continue       ! end of loop over elements

c-----------------------------------
c average value the normal vector
c at the nodes
c-----------------------------------

      Do i=1,npts
        par = float(itally(i))
        vna(i,1) = vna(i,1)/par
        vna(i,2) = vna(i,2)/par
        vna(i,3) = vna(i,3)/par
      End Do

c-----------------------
c account for the images
c-----------------------

c---
c save the points that have images
c---

      Do 76 i=1,npts

       If(Iedge(i,1).eq.0) Go to 76  ! no images

       vna_sv(i,1) = vna(i,1)
       vna_sv(i,2) = vna(i,2)
       vna_sv(i,3) = vna(i,3)

   76 Continue

c---
c average over the images
c---

      Do 77 i=1,npts

       noim = Iedge(i,1)   ! number of images

       If(noim.eq.0) Go to 77  ! no images

        Do j=2,noim+1     ! Iedge(i,1) images
         image = Iedge(i,j)
         vna(i,1) = vna(i,1) + vna_sv(image,1)
         vna(i,2) = vna(i,2) + vna_sv(image,2)
         vna(i,3) = vna(i,3) + vna_sv(image,3)
        End Do

        noim1 = noim+1

        vna(i,1) = vna(i,1)/noim1
        vna(i,2) = vna(i,2)/noim1
        vna(i,3) = vna(i,3)/noim1

  77  Continue

c-------------------------------------
c Normalize the averaged normal vector
c to make its length equal to unity
c-------------------------------------

      Do i=1,npts
        par = sqrt(vna(i,1)**2+vna(i,2)**2+vna(i,3)**2)
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

c---
c Done
c---

      return
      end

c========================================================

      subroutine interp_p
     +
     +    (x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +    ,x4,y4,z4
     +    ,x5,y5,z5
     +    ,x6,y6,z6
     +    ,al,be,ga
     +    ,xi,eta
     +    ,x,y,z
     +    ,DxDxi,DyDxi,DzDxi
     +    ,DxDet,DyDet,DzDet
     +    ,vnx,vny,vnz
     +    ,hxi,het,hs
     +    ,Ichoose
     +    )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement
c=========================================

c-------------------------------------------------
c Interpolates over an element
c to compute following geometrical variables:
c
c   position vector
c   tangential vectors in the xi and eta directions
c   unit normal vector
c   the directional and surface metrics
c
c   Ichoose = 1 only the position vector
c             2 position vector and rest of variables
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c prepare
c---

      al1 = 1.0D0-al
      be1 = 1.0D0-be
      ga1 = 1.0D0-ga

      alal1 = al*al1
      bebe1 = be*be1
      gaga1 = ga*ga1

c---
c  interpolation functions
c---

      ph2 = xi *(xi -al+eta*(al-ga)/ga1)   /al1
      ph3 = eta*(eta-be+xi*(be+ga-1.0)/ga)/be1
      ph4 = xi *(1.0-xi-eta)/alal1
      ph5 = xi*eta/gaga1
      ph6 = eta*(1.0-xi-eta)/bebe1
      ph1 = 1.0-ph2-ph3-ph4-ph5-ph6

c---
c  position vector (x, y, z)
c---

      x = x1*ph1 + x2*ph2 + x3*ph3 + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3 + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3 + z4*ph4 + z5*ph5 + z6*ph6

      If(Ichoose.eq.1) Go to 99

c---
c  xi derivatives of phi
c---

      dph2 =  (2.0*xi-al +eta*(al-ga)/ga1)/al1
      dph3 =  eta*(be+ga-1.0)/(ga*be1)
      dph4 =  (1.0-2.0*xi-eta)/alal1
      dph5 =  eta/gaga1
      dph6 = -eta/bebe1
      dph1 = -dph2-dph3-dph4-dph5-dph6

c---
c  compute dx/dxi from xi derivatives of phi
c---

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c---
c  eta derivatives of phi
c---

      pph2 =  xi*(al-ga)/(al1*ga1)
      pph3 =  (2.0*eta-be+xi*(be+ga-1.0)/ga)/be1
      pph4 =  -xi/alal1
      pph5 =   xi/gaga1
      pph6 =  (1.0-xi-2.0*eta)/bebe1
      pph1 = - pph2 - pph3 - pph4 - pph5 - pph6

c---
c  compute dx/deta from eta derivatives of phi
c---

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6
      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6
      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

c---
c  normal vector    vn = (DxDxi)x(DxDeta)
c  surface metric   hs = norm(vn)
c---

      vnx = DyDxi * DzDet - DyDet * DzDxi
      vny = DzDxi * DxDet - DzDet * DxDxi
      vnz = DxDxi * DyDet - DxDet * DyDxi

      hs  = sqrt( vnx*vnx + vny*vny + vnz*vnz )

c---
c  normalizations
c---

c     hxi = sqrt( DxDxi**2+DyDxi**2+DzDxi**2 )
c     DxDxi = DxDxi/hxi
c     DyDxi = DyDxi/hxi
c     DzDxi = DzDxi/hxi

c     het = sqrt( DxDet**2+DyDet**2+DzDet**2 )
c     DxDet = DxDet/het
c     DyDet = DyDet/het
c     DzDet = DzDet/het

      vnx = vnx/hs
      vny = vny/hs
      vnz = vnz/hs

c---
c Done
c---

  99  Continue

      Return
      End
