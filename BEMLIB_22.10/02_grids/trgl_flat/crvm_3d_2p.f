      subroutine crvm_3d_2p
     +
     +   (nelm,npts
     +   )

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c----------------------------------------
c Compute the mean curvature at the nodes
c of a periodic grid
c
c  SYMBOLS:
c  --------
c
c  itally(i): number of elements sharing a node
c             for averaging
c
c  crvm(i): mean curvature at the ith node
c           averaged over host elements
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension       p(1090,3)
      Dimension      ne(1090,9)
      Dimension     vna(1090,3)
      Dimension    crvm(1090)
      Dimension crvm_sv(1090)
      Dimension  Itally(1090)     ! internal

      Dimension  Iedge(1090,4)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension  xxi(6), eet(6)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/edgepoints/Iedge
      common/albega/alpha,beta,gamma
      common/geo2/vna
      common/geo3/crvm

c-----------
c initialize
c-----------

      Do i=1,npts
          crvm(i) = 0.0D0
        itally(i) = 0    ! for averaging at a node
      End Do

c------------------

      Do 1 k=1,nelm    ! loop over elements

       i1 = n(k,1)     ! global node labels
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be = beta (k)
       ga = gamma(k)

c---
c triangle coordinates
c of the nodes
c---

       xxi(1) = 0.0D0
       eet(1) = 0.0D0
       xxi(2) = 1.0D0
       eet(2) = 0.0D0
       xxi(3) = 0.0D0
       eet(3) = 1.0D0
       xxi(4) = al
       eet(4) = 0.0D0
       xxi(5) = ga
       eet(5) = 1.0D0-ga
       xxi(6) = 0.0D0
       eet(6) = be

      Do i=1,6

        xi  = xxi(i)
        eta = eet(i)

        call interp_crvm_3d_2p
     +
     +     (p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +
     +     ,vna(i1,1),vna(i1,2),vna(i1,3)
     +     ,vna(i2,1),vna(i2,2),vna(i2,3)
     +     ,vna(i3,1),vna(i3,2),vna(i3,3)
     +     ,vna(i4,1),vna(i4,2),vna(i4,3)
     +     ,vna(i5,1),vna(i5,2),vna(i5,3)
     +     ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +     ,al,be,ga
     +     ,xi,eta
     +
     +     ,crvmn
     +     )

        m = n(k,i)   ! global label of local node i
                     ! on the kth element

        crvm(m) = crvm(m) + crvmn

        itally(m) = itally(m)+1

      End Do

  1   Continue       ! end of loop over elements

c---------------------------
c average the mean curvature
c at the nodes
c---------------------------

      Do i=1,npts
        par = float(itally(i))
        crvm(i) = crvm(i)/par
      End Do

c-----------------------
c account for the images
c-----------------------

c---
c save the points that have periodic images
c---

      Do i=1,npts
       if(    Iedge(i,1).eq.1
     +    .or.Iedge(i,1).eq.3
     +   ) then
         crvm_sv(i) = crvm(i)
       end if
      End Do

c---
c average over images
c---

      Do i=1,npts

       noim = Iedge(i,1)   ! number of images

       if(   noim.eq.1
     +   .or.noim.eq.3
     +   ) then 
        Do j=2,noim+1      ! run over images
          img = Iedge(i,j)
          crvm(i) = crvm(i) + crvm_sv(img)
        End Do
        crvm(i) = crvm(i)/(noim+1.0D0)
       end if

      End Do

c---
c Done
c---

      return
      end

c===================================================

      subroutine interp_crvm_3d_2p
     +
     +    (x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +    ,x4,y4,z4
     +    ,x5,y5,z5
     +    ,x6,y6,z6
     +
     +    ,vx1,vy1,vz1
     +    ,vx2,vy2,vz2
     +    ,vx3,vy3,vz3
     +    ,vx4,vy4,vz4
     +    ,vx5,vy5,vz5
     +    ,vx6,vy6,vz6
     +
     +    ,al,be,ga
     +    ,xi,eta
     +    ,crvm
     +    )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c-------------------------------------------
c  Computation of the mean curvature at the
c  element nodes
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c--------
c prepare
c--------

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c------------------------------------------
c compute xi derivatives of basis functions
c------------------------------------------

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c------------------------------------------
c compute eta derivatives of basis functions
c------------------------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 = -xi/alalc
      pph5 =  xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c----------------------------------
c compute xi derivatives of x and n
c----------------------------------

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6
      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6
      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

c------------------------------------
c compute xi and eta derivatives of n
c------------------------------------

      DvxDxi = vx1*dph1 + vx2*dph2 + vx3*dph3 + vx4*dph4
     +       + vx5*dph5 + vx6*dph6
      DvyDxi = vy1*dph1 + vy2*dph2 + vy3*dph3 + vy4*dph4
     +       + vy5*dph5 + vy6*dph6
      DvzDxi = vz1*dph1 + vz2*dph2 + vz3*dph3 + vz4*dph4
     +       + vz5*dph5 + vz6*dph6

      DvxDet = vx1*pph1 + vx2*pph2 + vx3*pph3 + vx4*pph4
     +       + vx5*pph5 + vx6*pph6
      DvyDet = vy1*pph1 + vy2*pph2 + vy3*pph3 + vy4*pph4
     +       + vy5*pph5 + vy6*pph6
      DvzDet = vz1*pph1 + vz2*pph2 + vz3*pph3 + vz4*pph4
     +       + vz5*pph5 + vz6*pph6

c-----------------------------------------------
c compute the first and second fundamental forms
c of the surface and the mean curvature
c-----------------------------------------------

      gxx = DxDxi**2 + DyDxi**2 + DzDxi**2
      gee = DxDet**2 + DyDet**2 + DzDet**2
      gxe = DxDxi*DxDet + DyDxi*DyDet + DzDxi*DzDet

      fxx = -(DxDxi*DvxDxi + DyDxi*DvyDxi + DzDxi*DvzDxi)
      fee = -(DxDet*DvxDet + DyDet*DvyDet + DzDet*DvzDet)
      fxe = -(DxDxi*DvxDet + DyDxi*DvyDet + DzDxi*DvzDet)

c-----------------------------------------------
c mean curvature is computed as the ratio of the
c first and the second fundamental form
c-----------------------------------------------

      crvm = -0.5 * (gxx*fee - 2.0*gxe*fxe + gee*fxx)
     +            / (gxx*gee - gxe*gxe)

c---
c Done
c---

  99  Continue

      return
      end
