      subroutine surfact_fem_interp
     +
     +  (x1,y1,z1
     +  ,x2,y2,z2
     +  ,x3,y3,z3
     +  ,x4,y4,z4
     +  ,x5,y5,z5
     +  ,x6,y6,z6
     +
     +  ,dilt1,dilt2,dilt3
     +  ,dilt4,dilt5,dilt6
     +
     +  ,crvm1,crvm2,crvm3
     +  ,crvm4,crvm5,crvm6
     +
     +  ,vnx1,vny1,vnz1
     +  ,vnx2,vny2,vnz2
     +  ,vnx3,vny3,vnz3
     +  ,vnx4,vny4,vnz4
     +  ,vnx5,vny5,vnz5
     +  ,vnx6,vny6,vnz6
     +
     +  ,ux1,uy1,uz1
     +  ,ux2,uy2,uz2
     +  ,ux3,uy3,uz3
     +  ,ux4,uy4,uz4
     +  ,ux5,uy5,uz5
     +  ,ux6,uy6,uz6
     +
     +  ,al,be,ga
     +  ,xi,eta
     +
     +  ,dilt
     +  ,crvm
     +  ,vnx,vny,vnz
     +  ,ux,uy,uz
     +  ,ph,gph
     +  ,hs
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------------------
c Utility of the fem
c
c Interpolate over an element for:
c
c (a) the element functions ph
c (b) the surface gradient
c     of the element functions gph
c (c) the surface metric: hs
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ph(6),gph(6,3)

c--------
c prepare
c--------

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c----------------------------
c compute the basis functions
c----------------------------

      ph2 = xi*(xi-al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi*(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

      ph(1) = ph1
      ph(2) = ph2
      ph(3) = ph3
      ph(4) = ph4
      ph(5) = ph5
      ph(6) = ph6

c-----------------------
c interpolate dilatation
c-----------------------

      dilt = dilt1*ph1 + dilt2*ph2 + dilt3*ph3
     +     + dilt4*ph4 + dilt5*ph5 + dilt6*ph6

c-----------------------
c interpolate curvature
c-----------------------

      crvm = crvm1*ph1 + crvm2*ph2 + crvm3*ph3
     +     + crvm4*ph4 + crvm5*ph5 + crvm6*ph6

c--------------------------
c interpolate normal vector
c--------------------------

      vnx = vnx1*ph1 + vnx2*ph2 + vnx3*ph3
     +    + vnx4*ph4 + vnx5*ph5 + vnx6*ph6

      vny = vny1*ph1 + vny2*ph2 + vny3*ph3
     +    + vny4*ph4 + vny5*ph5 + vny6*ph6

      vnz = vnz1*ph1 + vnz2*ph2 + vnz3*ph3
     +    + vnz4*ph4 + vnz5*ph5 + vnz6*ph6

c---------------------
c interpolate velocity
c---------------------

      ux = ux1*ph1 + ux2*ph2 + ux3*ph3
     +   + ux4*ph4 + ux5*ph5 + ux6*ph6

      uy = uy1*ph1 + uy2*ph2 + uy3*ph3
     +   + uy4*ph4 + uy5*ph5 + uy6*ph6

      uz = uz1*ph1 + uz2*ph2 + uz3*ph3
     +   + uz4*ph4 + uz5*ph5 + uz6*ph6

c--------------------------------------------------
c compute the xi derivatives of the basis functions
c--------------------------------------------------

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c---------------------------------------------------
c compute the eta derivatives of the basis functions
c---------------------------------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 = -xi/alalc
      pph5 =  xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c----------------------------------------
c compute the xi and eta derivatives of x
c----------------------------------------

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

c--------------------------------------
c compute the raw normal vector and the
c surface metric hs
c--------------------------------------

      vnxr = DyDxi * DzDet - DyDet * DzDxi
      vnyr = DzDxi * DxDet - DzDet * DxDxi
      vnzr = DxDxi * DyDet - DxDet * DyDxi

      hs = Dsqrt(vnxr**2+vnyr**2+vnzr**2)

c-----------------------------------
c compute the tangential derivatives
c of the interpolation functions
c by solving three linear equations:
c
c dx/dxi . grad = df/dxi
c dx/det . grad = df/det
c vn     . grad = 0.0
c-----------------------------------

      A11 = DxDxi
      A12 = DyDxi
      A13 = DzDxi

      A21 = DxDet
      A22 = DyDet
      A23 = DzDet

      A31 = vnx
      A32 = vny
      A33 = vnz

      B1  = dph1
      B2  = pph1
      B3  = 0.0D0

      call cramer_33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,gph(1,1),gph(1,2),gph(1,3)
     +  )

      B1  = dph2
      B2  = pph2
      B3  = 0.0D0

      call cramer_33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,gph(2,1),gph(2,2),gph(2,3)
     +  )

      B1  = dph3
      B2  = pph3
      B3  = 0.0D0

      call cramer_33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,gph(3,1),gph(3,2),gph(3,3)
     +  )

      B1  = dph4
      B2  = pph4
      B3  = 0.0D0

      call cramer_33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,gph(4,1),gph(4,2),gph(4,3)
     +  )

      B1  = dph5
      B2  = pph5
      B3  = 0.0D0

      call cramer_33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,gph(5,1),gph(5,2),gph(5,3)
     +  )

      B1  = dph6
      B2  = pph6
      B3  = 0.0D0

      call cramer_33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,gph(6,1),gph(6,2),gph(6,3)
     +  )

c-----
c rone
c-----

      return
      end
