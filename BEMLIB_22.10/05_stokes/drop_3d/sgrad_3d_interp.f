      subroutine sgrad_3d_interp
     +
     +    (x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +    ,x4,y4,z4
     +    ,x5,y5,z5
     +    ,x6,y6,z6
     +
     +    ,f1,f2,f3,f4,f5,f6
     +
     +    ,al,be,ga
     +    ,xi,eta
     +    ,vnxa,vnya,vnza
     +    ,dfx,dfy,dfz
     +    )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------
c Compute the surface gradient
c of the scalar function f
c by interpolation
c------------------------------

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

c------------------------
c interpolation functions
c------------------------

      ph2 = xi *(xi -al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi *(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c---
c  xi derivatives of phi
c---

      dph2 =  (2.0D0*xi-al +eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c------------------------------------------
c  compute d/dxi from xi derivatives of phi
c------------------------------------------

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6

      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6

      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

      dfDxi = f1*dph1 + f2*dph2 + f3*dph3
     +      + f4*dph4 + f5*dph5 + f6*dph6

c------------------------
c  eta derivatives of phi
c------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 =  -xi/alalc
      pph5 =   xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = - pph2 - pph3 - pph4 - pph5 - pph6

c---------------------------------------------
c  compute dx/deta from eta derivatives of phi
c---------------------------------------------

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6

      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6

      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

      dfDet = f1*pph1 + f2*pph2 + f3*pph3
     +      + f4*pph4 + f5*pph5 + f6*pph6

c------------------------------
c  raw normal vector (non-unit):
c  vn = (DxDxi)x(DxDeta)
c------------------------------

      vnxr = DyDxi * DzDet - DyDet * DzDxi
      vnyr = DzDxi * DxDet - DzDet * DxDxi
      vnzr = DxDxi * DyDet - DxDet * DyDxi

c-----------------------------------
c compute the tangential derivatives
c by solving three linear equations:
c
c dx/dxi . grad = df/dxi
c dx/det . grad = df/det
c vn     . grad = 0.0
c-----------------------------------

      A11 = DxDxi
      A12 = DyDxi
      A13 = DzDxi
      B1  = dfDxi

      A21 = DxDet
      A22 = DyDet
      A23 = DzDet
      B2  = dfDet

      A31 = vnxa   ! averaged normal
      A32 = vnya
      A33 = vnza

c     A31 = vnxr   ! raw normal
c     A32 = vnyr
c     A33 = vnzr

      B3  = 0.0D0

      call cramer_33
     +
     +   (A11,A12,A13
     +   ,A21,A22,A23
     +   ,A31,A32,A33
     +   ,B1,B2,B3
     +   ,dfx,dfy,dfz
     +   )

c-----
c done
c-----

      return
      end
