      subroutine vs_3d_2p_curl_interp
     +
     +   (x1,y1,z1
     +   ,x2,y2,z2
     +   ,x3,y3,z3
     +   ,x4,y4,z4
     +   ,x5,y5,z5
     +   ,x6,y6,z6
     +
     +   ,ptlx1,ptly1,ptlz1
     +   ,ptlx2,ptly2,ptlz2
     +   ,ptlx3,ptly3,ptlz3
     +   ,ptlx4,ptly4,ptlz4
     +   ,ptlx5,ptly5,ptlz5
     +   ,ptlx6,ptly6,ptlz6
     +
     +   ,al,be,ga
     +   ,xi,eta
     +
     +   ,Dptlxx,Dptlyx,Dptlzx
     +   ,Dptlxy,Dptlyy,Dptlzy
     +   ,Dptlxz,Dptlyz,Dptlzz
     +   )

c======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the lisencing agreement
c======================================


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
      ph4 = xi *(1.0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0-xi-eta)/bebec
      ph1 = 1.0-ph2-ph3-ph4-ph5-ph6

c-----------------------
c  xi derivatives of phi
c-----------------------

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
      
      DptlxDxi = ptlx1*dph1 + ptlx2*dph2 + ptlx3*dph3
     +         + ptlx4*dph4 + ptlx5*dph5 + ptlx6*dph6

      DptlyDxi = ptly1*dph1 + ptly2*dph2 + ptly3*dph3
     +         + ptly4*dph4 + ptly5*dph5 + ptly6*dph6

      DptlzDxi = ptlz1*dph1 + ptlz2*dph2 + ptlz3*dph3
     +         + ptlz4*dph4 + ptlz5*dph5 + ptlz6*dph6

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

      DptlxDet = ptlx1*pph1 + ptlx2*pph2 + ptlx3*pph3
     +         + ptlx4*pph4 + ptlx5*pph5 + ptlx6*pph6

      DptlyDet = ptly1*pph1 + ptly2*pph2 + ptly3*pph3
     +         + ptly4*pph4 + ptly5*pph5 + ptly6*pph6

      DptlzDet = ptlz1*pph1 + ptlz2*pph2 + ptlz3*pph3
     +         + ptlz4*pph4 + ptlz5*pph5 + ptlz6*pph6

c----------------------------------------
c  normal vector    vn = (DxDxi)x(DxDeta)
c  un-normalized
c---------------------------------------

      vnx = DyDxi * DzDet - DyDet * DzDxi
      vny = DzDxi * DxDet - DzDet * DxDxi
      vnz = DxDxi * DyDet - DxDet * DyDxi

c----------------------------
c compute frechet derivatives
c----------------------------

      A31 = vnx
      A32 = vny
      A33 = vnz

      A11 = DxDxi
      A12 = DyDxi
      A13 = DzDxi

      A21 = DxDet
      A22 = DyDet
      A23 = DzDet

      B1 = DptlxDxi
      B2 = DptlxDet
      B3 = 0.0

      call cramer_33
     +
     +   (A11,A12,A13
     +   ,A21,A22,A23
     +   ,A31,A32,A33
     +   ,B1,B2,B3
     +   ,Dptlxx,Dptlyx,Dptlzx
     +   )

      B1  = DptlyDxi
      B2  = DptlyDet

      call cramer_33
     +
     +    (A11,A12,A13
     +    ,A21,A22,A23
     +    ,A31,A32,A33
     +    ,B1,B2,B3
     +    ,Dptlxy,Dptlyy,Dptlzy
     +             )

      B1  = DptlzDxi
      B2  = DptlzDet

      call cramer_33
     +
     +    (A11,A12,A13
     +    ,A21,A22,A23
     +    ,A31,A32,A33
     +    ,B1,B2,B3
     +    ,Dptlxz,Dptlyz,Dptlzz
     +    )

c---
c done
c---

      return
      end
