      subroutine rheology_interp 
     +
     +  (x1,y1,z1
     +  ,x2,y2,z2
     +  ,x3,y3,z3
     +  ,x4,y4,z4
     +  ,x5,y5,z5
     +  ,x6,y6,z6
     +
     +  ,dfx1,dfy1,dfz1
     +  ,dfx2,dfy2,dfz2
     +  ,dfx3,dfy3,dfz3
     +  ,dfx4,dfy4,dfz4
     +  ,dfx5,dfy5,dfz5
     +  ,dfx6,dfy6,dfz6
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
     +  ,x,y,z
     +  ,dfx,dfy,dfz
     +  ,vnx,vny,vnz
     +  ,ux,uy,uz
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

c-----------------------------------------
c Utility of the effective stress tensor 
c
c  Interpolates over an element to find the
c
c  (a) position vector
c  (b) traction discontinuity
c  (c) mean curvature (see TCFD, p. 169)
c  (d) velocity
c
c   iopt_int = 1 only the position vector
c              2 position vector and rest
c-----------------------------------------

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
c compute basis functions
c------------------------

      ph2 = xi *(xi -al+eta*(al-ga)/gac)   /alc
      ph3 = eta*(eta-be+xi *(be+ga-1.0D0)/ga)/bec
      ph4 = xi *(1.0D0-xi-eta)/alalc
      ph5 = xi*eta          /gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c------------------------
c interpolate x, df, n, u
c------------------------

      x = x1*ph1 + x2*ph2 + x3*ph3
     +  + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3
     +  + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3
     +  + z4*ph4 + z5*ph5 + z6*ph6

      dfx = dfx1*ph1 + dfx2*ph2 + dfx3*ph3
     +    + dfx4*ph4 + dfx5*ph5 + dfx6*ph6
      dfy = dfy1*ph1 + dfy2*ph2 + dfy3*ph3
     +    + dfy4*ph4 + dfy5*ph5 + dfy6*ph6
      dfz = dfz1*ph1 + dfz2*ph2 + dfz3*ph3
     +    + dfz4*ph4 + dfz5*ph5 + dfz6*ph6

      vnx = vnx1*ph1 + vnx2*ph2 + vnx3*ph3
     +    + vnx4*ph4 + vnx5*ph5 + vnx6*ph6
      vny = vny1*ph1 + vny2*ph2 + vny3*ph3
     +    + vny4*ph4 + vny5*ph5 + vny6*ph6
      vnz = vnz1*ph1 + vnz2*ph2 + vnz3*ph3
     +    + vnz4*ph4 + vnz5*ph5 + vnz6*ph6

      ux = ux1*ph1 + ux2*ph2 + ux3*ph3
     +   + ux4*ph4 + ux5*ph5 + ux6*ph6
      uy = uy1*ph1 + uy2*ph2 + uy3*ph3
     +   + uy4*ph4 + uy5*ph5 + uy6*ph6
      uz = uz1*ph1 + uz2*ph2 + uz3*ph3
     +   + uz4*ph4 + uz5*ph5 + uz6*ph6

c------------------------------------------
c compute xi derivatives of basis functions
c------------------------------------------

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c----------------------------
c compute xi derivatives of x
c----------------------------

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
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 =  -xi/alalc
      pph5 =   xi/gagac
      pph6 =  (1.0D0-xi-2.0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c-----------------------------
c compute eta derivatives of x
c-----------------------------

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

      rvnx = DyDxi * DzDet - DyDet * DzDxi
      rvny = DzDxi * DxDet - DzDet * DxDxi
      rvnz = DxDxi * DyDet - DxDet * DyDxi

      hs  = sqrt(rvnx**2 + rvny**2 + rvnz**2)

c-----
c done
c-----

      return
      end
