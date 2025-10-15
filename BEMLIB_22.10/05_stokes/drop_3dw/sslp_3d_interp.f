      subroutine sslp_3d_interp 
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
     +  ,al,be,ga
     +  ,xi,eta
     +
     +  ,x,y,z
     +  ,dfx,dfy,dfz
     +  ,hs
     +  )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c-------------------------------------------
c Utility of the slp integrator:
c
c Interpolate over an element for:
c
c (a) position: x, y, z
c (b) traction jump: Df
c (c) surface metric: hs
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

c----------------------------
c compute the basis functions
c----------------------------

      ph2 = xi*(xi-al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi*(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c-----------------------------
c Interpolate for the position
c-----------------------------

      x = x1*ph1 + x2*ph2 + x3*ph3 
     +  + x4*ph4 + x5*ph5 + x6*ph6

      y = y1*ph1 + y2*ph2 + y3*ph3
     +  + y4*ph4 + y5*ph5 + y6*ph6

      z = z1*ph1 + z2*ph2 + z3*ph3
     +  + z4*ph4 + z5*ph5 + z6*ph6

c-------------------
c Interpolate for df
c-------------------

      dfx = dfx1*ph1 + dfx2*ph2 + dfx3*ph3 
     +    + dfx4*ph4 + dfx5*ph5 + dfx6*ph6

      dfy = dfy1*ph1 + dfy2*ph2 + dfy3*ph3
     +    + dfy4*ph4 + dfy5*ph5 + dfy6*ph6

      dfz = dfz1*ph1 + dfz2*ph2 + dfz3*ph3
     +    + dfz4*ph4 + dfz5*ph5 + dfz6*ph6

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

      vnx = DyDxi * DzDet - DyDet * DzDxi
      vny = DzDxi * DxDet - DzDet * DxDxi
      vnz = DxDxi * DyDet - DxDet * DyDxi

      hs = Dsqrt(vnx*vnx +vny*vny +vnz*vnz)

c-----
c done
c-----

      return
      end
