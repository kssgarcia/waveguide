      subroutine interp_p 
     +
     +  (Iopt
     +
     +  ,x1,y1,z1
     +  ,x2,y2,z2
     +  ,x3,y3,z3
     +  ,x4,y4,z4
     +  ,x5,y5,z5
     +  ,x6,y6,z6
     +
     +  ,al,be,ga
     +  ,xi,eta
     +
     +  ,x,y,z
     +  ,DxDxi,DyDxi,DzDxi
     +  ,DxDet,DyDet,DzDet
     +  ,vnx,vny,vnz
     +  ,hxi,het,hs
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

c-------------------------------------------------
c interpolate over an element
c to compute following geometrical variables:
c
c (a) position vector
c (b) tangential vectors in the xi and eta directions
c (c) unit normal vector
c (d) directional and surface metrics
c
c Iopt = 1 only the position vector
c        2 position vector and rest
c--------------------------------------------------

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

c----------------
c basis functions
c----------------

      ph2 = xi *(xi -al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi *(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c--------------------------
c position vector (x, y, z)
c--------------------------

      x = x1*ph1 + x2*ph2 + x3*ph3 + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3 + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3 + z4*ph4 + z5*ph5 + z6*ph6

      if(Iopt.eq.1) Go to 99

c----------------------
c xi derivatives of phi
c----------------------

      dph2 =  (2.0D0*xi-al +eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c------------------------------------------
c compute dx/dxi from xi derivatives of phi 
c------------------------------------------

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c-----------------------
c eta derivatives of phi
c-----------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 =  -xi/alalc
      pph5 =   xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c--------------------------------------------
c compute dx/deta from eta derivatives of phi 
c--------------------------------------------

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6
      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6
      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

c----------------------------------------
c Compute:
c
c normal vector:    vn = (DxDxi)x(DxDeta) 
c surface metric:   hs = norm(vn)
c----------------------------------------

      vnx = DyDxi * DzDet - DyDet * DzDxi
      vny = DzDxi * DxDet - DzDet * DxDxi
      vnz = DxDxi * DyDet - DxDet * DyDxi

      hs  = sqrt( vnx**2 + vny**2 + vnz**2 )

c----------------
c  normalizations
c----------------

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

c-----
c done
c-----

  99  Continue

      return
      end
