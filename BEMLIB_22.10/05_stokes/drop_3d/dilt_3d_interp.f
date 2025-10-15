      subroutine dlt_3d_interp
     +
     +  (x1,y1,z1
     +  ,x2,y2,z2
     +  ,x3,y3,z3
     +  ,x4,y4,z4
     +  ,x5,y5,z5
     +  ,x6,y6,z6
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
     +  ,dilt
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c-------------------------------------------
c Computation of the mean curvature at the
c element nodes
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

c-------------------------------------------
c compute xi derivatives of basis functions
c-------------------------------------------

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c-------------------------------------------
c compute eta derivatives of basis functions
c-------------------------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 = -xi/alalc
      pph5 =  xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c------------------------------------
c compute xi and eta derivatives of x 
c------------------------------------

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
c compute xi and eta derivatives of u
c------------------------------------

      DuxDxi = ux1*dph1 + ux2*dph2 + ux3*dph3 + ux4*dph4
     +       + ux5*dph5 + ux6*dph6
      DuyDxi = uy1*dph1 + uy2*dph2 + uy3*dph3 + uy4*dph4
     +       + uy5*dph5 + uy6*dph6
      DuzDxi = uz1*dph1 + uz2*dph2 + uz3*dph3 + uz4*dph4
     +       + uz5*dph5 + uz6*dph6

      DuxDet = ux1*pph1 + ux2*pph2 + ux3*pph3 + ux4*pph4
     +       + ux5*pph5 + ux6*pph6
      DuyDet = uy1*pph1 + uy2*pph2 + uy3*pph3 + uy4*pph4
     +       + uy5*pph5 + uy6*pph6
      DuzDet = uz1*pph1 + uz2*pph2 + uz3*pph3 + uz4*pph4
     +       + uz5*pph5 + uz6*pph6

c-------------------------------------------
c compute  du/dxi x dx/det + dx/dxi x du/det
c-------------------------------------------

      crx = DuyDxi *  DzDet - DuzDxi *  DyDet
     +    +  DyDxi * DuzDet -  DzDxi * DuyDet

      cry = DuzDxi *  DxDet - DuxDxi *  DzDet
     +    +  DzDxi * DuxDet -  DxDxi * DuzDet

      crz = DuxDxi *  DyDet - DuyDxi *  DxDet
     +    +  DxDxi * DuyDet -  DyDxi * DuxDet

c----------------------------------------
c Compute:
c
c normal vector:    vn = (DxDxi)x(DxDeta)
c surface metric:   hs = norm(vn)
c                   hss = hs**2
c----------------------------------------

      vnx = DyDxi * DzDet - DyDet * DzDxi
      vny = DzDxi * DxDet - DzDet * DxDxi
      vnz = DxDxi * DyDet - DxDet * DyDxi

      hss = vnx**2 + vny**2 + vnz**2

      dilt = (crx*vnx + cry*vny + crz*vnz)/hss

c     write (6,*) " dilt_3d: ",crx,cry,crz,hss,dilt

c-----
c done
c-----

      return
      end
