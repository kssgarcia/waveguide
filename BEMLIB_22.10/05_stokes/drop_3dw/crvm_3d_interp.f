      subroutine crvm_3d_interp
     +
     +  (x1,y1,z1
     +  ,x2,y2,z2
     +  ,x3,y3,z3
     +  ,x4,y4,z4
     +  ,x5,y5,z5
     +  ,x6,y6,z6
     +
     +  ,vx1,vy1,vz1
     +  ,vx2,vy2,vz2
     +  ,vx3,vy3,vz3
     +  ,vx4,vy4,vz4
     +  ,vx5,vy5,vz5
     +  ,vx6,vy6,vz6
     +
     +  ,al,be,ga
     +  ,xi,eta
     +  ,crvm
     +  )

c=========================================
c FDLIB
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

      fxx = DxDxi*DvxDxi + DyDxi*DvyDxi + DzDxi*DvzDxi
      fee = DxDet*DvxDet + DyDet*DvyDet + DzDet*DvzDet
      fxe = DxDxi*DvxDet + DyDxi*DvyDet + DzDxi*DvzDet

      fxx = - fxx
      fee = - fee
      fxe = - fxe

      crvm = -0.5D0 * (gxx*fee - 2.0D0*gxe*fxe + gee*fxx)
     +              / (gxx*gee - gxe**2)

c-----
c done
c-----

      return
      end
