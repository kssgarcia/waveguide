      subroutine crvmten6_interp
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
     +
     +  ,Kxx,Kxe
     +  ,Kex,Kee
     +
     +  ,crvten
     +  ,crvm
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------
c Compute the curvature tensor at the nodes
c of a 6-node triangle
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double precision Kxx,Kxe,Kex,Kee,Jac

      Dimension convx(3),conve(3),crvten(3,3)

c--------
c prepare
c--------

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c-------------------------------
c  interpolating basis functions
c-------------------------------

      ph2 = xi *(xi -al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0)/ga)/bec
      ph4 = xi *(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c----------------------------------
c xi derivatives of basis functions
c----------------------------------

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c-----------------------------------
c eta derivatives of basis functions
c-----------------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 = -xi/alalc
      pph5 =  xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c--------------------------
c position vector (x, y, z)
c--------------------------

      x = x1*ph1 + x2*ph2 + x3*ph3 + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3 + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3 + z4*ph4 + z5*ph5 + z6*ph6

c------------------------
c normal vector (x, y, z)
c------------------------

      vx = vx1*ph1 + vx2*ph2 + vx3*ph3 + vx4*ph4 + vx5*ph5 + vx6*ph6
      vy = vy1*ph1 + vy2*ph2 + vy3*ph3 + vy4*ph4 + vy5*ph5 + vy6*ph6
      vz = vz1*ph1 + vz2*ph2 + vz3*ph3 + vz4*ph4 + vz5*ph5 + vz6*ph6

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
      gex = gxe

      fxx = DxDxi*DvxDxi + DyDxi*DvyDxi + DzDxi*DvzDxi
      fee = DxDet*DvxDet + DyDet*DvyDet + DzDet*DvzDet
      fxe = DxDxi*DvxDet + DyDxi*DvyDet + DzDxi*DvzDet
      fex = DxDet*DvxDxi + DyDet*DvyDxi + DzDet*DvzDxi

      fxx = -fxx
      fee = -fee
      fxe = -fxe
      fex = -fex

c--
c triple scalar product
c--

      Jac = (DyDxi*DzDet - DzDxi*DyDet)*vx
     +    + (DzDxi*DxDet - DxDxi*DzDet)*vy
     +    + (DxDxi*DyDet - DyDxi*DxDet)*vz

c---
c contravariant base vectors
c---

      convx(1) = (DyDet*vz - DzDet*vy)/Jac
      convx(2) = (DzDet*vx - DxDet*vz)/Jac
      convx(3) = (DxDet*vy - DyDet*vx)/Jac

      conve(1) = (-DyDxi*vz + DzDxi*vy)/Jac
      conve(2) = (-DzDxi*vx + DxDxi*vz)/Jac
      conve(3) = (-DxDxi*vy + DyDxi*vx)/Jac

c---
c covariant curvature components
c---

      Kxx = DvxDxi*DxDxi + DvyDxi*DyDxi + DvzDxi*DzDxi
      Kxe = DvxDxi*DxDet + DvyDxi*DyDet + DvzDxi*DzDet
      Kex = DvxDet*DxDxi + DvyDet*DyDxi + DvzDet*DzDxi
      Kee = DvxDet*DxDet + DvyDet*DyDet + DvzDet*DzDet

      Do i=1,3
       Do j=1,3
        crvten(i,j) = Kxx*convx(i)*convx(j) + Kxe*convx(i)*conve(j)
     +              + Kex*conve(i)*convx(j) + Kee*conve(i)*conve(j)
       End Do
      End Do

c---
c mean curvature as the trace of the curvature tensor
c---

      crvm = 0.5*(crvten(1,1)+crvten(2,2)+crvten(3,3))

c---
c mean curvature otherwise
c---
 
      crvmother = -0.5D0*(gxx*fee - 2.0*gxe*fxe + gee*fxx)
     +               /(gxx*gee - gxe*gxe)

c-----
c done
c-----

  99  Continue

      return
      end
