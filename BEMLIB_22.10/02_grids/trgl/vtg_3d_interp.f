      subroutine vtg_3d_interp
     +
     +  (x1,y1,z1
     +  ,x2,y2,z2
     +  ,x3,y3,z3
     +  ,x4,y4,z4
     +  ,x5,y5,z5
     +  ,x6,y6,z6
     +
     +  ,vnx1,vny1,vnz1
     +  ,vnx2,vny2,vnz2
     +  ,vnx3,vny3,vnz3
     +  ,vnx4,vny4,vnz4
     +  ,vnx5,vny5,vnz5
     +  ,vnx6,vny6,vnz6
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
     +  ,vtg
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c---------------------------------------
c Facilitry for computing the tangential gradient
c at the element nodes
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension vtg(3,3)

c--------
c prepare
c--------

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c-----------------------------
c  compute the basis functions
c-----------------------------

      ph2 = xi*(xi-al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi*(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c--------------------------
c compute the normal vector
c--------------------------

      vnx = vnx1*ph1 + vnx2*ph2 + vnx3*ph3 + vnx4*ph4
     +    + vnx5*ph5 + vnx6*ph6

      vny = vny1*ph1 + vny2*ph2 + vny3*ph3 + vny4*ph4
     +    + vny5*ph5 + vny6*ph6

      vnz = vnz1*ph1 + vnz2*ph2 + vnz3*ph3 + vnz4*ph4
     +    + vnz5*ph5 + vnz6*ph6

c------------------------------------------
c compute xi derivatives of basis functions
c------------------------------------------

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c----------------------------------
c compute xi derivatives of x and v
c----------------------------------

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

      DvxDxi = vx1*dph1 + vx2*dph2 + vx3*dph3 + vx4*dph4
     +       + vx5*dph5 + vx6*dph6
      DvyDxi = vy1*dph1 + vy2*dph2 + vy3*dph3 + vy4*dph4
     +       + vy5*dph5 + vy6*dph6
      DvzDxi = vz1*dph1 + vz2*dph2 + vz3*dph3 + vz4*dph4
     +       + vz5*dph5 + vz6*dph6

c-------------------------------------------
c compute eta derivatives of basis functions
c-------------------------------------------

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 = -xi/alalc
      pph5 =  xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

c-----------------------------------
c compute eta derivatives of x and n
c-----------------------------------

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6
      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6
      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

      DvxDet = vx1*pph1 + vx2*pph2 + vx3*pph3 + vx4*pph4
     +       + vx5*pph5 + vx6*pph6
      DvyDet = vy1*pph1 + vy2*pph2 + vy3*pph3 + vy4*pph4
     +       + vy5*pph5 + vy6*pph6
      DvzDet = vz1*pph1 + vz2*pph2 + vz3*pph3 + vz4*pph4
     +       + vz5*pph5 + vz6*pph6

c--------------------------
c solve three systems
c of three linear equations
c--------------------------

      a11 = DxDxi
      a12 = DyDxi
      a13 = DzDxi

      a21 = DxDet
      a22 = DyDet
      a23 = DzDet

      a31 = vnx
      a32 = vny
      a33 = vnz

      b1 = DvxDxi
      b2 = DvxDet
      b3 = 0.0D0

      call cramer33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,vtg(1,1)
     +  ,vtg(2,1)
     +  ,vtg(3,1)
     +  )

      b1 = DvyDxi
      b2 = DvyDet
      b3 = 0.0D0

      call cramer33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,vtg(1,2)
     +  ,vtg(2,2)
     +  ,vtg(3,2)
     +  )

      b1 = DvzDxi
      b2 = DvzDet
      b3 = 0.0D0

      call cramer33
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +  ,vtg(1,3)
     +  ,vtg(2,3)
     +  ,vtg(3,3)
     +  )

c-----
c done
c-----

  100 Format (20(1x,f10.5))

  99  Continue

      return
      end
