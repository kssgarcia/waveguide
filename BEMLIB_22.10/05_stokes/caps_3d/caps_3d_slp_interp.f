      subroutine caps_3d_slp_interp
     +
     +    (Iopt_int
     +
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +    ,x4,y4,z4
     +    ,x5,y5,z5
     +    ,x6,y6,z6
     +
     +    ,dfx1,dfy1,dfz1
     +    ,dfx2,dfy2,dfz2
     +    ,dfx3,dfy3,dfz3
     +    ,dfx4,dfy4,dfz4
     +    ,dfx5,dfy5,dfz5
     +    ,dfx6,dfy6,dfz6
     +
     +    ,al,be,ga
     +    ,xi,eta
     +    ,x,y,z
     +    ,hs
     +    ,dfx,dfy,dfz
     +
     +    ,vx1,vy1,vz1
     +    ,vx2,vy2,vz2
     +    ,vx3,vy3,vz3
     +    ,vx4,vy4,vz4
     +    ,vx5,vy5,vz5
     +    ,vx6,vy6,vz6
     +
     +    ,tn1,tn2,tn3
     +    ,tn4,tn5,tn6
     +
     +    ,vx,vy,vz
     +    ,crvm
     +    ,tn
     +    ,gradtnx,gradtny,gradtnz
     +    )

c=============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=============================================

c-------------------------------------------
c   Utility of the slp integrator:
c
c   Interpolates over an element for:
c
c   position vector
c   normal vector
c   mean curvature (see Pozrikidis, 1997, TCFD, p. 169)
c   surface metric
c   surface tension
c   surface tension gradient
c  
c   Iopt_int = 1 only the position vector
c              2 position vector etc
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c prepare
c---

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c---
c compute the basis functions
c---

      ph2 = xi*(xi-al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi*(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c------------------------------------
c interpolate for the position vector
c------------------------------------

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

c----------------------------------------------
c compute the xi derivatives of basis functions
c----------------------------------------------

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

c-----------------------------------------------
c compute the eta derivatives of basis functions
c-----------------------------------------------

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

      rvnx = DyDxi * DzDet - DyDet * DzDxi
      rvny = DzDxi * DxDet - DzDet * DxDxi
      rvnz = DxDxi * DyDet - DxDet * DyDxi

      hs = sqrt(rvnx**2+rvny**2+rvnz**2 )

      if(Iopt_int.eq.1) Go to 99

c------------------------------
c interpolate the normal vector
c------------------------------

      vx = vx1*ph1 +vx2*ph2 +vx3*ph3 +vx4*ph4 +vx5*ph5 +vx6*ph6
      vy = vy1*ph1 +vy2*ph2 +vy3*ph3 +vy4*ph4 +vy5*ph5 +vy6*ph6
      vz = vz1*ph1 +vz2*ph2 +vz3*ph3 +vz4*ph4 +vz5*ph5 +vz6*ph6

c------------------------
c interpolate the tension
c------------------------

      tn = tn1*ph1 +tn2*ph2 +tn3*ph3 +tn4*ph4
     +   + tn5*ph5 +tn6*ph6

c-----------------------------------------
c compute the xi and eta derivatives of n
c-----------------------------------------

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

c--------------------------------------------------
c compute the xi and eta derivatives of the tension
c--------------------------------------------------

      DtnDxi = tn1*dph1 + tn2*dph2 + tn3*dph3 + tn4*dph4
     +       + tn5*dph5 + tn6*dph6

      DtnDet = tn1*pph1 + tn2*pph2 + tn3*pph3 + tn4*pph4
     +       + tn5*pph5 + tn6*pph6

c-----------------------------------------
c compute the coefficients of the first and
c second fundamental
c form of the surface and the mean curvature
c-------------------------------------------

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
     +              / (gxx*gee-gxe**2)

c--------------------------------------------
c compute the surface gradient of the tension
c by solving a 3x3 system
c--------------------------------------------

      A11 = DxDxi
      A12 = DyDxi
      A13 = DzDxi
      A21 = DxDet
      A22 = DyDet
      A23 = DzDet

      A31 = vx
      A32 = vy
      A33 = vz

      B1 = DtnDxi
      B2 = DtnDet
      B3 = 0.0

      Det =  A11*( A22*A33-A23*A32 )
     +     - A12*( A21*A33-A23*A31 )
     +     + A13*( A21*A32-A22*A31 )

      Det1 =  B1*( A22*A33-A23*A32 )
     +     - A12*(  B2*A33-A23*B3  )
     +     + A13*(  B2*A32-A22*B3  )

      Det2 = A11*( B2 *A33-A23*B3  )
     +     -  B1*( A21*A33-A23*A31 )
     +     + A13*( A21* B3-B2 *A31 )

      Det3 = A11*( A22* B3-A32* B2 )
     +     - A12*( A21* B3-A31* B2 )
     +     +  B1*( A21*A32-A22*A31 )

      gradtnx = Det1/Det
      gradtny = Det2/Det
      gradtnz = Det3/Det

c-----------
c compute Df
c-----------

      cf = 2.0*crvm*tn

      dfx = cf*vnx - gradtnx
      dfy = cf*vny - gradtny
      dfz = cf*vnz - gradtnz

c---
c done
c---

  99  Continue

      return
      end
