      subroutine prtcl_sw_slp
     +
     +   (X0,Y0,T0
     +   ,X1,Y1,T1
     +   ,X2,Y2,T2
     +   ,NGL
     +   ,Ising
     +   ,Itype
     +   ,Rad,xcnt,ycnt
     +   ,QQQ
     +            )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------
c Compute the single-layer potential over
c a straight segment or a circular arc
c
c SYMBOLS:
c -------
c
c QQQ: Single-layer potential
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

      common/piii/pi,pih,pi2,pi4

      common/ZZWW/ZZ,WW

c-----------
c initialize
c-----------

      QQQ = 0.0D0

c--------------------------
c prepare for the quadrture
c--------------------------

      If(Itype.eq.1) then     ! straight segments

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)

      Else                     ! circular arcs

        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = Rad*abs(TD)
        ornt = 1.0D0             ! orientation index
        If(TD.lt.0) ornt = -1.0D0

      End If

c--------------------------
c loop over Gaussian points
c--------------------------

      Do 1 i=1,NGL

        If(Itype.eq.1) then

          X = XM + XD*ZZ(i)
          Y = YM + YD*ZZ(i)

        Else

          T  = TM + TD*ZZ(i)
          cs = cos(t)
          sn = sin(t)
          X  = xcnt + Rad*cs
          Y  = ycnt + Rad*sn

        End If

c---
c compute the Green function
c---

        Dx   = x-x0
        Dxs  = Dx*Dx
        yy0s = (y+y0)**2

        rks = 4.0D0*y*y0/(Dxs+yy0s)

        call ell_int (rks,F,E)

        RJ11 = ((2.0D0-rks)*F-2.0D0*E)/rks
        tmp1 = Dsqrt(Dxs+yy0s)
        RI11 = 4.0D0*RJ11/tmp1

        G = RI11*Y/pi4

c--- 
c treat slp singularity
c--- 

      If(Ising.eq.1) then
         If(Itype.eq.1) Dists = (X-X0)**2+(Y-Y0)**2
         If(Itype.eq.2) Dists = ( Rad*(T0-T) )**2
         G  = G + log(Dists)/pi4
      End If

      WI  = WW(I)
      QQQ = QQQ + G*WI

   1  Continue

c-------------------------
c finish up the quadrature
c-------------------------

      QQQ = QQQ * DR

c------------------------------------
c add slp singularity back to the slp
c------------------------------------

      If(Ising.eq.1) then
        QQQ = QQQ - 4.0*DR*(LOG(DR)-1.0)/pi4
      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
