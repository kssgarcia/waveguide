      subroutine body_ax_sdlp
     +
     +   (X0,Y0,T0
     +   ,X1,Y1,T1
     +   ,X2,Y2,T2
     +   ,NGL
     +   ,Ising
     +   ,Itype
     +   ,Rad,xcnt,ycnt
     +   ,QQQ
     +   ,WWW
     +   )

c=========================================
c FDLIB, BEMLIB, CFDLAB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=========================================

c----------------------------------------------------------
c Compute the single-layer and double-layer potential over
c a straight segment or circular arc
c
c LEGEND:
c -------
c
c QQQ: single-layer potential
c WWW: double-layer potential
c----------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

      HF = 0.5D0

c-----------
c initialize
c-----------

      Iopt = 2   ! for the Green's function

      QQQ = 0.0D0
      WWW = 0.0D0

c---------------------------
c prepare for the quadrature
c---------------------------

      If(Itype.eq.1) then     ! straight segments

        XM = HF*(X2+X1)
        XD = HF*(X2-X1)
        YM = HF*(Y2+Y1)
        YD = HF*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)

        vnx =  YD/DR           ! unit normal vector
        vny = -XD/DR

      Else                     ! circular arcs

        TM = HF*(T2+T1)
        TD = HF*(T2-T1)
        DR = Rad*abs(TD)
        ornt = 1.0D0             ! orientation index
        If(TD.lt.0) ornt = -1.0D0

      End If

c---
c loop over Gaussian points
c---

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
          vnx = cs * ornt  ! unit normal vector
          vny = sn * ornt  ! when arc is counter-clockwise, 
                           ! normal vector points away from center
        End If

        call lgf_ax_fs 
     +
     +    (Iopt
     +    ,X,Y
     +    ,X0,Y0
     +    ,G
     +    ,dGdx
     +    ,dGdy
     +    )

c--------------------------------------------------
c  treat the slp singularity
c
c  Subtract off
c  the logarithmic singularity corresponding to the
c  free-space Green's function
c
c  NOTE: The double-layer singularity is not treated
c
c--------------------------------------------------

      If(Ising.eq.1) then
        If(Itype.eq.1) Dists = (X-X0)**2+(Y-Y0)**2
        If(Itype.eq.2) Dists = ( Rad*(T0-T) )**2
        G  = G + Dlog(Dists)/pi4
      End If

      QQQ = QQQ + G * Y * WW(i)

      WWW = WWW + (dGdx*vnx+dGdy*vny) * Y * WW(i)

   1  Continue

c-------------------------
c finish up the quadrature
c-------------------------

      QQQ = QQQ * DR
      WWW = WWW * DR

c------------------------------------
c add slp singularity back to the slp
c------------------------------------

      If(Ising.eq.1) then
        QQQ = QQQ - 2.0D0*DR*(LOG(DR)-1.0D0)/pi2
      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
