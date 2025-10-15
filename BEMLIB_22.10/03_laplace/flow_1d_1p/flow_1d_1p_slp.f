      subroutine flow_1d_1p_slp
     +
     +   (Iflow
     +   ,RL
     +   ,X0,Y0,T0
     +   ,X1,Y1,T1
     +   ,X2,Y2,T2
     +   ,NGL
     +   ,Ising
     +   ,Itype
     +   ,Rad,xcnt,ycnt
     +   ,QQQ            ! single-layer potential
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------
c Computation of the single-layer potential
c over a straight segment or circular arc
c
c SYMBOLS:
c -------
c
c QQQ: single-layer potential
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi = 3.1415 92653 58979 32384 D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c-----------
c initialize
c-----------

      Iopt = 1     ! need only G, not grad(G)

      QQQ = 0.0D0

c---
c prepare for the quadrture
c---

      If(Itype.eq.1) then     ! straight segments
        XM  = 0.5D0*(X2+X1)
        XD  = 0.5D0*(X2-X1)
        YM  = 0.5D0*(Y2+Y1)
        YD  = 0.5D0*(Y2-Y1)
        DR  = Dsqrt(XD**2+YD**2)
        vnx =  YD/DR           ! unit normal vector outward
        vny = -XD/DR
      Else                     ! circular arcs
        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = Rad*Dabs(TD)
        ornt = 1.0D0
        If(TD.lt.0D0) ornt = -1.0D0
      End If

c--------------------------
c loop over Gaussian points
c--------------------------

      Do 1 I=1,NGL

        If(Itype.eq.1) then
          X = XM + XD*ZZ(I)
          Y = YM + YD*ZZ(I)
        Else
          t   = TM + TD*ZZ(I)
          cs  = Dcos(t)
          sn  = Dsin(t)
          X   = xcnt + Rad*cs
          Y   = ycnt + Rad*sn
          vnx = cs * ornt    ! unit normal vector
          vny = sn * ornt    ! when arc is clockwise, 
                             ! normal vector points toward arc center
        End If

       call lgf_2d_1p
     +
     +     (Iopt
     +     ,x,y
     +     ,x0,y0
     +     ,RL
     +     ,G
     +     ,Gx,Gy
     +     )

        G = G-(y-y0)/(2.0D0*RL)     ! bias the Green's function upward

c---
c subtract off the slp singularity
c---

      If(Ising.eq.1) then
         If(Itype.eq.1) Dists = (X-X0)**2+(Y-Y0)**2
         IF(Itype.eq.2) Dists = ( Rad*(T0-T) )**2
         G = G + log(Dists)/pi4
      End If

      QQQ = QQQ + G * WW(I)

   1  Continue

c-------------------------
c finish up the quadrature
c-------------------------

      QQQ = QQQ * DR

c-----------------------------
c add back the slp singularity
c-----------------------------

      If(Ising.eq.1) then
        QQQ = QQQ - 2.0D0*DR*(LOG(DR)-1.0)/pi2
      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
