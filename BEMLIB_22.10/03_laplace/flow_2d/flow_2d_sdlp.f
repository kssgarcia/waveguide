      subroutine flow_2d_sdlp
     +
     +  (Iflow
     +  ,X0,Y0,T0
     +  ,X1,Y1,T1
     +  ,X2,Y2,T2
     +  ,NGL
     +  ,Ising
     +  ,Itype
     +  ,rad,xcnt,ycnt
     +  ,QQQ
     +  ,WWW
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c----------------------------------------------------------
c Computes the single-layer and double-layer potential over
c a straight segment or a circular arc
c
c SYMBOLS:
c -------
c
c QQQ: single-layer
c WWW: double-layer
c
c vni:	normal vector index to indicate side of the flow
c----------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)
      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c-----------
c initialize
c-----------

      Iopt = 2   ! for the Green's function

      QQQ = 0.0D0
      WWW = 0.0D0

      vni = 1.0D0

c---------------------------
c prepare for the quadrature
c---------------------------

      If(Itype.eq.1) then     ! straight segments

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)

        vnx = -YD/DR        ! unit normal vector
        vny =  XD/DR

      Else                     ! circular arcs

        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = rad*abs(TD)
        ornt = 1.0             ! orientation index
        If(TD.lt.0) ornt = -1.0

      End If

c---
c loop over Gaussian points
c---

      Do 1 i=1,NGL

        If(Itype.eq.1) then   ! straight segments
          X = XM + XD*ZZ(i)
          Y = YM + YD*ZZ(i)
        Else                  ! circular arcs
          T  = TM + TD*ZZ(i)
          cs = Dcos(t)
          sn = Dsin(t)
          X  = xcnt + rad*cs
          Y  = ycnt + rad*sn
          vnx = -cs * ornt  ! unit normal vector
          vny = -sn * ornt  ! when arc is clockwise, 
                            ! normal vector points toward the center
        End If

          call lgf_2d_fs 
     +
     +        (Iopt
     +        ,X,Y
     +        ,X0,Y0
     +        ,G
     +        ,dGdx
     +        ,dGdy
     +        )

c---
c treat the slp singularity
c---

      If(Ising.eq.1) then
         If(Itype.eq.1) Dists = (X-X0)**2+(Y-Y0)**2
         If(Itype.eq.2) Dists = ( rad*(T0-T) )**2
         G  = G + log(Dists)/pi4
      End If

      WI = WW(I)

      QQQ = QQQ + G * WI

      WWW = WWW + vni*(dGdx*vnx+dGdy*vny) * WI

   1  Continue

c---
c finish up
c---

      QQQ = QQQ * DR
      WWW = WWW * DR

c---
c add the slp singularity back to the slp
c---

      If(Ising.eq.1) then
        QQQ = QQQ - 2.0D0*DR*(LOG(DR)-1.0D0)/pi2
      End If

c----------------------------------
c Analytical integration of the dlp
c for the free space GF
c over the singular elements
c
c Note that the dlp of the free-space
c Green's function vanishes over
c singular straight segments
c----------------------------------

      If(Ising.eq.1) then

        If(Itype.eq.1) then        ! straight segments
         WWW = 0.0D0
        Else If(Itype.eq.2) then   ! circular arcs
         WWW = (T2-T1)/pi4         ! independent of T0 !
        End If

      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
