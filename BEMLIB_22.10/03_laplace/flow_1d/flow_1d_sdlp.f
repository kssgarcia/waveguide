      subroutine flow_1d_sdlp
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

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------------------
c Computate the single-layer 
c and double-layer potential
c over a straight segment or a circular arc
c
c SYMBOLS:
c -------
c
c QQQ: single-layer potential
c WWW: double-layer potential
c----------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension ZZ(20),WW(20)

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi = 3.1415 92653 58979 32384 D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c---
c initialize
c---

      Iopt = 2   ! for the Green's function

      vni = 1.0D0  ! indicates the orientation of the normal vector

      QQQ = 0.0D0
      WWW = 0.0D0

c---
c prepare for the quadrture
c---

      if(Itype.eq.1) then     ! straight segments
        XM  = 0.5D0*(X2+X1)
        XD  = 0.5D0*(X2-X1)
        YM  = 0.5D0*(Y2+Y1)
        YD  = 0.5D0*(Y2-Y1)
        DR  = Dsqrt(XD**2+YD**2)
        vnx = -YD/DR           ! unit normal vector
        vny =  XD/DR
      else                     ! circular arcs
        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = rad*abs(TD)
        ornt = 1.0D0
        if(TD.lt.0) ornt = -1.0D0
      end if

c---
c loop over Gaussian points
c---

      Do 1 i=1,NGL

        if(Itype.eq.1) then
          X = XM + XD*ZZ(i)
          Y = YM + YD*ZZ(i)
        Else
          t  = TM + TD*ZZ(i)
          cs = cos(t)
          sn = sin(t)
          X  = xcnt + rad*cs
          Y  = ycnt + rad*sn
          vnx = -cs * ornt    ! unit normal vector
          vny = -sn * ornt    ! when arc is clockwise, 
c                             ! normal vector points toward center
        end if

          call lgf_2d_fs 
     +
     +       (Iopt
     +       ,X,Y
     +       ,X0,Y0
     +       ,G
     +       ,dGdx
     +       ,dGdy
     +       )

c---
c  treat slp singularity
c---

      if(Ising.eq.1) then
         If(Itype.eq.1) Dists = (X-X0)**2+(Y-Y0)**2
         IF(Itype.eq.2) Dists = ( rad*(T0-T) )**2
         G = G + log(Dists)/pi4
      end if

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
c add slp singularity back to the slp
c---

      if(Ising.eq.1) then
        QQQ = QQQ - 2.0D0*DR*(LOG(DR)-1.0D0)/pi2
      end if

c----------------------------
c Analytical integration of the dlp
c for the free space GF
c over the singular elements
c----------------------------

      if(Ising.eq.1) then

        if(Itype.eq.1) then        ! straight segments
         WWW = 0.0D0
        else if(Itype.eq.2) then   ! circular arcs
         WWW = vni*(T2-T1)/pi4     ! integral is independent of T0 !
        end if

      end if

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      return
      end
