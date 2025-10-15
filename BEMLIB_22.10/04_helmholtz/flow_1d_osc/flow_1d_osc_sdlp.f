      subroutine flow_1d_osc_sdlp
     +
     +   (X0,Y0,T0
     +   ,X1,Y1,T1
     +   ,X2,Y2,T2
     +   ,NGL
     +   ,Ising
     +   ,Itype
     +   ,delta
     +   ,Rad,xcnt,ycnt
     +   ,QQQr,QQQi        ! single layer
     +   ,WWWr,WWWi        ! double layer
     +   )

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------------------
c Computes the single-layer and double-layer potential 
c of the 2d complex Helmholtz operator 
c over a straight segment or a circular arc
c
c SYMBOLS:
c -------
c
c QQQr QQQi: real and imaginary parts
c            of the single-layer potential
c
c WWWr WWWi: real and imaginary parts
c            of the double-layer potential
c
c-----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi = 3.1415 92653 58979 32384  D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c-----------
c initialize
c-----------

      Iopt = 2   ! for the Green's function

c--------
c prepare
c--------

      QQQr = 0.0D0
      QQQi = 0.0D0
      WWWr = 0.0D0
      WWWi = 0.0D0

c----------------------------
c prepare for the quadrature
c---------------------------

      If(Itype.eq.1) then     ! straight segments

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)

        vnx = -YD/DR           ! unit normal vector
        vny =  XD/DR           ! pointing inward

      Else                     ! circular arcs

        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = Rad*abs(TD)
        ornt = 1.0D0
        If(TD.lt.0) ornt = -1.0D0

      End If

c---
c loop over Gaussian points
c---

      Do 1 I=1,NGL

        If(Itype.eq.1) then
          X = XM + XD*ZZ(I)
          Y = YM + YD*ZZ(I)
        Else
          t  = TM + TD*ZZ(I)
          cs = Dcos(t)
          sn = Dsin(t)
          X  = xcnt + Rad*cs
          Y  = ycnt + Rad*sn
          vnx = -cs * ornt    ! unit normal vector
          vny = -sn * ornt    ! when arc is clockwise, 
                              ! normal vector points toward center
        End If

        call hgf_2dc_fs
     +
     +    (Iopt
     +    ,delta
     +    ,X,Y
     +    ,X0,Y0
     +    ,Gr,Gi
     +    ,dGrdx,dGrdy
     +    ,dGidx,dGidy
     +    )

c----------------------- 
c  treat slp singularity: 
c  real part behaves like -ln(r)/(2 pi)
c  imaginary part is regular
c
c  treat dlp singularity
c-----------------------

      if(Ising.eq.1) then

         If(Itype.eq.1) then
            Dists = (X-X0)**2+(Y-Y0)**2
         Else If(Itype.eq.2) then 
            Dists = ( Rad*(T0-T) )**2
         End If

         Gr = Gr + log(Dists)/pi4

         dx = x-x0
         dy = y-y0
         rs = dx*dx+dy*dy
         dGrdx = dGrdx + dx/(pi2*rs)
         dGrdy = dGrdy + dy/(pi2*rs)

      end if

      WI = WW(I)

      QQQr = QQQr + Gr * WI
      QQQi = QQQi + Gi * WI

      WWWr = WWWr + (dGrdx*vnx+dGrdy*vny) * WI
      WWWi = WWWi + (dGidx*vnx+dGidy*vny) * WI

   1  Continue

c-------------------------
c finish up the quadrature
c-------------------------

      QQQr = QQQr * DR
      QQQi = QQQi * DR
      WWWr = WWWr * DR
      WWWi = WWWi * DR

c------------------------------------
c add slp singularity back to the slp
c------------------------------------

      if(Ising.eq.1) then
        QQQr = QQQr - 2.0D0*DR*(LOG(DR)-1.0)/pi2
      end if

c------------------------
c analytical integration
c of the laplace dlp
c for the free space GF
c over singular elements
c------------------------

      if(Ising.eq.1) then
        if(Itype.eq.2) then       ! circular arcs
         WWWr = WWWr+(T2-T1)/pi4  ! integral is independent of T0 !
        end if
      end if

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
