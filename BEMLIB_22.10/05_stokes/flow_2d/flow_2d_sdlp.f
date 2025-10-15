      subroutine flow_2d_sdlp
     +
     +  (Iflow
     +  ,X0,Y0,T0
     +  ,X1,Y1,T1
     +  ,X2,Y2,T2
     +  ,NGL
     +  ,Ising
     +  ,Itype
     +  ,Rad,xcnt,ycnt
     +  ,Qxx,Qxy
     +  ,Qyx,Qyy
     +  ,Wxx,Wyx
     +  ,Wxy,Wyy
     +  )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------------------
c Computes the single-layer and double-layer potential
c over a straight segment or a circular arc
c
c SYMBOLS:
c -------
c
c Qij:	components of the slp
c Wij:	components of the dlp
c
c If Ising = 1 will desingularize the integrals
c
c RL: cylinder separation for Iflow = 91, 92
c
c-----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

      common/flow_91/RL,Uslip

      common/ZZWW/ZZ,WW

      common/ppii/pi,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      Null = 0

c----------
c set flags
c----------

      Iopt = 2      ! default for sgf

      If(Iflow.eq.71) Iopt = 1   ! only GF for velocity needed
      If(Iflow.eq.81) Iopt = 1
      If(Iflow.eq.91) Iopt = 1 
      If(Iflow.eq.92) Iopt = 1 

      walll = 0.0D0           ! for sgf_2d_w

c-----------
c initialize
c-----------

      Qxx = 0.0D0
      Qxy = 0.0D0
      Qyx = 0.0D0
      Qyy = 0.0D0

      Wxx = 0.0D0
      Wxy = 0.0D0
      Wyx = 0.0D0
      Wyy = 0.0D0

c---
c prepare for the quadrture
c---

      If(Itype.eq.1) then

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = sqrt(XD*XD+YD*YD)
        vnx = -YD/DR           ! unit normal vector
        vny =  XD/DR           ! points into the flow

      Else

        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = Rad*Dabs(TD)

        ornt = 1.0D0            ! unit normal vector index
        If(TD.lt.0) ornt = -1.0D0

      End If

c--------------------------
c loop over Gaussian points
c--------------------------

      Do 1 I=1,NGL

        If(Itype.eq.1) then
          X = XM + XD*ZZ(I)
          Y = YM + YD*ZZ(I)
        Else
          T  = TM + TD*ZZ(I)
          cs = Dcos(t)
          sn = Dsin(t)
          X  = xcnt + Rad*cs
          Y  = ycnt + Rad*sn
          vnx = -cs * ornt    ! unit normal vector
          vny = -sn * ornt    ! points into the flow
        End If

c-------------------------------------
c call the free-space Green's function
c-------------------------------------

        If   (Iflow.eq.1 
     +    .or.Iflow.eq.11
     +    .or.Iflow.eq.41
     +    .or.Iflow.eq.51
     +       ) then

        call sgf_2d_fs 
     +
     +      (Iopt
     +      ,X,Y
     +      ,X0,Y0
     +      ,Gxx,Gxy
     +      ,Gyx,Gyy
     +      ,px,py
     +      ,Txxx,Txxy,Tyxx,Tyxy
     +      ,Txyx,Txyy,Tyyx,Tyyy
     +      )

c--------------------------------
        Else If(Iflow.eq.71.or.Iflow.eq.81) then
c--------------------------------

         call sgf_2d_w
     +
     +       (Iopt
     +       ,x,y
     +       ,x0,y0
     +       ,walll
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

c-----------------------------------------------
        Else If(Iflow.eq.91.or.Iflow.eq.92) then
c-----------------------------------------------

         call sgf_2d_1p
     +
     +      (Iopt
     +      ,X,Y
     +      ,X0,Y0
     +      ,RL
     +      ,Null
     +      ,Gxx,Gxy
     +      ,Gyx,Gyy
     +      ,px,py
     +      ,Txxx,Txxy,Tyxx,Tyxy
     +      ,Txyx,Txyy,Tyyx,Tyyy
     +      )

        Gxx = Gxx - (y-y0)*pi2/RL

c-------------
        Else
c-------------

        write (6,*) " flow_2d_sdlp: sorry this option not available"
        stop

c-------------
        End If
c-------------

c---------------------------------------------
c subtract off the logarithmic singularity
c of the single-layer potential
c
c The double-layer potential will be computed
c analytically at the end
c--------------------------------------------

       If(Ising.eq.1) then

        If(Itype.eq.1) then
          Dists = (X-X0)**2+(Y-Y0)**2   ! straight segments
        Else If(Itype.eq.2)  then
          Dists = (Rad*(T0-T))**2       ! circular arcs
        End If

        DD = 0.5D0*Dlog(Dists)

        Gxx = Gxx + DD
        Gyy = Gyy + DD

       End If

c---
c proceed with the quadrature
c---

       cf = WW(I)

       Qxx = Qxx + Gxx*cf
       Qxy = Qxy + Gxy*cf
       Qyx = Qyx + Gyx*cf
       Qyy = Qyy + Gyy*cf

        If(Iopt.eq.2) then

          Wxx = Wxx + (Txxx*vnx+Txxy*vny) * cf
          Wxy = Wxy + (Txyx*vnx+Txyy*vny) * cf
          Wyx = Wyx + (Tyxx*vnx+Tyxy*vny) * cf
          Wyy = Wyy + (Tyyx*vnx+Tyyy*vny) * cf

        End If

   1  Continue     ! loop over Gaussian points

c---
c finish up
c---

      Qxx = Qxx * DR
      Qxy = Qxy * DR 
      Qyx = Qyx * DR
      Qyy = Qyy * DR

      Wxx = Wxx * DR
      Wyx = Wyx * DR
      Wxy = Wxy * DR
      Wyy = Wyy * DR

c----------------------------
c add singularity back to the
c single-layer integral
c
c analytical integration of the 
c double-layer potential
c for the free space Green's function
c----------------------------

      If(Ising.eq.1) then

       corr = - 2.0D0*DR*(DLOG(DR)-1.0D0)

       Qxx = Qxx + corr
       Qyy = Qyy + corr

       If(Iopt.eq.2) then

        If(Itype.eq.1) then  ! vanishes over straight segments
          Wxx = 0.0D0
          Wyx = 0.0D0
          Wxy = 0.0D0
          Wyy = 0.0D0
        Else If(Itype.eq.2) then
          Wxx = -sin(T2+T0)+T2 +sin(T1+T0)-T1
          Wyx = cos(T2+T0) - cos(T1+T0)
          Wxy = Wyx
          Wyy = sin(T2+T0)+T2 - sin(T1+T0)-T1
        End If

       End If

      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
