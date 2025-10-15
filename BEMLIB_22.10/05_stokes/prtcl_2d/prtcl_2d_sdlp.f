      subroutine prtcl_2d_sdlp
     +
     +   (Iopt
     +   ,Iflow
     +   ,X0,Y0,T0     ! collocation point
     +   ,X1,Y1,T1
     +   ,X2,Y2,T2
     +   ,NGL
     +   ,Ising
     +   ,Itype
     +   ,axis1,axis2
     +   ,xcntr,ycntr,tilt
     +   ,Qxx,Qxy
     +   ,Qyx,Qyy
     +   ,Wxx,Wyx
     +   ,Wxy,Wyy
     +   )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------------------
c Computes the single- and double-layer potential
c at the collocation point X0, Y0
c over a straight segment 
c or a native element of an ellipse
c
c SYMBOLS:
c -------
c
c Qij:	components of the slp
c Wij:	components of the dlp
c
c If Iopt =  1 compute Qij but not Wij
c If Iopt ne 1 compute Qij and Wij
c
c-----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/REAL1/visc,shrt,delta,Uinf,wall
      common/REAL2/Uprtcl,Vprtcl,Aprtcl

      common/CHANNELI/IQPD,NGFWW
      common/CHANNELR/U1,U2,RL,h

      common/aaaa/a11,a12,a21,a22
      common/bbbb/b11,b12,b21,b22
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

      common/ZZWW/ZZ,WW

c-----------
c initialize
c-----------

      Qxx = 0.0D0   ! single-layer
      Qxy = 0.0D0
      Qyx = 0.0D0
      Qyy = 0.0D0

      Wxx = 0.0D0   ! double-layer
      Wxy = 0.0D0
      Wyx = 0.0D0
      Wyy = 0.0D0

c---------------------------
c prepare for the quadrature
c---------------------------

      If(Itype.eq.1) then      ! straight segments

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)

        vnx =  YD/DR           ! unit normal vector
        vny = -XD/DR           ! points into the flow
        fcc = DR               ! factor for numerical integration

      Else                     ! elliptical elements

        cs = cos(tilt)
        sn = sin(tilt)

        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        ornt = 1.0D0
        If(TD.lt.0) ornt = -1.0D0

        css  = Dcos(T0)           ! for desingularization
        snn  = Dsin(T0)
        tmpx = axis2*css
        tmpy = axis1*snn
        alm0 = Dsqrt(tmpx**2+tmpy**2)   ! arc length metric
        DR   = alm0*abs(TD)

      End If

c--------------------------
c loop over Gaussian points
c--------------------------

      Do 1 I=1,NGL

        If(Itype.eq.1) then
          X = XM + XD*ZZ(I)
          Y = YM + YD*ZZ(I)
        Else
          T    = TM + TD*ZZ(I)
          css  = cos(t)
          snn  = sin(t)
          tmpx = axis1*css
          tmpy = axis2*snn
          X    = xcntr + tmpx*cs-tmpy*sn
          Y    = ycntr + tmpx*sn+tmpy*cs
          tmpx = axis2*css
          tmpy = axis1*snn
          alm  = dsqrt(tmpx**2+tmpy**2)    ! arc length metric
          tmpx = tmpx/alm
          tmpy = tmpy/alm
          vnx  = tmpx*cs-tmpy*sn
          vny  = tmpx*sn+tmpy*cs
          vnx  = vnx * ornt    ! unit normal vector
          vny  = vny * ornt    ! points into the flow
          fcc  = alm*TD        ! factor for numerical integration
        End If

c---
        If(Iflow.eq.1) then
c---

        call sgf_2d_w
     +
     +    (Iopt
     +    ,x,y
     +    ,x0,y0
     +    ,wall
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

c---
        Else If(Iflow.eq.2) then
c---
        Ising_gf = 0

        call sgf_2d_1p_w
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,Wall
     +   ,RL
     +   ,Ising_gf
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

c---
        Else If(Iflow.eq.3) then
c---
         call sgf_2d_1p_ww
     +
     +   (Iopt
     +   ,IQPD
     +   ,X,Y
     +   ,X0,Y0
     +   ,RL,NGFWW,h
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

c---
        Else If(Iflow.eq.10) then
c---
         call sgf_2d_2p
     +
     +    (Iopt
     +    ,x,y
     +    ,x0,y0
     +    ,a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,tau
     +    ,Max1,Max2
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    ,px,py
     +    ,txxx,txxy,tyxx,tyxy
     +    ,txyx,txyy,tyyx,tyyy
     +    )

c------------
       End If
c------------

c------
c Subtract off the singularity
c------

      If(Ising.eq.1) then

         If(Itype.eq.1) then
           Dists = (X-X0)**2+(Y-Y0)**2
           DD    = 0.5D0*dlog(Dists)
         Else If(Itype.eq.2) then
           Dists = ( alm0*(T0-T) )**2
           DD    = 0.5D0*alm0/alm*dlog(Dists)
         End If

         Gxx = Gxx + DD
         Gyy = Gyy + DD

         If(Iopt.eq.2) then   ! subtract off the free-space GF

           call sgf_2d_fs
     +
     +       (Iopt
     +       ,x,y
     +       ,x0,y0
     +       ,Exx,Exy
     +       ,Eyx,Eyy
     +       ,zx,zy
     +       ,Rxxx,Rxxy,Ryxx,Ryxy
     +       ,Rxyx,Rxyy,Ryyx,Ryyy
     +       )

          Txxx = Txxx - Rxxx
          Txxy = Txxy - Rxxy
          Tyxx = Tyxx - Ryxx
          Tyxy = Tyxy - Ryxy
          Txyx = Txyx - Rxxx
          Txyy = Txyy - Rxyy
          Tyyx = Tyyx - Ryyx
          Tyyy = Tyyy - Ryyy

         End If

      End If

c-----------

      WI = WW(I)*fcc

      Qxx = Qxx + Gxx * WI
      Qxy = Qxy + Gxy * WI
      Qyx = Qyx + Gyx * WI
      Qyy = Qyy + Gyy * WI

      If(Iopt.eq.2) then
       Wxx = Wxx + (Txxx*vnx+Txxy*vny) * WI
       Wxy = Wxy + (Txyx*vnx+Txyy*vny) * WI
       Wyx = Wyx + (Tyxx*vnx+Tyxy*vny) * WI
       Wyy = Wyy + (Tyyx*vnx+Tyyy*vny) * WI
      End If

   1  Continue

c--------------------------------------
c Add the singularity back to the slp 
c
c Add the free-space dlp 
c noting that it vanishes over straight segments
c-----------------------------------------

      If(Ising.eq.1) then

        sing_cont =  2.0D0*DR*(LOG(DR)-1.0D0)

        Qxx = Qxx - sing_cont
        Qyy = Qyy - sing_cont

        If(Iopt.eq.2) then

         If(Itype.eq.2) then    ! caution: this only 
                                !           works for circular arcs

         Wxx = Wxx-sin(T2+T0)+T2+sin(T1+T0)-T1
         Wyx = Wyx+cos(T2+T0) - cos(T1+T0)
         Wxy = Wyx
         Wyy = Wyy+sin(T2+T0)+T2-sin(T1+T0)-T1

         Wxx = -Wxx 
         Wxy = -Wxy  ! invert the sign
         Wyx = -Wyx  ! because of the orientation
         Wyy = -Wyy  ! of the normal vector

         End If

        End If

      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
