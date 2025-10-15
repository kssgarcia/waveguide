      subroutine prtcl_ax_slp
     +
     +   (Iflow
     +   ,X0,Y0,T0
     +   ,X1,Y1,T1
     +   ,X2,Y2,T2
     +   ,NGL
     +   ,Ising
     +   ,Itype
     +   ,xcntr,ycntr
     +   ,amaj,amin,tilt
     +   ,Qxx,Qxy
     +   ,Qyx,Qyy
     +   )

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
c Computation of the single-layer potential
c at collocation points located at the mid-point of
c a straight element or a native segment of an ellipse
c
c SYMBOLS:
c -------
c
c Qij:	components of the slp
c-----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/REAL1/visc,Uinf,wall,pg,RL,sc
      common/sgfaxct/Nsum,Np
      common/ZZWW/ZZ,WW

c--------
c prepare
c--------

      Iopt = 1

c-----------
c initialize
c-----------

      Qxx = 0.0D0
      Qxy = 0.0D0
      Qyx = 0.0D0
      Qyy = 0.0D0

c---
c prepare for the quadrture
c---

      if(Itype.eq.1) then             ! straight segments

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)
        vnx =  YD/DR           ! unit normal vector
        vny = -XD/DR           ! points into the flow
        fcc = DR               ! factor for numerical integration

      else                     ! elliptical elements

        cs = Dcos(tilt)
        sn = Dsin(tilt)
        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        ornt = 1.0D0
        if(TD.lt.0) ornt = -1.0D0

        css  = Dcos(T0)           ! for desingularization
        snn  = Dsin(T0)
        tmpx = amin*css
        tmpy = amaj*snn
        alm0 = Dsqrt(tmpx*tmpx+tmpy*tmpy)   ! arc length metric
        DR   = alm0*abs(TD)

      end if

c---
c loop over Gaussian points
c---

      Do 1 i=1,NGL

        if(Itype.eq.1) then
          X = XM + XD*ZZ(i)
          Y = YM + YD*ZZ(i)
        else
          T    = TM + TD*ZZ(i)
          css  = Dcos(T)
          snn  = Dsin(T)
          tmpx = amaj*css
          tmpy = amin*snn
          X    = xcntr + tmpx*cs-tmpy*sn
          Y    = ycntr + tmpx*sn+tmpy*cs
          tmpx = amin*css
          tmpy = amaj*snn
          alm  = Dsqrt(tmpx*tmpx+tmpy*tmpy)    ! arc length metric
          tmpx = tmpx/alm
          tmpy = tmpy/alm
          vnx  = tmpx*cs-tmpy*sn
          vny  = tmpx*sn+tmpy*cs
          vnx  = vnx * ornt    ! unit normal vector
          vny  = vny * ornt    ! points into the flow
          fcc  = alm*TD        ! factor for numerical integration
        end if

c--------------------------
        if(Iflow.eq.1) then      ! flow in free space
c--------------------------

        call sgf_ax_fs
     +
     +    (Iopt
     +    ,X,Y
     +    ,X0,Y0
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    ,QXXX,QXXY,QXYX,QXYY
     +    ,QYXX,QYXY,QYYX,QYYY
     +    ,PXX,PXY,PYX,PYY
     +    ,Iaxis
     +    )

c-------------------------------
        elseif(Iflow.eq.3) then     ! toward a wall
c-------------------------------

        call sgf_ax_w
     +
     +     (Iopt
     +     ,X,Y
     +     ,X0,Y0
     +     ,wall
     +     ,Gxx,Gxy
     +     ,Gyx,Gyy
     +     ,QXXX,QXXY,QXYX,QXYY
     +     ,QYXX,QYXY,QYYX,QYYY
     +     ,PXX,PXY,PYX,PYY
     +     ,Iaxis
     +     )

c-------------------------------
        elseif(Iflow.eq.5) then     ! flow in a circular tube
c-------------------------------

        call sgf_ax_1p_ct
     +
     +    (X0,Y0
     +    ,X,Y
     +    ,sc,RL
     +    ,Nsum,Np
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    )

c------------
       End If
c------------

c-------
c subtract off the singularity
c-------

      if(Ising.eq.1) then

         if(Itype.eq.1) then
           Dists = (X-X0)**2+(Y-Y0)**2
           DD    = 0.5D0*Dlog(Dists)
         else if(Itype.eq.2) then
           Dists = ( alm0*(T0-T) )**2
           DD    = 0.5D0*alm0/alm * Dlog(Dists)
         end if

         Gxx = Gxx + 2.0D0*DD
         Gyy = Gyy + 2.0D0*DD

      end if

c----------
c end of subtract off the singularity
c----------

      WI = WW(i)*fcc

      Qxx = Qxx + Gxx * WI
      Qxy = Qxy + Gxy * WI
      Qyx = Qyx + Gyx * WI
      Qyy = Qyy + Gyy * WI

   1  Continue

c--------------------------------
c add singularity back to the slp 
c--------------------------------

      If(Ising.eq.1) then

        Qxx = Qxx - 4.0D0* DR*(Dlog(DR)-1.0D0)
        Qyy = Qyy - 4.0D0* DR*(Dlog(DR)-1.0D0)

      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      return
      end
