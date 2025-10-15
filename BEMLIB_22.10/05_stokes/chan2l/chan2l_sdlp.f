      subroutine chan2l_sdlp
     +
     +  (RL
     +  ,h
     +  ,X0,Y0  ! evaluation point
     +  ,X1,Y1
     +  ,X2,Y2
     +  ,NGL
     +  ,Ising
     +  ,Qxx,Qxy
     +  ,Qyx,Qyy
     +  ,Wxx,Wyx
     +  ,Wxy,Wyy
     +  )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c-----------------------------------------------------
c Compute the single-layer and double-layer potentials
c over a straight segment
c
c SYMBOLS:
c -------
c
c Qij:	components of the slp
c Wij:	components of the dlp
c
c If Ising = 1 desingularize the integrals
c----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

      common/ZZWW/ZZ,WW

      common/ppii/pi,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      Null = 0

c----------
c for the Green's function 
c----------

      Iopt = 2      ! default for sgf
      IQPD = 0
      NP = 10

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
c prepare for the quadrature
c---

      XM = 0.5D0*(X2+X1)
      XD = 0.5D0*(X2-X1)
      YM = 0.5D0*(Y2+Y1)
      YD = 0.5D0*(Y2-Y1)
      DR = SQRT(XD**2+YD**2)

      vnx = -YD/DR     ! unit normal vector
      vny =  XD/DR     ! points upward

c--------------------------
c loop over Gaussian points
c--------------------------

      Do i=1,NGL

         X = XM + XD*ZZ(i)
         Y = YM + YD*ZZ(i)

         call sgf_2d_1p_ww
     +
     +    (Iopt
     +    ,IQPD
     +    ,X,Y
     +    ,X0,Y0
     +    ,RL
     +    ,NP
     +    ,h
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

c------------------------------------------
c subtract off the logarithmic singularity
c of the single- and double-layer potential
c------------------------------------------

       if(Ising.eq.1) then

        Dists = (X-X0)**2+(Y-Y0)**2 

        DD = 0.5*Dlog(Dists)

        Gxx = Gxx + DD
        Gyy = Gyy + DD

       end if

c---
c proceed with the quadrature
c---

       cf = WW(I)

       Qxx = Qxx + Gxx * cf
       Qxy = Qxy + Gxy * cf
       Qyx = Qyx + Gyx * cf
       Qyy = Qyy + Gyy * cf

        if(Iopt.eq.2) then

          if(Ising.eq.1) then

          call sgf_2d_fs
     +
     +     (Iopt
     +     ,X,Y
     +     ,X0,Y0
     +     ,Gxx,Gxy
     +     ,Gyx,Gyy
     +     ,px,py
     +     ,Txxx_fs,Txxy_fs,Tyxx_fs,Tyxy_fs
     +     ,Txyx_fs,Txyy_fs,Tyyx_fs,Tyyy_fs
     +     )

          Txxx = Txxx-Txxx_fs
          Tyxx = Tyxx-Tyxx_fs
          Txyx = Txyx-Txyx_fs
          Txxy = Txxy-Txxy_fs
          Tyyx = Tyyx-Tyyx_fs
          Tyxy = Tyxy-Tyxy_fs
          Txyy = Txyy-Txyy_fs
          Tyyy = Tyyy-Tyyy_fs

          end if

          Wxx = Wxx + (Txxx*vnx+Txxy*vny) * cf
          Wxy = Wxy + (Txyx*vnx+Txyy*vny) * cf
          Wyx = Wyx + (Tyxx*vnx+Tyxy*vny) * cf
          Wyy = Wyy + (Tyyx*vnx+Tyyy*vny) * cf

        end if

       End Do     ! loop over Gaussian points

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
c add the singularity back to the
c single-layer integral
c
c analytical integration of the 
c double-layer potential
c for the free-space Green's function
c----------------------------

      if(Ising.eq.1) then

       corr = - 2.0D0*DR*(Dlog(DR)-1.0D0)

       Qxx = Qxx + corr
       Qyy = Qyy + corr

       if(Iopt.eq.2) then

          Wxx = Wxx + 0.0D0
          Wyx = Wyx + 0.0D0
          Wxy = Wxy + 0.0D0
          Wyy = Wyy + 0.0D0

       end if

      end if

c-----
c done
c-----

 100  Format (3(1x,f15.10))

      return
      end
