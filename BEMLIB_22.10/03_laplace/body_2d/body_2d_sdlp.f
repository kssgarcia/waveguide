      subroutine body_2d_sdlp
     +
     +   (Iwall
     +   ,wall
     +   ,X0,Y0,T0
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

c-------------------------------------------
c Compute the single-layer and double-layer
c potential over a straight segment
c or circular arc
c 
c
c SYMBOLS:
c --------
c
c QQQ: single-layer
c WWW: double-layer
c-------------------------------------------

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

c-----------
c initialize
c-----------

      Iopt = 2                  ! for the Green's function
      if(Iwall.eq.1) Ign = 2    ! need the Neumann function

      QQQ = 0.0D0
      WWW = 0.0D0

c---
c prepare for the quadrature
c---

      if(Itype.eq.1) then     ! straight segments

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)
        vnx =  YD/DR           ! unit normal vector
        vny = -XD/DR

      else                     ! circular arcs

        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = Rad*abs(TD)
        ornt = 1.0D0             ! orientation index
        if(TD.lt.0) ornt = -1.0D0

      end if

c---
c loop over Gaussian points
c---

      Do 1 i=1,NGL

        if(Itype.eq.1) then     ! straight segments

          X = XM + XD*ZZ(i)
          Y = YM + YD*ZZ(i)

        else                    ! circular arcs

          T  = TM + TD*ZZ(i)
          cs = Dcos(t)
          sn = Dsin(t)
          X  = xcnt + Rad*cs
          Y  = ycnt + Rad*sn
          vnx = cs * ornt  ! unit normal vector
          vny = sn * ornt  ! when arc is counter-clockwise, 
                           ! normal vector points away from center
        end if

c-------
        if(Iwall.eq.0) then
c-------

          call lgf_2d_fs 
     +
     +       (Iopt
     +       ,X,Y
     +       ,X0,Y0
     +       ,G
     +       ,dGdx
     +       ,dGdy
     +       )

c-------
        else
c-------
          call lgf_2d_w 
     +
     +       (Iopt
     +       ,Ign
     +       ,X,Y
     +       ,X0,Y0
     +       ,wall
     +       ,G
     +       ,dGdx
     +       ,dGdy
     +       )

c-------
        end if
c-------

c---
c treat the slp singularity
c---

      if(Ising.eq.1) then

         if(Itype.eq.1) Dists = (X-X0)**2+(Y-Y0)**2
         if(Itype.eq.2) Dists = ( Rad*(T0-T) )**2

         G  = G + log(Dists)/pi4

         call lgf_2d_fs         ! subtract off the free-space kernel
     +
     +       (Iopt
     +       ,X,Y
     +       ,X0,Y0
     +       ,Gfs
     +       ,dGdxfs
     +       ,dGdyfs
     +       )

         dGdx = dGdx - dGdxfs
         dGdy = dGdy - dGdyfs

      End If

      WI = WW(I)

      QQQ = QQQ + G*WI

      WWW = WWW + (dGdx*vnx+dGdy*vny)*WI

   1  Continue

c----------
c finish up
c----------

      QQQ = QQQ * DR
      WWW = WWW * DR

c----------------------------------
c add the integrated sl singularity
c of the free-space GF
c back to the slp
c----------------------------------

      If(Ising.eq.1) then
        QQQ = QQQ - 2.0D0*DR*(LOG(DR)-1.0D0)/pi2
      End If

c--------------------------------------------
c analytical integration of the dlp
c for the free space GF over singular elements
c---------------------------------------------

      if(Ising.eq.1) then
        if(Itype.eq.1) then        ! straight segments
         WWW = WWW
        else if(Itype.eq.2) then   ! circular arcs
         WWW = WWW-(T2-T1)/pi4     ! independent of T0 !
        end if
      end if

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
