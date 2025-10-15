      subroutine prtcl_ax_geo (Iflow)

c-----------------------------------------
c BEMLIB, FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------
c Collocation points:
c
c Compute:
c
c 1) coordinates, normal and tangential vector,
c 2) element surface areas: elar
c 3) host particle of collocation point
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(25),  Itp(25)
      Dimension xcntr(25),ycntr(25)
      Dimension axis1(25),axis2(25),tilt(25)

      Dimension xw(25,65),yw(25,65),tw(25,65)
      Dimension fx(25,64),fy(25,64)

      Dimension   X0(3200), Y0(3200),T0(3200)
      Dimension  ux0(3200), uy0(3200)
      Dimension vnx0(3200),vny0(3200),elar(3200)
      Dimension tnx0(3200),tny0(3200)
      Dimension Iprt(3200)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NGL,NE,Itp,Ncl
      common/points/xw,yw,tw
      common/particles/xcntr,ycntr,axis1,axis2,tilt

      common/colloc1/fx,fy,ux0,uy0
      common/colloc2/vnx0,vny0,elar
      common/colloc3/tnx0,tny0
      common/colloc4/x0,y0,t0
      common/colloc5/Iprt

      common/REAL1/visc,Uinf,wall,pg,RL,sc
      common/REAL2/Uprtcl

      common/piii/pi,pi2,pi4,pi8

c-------------------
c start construction
c-------------------

      Ic = 0     ! counter

      Do i=1,Nprtcl

c-------------------------
       If(Itp(i).eq.1) then                  ! straight segments
c-------------------------

        Do j=1,NE(i)

         Ic = Ic+1

         X0(Ic) = 0.5D0*(XW(i,j)+XW(i,j+1))
         Y0(Ic) = 0.5D0*(YW(i,j)+YW(i,j+1))
         XD = XW(i,J+1)-XW(i,J)
         YD = YW(i,J+1)-YW(i,J)
         ddl = dsqrt(XD*XD+YD*YD)
         elar(Ic) = pi2*ddl*Y0(Ic)     ! area of surface of revolution
         vnx0(Ic) = YD/elar(Ic)        ! normal vector into the flow
         vny0(Ic) =-XD/elar(Ic)
         tnx0(Ic) =-vny0(Ic)
         tny0(Ic) = vnx0(Ic)

         Iprt(Ic) = i         ! particle number
                              ! hosting Ic collocation point
        End Do

c----------
       Else               ! native elements of an ellipse
c----------

        cs = Dcos(tilt(i))
        sn = Dsin(tilt(i))

        Do j=1,NE(i)

         Ic = Ic+1

         t0(Ic) = 0.5D0*(TW(i,j)+TW(i,j+1))
         css    = Dcos(t0(Ic))
         snn    = Dsin(t0(Ic))
         tmpx   = axis1(i)*css
         tmpy   = axis2(i)*snn
         X0(Ic) = xcntr(i)+tmpx*cs-tmpy*sn    ! rotate to tilt
         Y0(Ic) = ycntr(i)+tmpx*sn+tmpy*cs

         tmpx     = axis2(i)*css
         tmpy     = axis1(i)*snn
         alm      = Dsqrt(tmpx**2+tmpy**2)  ! arc length metric
         elar(Ic) = alm*Dabs(TW(i,j+1)-TW(i,j))*pi2*Y0(Ic)

         tmpx = tmpx/alm
         tmpy = tmpy/alm

         vnx0(Ic) = tmpx*cs-tmpy*sn    ! rotate to tilt
         vny0(Ic) = tmpx*sn+tmpy*cs

         tnx0(Ic) = -vny0(Ic)
         tny0(Ic) =  vnx0(Ic)

         Iprt(Ic) = i         ! particle number
                              ! hosting Ic collocation point
        End Do

c------------
       End If
c------------

      End Do              ! over particles

c-----------------------------
c Compute disturbance velocity
c        at collocation points
c------------------------------

      Do i=1,Ncl

c--------------------------
        If(Iflow.eq.1) then   ! infinite uniform flow
c--------------------------

         ux0(i) = - Uinf + Uprtcl   ! incident and particle velocity
         uy0(i) =  0.0D0

c-------------------------------
        Else If(Iflow.eq.3) then  ! flow toward a wall
c-------------------------------

         ux0(i) = Uprtcl          ! particle velocity
         uy0(i) = 0.0D0

c-------------------------------
        Else If(Iflow.eq.5) then  ! flow in a circular tube
c-------------------------------

         ux0(i) = Uprtcl - 0.25D0*pg/visc * (sc*sc-Y0(i)**2)
         uy0(i) = 0.0D0

c-------------
        End If
c-------------

      End Do

c-----------------
c printing session
c-----------------
c
c     write (6,*)
c     write (6,*) " Position and disturbance velocity"
c     write (6,*) " at collocation points"
c     write (6,*)
c
c     Do i=1,Ncl
c      write (6,200) i,Iprt(i),X0(i),Y0(i),elar(i),vnx0(i),vny0(i)
c    +                ,tnx0(i),tny0(i)
c     End Do
c
c---
c     write (6,*)
c     write (6,*) " Disturbance velocity"
c     write (6,*) " at collocation points"
c     write (6,*)
c
c     Do i=1,Ncl
c      write (6,201) i,ux0(i),uy0(i)
c     End Do
c
c-----
c Done
c-----

  200 Format (1x,i3,1x,i3,10(f10.5))
  201 Format (1x,i3,10(f10.5))

      Return
      End
