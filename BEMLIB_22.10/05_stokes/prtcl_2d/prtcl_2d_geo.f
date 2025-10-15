      subroutine coll_points (Iflow)

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------
c Discretization of the particle contour
c and definition of collocation points
c
c Will compute:
c
c (1) Collocation points
c (2) Element arc lengths
c (3) Normal and tangential vector
c (4) Host particle of collocation points (Iprt)
c
c Capacity:
c --------
c
c 72 particles
c 64 elements around each particle
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(72),  Itp(72)
      Dimension axis1(72),axis2(72)
      Dimension xcntr(72),ycntr(72),tilt(72)

      Dimension xw(72,65),yw(72,65),tw(72,65)

      Dimension   X0(4608),  Y0(4608),T0(4608)
      Dimension  ux0(4608), uy0(4608)
      Dimension vnx0(4608),vny0(4608),elml(4608)
      Dimension tnx0(4608),tny0(4608)
      Dimension Iprt(4608)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl

      common/REAL1/visc,shrt,delta,Uinf,wall
      common/REAL2/Uprtcl,Vprtcl,Aprtcl

      common/CHANNELR/U1,U2,RL,h

      common/pax/axis1,axis2
      common/pap/xcntr,ycntr,tilt

      common/points/xw,yw,tw

      common/colloc1/ux0,uy0
      common/colloc2/vnx0,vny0,elml
      common/colloc3/tnx0,tny0
      common/colloc4/x0,y0,t0
      common/colloc5/Iprt

      common/piii/pi,pi2,pi4

c----------
c constants
c----------

      HF = 0.50D0

c-------------------------
c start the discretization
c-------------------------

      Ic = 0     ! counter

      Do i=1,Nprtcl

c-------------------------
       If(Itp(i).eq.1) then                  ! straight segments
c-------------------------

        Do j=1,NE(i)

         Ic = Ic+1

         X0(Ic) = HF*(XW(i,j)+XW(i,j+1))
         Y0(Ic) = HF*(YW(i,j)+YW(i,j+1))
         XD     = XW(i,J+1)-XW(i,J)
         YD     = YW(i,J+1)-YW(i,J)

         elml(Ic) = Dsqrt(XD*XD+YD*YD)  ! element length

         vnx0(Ic) = YD/elml(Ic)        ! normal vector into the flow
         vny0(Ic) =-XD/elml(Ic)

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

         t0(Ic) = HF*(TW(i,j)+TW(i,j+1))
         css    = Dcos(t0(Ic))
         snn    = Dsin(t0(Ic))
         tmpx   = axis1(i)*css
         tmpy   = axis2(i)*snn
         X0(Ic) = xcntr(i)+tmpx*cs-tmpy*sn    ! rotate to tilt
         Y0(Ic) = ycntr(i)+tmpx*sn+tmpy*cs

         tmpx     = axis2(i)*css
         tmpy     = axis1(i)*snn
         alm      = dsqrt(tmpx**2+tmpy**2)  ! arc length metric
         elml(Ic) = alm*abs(TW(i,j+1)-TW(i,j))

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

c------
       If(Iflow.eq.1.or.Iflow.eq.2) then   ! shear flow above a wall
c------

        yref   = Y0(i)-wall
        ux0(i) = Uprtcl
     +         - Aprtcl*(Y0(i)-ycntr(Iprt(i)))
     +         - shrt*yref - HF*delta*yref**2/visc

        uy0(i) = Vprtcl
     +         + Aprtcl*(X0(i)-xcntr(Iprt(i)))

c------
       Else If(Iflow.eq.3) then      ! channel flow
c------

         ux0(i) = Uprtcl
     +          - U1-(U2-U1)*(Y0(i)+h)/(2.0*h)
     +           -HF*delta*Y0(i)**2/visc
         uy0(i) = Vprtcl

c------
       Else If(Iflow.eq.10) then     ! uniform (streaming) flow
c------

         ux0(i) = -Uinf
         uy0(i) =  0.0D0

c------
       End If
c------

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
c      write (6,200) i,Iprt(i),X0(i),Y0(i),elml(i),vnx0(i),vny0(i)
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
