      subroutine flow_1d_osc_geo (Ncl)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------
c Boundary distribution
c
c Important:
c ----------
c
c elements must be distributed
c in the countercloskwise sense
c
c LEGEND:
c -------
c
c velamp: amplitude of the boundary velocity
c------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(129),Ye(129),Te(129),Se(129)
      Dimension Xm(128),Ym(128),Tm(128),Sm(128)

      Dimension NE(10),RT(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)

      Dimension   X0(500),  Y0(500),T0(500),S0(0:500),w0(500)
      Dimension tnX0(500),tnY0(500)
      Dimension vnX0(500),vnY0(500)
      Dimension elml(500)

c--------------
c common blocks
c--------------

      common/VEL00/Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw,w0
      common/VEL03/actis,xcntr,ycntr

      common/xxx04/visc,den,omega
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx07/elml
      common/xxx08/xwmin,ywmin,xwmax,ywmax
      common/xxx10/xcnt,ycnt

c----------
c constants
c----------

      pi = 3.1415 92653 58979 32384 D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

c------------------
c CIRCULAR TUBE
c
c Important:
c elements must be arranged in the countercloskwise sense
c-------------------

      If(Iflow.eq.1) then

      open (4,file='circle.dat')

        read (4,*) visc        ! viscosity
        read (4,*) den         ! density
        read (4,*) NGL         ! number of Gauss-Legendre points
        read (4,*) radius      ! tube radius
        read (4,*)
        read (4,*) velamp      ! amplitude of tube velocity
        read (4,*) omega       ! angular frequency
        read (4,*)
        read (4,*) NE(1),RT(1)

      close (4)

      NSG = 1      ! one segment only

      Ic    = 0       ! collocation points counter
      sinit = 0.0D0   ! initialize arc length

c---
c element distribution
c---

      Itp(1) = 2    ! elements are arcs

      actis(1) = radius

      angle1 = 0.0
      angle2 = pi2

      xstart = radius
      ystart = 0.0

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(1) = xcnt
      ycntr(1) = ycnt

      Isym = 0

      call elm_arc        ! circle
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,xcnt,ycnt
     +   ,radius
     +   ,angle1,angle2
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,Te,se
     +   ,Xm,Ym,Tm,sm
     +   )

      Do i=1,NE(1)+1
        XW(1,i) = Xe(i)
        YW(1,i) = Ye(i)
        TW(1,i) = Te(i)
      End Do

      Do i=1,NE(1)

        Ic = Ic + 1

        x0  (Ic) = xm(i)   ! collocation points
        y0  (Ic) = ym(i)
        t0  (Ic) = tm(i)
        S0  (Ic) = sm(i)

        elml(Ic) = se(i+1)-se(i)

        tnx0(Ic) =-sin(t0(Ic))    ! tangential vector
        tny0(Ic) = cos(t0(Ic))
        vnx0(Ic) =-tny0(Ic)       ! normal vector
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp

      End Do

      Ncl = Ic          ! number of collocation points

c------------------
c Elliptical tube
c
c Important:
c elements must be arranged in the countercloskwise sense
c-------------------

      Else If(Iflow.eq.2) then

      open (4,file='ellipse.dat',status='unknown')

        read (4,*) visc        ! viscosity
        read (4,*) den         ! density
        read (4,*) NGL
        read (4,*) aa
        read (4,*) bb
        read (4,*) 
        read (4,*) velamp      ! amplitude of tube velocity
        read (4,*) omega       ! angular frequency
        read (4,*) 
        read (4,*) NE(1)

      close (4)

c---
c preparations
c---

      NSG = 1       ! one segment

      Ic    = 0       ! collocation points counter
      sinit = 0.0D0   ! initialize arc length

      Itp(1) = 1      ! straight segments
      Isym   = 1

c---
c generate the boundary elements (straight segments)
c---

      step = pi2/NE(1)

      Do i=1,NE(1)+1
        tt = (i-1.0D0)*step
        xw(1,i) = aa*cos(tt)
        yw(1,i) = bb*sin(tt)
      End Do

      arcl = 0.0D0

      Do i=1,NE(1)

        Ic = Ic + 1

        x0(Ic) = 0.5*(xw(1,i)+xw(1,i+1))
        y0(Ic) = 0.5*(yw(1,i)+yw(1,i+1))
        ddx    = xw(1,i+1)-xw(1,i)
        ddy    = yw(1,i+1)-yw(1,i)
        ddl    = dsqrt(ddx*ddx+ddy*ddy)
        elml(Ic) = ddl     ! element arce length

        arcl   = arcl + 0.5*ddl
        S0(Ic) = arcl
        arcl   = arcl + 0.5*ddl

        tnx0(Ic) = ddx/ddl     ! tangential vector
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)    ! normal vector
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp      ! boundary-driven flow

      End Do

      Ncl = Ic          ! number of collocation points

c------------------
c Rectangular tube
c
c Important:
c elements must be arranged
c in the countercloskwise sense
c-------------------

      Else If(Iflow.eq.3) then

      open (4,file='rectangle.dat',status='unknown')

        read (4,*) visc        ! viscosity
        read (4,*) den         ! density
        read (4,*) NGL
        read (4,*) sizex
        read (4,*) sizey
        read (4,*)
        read (4,*) velamp      ! amplitude of tube velocity
        read (4,*) omega       ! angular frequency
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)
        read (4,*) NE(4),RT(4)

      close (4)

c---
c preparations
c---

      NSG = 4       ! four walls

      Ic    = 0        ! collocation points counter
      sinit = 0.0D0    ! initialize arc length

      sizexm = -sizex
      sizeym = -sizey

c---
c  side # 1 (right side)
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,sizex,sizeym
     +   ,sizex,sizey
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(1)+1
        xw(1,i) = Xe(i)
        yw(1,i) = Ye(i)
      End Do

      Do i=1,NE(1)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx*ddx+ddy*ddy)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp    ! boundary-driven flow

      End Do

      sinit = se(NE(1)+1)

c---
c  side # 2 (top side)
c---

      Itp(2) = 1    ! straight segment
      Isym   = 1

      call elm_line 
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,sizex ,sizey
     +   ,sizexm,sizey
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(2)+1
       xw(2,i) = Xe(i)
       yw(2,i) = Ye(i)
      End Do

      Do i=1,NE(2)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx*ddx+ddy*ddy)
        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp    ! boundary-driven flow

      End Do

      sinit = se(NE(2)+1)

c---
c  side # 3  (left side)
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,sizexm,sizey
     +   ,sizexm,sizeym
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
      End Do

      Do i=1,NE(3)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp    ! boundary-driven flow

      End Do

      sinit = se(NE(3)+1)

c---
c  side # 4   (bottom side)
c---

      Itp(4) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +   (NE(4)
     +   ,RT(4)
     +   ,sizexm,sizeym
     +   ,sizex ,sizeym
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(4)+1
        XW(4,i) = Xe(i)
        YW(4,i) = Ye(i)
      End Do

      Do i=1,NE(4)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx*ddx+ddy*ddy)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp    ! boundary-driven flow

      End Do

      Ncl = Ic          ! number of collocation points

c----------------
c Triangular tube
c
c Important:
c elements arranged in the counterclockwise sense
c-----------------

      Else If(Iflow.eq.4) then

      open (4,file='triangle.dat',status='unknown')

        read (4,*) visc        ! viscosity
        read (4,*) den         ! density
        read (4,*) NGL
        read (4,*) xfirst,yfirst
        read (4,*) xsecond,ysecond
        read (4,*) xthird,ythird
        read (4,*)
        read (4,*) velamp      ! amplitude of tube velocity
        read (4,*) omega       ! angular frequency
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)

      close (4)

c---
c preparations
c---

      NSG   = 3    ! three sides

      Ic    = 0       ! collocation point counter
      sinit = 0.0D0   ! initialize arc length

c---
c  side # 1
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,xfirst,yfirst
     +   ,xsecond,ysecond
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(1)+1
        XW(1,i) = Xe(i)
        YW(1,i) = Ye(i)
      End Do

      Do i=1,NE(1)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx*ddx+ddy*ddy)
        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp    ! boundary-driven flow

      End Do

      sinit = se(NE(1)+1)

c---
c  side # 2
c---

      Itp(2) = 1    ! straight segment
      Isym   = 1

      call elm_line 
     +
     +  (NE(2)
     +  ,RT(2)
     +  ,xsecond,ysecond
     +  ,xthird,ythird
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,se
     +  ,Xm,Ym,sm
     +  )

      Do i=1,NE(2)+1
        XW(2,i) = Xe(i)
        YW(2,i) = Ye(i)
      End Do

      Do i=1,NE(2)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx*ddx+ddy*ddy)
        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp    ! boundary-driven flow

      End Do

      sinit = se(NE(2)+1)

c---
c  side # 3
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1

      call elm_line 
     +
     +  (NE(3)
     +  ,RT(3)
     +  ,xthird,ythird
     +  ,xfirst,yfirst
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,se
     +  ,Xm,Ym,sm
     +  )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
      End Do

      Do i=1,NE(3)
        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)
        ddx    = xe(i+1)-xe(i)
        ddy    = ye(i+1)-ye(i)
        ddl    = dsqrt(ddx*ddx+ddy*ddy)
        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          w0(Ic) = velamp    ! boundary-driven flow

      End Do

      Ncl = Ic          ! number of collocation points

c-----------
      End If           ! End of geometry module
c-----------

c------------------------------------
c  End of boundary element generation 
c------------------------------------


c-----
c Done
c-----

 100  Format (1x,i3,20(1x,f15.10))
 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 107  Format (" Velocity: ",f10.5)
 150  Format (1X,10(1X,f10.5))

      Return
      End
