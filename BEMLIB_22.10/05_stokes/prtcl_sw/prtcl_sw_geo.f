      subroutine prtcl_sw_geo 
     +                       (Iflow
     +                       ,Ncl)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------
c Boundary element generator
c
c LEGEND
c ------
c
c Ncl:    number of collocation points
c X0, Y0: collocation points
c
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(100),Ye(100),Se(100)
      Dimension Te(100)
      Dimension Xm(100),Ym(100),Sm(100)
c     Dimension Tm(100)

      Dimension NE(10),RT(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)

      Dimension    X0(500),Y0(500),T0(500),s0(0:500)
      Dimension arel(500)

c--------------
c common blocks
c--------------

      common/piii/pi,pih,pi2,pi4

      common/XF1/visc,Omega
      common/XF2/NSG,NE,Itp,NGL
      common/XF3/xw,yw,tw
      common/XF4/X0,Y0,T0,S0
      common/XF5/arel
      common/XF5/actis,xcntr,ycntr

c-------
c Sphere
c-------

      If(Iflow.eq.1) then

      open (4,file='sphere.dat')

        read (4,*) visc
        read (4,*) NGL
        read (4,*) rad
        read (4,*) xcenter
        read (4,*) Omega
        read (4,*)
        read (4,*) NE(1)    ! one segment consisting of arc elms

      close (4)

      ycenter = 0.0D0        ! sphere center is on the axis

      dth = pi/NE(1)

      Do i=1,NE(1)+1       ! points arranged in the countecl sense
        angle = (i-1.0D0)*dth
        te(i) = angle
        se(i) = angle*rad
        xe(i) = xcenter+rad*cos(angle)
        ye(i) =         rad*sin(angle)
      End Do

c---
c prepare
c---

      NSG = 1
      Ncl = NE(1)

      actis(1) = rad
      xcntr(1) = xcenter
      ycntr(1) = ycenter

c---
c  semi-circular contour
c---

      Itp(1) = 2    ! circular arcs

      Do i=1,NE(1)+1        ! boundary points on segment #1
       tw(1,i) = te(i)
       xw(1,i) = Xe(i)
       yw(1,i) = Ye(i)
      End Do

      Do i=1,NE(1)        ! collocation points

        t0(i) = 0.5D0*(te(i)+te(i+1))
        x0(i) = xcenter+rad*cos(t0(i))
        y0(i) = ycenter+rad*sin(t0(i))
        s0(i) = 0.5*(se(i)+se(i+1))

        arel(i) = dth*rad*pi2*y0(i)

      End Do

c---------
c Spheroid
c---------

      Else If(Iflow.eq.2) then

      open (4,file='spheroid.dat')

        read (4,*) visc
        read (4,*) NGL
        read (4,*) amaj
        read (4,*) amin
        read (4,*) xcenter
        read (4,*) Omega
        read (4,*)
        read (4,*) NE(1)    ! one segment consisting of linear elms

      close (4)

      ycenter = 0.0D0        ! spheroid center is on the axis

      dth = pi/NE(1)

      Do i=1,NE(1)+1       ! points arranged in the countecl sense
        angle = (i-1.0D0)*dth
        xe(i) = xcenter+amaj*cos(angle)
        ye(i) =         amin*sin(angle)
      End Do

c---
c Boundary and collocation points
c---

      Ncl = NE(1)
      NSG = 1

      Itp(1) = 1    ! straight segments

      Do i=1,NE(1)+1
       XW(1,i) = Xe(i)
       YW(1,i) = Ye(i)
      End Do

      Do i=1,NE(1)

       x0(i) = 0.5D0*(xe(i)+xe(i+1))
       y0(i) = 0.5D0*(ye(i)+ye(i+1))
       ddx = xe(i+1)-xe(i)
       ddy = ye(i+1)-ye(i)
       den = dsqrt(ddx**2+ddy**2)
       arel(i) = den*pi2*y0(i)

      End Do

c-----------------
c Triangular torus
c-----------------

c---------------
c Important:
c Segments must be arranged in the counterclockwise sense
c
c Normal vector points into the fluid
c-----------------

      Else If(Iflow.eq.3) then

      open (4,file='torus_trgl.dat')

        read (4,*) visc
        read (4,*) NGL
        read (4,*) xfirst,yfirst
        read (4,*) xsecond,ysecond
        read (4,*) xthird,ythird
        read (4,*) Omega
        read (4,*)
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)

      close (4)

c---
c preparations
c---

      NSG   = 3

      Ic    = 0       ! collocation point counter
      sinit = 0.0D0   ! initialize arc length

c--- side # 1

      Itp(1) = 1    ! straight segment
      Isym   = 1

      call elm_line
     +
     +         (NE(1)
     +         ,RT(1)
     +         ,xfirst,yfirst
     +         ,xsecond,ysecond
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

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
        den = sqrt(ddx**2+ddy**2)

        arel(Ic) = den*pi2*y0(Ic)

      End Do

      sinit = se(NE(1)+1)

c--- side # 2

      Itp(2) = 1    ! straight segment
      Isym   = 1

      call elm_line
     +
     +         (NE(2)
     +         ,RT(2)
     +         ,xsecond,ysecond
     +         ,xthird,ythird
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

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
        den = sqrt(ddx**2+ddy**2)

        arel(Ic) = den*pi2*y0(Ic)

      End Do

c--- side # 3

      sinit = se(NE(2)+1)

      Itp(3) = 1    ! straight segment
      Isym   = 1

      call elm_line
     +
     +         (NE(3)
     +         ,RT(3)
     +         ,xthird,ythird
     +         ,xfirst,yfirst
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

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
        den = sqrt(ddx**2+ddy**2)

        arel(Ic) = den*pi2*y0(Ic)

      End Do

      Ncl = Ic          ! number of collocation points

c-----------
      End If       ! End of geometry module
c-----------

c-----
c Done
c-----

      Return
      End
