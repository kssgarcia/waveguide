      subroutine flow_1d_geo
     +
     +    (Ncl
     +    )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------
c Boundary distribution
c
c Important:
c --------
c 
c elements must be distributed 
c in the countercloskwise sense
c------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(300),Ye(300),Te(300),Se(300)
      Dimension Xm(300),Ym(300),Tm(300),Sm(300)

      Dimension NE(10),RT(10),Itp(10),Actis(10)
      Dimension xcntr(10),ycntr(10)

      Dimension xw(10,300),yw(10,300),tw(10,300)

      Dimension   x0(900),  y0(900),t0(900),s0(0:900),u0(900)
      Dimension tnX0(900),tnY0(900)
      Dimension vnX0(900),vnY0(900)
      Dimension elml(900)

c--------------
c common blocks
c--------------

      common/VEL00/Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw,u0
      common/VEL03/actis,xcntr,ycntr

      common/xxx04/visc,pg,orgx,orgy
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx07/elml
      common/xxx08/xwmin,ywmin,xwmax,ywmax
      common/xxx10/xcnt,ycnt

      common/rectangle/sizex,sizey

c----------
c constants
c----------

      pi = 3.1415926535897932384 D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

c--------------
c Circular tube
c--------------

      if(Iflow.eq.1) then

      open (4,file='circle.dat',status='unknown')

        read (4,*) pg          ! negative of the pressure gradient
        read (4,*) visc        ! fluid viscosity
        read (4,*) NGL         ! number of Gauss Legendre points
        read (4,*) radius      ! tube radius
        read (4,*)
        read (4,*) orgx,orgy   ! shifted origin of particular solution
        read (4,*)
        read (4,*) NE(1),RT(1) ! number of elements, strech ratio

      close (4)

      if(radius.eq.1.0) then
        write (6,*) 
        write (6,*) " flow_1d_geo: the solution of the integral"
        write (6,*) "               is not unique"
        write (6,*) 
        write (6,*) "  Please rescale"
        stop
      end if

c---
c prepare
c---

      fc = 0.25D0*pg/visc

      NSG = 1      ! one segment only

      Ic    = 0    ! collocation points counter
      sinit = 0.0  ! initialize arc length

      Itp(1) = 2   ! circular elements

      actis(1) = radius

      angle1 = 0.0D0
      angle2 = pi2

      xstart = radius    ! first point of the circle
      ystart = 0.0D0

      xcnt = xstart-radius*cos(angle1)    ! circle center
      ycnt = ystart-radius*sin(angle1)

      xcntr(1) = xcnt
      ycntr(1) = ycnt

      Isym = 0

c---
c element distribution
c---

      call elm_arc        ! circle
     +
     +    (NE(1)
     +    ,RT(1)
     +    ,xcnt,ycnt
     +    ,radius
     +    ,angle1,angle2
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,Te,se
     +    ,Xm,Ym,Tm,sm
     +    )

      Do i=1,NE(1)+1
        XW(1,i) = Xe(i)
        YW(1,i) = Ye(i)
        TW(1,i) = Te(i)
      End Do

      Do i=1,NE(1)

        Ic = Ic + 1

        x0(Ic) = xm(i)   ! collocation points
        y0(Ic) = ym(i)
        t0(Ic) = tm(i)
        s0(Ic) = sm(i)

        elml(Ic) = se(i+1)-se(i)    ! element length

        tnx0(Ic) =-sin(t0(Ic))      ! tangential vector
        tny0(Ic) = cos(t0(Ic))
        vnx0(Ic) =-tny0(Ic)         ! normal vector (inward)
        vny0(Ic) = tnx0(Ic)
          U0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)
      End Do

      Ncl = Ic          ! number of collocation points

c----------------
c Elliptical tube
c
c discretized into straight segments
c----------------

      else if(Iflow.eq.2) then

      open(4,file='ellipse.dat',status='unknown')

        read (4,*) pg       ! negative of the pressure gradient
        read (4,*) visc     ! viscosity
        read (4,*) NGL      ! number of Gauss-Legendre points
        read (4,*) aa       ! major semi-axis
        read (4,*) bb       ! minor semi-axis
        read (4,*) 
        read (4,*) NE(1)

      close (4)

c---
c prepare
c---

      fc = 0.25D0*pg/visc

      orgx = 0.0D0    ! shifted origin of particular solution
      orgy = 0.0D0

      NSG = 1       ! one segment

      Ic = 0        ! collocation points counter
      sinit = 0.D0  ! initialize arc length

      Itp(1) = 1    ! straight segments

      step = pi2/NE(1)

c---
c generate the boundary elements (straight segments)
c---

      Do i=1,NE(1)+1
        tt = (i-1.0D0)*step
        xw(1,i) = aa*cos(tt)
        yw(1,i) = bb*sin(tt)
      End Do

      arcl = 0.0D0

      Do i=1,NE(1)
        Ic = Ic+1

        x0(Ic)   = 0.5D0*(xw(1,i)+xw(1,i+1))  ! collocation points
        y0(Ic)   = 0.5D0*(yw(1,i)+yw(1,i+1))
        ddx      = xw(1,i+1)-xw(1,i)
        ddy      = yw(1,i+1)-yw(1,i)
        ddl      = sqrt(ddx**2+ddy**2)
        elml(Ic) = ddl

        arcl   = arcl + 0.5D0*ddl
        s0(Ic) = arcl
        arcl   = arcl + 0.5D0*ddl

        tnx0(Ic) = ddx/ddl    ! tangential vector
        tny0(Ic) = ddy/ddl

        vnx0(Ic) =-tny0(Ic)   ! normal vector
        vny0(Ic) = tnx0(Ic)
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      Ncl = Ic          ! number of collocation points

c-----------------
c rectangular tube
c-----------------

      else if(Iflow.eq.3) then

      open (4,file='rectangle.dat',status='unknown')

        read (4,*) pg
        read (4,*) visc
        read (4,*) NGL
        read (4,*) sizex
        read (4,*) sizey
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)
        read (4,*) NE(4),RT(4)

      close (4)

c---
c prepare
c---

      fc = 0.25D0*pg/visc

      orgx = 0.0D0  ! shifted origin of particular solution
      orgy = 0.0D0

      NSG = 4       ! four walls

      Ic = 0        ! collocation points counter
      sinit = 0.0D0 ! initialize arc length

      sizexm = -sizex
      sizeym = -sizey

c---
c  side # 1 (right side)
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line
     +
     +    (NE(1)
     +    ,RT(1)
     +    ,sizex,sizeym
     +    ,sizex,sizey
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      sinit = se(NE(1)+1)

c---
c  side # 2 (top side)
c---

      Itp(2) = 1    ! straight segment
      Isym   = 1

      call elm_line 
     +
     +    (NE(2)
     +    ,RT(2)
     +    ,sizex ,sizey
     +    ,sizexm,sizey
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      sinit = se(NE(2)+1)

c---
c  side # 3  (left side)
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +    (NE(3)
     +    ,RT(3)
     +    ,sizexm,sizey
     +    ,sizexm,sizeym
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      Ncl = Ic          ! number of collocation points

c----------------
c Triangular tube
c----------------

      else if(Iflow.eq.4) then

      open (4,file='triangle.dat',status='unknown')

        read (4,*) pg
        read (4,*) visc
        read (4,*) NGL
        read (4,*) 
        read (4,*) xfirst,yfirst
        read (4,*) xsecond,ysecond
        read (4,*) xthird,ythird
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)

      close (4)

c---
c prepare
c---

      fc = 0.25D0*pg/visc

      NSG   = 3    ! three sides

      Ic    = 0      ! collocation point counter
      sinit = 0.D0   ! initialize arc length

c---
c shifted origin of particular solution
c placed at the centroid
c---

      orgx = (xfirst+xsecond+xthird)/3.0
      orgy = (yfirst+ysecond+ythird)/3.0

c---
c  side # 1
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line
     +
     +  (NE(1)
     +  ,RT(1)
     +  ,xfirst,yfirst
     +  ,xsecond,ysecond
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,se
     +  ,Xm,Ym,sm
     +  )

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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      sinit = se(NE(1)+1)

c---
c  side # 2
c---

      Itp(2) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line 
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,xsecond,ysecond
     +   ,xthird,ythird
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      sinit = se(NE(2)+1)

c---
c  side # 3
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line 
     +
     +    (NE(3)
     +    ,RT(3)
     +    ,xthird,ythird
     +    ,xfirst,yfirst
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          u0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      Ncl = Ic          ! number of collocation points

c----------------
c Read end-points from file: contour.dat
c
c discretize into straight segments
c----------------

      else if(Iflow.eq.10) then

      open (4,file='contour.dat',status='unknown')

c-------
c       read (4,*) pg       ! negative of the pressure gradient
c       read (4,*) visc     ! viscosity
c       read (4,*) NGL      ! number of Gauss-Legendre points
c-------

        pg   = 1.0D0   ! negative of the pressure gradient
        visc = 1.0D0   ! viscosity
        NGL  = 6       ! number of Gauss-Legendre points

        read (4,*) NE1

        Do i=1,NE1
          read (4,*) idle,xw(1,i),yw(1,i)
        End Do

        NE(1) = NE1-1

      close (4)

c---
c prepare
c---

      fc = 0.25D0*pg/visc

      orgx = 0.0D0    ! shifted origin of particular solution
      orgy = 0.0D0

      NSG = 1       ! one segment

      Ic = 0         ! collocation point counter
      sinit = 0.0D0  ! initialize arc length

      Itp(1) = 1    ! straight segments

c---
c generate the boundary elements (straight segments)
c---

      arcl = 0.0D0

      Do i=1,NE(1)

        Ic = Ic + 1
        x0(Ic) = 0.5D0*(xw(1,i)+xw(1,i+1))  ! collocation points
        y0(Ic) = 0.5D0*(yw(1,i)+yw(1,i+1))
        ddx    = xw(1,i+1)-xw(1,i)
        ddy    = yw(1,i+1)-yw(1,i)
        ddl    = Dsqrt(ddx**2+ddy**2)
        elml(Ic) = ddl

        arcl   = arcl + 0.5D0*ddl
        s0(Ic) = arcl
        arcl   = arcl + 0.5D0*ddl

        tnx0(Ic) = ddx/ddl    ! tangential vector
        tny0(Ic) = ddy/ddl

        vnx0(Ic) =-tny0(Ic)   ! normal vector
        vny0(Ic) = tnx0(Ic)

        U0(Ic) = fc*((x0(Ic)-orgx)**2+(y0(Ic)-orgy)**2)

      End Do

      Ncl = Ic          ! number of collocation points

c-----------
      end if           ! End of geometry options
c-----------

c-----
c Done
c-----

 100  Format (1x,i3,20(1x,f15.10))
 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 107  Format (" Velocity: ",f10.5)
 150  Format (1X,10(1X,f10.5))

      return
      end
