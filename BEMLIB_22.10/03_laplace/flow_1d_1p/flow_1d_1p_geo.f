      subroutine flow_1d_1p_geo
     +
     +   (shrt
     +   ,Ncl
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c Steady unidirectional shear flow over a
c periodic array of cylinders.
c
c SYMBOLS:
c --------
c
c NSG: Number of segments
c
c NE(i): Number of elements on ith segment 
c
c RT(i): Stretch ratio of elements on the ith segment 
c
c elml(i): length of the ith element
c
c Itp(i): Index for shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c Xe, Ye:  end-nodes of elements on a segment
c Xm, Ym:  mid-nodes of elements on a segment
c Xw, Yw:  end-nodes of elements on all segments
c
c X0, Y0:  coordinates of collocation points
c
c t0: angle measured around the center a circular element
c
c shrt: shear rate
c
c Ncl: number of collocation points
c
c Note:
c -----
c
c The normal vector points into the fluid
c
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(129),Ye(129),Te(129),Se(129)
      Dimension Xm(128),Ym(128),Tm(128),Sm(128)

      Dimension    NE(10),   RT(10),Itp(10),Actis(10)
      Dimension xcntr(10),ycntr(10)

      Dimension xw(10,129),yw(10,129),tw(10,129)

      Dimension   x0(1280),  y0(1280),t0(1280)
      Dimension   s0(1280),  u0(1280)
      Dimension tnX0(1280),tnY0(1280)
      Dimension vnX0(1280),vnY0(1280)
      Dimension elml(1280)

c---
c common blocks
c---

      common/VEL00/RL,visc,Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw
      common/VEL03/actis,xcntr,ycntr
      common/VEL04/x0,y0,t0,s0,tnx0,tny0,vnx0,vny0,elml
      common/VEL05/u0

      common/REC01/recx,recy

c----------
c constants
c----------

      pi = 3.1415 92653 58979 32384 D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

      zero = 0.0D0

c------------------
c Circular cylinders
c
c Important:
c
c elements must be distributed 
c in the countercloskwise sense
c-------------------

      If(Iflow.eq.1) then

      open (4,file='circle.dat',status='unknown')

        read (4,*) shrt        ! shear rate
        read (4,*) visc        ! fluid viscosity
        read (4,*) NGL         ! number of Gauss Legendre points
        read (4,*) radius      ! tube radius
        read (4,*) RL          ! period
        read (4,*)
        read (4,*) NE(1),RT(1) ! number of elements, strech ratio

      close (4)

      If(radius.eq.1.0) then
        write (6,*) 
        write (6,*) " flow_1d_1p: the solution of the integral"
        write (6,*) "             equation is not unique"
        write (6,*) 
        write (6,*) "  Please rescale"
        write (6,*) 
        stop
      End If

c---
c initialize
c---

      NSG    = 1      ! one segment only
      Ic     = 0      ! collocation point counter
      sinit  = 0.0D0  ! initialize arc length
      Itp(1) = 2      ! circular elements

c---
c prepare
c---

      actis(1) = radius

      angle1 = 0.0D0
      angle2 = pi2

      xstart = radius    ! first point
      ystart = 0.0       ! first point

      xcnt = xstart-radius*Dcos(angle1)    ! circle center
      ycnt = ystart-radius*Dsin(angle1)

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

c---
c collocation points
c---

      Do i=1,NE(1)

        Ic = Ic + 1

        x0(Ic) = xm(i) 
        y0(Ic) = ym(i)
        t0(Ic) = tm(i)
        s0(Ic) = sm(i)
        u0(Ic) = shrt*y0(Ic)

        elml(Ic) = se(i+1)-se(i)    ! element length
        tnx0(Ic) =-sin(t0(Ic))      ! tangential vector
        tny0(Ic) = cos(t0(Ic))
        vnx0(Ic) = tny0(Ic)         ! normal vector outward
        vny0(Ic) =-tnx0(Ic)

      End Do

      Ncl = Ic          ! number of collocation points

c------------------
c Elliptical cylinders
c
c Important:
c elements must be arranged 
c in the countercloskwise sense
c-------------------

      Else If(Iflow.eq.2) then

      open (4,file='ellipse.dat',status='unknown')

        read (4,*) shrt     ! shear rate
        read (4,*) visc     ! viscosity
        read (4,*) NGL      ! number of Gauss-Legendre points
        read (4,*) ael
        read (4,*) bel
        read (4,*) RL       ! period
        read (4,*) 
        read (4,*) NE(1)

      close (4)

c---
c prepare
c---

      NSG    = 1    ! one segment
      Ic     = 0    ! collocation points counter
      sinit  = 0.0  ! initialize arc length
      Itp(1) = 1    ! straight segments

c---
c generate the boundary elements (straight segments)
c---

      step = pi2/NE(1)

      Do i=1,NE(1)+1
        tt = (i-1.0D0)*step
        xw(1,i) = ael*cos(tt)
        yw(1,i) = bel*sin(tt)
      End Do

c---
c collocation points
c---

      arcl = 0.0D0

      Do i=1,NE(1)

        Ic = Ic + 1

        x0(Ic)   = 0.5D0*(xw(1,i)+xw(1,i+1))  ! collocation points
        y0(Ic)   = 0.5D0*(yw(1,i)+yw(1,i+1))
        u0(Ic  ) = shrt*y0(Ic)
        ddx      = xw(1,i+1)-xw(1,i)
        ddy      = yw(1,i+1)-yw(1,i)
        ddl      = Dsqrt(ddx**2+ddy**2)
        elml(Ic) = ddl

        arcl   = arcl + 0.5D0*ddl
        s0(Ic) = arcl
        arcl   = arcl + 0.5D0*ddl

        tnx0(Ic) = ddx/ddl    ! tangential vector
        tny0(Ic) = ddy/ddl

        vnx0(Ic) = tny0(Ic)   ! normal vector outward
        vny0(Ic) =-tnx0(Ic)

      End Do

      Ncl = Ic          ! number of collocation points

c------------------
c Rectangular cylinders
c
c Important:
c elements must be arranged
c in the countercloskwise sense
c-------------------

      Else If(Iflow.eq.3) then

      open (4,file='rec.dat',status='unknown')

        read (4,*) shrt
        read (4,*) visc
        read (4,*) NGL
        read (4,*) recx
        read (4,*) recy
        read (4,*) RL          ! period
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)
        read (4,*) NE(4),RT(4)

      close (4)

c---
c prepare
c---

      NSG = 4       ! four segments

      Ic    = 0       ! collocation points counter
      sinit = 0.0D0   ! initialize arc length

      recxm = -recx
      recym = -recy

c---
c Default cylinder center is (0,0)
c Shift if desired
c---

      shiftx = 0.0D0
      shifty = 0.0D0
      shifty = recy

c---
c  side # 1 (right side)
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,recx-shiftx, recym-shifty
     +   ,recx-shiftx, recy -shifty
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(1)+1
        xw(1,i) = Xe(i)
        yw(1,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Do i=1,NE(1)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)
        u0(Ic) = shrt*y0(Ic)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

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
     +   ,recx -shiftx,recy-shifty
     +   ,recxm-shiftx,recy-shifty
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(2)+1
       xw(2,i) = Xe(i)
       yw(2,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Do i=1,NE(2)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)
        u0(Ic) = shrt*y0(Ic)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

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
     +   ,recxm-shiftx,recy -shifty
     +   ,recxm-shiftx,recym-shifty
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Do i=1,NE(3)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)
        u0(Ic) = shrt*y0(Ic)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

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
     +   ,recxm-shiftx,recym-shifty
     +   ,recx -shiftx,recym-shifty
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(4)+1
        XW(4,i) = Xe(i)
        YW(4,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Do i=1,NE(4)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)
        u0(Ic) = shrt*y0(Ic)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)


      End Do

      Ncl = Ic          ! number of collocation points

c----------------
c Triangular cylinders
c
c Important:
c
c Vertices and thus elements
c must be arranged in the counterclockwise sense
c-----------------

      Else If(Iflow.eq.4) then

      open (4,file='triangle.dat',status='unknown')

        read (4,*) shrt
        read (4,*) visc
        read (4,*) NGL
        read (4,*) 
        read (4,*) xfirst,yfirst
        read (4,*) xsecond,ysecond
        read (4,*) xthird,ythird
        read (4,*) RL            ! period
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)

      close (4)

c---
c preparations
c---

      NSG   = 3    ! three sides

      Ic    = 0    ! collocation point counter
      sinit = 0.   ! initialize arc length

c---
c  side # 1
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line
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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

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
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

      End Do

      sinit = se(NE(2)+1)

c---
c  side # 3
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line 
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,xthird,ythird
     +   ,xfirst,yfirst
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
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

      End Do

      Ncl = Ic          ! number of collocation points

c------------------
c Rectangular protrusion
c
c Important:
c elements must be arranged
c in the countercloskwise sense
c-------------------

      Else If(Iflow.eq.10) then

      open (4,file='prt_rec.dat',status='unknown')

        read (4,*) shrt
        read (4,*) visc
        read (4,*) NGL
        read (4,*) width
        read (4,*) height
        read (4,*) RL          ! period
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3)
        read (4,*) NE(4),RT(4)
        read (4,*) NE(5),RT(5)

      close (4)

c---
c preparations
c---

      NSG = 5       ! five segments

      Ic    = 0     ! collocation points counter
      sinit = 0.0   ! initialize arc length

      widthh  = 0.50D0*width
      widthhm = -widthh
      RLh     = 0.50D0*RL
      RLhm    = -RLh

c---
c  side # 1 (right side)
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1    ! symmetric distribution

      call elm_line
     +
     +  (NE(1)
     +  ,RT(1)
     +  ,RLh,zero
     +  ,widthh,zero
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,se
     +  ,Xm,Ym,sm
     +  )

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
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

      End Do

      sinit = se(NE(1)+1)

c---
c  side # 2
c---

      Itp(2) = 1    ! straight segment
      Isym   = 1

      call elm_line 
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,widthh,zero
     +   ,widthh,height
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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

      End Do

      sinit = se(NE(2)+1)

c---
c  side # 3 
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,widthh,height
     +   ,widthhm,height
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
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
        ddl = sqrt(ddx**2+ddy**2)

        elml(Ic) = ddl

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

      End Do

      sinit = se(NE(3)+1)

c---
c  side # 4
c---

      Itp(4) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +         (NE(4)
     +         ,RT(4)
     +         ,widthhm,height
     +         ,widthhm,zero
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

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
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

      End Do

      sinit = se(NE(4)+1)

c---
c  side # %
c---

      Itp(5) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +   (NE(5)
     +   ,RT(5)
     +   ,widthhm,zero
     +   ,RLhm,zero
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(5)+1
        XW(5,i) = Xe(i)
        YW(5,i) = Ye(i)
      End Do

      Do i=1,NE(5)

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
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        u0(Ic) = shrt*y0(Ic)

      End Do

      Ncl = Ic          ! number of collocation points

c-----------
      End If           ! End of geometry module
c-----------


      Return
      End
