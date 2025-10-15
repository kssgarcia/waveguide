      subroutine flow_2d_geo
     +
     +   (Ncl
     +   ,Iwall
     +   ,wall
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

c---------------------------------------------
c Element distribution,
c computation of collocation points
c and disturbance normal velocity
c
c SYMBOLS
c -------
c
c NSG: Number of segments
c
c NE(i): Number of elements on ith segment
c
c RT(i): Stretch ratio of elements on ith segment
c
c Itp(i): Index for shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c (Xe, Ye):  end-nodes of elements on a segment
c (Xm, Ym):  mid-nodes of elements on a segment
c (Xw, Yw):  end-nodes of elements on all segments
c
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(100),Ye(100),Te(100),Se(100)
      Dimension Xm(100),Ym(100),Tm(100),Sm(100)

      Dimension NE(10),RT(10),Itp(10),Actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)
      Dimension phi(10,200)

      Dimension X0(500),Y0(500),T0(500),S0(0:500),dphidn0(500)
      Dimension tnX0(500),tnY0(500)
      Dimension vnX0(500),vnY0(500)

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw,phi,dphidn0
      common/xxx03/actis,xcntr,ycntr
      common/xxx04/Vx,Vy
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx08/xwmin,ywmin,xwmax,ywmax

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

c--------------------------------
c Circular cavity on a plane wall
c or in a channel
c--------------------------------

      If(Iflow.eq.1) then

      open (4,file='cvt_crc.dat')

      read (4,*) NGL
      read (4,*) radius
      read (4,*) angle
      read (4,*) 
      read (4,*) NE(1),RT(1),trunc1
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*) NE(4),RT(4),trunc4
      read (4,*) 
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) NE(5),RT(5),trunc5
      read (4,*) NE(6),RT(6),trunc6
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c---
c preparations
c---

      angl     = angle*pi
      opening  = radius*sin(angl)
      openingm = - opening

      Vx = 1.0D0   ! uniform flow
      Vy = 0.0D0

      Ic    = 0       ! collocation point counter
      sinit = 0.0D0   ! initialize arc length

c---
c truncated wall
c---

      Itp(1) = 1 
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +    (NE(1)
     +    ,RT(1)
     +    ,trunc1,zero
     +    ,openingm,zero
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

      Do i=1,NE(1)+1
        XW(1,i) = Xe(i)
        YW(1,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Do i=1,NE(1)
        Ic = Ic+1
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(1)+1)

c---
c left part of the cavity
c---

      Itp(2) = 2
      Isym   = 0

      actis(2) = radius

      angle1 = 3.0*pih-angl
      angle2 = 3.0*pih

      xstart = XW(1,NE(1)+1)
      ystart = YW(1,NE(1)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(2) = xcnt
      ycntr(2) = ycnt

      call elm_arc        ! cavity
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,Xcnt,Ycnt
     +   ,radius
     +   ,angle1,angle2
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,Te,se
     +   ,Xm,Ym,Tm,sm
     +   )

      Do i=1,NE(2)+1
        XW(2,i) = Xe(i)
        YW(2,i) = Ye(i)
        TW(2,i) = Te(i)
      End Do

      Do i=1,NE(2)
        Ic = Ic + 1
        x0  (Ic) = xm(i)   ! collocation points
        y0  (Ic) = ym(i)
        t0  (Ic) = tm(i)
        s0  (Ic) = sm(i)
        tnx0(Ic) =-sin(t0(Ic))
        tny0(Ic) = cos(t0(Ic))
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(2)+1)

c---
c right part of the cavity
c---

      Itp(3) = 2
      Isym   = 0

      actis(3) = radius

      angle1 = 3.0*pih
      angle2 = 3.0*pih+angl

      xstart = XW(2,NE(2)+1)
      ystart = YW(2,NE(2)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(3) = xcnt
      ycntr(3) = ycnt

      call elm_arc        ! cavity
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,Xcnt,Ycnt
     +   ,radius
     +   ,angle1,angle2
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,Te,se
     +   ,Xm,Ym,Tm,sm
     +   )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
        TW(3,i) = Te(i)
      End Do

      Do i=1,NE(3)
        Ic = Ic + 1
        x0  (Ic) = xm(i)   ! collocation points
        y0  (Ic) = ym(i)
        t0  (Ic) = tm(i)
        s0  (Ic) = sm(i)
        tnx0(Ic) =-sin(t0(Ic))
        tny0(Ic) = cos(t0(Ic))
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(3)+1)

c---
c  right part of the wall
c---

      Itp(4) = 1    ! straight segment
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +   (NE(4)
     +   ,RT(4)
     +   ,opening,zero
     +   ,trunc4,zero
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
        den = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      NSG = 4               ! number of segments

c---
c  upper wall
c---

      If(Iwall.eq.1) then

        Itp(5) = 1
        Isym   = 0
        sinit  = 0.0D0

        call elm_line         ! truncated upper wall
     +
     +     (NE(5)
     +     ,RT(5)
     +     ,trunc5,wall
     +     ,zero,wall
     +     ,sinit
     +     ,Isym
     +     ,Xe,Ye,se
     +     ,Xm,Ym,sm
     +     )

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
          den = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
          vnx0(Ic) =-tny0(Ic)
          vny0(Ic) = tnx0(Ic)
          dphidn0(Ic) = - Vx*vnx0(Ic)
        End Do

        sinit = se(NE(5)+1)

        Itp(6) = 1
        Isym  = 0

        call elm_line         ! truncated upper wall
     +
     +         (NE(6)
     +         ,RT(6)
     +         ,zero,wall
     +         ,trunc6,wall
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

        Do i=1,NE(6)+1
          XW(6,i) = Xe(i)
          YW(6,i) = Ye(i)
        End Do

c---
c collocation points
c---

        Do i=1,NE(6)
          Ic = Ic + 1
          x0(Ic) = xm(i)
          y0(Ic) = ym(i)
          s0(Ic) = sm(i)
          ddx = xe(i+1)-xe(i)
          ddy = ye(i+1)-ye(i)
          den = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
          vnx0(Ic) =-tny0(Ic)
          vny0(Ic) = tnx0(Ic)
          dphidn0(Ic) = - Vx*vnx0(Ic)
        End Do

        NSG = 6

      End If

      Ncl = Ic          ! number of collocation points

c-----------------------------------
c Rectangular cavity on a plane wall
c or in a channel
c-----------------------------------

      Else If(Iflow.eq.11) then

      open (4,file='cvt_rec.dat')

      read (4,*) NGL
      read (4,*) depth
      read (4,*) opening
      read (4,*) 
      read (4,*) NE(1),RT(1),trunc1
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*) NE(4),RT(4)
      read (4,*) NE(5),RT(5),trunc5
      read (4,*) 
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) NE(6),RT(6),trunc6
      read (4,*) NE(7),RT(7),trunc7
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c---
c preparations
c---

      Ic    = 0       ! collocation point counter
      sinit = 0.0D0   ! initialize arc length

      depthm   = -depth
      openingm = -opening

      Vx = 1.0D0
      Vy = 0.0D0

c---
c  left part of the wall
c---

      Itp(1) = 1    ! straight segment
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,trunc1,zero
     +   ,openingm,zero
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
        den = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(1)+1)

c---
c  left side of the cavity
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +         (NE(2)
     +         ,RT(2)
     +         ,openingm,zero
     +         ,openingm,depthm
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
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(2)+1)

c---
c  Bottom part of the cavity
c---

      Itp(3) = 1
      Isym   = 1

      call elm_line            ! bottom of cavity
     +
     +         (NE(3)
     +         ,RT(3)
     +         ,openingm,depthm
     +         ,opening ,depthm
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
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(3)+1)

c---
c  right side of the cavity
c---

      Itp(4) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +         (NE(4)
     +         ,RT(4)
     +         ,opening ,depthm
     +         ,opening ,zero
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
        den = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(4)+1)

c---
c  right side of the wall
c---

      Itp(5) = 1
      Isym   = 0

      call elm_line           ! truncated wall
     +
     +         (NE(5)
     +         ,RT(5)
     +         ,opening,zero
     +         ,trunc5,zero
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

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
        den = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      NSG = 5               ! number of segments

c---
c  upper wall
c---

      If(Iwall.eq.1) then

        Itp(6) = 1
        Isym   = 0
        sinit  = 0.0

        call elm_line         ! truncated upper wall
     +
     +         (NE(6)
     +         ,RT(6)
     +         ,trunc6,wall
     +         ,zero,wall
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

        Do i=1,NE(6)+1
          XW(6,i) = Xe(i)
          YW(6,i) = Ye(i)
        End Do

        Do i=1,NE(6)
          Ic = Ic + 1
          x0(Ic) = xm(i)
          y0(Ic) = ym(i)
          s0(Ic) = sm(i)
          ddx = xe(i+1)-xe(i)
          ddy = ye(i+1)-ye(i)
          den = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
          vnx0(Ic) =-tny0(Ic)
          vny0(Ic) = tnx0(Ic)
          dphidn0(Ic) = - Vx*vnx0(Ic)
        End Do

        sinit = se(NE(6)+1)

        Itp(7) = 1
        Isym  = 0

        call elm_line         ! truncated upper wall
     +
     +         (NE(7)
     +         ,RT(7)
     +         ,zero,wall
     +         ,trunc7,wall
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

        Do i=1,NE(7)+1
          XW(7,i) = Xe(i)
          YW(7,i) = Ye(i)
        End Do

        Do i=1,NE(7)
          Ic = Ic + 1
          x0(Ic) = xm(i)
          y0(Ic) = ym(i)
          s0(Ic) = sm(i)
          ddx = xe(i+1)-xe(i)
          ddy = ye(i+1)-ye(i)
          den = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
          vnx0(Ic) =-tny0(Ic)
          vny0(Ic) = tnx0(Ic)
          dphidn0(Ic) = - Vx*vnx0(Ic)
        End Do

        NSG = 7

      End If

      Ncl = Ic          ! number of collocation points

c------------------------------------
c Circular protrusion on a plane wall
c or in a channel
c------------------------------------

      Else If(Iflow.eq.41) then

      open (4,file='prt_crc.dat')

      read (4,*) NGL
      read (4,*) radius
      read (4,*) angle
      read (4,*) 
      read (4,*) NE(1),RT(1),trunc1
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*) NE(4),RT(4),trunc4
      read (4,*) 
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) NE(5),RT(5),trunc5
      read (4,*) NE(6),RT(6),trunc6
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c---
c preparations
c---

      angl     = angle*pi
      opening  = radius*sin(angl)
      openingm = - opening

      Ic    = 0    ! collocation point counter
      sinit = 0.   ! initialize arc length

      Vx = 1.0
      Vy = 0.0
c---

      Itp(1) = 1    ! straight segment
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +         (NE(1)
     +         ,RT(1)
     +         ,trunc1,zero
     +         ,openingm,zero
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
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(1)+1)

c---
c left part of the protrusion
c---

      Itp(2) = 2
      Isym   = 0

      actis(2) = radius

      angle1 = pih+angl
      angle2 = pih

      xstart = XW(1,NE(1)+1)
      ystart = YW(1,NE(1)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(2) = xcnt
      ycntr(2) = ycnt

      call elm_arc        ! cavity
     +
     +         (NE(2)
     +         ,RT(2)
     +         ,Xcnt,Ycnt
     +         ,radius
     +         ,angle1,angle2
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,Te,se
     +         ,Xm,Ym,Tm,sm
     +         )

      Do i=1,NE(2)+1
        XW(2,i) = Xe(i)
        YW(2,i) = Ye(i)
        TW(2,i) = Te(i)
      End Do

      Do i=1,NE(2)
        Ic = Ic + 1
        x0  (Ic) = xm(i)   ! collocation points
        y0  (Ic) = ym(i)
        t0  (Ic) = tm(i)
        s0  (Ic) = sm(i)
        tnx0(Ic) = sin(t0(Ic))
        tny0(Ic) =-cos(t0(Ic))
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(2)+1)

c------
c right part of the protrusion
c------

      Itp(3) = 2
      Isym   = 0

      actis(3) = radius

      angle1 = pih
      angle2 = pih-angl

      xstart = XW(2,NE(2)+1)
      ystart = YW(2,NE(2)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(3) = xcnt
      ycntr(3) = ycnt

      call elm_arc        ! cavity
     +
     +         (NE(3)
     +         ,RT(3)
     +         ,Xcnt,Ycnt
     +         ,radius
     +         ,angle1,angle2
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,Te,se
     +         ,Xm,Ym,Tm,sm
     +         )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
        TW(3,i) = Te(i)
      End Do

      Do i=1,NE(3)
        Ic = Ic + 1
        x0  (Ic) = xm(i)   ! collocation points
        y0  (Ic) = ym(i)
        t0  (Ic) = tm(i)
        s0  (Ic) = sm(i)
        tnx0(Ic) = sin(t0(Ic))
        tny0(Ic) =-cos(t0(Ic))
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      sinit = se(NE(3)+1)

c---
c   right part of the wall
c---

      Itp(4) = 1    ! straight segment
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +         (NE(4)
     +         ,RT(4)
     +         ,opening,zero
     +         ,trunc4,zero
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
        den = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        vnx0(Ic) =-tny0(Ic)
        vny0(Ic) = tnx0(Ic)
        dphidn0(Ic) = - Vx*vnx0(Ic)
      End Do

      NSG = 4               ! number of segments

c--- 
c  upper wall
c--- 

      If(Iwall.eq.1) then

        sinit = 0.0

        Itp(5) = 1
        Isym   = 0
        sinit  = 0.0

        call elm_line         ! truncated upper wall
     +
     +         (NE(5)
     +         ,RT(5)
     +         ,trunc5,wall
     +         ,zero,wall
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

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
          den = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
          vnx0(Ic) =-tny0(Ic)
          vny0(Ic) = tnx0(Ic)
          dphidn0(Ic) = - Vx*vnx0(Ic)
        End Do

        sinit = se(NE(5)+1)

        Itp(6) = 1
        Isym  = 0

        call elm_line         ! truncated upper wall
     +
     +         (NE(6)
     +         ,RT(6)
     +         ,zero,wall
     +         ,trunc6,wall
     +         ,sinit
     +         ,Isym
     +         ,Xe,Ye,se
     +         ,Xm,Ym,sm
     +         )

        Do i=1,NE(6)+1
          XW(6,i) = Xe(i)
          YW(6,i) = Ye(i)
        End Do

        Do i=1,NE(6)
          Ic = Ic + 1
          x0(Ic) = xm(i)
          y0(Ic) = ym(i)
          s0(Ic) = sm(i)
          ddx = xe(i+1)-xe(i)
          ddy = ye(i+1)-ye(i)
          den = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
          vnx0(Ic) =-tny0(Ic)
          vny0(Ic) = tnx0(Ic)
          dphidn0(Ic) = - Vx*vnx0(Ic)
        End Do

        NSG = 6

      End If

      Ncl = Ic          ! number of collocation points

c-----------
      End If       ! End of geometry module
c-----------

c-----
c Done
c-----

      Return
      End
