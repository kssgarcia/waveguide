      subroutine flow_2d_geo
     +
     +  (Ncl
     +  ,Iwall
     +  ,Itry
     +  )

c-----------------------------------------
c FDLIB - BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------
c Element distribution
c and computation of collocation points
c and disturbance velocity
c at collocation points
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
c Itp(i): Index for the shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c (Xe, Ye):  end-nodes of elements
c (Xm, Ym):  mid-nodes of elements
c
c (X0, Y0):  coordinates of collocation points
c
c t0:	     angle subtended from a circular segment center
c
c (ux0, uy0): perturbation velocity components at collocation points
c
c fx, fy, Cartesian components of the traction
c
c NGL: Number of Gaussian points for integration over each element
c
c shrt:	shear rate of incident flow at the wall
c
c delta: pressure gradient of incident flow
c
c ftn: tangential component of the traction
c
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(128),Ye(128),Te(128),Se(128)
      Dimension Xm(128),Ym(128),Tm(128),Sm(128)

      Dimension NE(10),RT(10),Itp(10),actis(10)
      Dimension     xcntr(10),        ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)
      Dimension fx(10,200),fy(10,200)

      Dimension   X0(900),  Y0(900),T0(900),S0(900)
      Dimension tnX0(900),tnY0(900)
      Dimension  ux0(900), uy0(900)
      Dimension elml(1280)

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw,fx,fy,ux0,uy0
      common/xxx03/actis,xcntr,ycntr
      common/xxx04/visc,shrt,delta
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx07/elml
      common/xxx08/xwmin,ywmin,xwmax,ywmax
      common/xxx09/wallnew,rotnew,ycntnew
      common/xxx10/xcnt,ycnt

      common/flow_91/RL,Uslip
      common/flow_81/Utrans,Orot

      common/wwww/wall,rotation

c----------
c constants
c----------

      pi = 3.14159 265358D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pi6 = 6.0D0*pi
      pi8 = 8.0D0*pi

      zero = 0.0D0

c--------------------------------
c circular cavity on a plane wall
c--------------------------------

      If(Iflow.eq.1) then

      open (4,file='cvt_crc.dat')

      read (4,*) visc
      read (4,*) shrt
      read (4,*) delta
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

      if(Iwall.eq.1) delta = -2.0D0*shrt*visc/wall

c     if(Itry.gt.1)  then
c       wall = wallnew
c     end if

      gwnia    = angle*pi
      opening  = radius*sin(gwnia)
      openingm = - opening

      Ic    = 0      !   collocation point counter
      sinit = 0.0D0  !   initialize arc length

      cf = 0.5D0*delta/visc

c---
c left part of the truncated wall
c---

      Itp(1) = 1     ! straight segments
      Isym   = 0     ! element distribution non-symmetric

      call elm_line
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
c       write (6,*) i,Xe(i),Ye(i)
      End Do

      Do i=1,NE(1)     ! collocation points
        Ic = Ic + 1
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        ddl      = Dsqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(1)+1)

c--- 
c left part of the cavity
c--- 

      Itp(2) = 2     ! circular arcs
      Isym   = 0

      angle1 = 3.0D0*pih-gwnia
      angle2 = 3.0D0*pih

      xstart = XW(1,NE(1)+1)
      ystart = YW(1,NE(1)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      actis(2) = radius
      xcntr(2) = xcnt
      ycntr(2) = ycnt

      call elm_arc
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,xcnt,ycnt
     +   ,radius
     +   ,angle1,angle2
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,Te,se
     +   ,Xm,Ym,Tm,sm
     +   )

      Do i=1,NE(2)+1
        xw(2,i) = Xe(i)
        yw(2,i) = Ye(i)
        tw(2,i) = Te(i)
c       write (6,*) i,Xe(i),Ye(i)
      End Do

      Do i=1,NE(2)
        Ic = Ic + 1
          x0(Ic) = xm(i)   ! collocation points
          y0(Ic) = ym(i)
          t0(Ic) = tm(i)
          s0(Ic) = sm(i)
        tnx0(Ic) =-sin(t0(Ic))
        tny0(Ic) = cos(t0(Ic))
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(2)+1)

c------
c  right part of the cavity
c------

      Itp(3) = 2
      Isym   = 0

      angle1 = 3.0*pih
      angle2 = 3.0*pih+gwnia

      xstart = XW(2,NE(2)+1)
      ystart = YW(2,NE(2)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      actis(3) = radius
      xcntr(3) = xcnt
      ycntr(3) = ycnt

      call elm_arc
     +
     +  (NE(3)
     +  ,RT(3)
     +  ,xcnt,ycnt
     +  ,radius
     +  ,angle1,angle2
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,Te,se
     +  ,Xm,Ym,Tm,sm
     +  )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
        TW(3,i) = Te(i)
      End Do

      Do i=1,NE(3)           ! collocation points
        Ic = Ic + 1
        x0  (Ic) = xm(i) 
        y0  (Ic) = ym(i)
        t0  (Ic) = tm(i)
        s0  (Ic) = sm(i)
        tnx0(Ic) =-Dsin(t0(Ic))
        tny0(Ic) = Dcos(t0(Ic))
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      sinit = se(NE(3)+1)

c---
c right part of the truncated wall
c---

      Itp(4) = 1    ! straight segment
      Isym   = 0

      call elm_line
     +
     +  (NE(4)
     +  ,RT(4)
     +  ,opening,zero
     +  ,trunc4,zero
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,se
     +  ,Xm,Ym,sm
     +  )

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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      NSG = 4               ! number of segments

c---
c upper wall
c---

      If(Iwall.eq.1) then

        Itp(5) = 1
        Isym   = 0

        sinit  = 0.0D0   ! re-intialize arc length

        call elm_line
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
          ddl = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0D0
        End Do

        sinit = se(NE(5)+1)

c---
c left part of truncated upper wall
c---

        Itp(6) = 1
        Isym  = 0

        call elm_line
     +
     +    (NE(6)
     +    ,RT(6)
     +    ,zero,wall
     +    ,trunc6,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          ddl = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0D0
        End Do

        NSG = 6

      End If

      Ncl = Ic          ! number of collocation points

c------------------
c rectangular cavity
c on a plane wall
c------------------

      else if(Iflow.eq.11) then

      open (4,file='cvt_rec.dat')

      read (4,*) visc
      read (4,*) shrt
      read (4,*) delta
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
      If(Iwall.eq.1) delta = -2.0*shrt*visc/wall

      If(Itry.gt.1)  then
         wall = wallnew
      End If

      Ic    = 0     ! collocation point counter
      sinit = 0.0   ! initialize arc length

      depthm   = -depth
      openingm = -opening

      cf = 0.5*delta/visc

c---
c  left part of the truncated wall
c---

      Itp(1) = 1    ! straight segment
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +  (NE(1)
     +  ,RT(1)
     +  ,trunc1,zero
     +  ,openingm,zero
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
        den = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      sinit = se(NE(1)+1)

c---
c  left side of the cavity
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +  (NE(2)
     +  ,RT(2)
     +  ,openingm,zero
     +  ,openingm,depthm
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
        ddl = Dsqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      sinit = se(NE(2)+1)

c---
c  Bottom of the cavity
c---

      Itp(3) = 1
      Isym   = 1

      call elm_line            ! bottom of cavity
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,openingm,depthm
     +   ,opening ,depthm
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
        ddl = Dsqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0
      End Do

      sinit = se(NE(3)+1)

c---
c  right side of the cavity
c---

      Itp(4) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +    (NE(4)
     +    ,RT(4)
     +    ,opening ,depthm
     +    ,opening ,zero
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf * y0(Ic)**2
        uy0 (Ic) = 0.0
      End Do

      sinit = se(NE(4)+1)

c---
c  right side of the truncated wall
c---

      Itp(5) = 1
      Isym   = 0

      call elm_line           ! truncated wall
     +
     +    (NE(5)
     +    ,RT(5)
     +    ,opening,zero
     +    ,trunc5,zero
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      NSG = 5               ! number of segments

c------------
c  upper wall
c------------

      If(Iwall.eq.1) then

c---
c  right part of the upper wall
c---

        Itp(6) = 1
        Isym   = 0

        sinit  = 0.0          ! re-initialize arc length

        call elm_line         ! truncated upper wall
     +
     +    (NE(6)
     +    ,RT(6)
     +    ,trunc6,wall
     +    ,zero,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0
        End Do

        sinit = se(NE(6)+1)

c---
c  left part of the upper wall
c---

        Itp(7) = 1
        Isym   = 0

        call elm_line         ! truncated upper wall
     +
     +    (NE(7)
     +    ,RT(7)
     +    ,zero,wall
     +    ,trunc7,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0D0
        End Do

        NSG = 7

      End If

      Ncl = Ic          ! number of collocation points

c--------------------
c Circular protrusion
c on a plane wall
c--------------------

      elseif(Iflow.eq.41) then

      open (4,file='prt_crc.dat')

      read (4,*) visc
      read (4,*) shrt
      read (4,*) delta
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
      if(Iwall.eq.1) delta = -2.0*shrt*visc/wall

      if(Itry.gt.1)  then
         wall = wallnew
      endif

      gwnia    = angle*pi
      opening  = radius*sin(gwnia)
      openingm = -opening

      cf = 0.5*delta/visc

      Ic    = 0     ! collocation point counter
      sinit = 0.0   ! initialize arc length

c---
c left part of truncated wall
c---

      Itp(1) = 1    ! straight segment
      Isym   = 0

      call elm_line
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

      Do i=1,NE(1)
        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)
        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = -shrt*y0(Ic)-cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      sinit = se(NE(1)+1)

c---
c  left part of the protrusion
c---

      Itp(2) = 2
      Isym   = 0

      actis(2) = radius

      angle1 = pih+gwnia
      angle2 = pih

      xstart = XW(1,NE(1)+1)
      ystart = YW(1,NE(1)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(2) = xcnt
      ycntr(2) = ycnt

      call elm_arc        ! protrusion
     +
     +  (NE(2)
     +  ,RT(2)
     +  ,xcnt,ycnt
     +  ,radius
     +  ,angle1,angle2
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,Te,se
     +  ,Xm,Ym,Tm,sm
     +  )

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
        tnx0(Ic) =-Dsin(t0(Ic))
        tny0(Ic) = Dcos(t0(Ic))
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = -shrt*y0(Ic)-cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      sinit = se(NE(2)+1)

c---
c  right part of the protrusion
c---

      Itp(3) = 2
      Isym   = 0

      actis(3) = radius

      angle1 = pih
      angle2 = pih-gwnia

      xstart = XW(2,NE(2)+1)
      ystart = YW(2,NE(2)+1)

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(3) = xcnt
      ycntr(3) = ycnt

      call elm_arc        ! protrusion
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,xcnt,ycnt
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
        tnx0(Ic) =-Dsin(t0(Ic))
        tny0(Ic) = Dcos(t0(Ic))
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = -shrt*y0(Ic)-cf*y0(Ic)**2
        uy0 (Ic) = 0.0
      End Do

      sinit = se(NE(3)+1)

c---
c   right part of the wall
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
        ddl = Dsqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = -shrt*y0(Ic)-cf*y0(Ic)**2
        uy0 (Ic) = 0.0
      End Do

      sinit = se(NE(4)+1)

      NSG = 4               ! number of segments

c---
c  upper wall
c---

      if(Iwall.eq.1) then

c---
c right part of the upper wall
c---

        Itp(5) = 1
        Isym   = 0

        sinit  = 0.0     ! re-initialize arce legnth

        call elm_line
     +
     +    (NE(5)
     +    ,RT(5)
     +    ,trunc5,wall
     +    ,zero,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          ddl = Dsqrt(ddx*ddx+ddy*ddy)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
           ux0(Ic) = -shrt*y0(Ic)-cf*y0(Ic)**2
           uy0(Ic) = 0.0
        End Do

        sinit = se(NE(5)+1)

c---
c left part of the upper wall
c---

        Itp(6) = 1
        Isym   = 0

        call elm_line 
     +
     +    (NE(6)
     +    ,RT(6)
     +    ,zero,wall
     +    ,trunc6,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          ddl = Dsqrt(ddx*ddx+ddy*ddy)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0
        End Do

        NSG = 6

      End If

      Ncl = Ic          ! number of collocation points

c------------------
c Rectangular protrusion
c on a plane wall
c------------------

      elseif(Iflow.eq.51) then

      open (4,file='prt_rec.dat')

      read (4,*) visc
      read (4,*) shrt
      read (4,*) delta
      read (4,*) NGL
      read (4,*) height
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
      If(Iwall.eq.1) delta = -2.0*shrt*visc/wall

      If(Itry.gt.1)  then
         wall = wallnew
      End If

      Ic    = 0     ! collocation point counter
      sinit = 0.0D0 ! initialize arc length

      openingm = -opening

      cf = 0.5*delta/visc

c---
c  left part of the truncated wall
c---

      Itp(1) = 1    ! straight segment
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +  (NE(1)
     +  ,RT(1)
     +  ,trunc1,zero
     +  ,openingm,zero
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      sinit = se(NE(1)+1)

c---
c  left side of the protrusion
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +    (NE(2)
     +    ,RT(2)
     +    ,openingm,zero
     +    ,openingm,height
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0
      End Do

      sinit = se(NE(2)+1)

c---
c top of the protrusion
c---

      Itp(3) = 1
      Isym   = 1

      call elm_line            ! bottom of cavity
     +
     +    (NE(3)
     +    ,RT(3)
     +    ,openingm,height
     +    ,opening ,height
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      sinit = se(NE(3)+1)

c---
c  right side of the protrusion
c---

      Itp(4) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +    (NE(4)
     +    ,RT(4)
     +    ,opening,height
     +    ,opening,zero
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf * y0(Ic)**2
        uy0 (Ic) = 0.0
      End Do

      sinit = se(NE(4)+1)

c---
c  right side of the truncated wall
c---

      Itp(5) = 1
      Isym   = 0

      call elm_line           ! truncated wall
     +
     +    (NE(5)
     +    ,RT(5)
     +    ,opening,zero
     +    ,trunc5,zero
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      NSG = 5               ! number of segments

c------------
c  upper wall
c------------

      If(Iwall.eq.1) then

c---
c  right part of the upper wall
c---

        Itp(6) = 1
        Isym   = 0

        sinit = 0.0D0          ! re-initialize arc length

        call elm_line         ! truncated upper wall
     +
     +     (NE(6)
     +     ,RT(6)
     +     ,trunc6,wall
     +     ,zero,wall
     +     ,sinit
     +     ,Isym
     +     ,Xe,Ye,se
     +     ,Xm,Ym,sm
     +     )

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
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0
        End Do

        sinit = se(NE(6)+1)

c---
c  left part of the upper wall
c---

        Itp(7) = 1
        Isym   = 0

        call elm_line         ! truncated upper wall
     +
     +     (NE(7)
     +     ,RT(7)
     +     ,zero,wall
     +     ,trunc7,wall
     +     ,sinit
     +     ,Isym
     +     ,Xe,Ye,se
     +     ,Xm,Ym,sm
     +     )

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
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0D0
        End Do

        NSG = 7

      End If

      Ncl = Ic          ! number of collocation points

c-----------------------------
c square cylinder above a wall
c-----------------------------

      Else If(Iflow.eq.71) then

      open (4,file='square.dat')

      read (4,*) visc
      read (4,*) shrt
      read (4,*) delta
      read (4,*) NGL
      read (4,*) width
      read (4,*) height
      read (4,*) xcnt,ycnt
      read (4,*) rotation
      read (4,*) 
      read (4,*) NE(1),RT(1)
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*) NE(4),RT(4)
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

      If(Iwall.eq.1) delta = -2.0D0*shrt*visc/wall

      If(Itry.eq.1)  then
       rotation = pi*rotation
      End If

      If(Itry.gt.1)  then
c       write (6,*) wallnew,ycntnew
        wall      = wallnew
        rotation  = rotnew
        ycnt      = ycntnew
      End If

      Ic    = 0       ! collocation point counter
      sinit = 0.0D0   ! initialize arc length

      cf = 0.5D0*delta/visc

      widh = 0.50D0*width
      heih = 0.50D0*height

      xx1 = xcnt-widh
      yy1 = ycnt+heih
      xx2 = xcnt+widh
      yy2 = ycnt+heih
      xx3 = xcnt+widh
      yy3 = ycnt-heih
      xx4 = xcnt-widh
      yy4 = ycnt-heih

      cs = dcos(rotation)
      sn = dsin(rotation)

      vx1 = (xx1-xcnt)*cs+(yy1-ycnt)*sn + xcnt
      vy1 =-(xx1-xcnt)*sn+(yy1-ycnt)*cs + ycnt
      vx2 = (xx2-xcnt)*cs+(yy2-ycnt)*sn + xcnt
      vy2 =-(xx2-xcnt)*sn+(yy2-ycnt)*cs + ycnt
      vx3 = (xx3-xcnt)*cs+(yy3-ycnt)*sn + xcnt
      vy3 =-(xx3-xcnt)*sn+(yy3-ycnt)*cs + ycnt
      vx4 = (xx4-xcnt)*cs+(yy4-ycnt)*sn + xcnt
      vy4 =-(xx4-xcnt)*sn+(yy4-ycnt)*cs + ycnt

c---
c  first side
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +    (NE(1)
     +    ,RT(1)
     +    ,vx1,vy1
     +    ,vx2,vy2
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(1)+1)

c---
c second side
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,vx2,vy2
     +   ,vx3,vy3
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(2)+1)

c---
c  third side
c---

      Itp(3) = 1
      Isym   = 1

      call elm_line            ! bottom of cavity
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,vx3,vy3
     +   ,vx4,vy4
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(3)+1)

c---
c  fourth side
c---

      Itp(4) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +   (NE(4)
     +   ,RT(4)
     +   ,vx4,vy4
     +   ,vx1,vy1
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic) - cf * y0(Ic)**2
        uy0 (Ic) = 0.0D0
      End Do

      NSG = 4               ! number of segments

c------------
c  upper wall
c------------

      If(Iwall.eq.1) then

c---
c  right part of the upper wall
c---

        Itp(5) = 1
        Isym   = 0

        sinit  = 0.0          ! re-initialize arc length

        call elm_line         ! truncated upper wall
     +
     +    (NE(5)
     +    ,RT(5)
     +    ,trunc5,wall
     +    ,zero,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

        Do i=1,NE(6)+1
          XW(5,i) = Xe(i)
          YW(5,i) = Ye(i)
        End Do

        Do i=1,NE(6)
          Ic = Ic + 1
          x0(Ic) = xm(i)
          y0(Ic) = ym(i)
          s0(Ic) = sm(i)
          ddx = xe(i+1)-xe(i)
          ddy = ye(i+1)-ye(i)
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0
        End Do

        sinit = se(NE(5)+1)

c---
c  left part of the upper wall
c---

        Itp(6) = 1
        Isym   = 0

        call elm_line         ! truncated upper wall
     +
     +    (NE(6)
     +    ,RT(6)
     +    ,zero,wall
     +    ,trunc6,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0D0
        End Do

        NSG = 6

      End If

      Ncl = Ic          ! number of collocation points

c-------------------------------------
c Circular cylinder above a plane wall
c-------------------------------------

      elseif(Iflow.eq.81) then

      open (4,file='circle.dat')

      read (4,*) visc
      read (4,*) shrt
      read (4,*) delta
      read (4,*) Utrans
      read (4,*) Orot
      read (4,*) NGL
      read (4,*) radius
      read (4,*) xcnt,ycnt
      read (4,*) 
      read (4,*) NE(1),RT(1)
      read (4,*)
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) NE(2),RT(2),trunc2
      read (4,*) NE(3),RT(3),trunc3
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c---
c preparations
c---

      If(Iwall.eq.1) delta = -2.0D0*shrt*visc/wall

      If(Itry.gt.1)  then
        wall = wallnew
        ycnt = ycntnew
      End If

      Ic    = 0      !   collocation point counter
      sinit = 0.0D0  !   initialize arc length

      cf = 0.5D0*delta/visc

c---
c circle
c---

      Itp(1) = 2     ! circular elements
      Isym   = 1     ! symmetric element distribution

      actis(1) = radius
      xcntr(1) = xcnt
      ycntr(1) = ycnt

      angle1 = 1.5D0*pi
      angle2 =-0.5D0*pi

      call elm_arc
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
        xw(1,i) = Xe(i)
        yw(1,i) = Ye(i)
        tw(1,i) = Te(i)
c       write (6,100) i,Xe(i),Ye(i),Te(i)
      End Do

      Do i=1,NE(1)
        Ic = Ic + 1
          x0(Ic) = xm(i)   ! collocation points
          y0(Ic) = ym(i)
          t0(Ic) = tm(i)
          s0(Ic) = sm(i)
        tnx0(Ic) =-Dsin(t0(Ic))
        tny0(Ic) = Dcos(t0(Ic))
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
     +           + Utrans 
     +           - Orot*(y0(Ic)-ycnt)
         uy0(Ic) = 0.0D0
     +           + Orot*(x0(Ic)-xcnt)
      End Do

      NSG = 1

c------------
c  upper wall
c------------

      If(Iwall.eq.1) then

c---
c  right part of the upper wall
c---

        Itp(2) = 1
        Isym   = 0

        sinit  = 0.0          ! re-initialize arc length

        call elm_line         ! truncated upper wall
     +
     +    (NE(2)
     +    ,RT(2)
     +    ,trunc2,wall
     +    ,zero,wall
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

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
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0
        End Do

        sinit = se(NE(5)+1)

c---
c  left part of the upper wall
c---

        Itp(3) = 1
        Isym   = 0

        call elm_line         ! truncated upper wall
     +
     +    (NE(3)
     +    ,RT(3)
     +    ,zero,wall
     +    ,trunc3,wall
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
          ddl = Dsqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/ddl
          tny0(Ic) = ddy/ddl
          elml(Ic) = se(i+1)-se(i)    ! element length
          ux0 (Ic) = - shrt*y0(Ic) - cf*y0(Ic)**2
          uy0 (Ic) = 0.0D0
        End Do

        NSG = 3

      End If

      Ncl = Ic          ! number of collocation points

c---------------------------------------
c A periodic array of circular cylinders
c---------------------------------------

      elseif(Iflow.eq.91) then

      open (4,file='circle_1p.dat')

      read (4,*) visc        ! fluid viscosity
      read (4,*) shrt        ! shear rate
      read (4,*) NGL         ! number of Gauss-Legendre base points
      read (4,*) radius      ! radius of a cylinder
      read (4,*) xcyl,ycyl   ! coordinates of the cylinder center
      read (4,*) RL          ! cylinder separation
      read (4,*) 
      read (4,*) NE(1)       ! number of elements around a cylinder
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

      Iwall = 0     ! no upper wall
      NSG   = 1     !  one circular segment

        Itp(1) = 2     ! circular arcs
      actis(1) = radius
      xcntr(1) = xcyl
      ycntr(1) = ycyl

      sinit = 0.0D0  !  initialize arc length
      Ic = 0         !  collocation point counter

c---
c generate end-points
c---

      dth = pi2/NE(1)

      Do i=1,NE(1)+1
       theta = (i-1.0)*Dth
       Te(i) = theta
       Xe(i) = xcyl + radius*Dcos(theta)
       Ye(i) = ycyl + radius*Dsin(theta)
       Se(i) = theta*radius
      End Do
      
c---
c generate mid-points
c---

      Do i=1,NE(1)
       theta = (i-0.5D0)*Dth
       Tm(i) = theta
       Xm(i) = xcyl + radius*Dcos(theta)
       Ym(i) = ycyl + radius*Dsin(theta)
       Sm(i) = theta*radius
      End Do

      Do i=1,NE(1)+1
        xw(1,i) = xe(i)
        yw(1,i) = ye(i)
        tw(1,i) = te(i)
      End Do

      Do i=1,NE(1)
        Ic = Ic + 1
          x0(Ic) = xm(i)   ! collocation points
          y0(Ic) = ym(i)
          t0(Ic) = tm(i)
          s0(Ic) = sm(i)
        tnx0(Ic) =-Dsin(t0(Ic))
        tny0(Ic) = Dcos(t0(Ic))
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = -shrt*Y0(Ic)
         uy0(Ic) = 0.0D0
      End Do

      Ncl = NE(1)

c------------------------------------------
c A periodic array of rectangular cylinders
c------------------------------------------

      Else If(Iflow.eq.92) then

      open (4,file='rec_1p.dat')

      read (4,*) visc
      read (4,*) shrt
      read (4,*) NGL
      read (4,*) recx
      read (4,*) recy
      read (4,*) RL          ! period
      read (4,*)
      read (4,*) NE(1),RT(1)
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*) NE(4),RT(4)
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

      Iwall = 0     ! no upper wall
      NSG   = 4     ! four straight segments

      sinit = 0.0D0  !  initialize arc length
      Ic = 0         !  collocation point counter

c     write (6,*) NE(1),NE(2),NE(3),NE(4)
c     write (6,*) RT(1),RT(2),RT(3),RT(4)
c     stop

c---
c Default cylinder center is (0,0)
c Shift if desired
c---

      shiftx = 0.0D0
      shifty = 0.0D0
      shifty = recy

c---
c  first side
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1

      call elm_line
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,recx-shiftx, -recy-shifty
     +   ,recx-shiftx,  recy -shifty
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = -shrt*Y0(Ic)
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(1)+1)

c---
c second side
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line
     +
     +   (NE(2)
     +   ,RT(2)
     +   , recx-shiftx,recy-shifty
     +   ,-recx-shiftx,recy-shifty
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic)
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(2)+1)

c---
c  third side
c---

      Itp(3) = 1
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,-recx-shiftx, recy-shifty
     +   ,-recx-shiftx,-recy-shifty
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
         ux0(Ic) = - shrt*y0(Ic)
         uy0(Ic) = 0.0D0
      End Do

      sinit = se(NE(3)+1)

c---
c  fourth side
c---

      Itp(4) = 1
      Isym   = 1

      call elm_line         ! truncated wall
     +
     +   (NE(4)
     +   ,RT(4)
     +   ,-recx-shiftx,-recy-shifty
     +   , recx-shiftx,-recy-shifty
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
        ddl = Dsqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        elml(Ic) = se(i+1)-se(i)    ! element length
        ux0 (Ic) = - shrt*y0(Ic)
        uy0 (Ic) = 0.0D0
      End Do

      Ncl = Ic          ! number of collocation points

c-----------
      end if  ! End of boundary-element generation
c-----------

c-----
c done
c-----

 100  Format (1x,i4,3(1x,f15.10))

      Return
      End
