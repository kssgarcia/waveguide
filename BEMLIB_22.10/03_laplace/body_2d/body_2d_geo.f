      subroutine body_2d_geo 
     +
     +  (Ncl
     +  ,Iwall
     +  ,wall
     +  )

c=========================================
c FDLIB, BEMLIB, CFDLAB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=========================================
 
c---------------------------------------------
c Element distribution
c and computation of collocation points
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
c CAPACITY
c --------
c
c   10 boundary segments
c  128 elements per segment
c
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(129),Ye(129),Se(129),Te(129)
      Dimension Xm(128),Ym(128),Sm(128)

      Dimension    NE(10),   RT(10),  Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension  xw(10,129),yw(10,129),tw(10,129)

      Dimension  X0(1280),Y0(1280),t0(1280),s0(1280)
      Dimension      dphidn0(1280)
      Dimension         elml(1280)
      Dimension         tnX0(1280),tnY0(1280)
      Dimension         vnX0(1280),vnY0(1280)

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw
      common/xxx03/actis,xcntr,ycntr

      common/xxx04/Vx,Vy,cr,Xpv,Ypv,theta

      common/xxx05/X0,Y0,T0,S0,dphidn0
      common/xxx06/tnx0,tny0,vnx0,vny0,elml
      common/xxx07/xcenter,ycenter
      common/xxx08/xwmin,ywmin,xwmax,ywmax ! streamline window

      common/gr1/cr_new,ycenter_new        ! graphics
      common/gr2/Itry                      ! graphics

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

c-------
c CIRCLE
c-------

      If(Iflow.eq.50) then

      open (4,file='circle.dat')

      read (4,*) NGL
      read (4,*) rad
      read (4,*) xcenter,ycenter
      read (4,*) Vx,Vy
      read (4,*) cr
      read (4,*) 
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) 
      read (4,*) NE(1)    ! one segment consisting of arc elms
      read (4,*) 
      read (4,*) xwmin,xwmax   ! streamline window
      read (4,*) ywmin,ywmax   ! streamline window

c----------
c over-ride
c----------

      If(Itry.gt.1) then               ! graphics
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      End If                           ! graphics

c-------------------------------------
c place the point vortex at the center
c of the cylinder
c-------------------------------------

      Xpv = xcenter
      Ypv = ycenter

c---
c distribute points around the circle
c in the counteclockwise direction
c---

      dth = pi2/NE(1)

      Do i=1,NE(1)+1
        angle = (i-1.0D0)*dth
        te(i) = angle
        se(i) = angle*rad
        xe(i) = xcenter+rad*Dcos(angle)
        ye(i) = ycenter+rad*Dsin(angle)
      End Do

c---
c prepare
c---

      NSG = 1            ! one circular segment
      actis(1) = rad
      xcntr(1) = xcenter
      ycntr(1) = ycenter

c---
c circle contour
c---

      Itp(1) = 2           ! circular arc

      Do i=1,NE(1)+1
       tw(1,i) = te(i)
       xw(1,i) = Xe(i)
       yw(1,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Ncl = NE(1)        ! number of collocation points

      Do i=1,NE(1)        ! collocation points

          t0(i) = 0.5D0*(te(i)+te(i+1))
          x0(i) = xcenter+rad*Dcos(t0(i)) 
          y0(i) = ycenter+rad*Dsin(t0(i)) 
          s0(i) = 0.5*(se(i)+se(i+1))
        elml(i) = dth*rad
        tnx0(i) =-Dsin(t0(i))
        tny0(i) = Dcos(t0(i))
        vnx0(i) = tny0(i)
        vny0(i) =-tnx0(i)

        dphidn0(i) = -Vx*vnx0(i)-Vy*vny0(i)

        xpvd  = x0(i)-xpv              ! point vortex
        ypvd  = y0(i)-ypv
        rpvds = xpvd**2+ypvd**2
        dphidn0(i) = dphidn0(i) 
     +              -cr*(-ypvd*vnx0(i)+xpvd*vny0(i))
     +             /(pi2*rpvds)

        If(Iwall.eq.1) then            ! image point vortex
         xpvdi = xpvd
         ypvdi = y0(i)+ypv-2.0*wall
         rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(i) = dphidn0(i) 
     +    +cr*(-ypvdi*vnx0(i)+xpvdi*vny0(i))/(pi2*rpvdsi)
        End If

      End Do

c------------------
c RECTANGULAR BLOCK 
c
c Important:
c segments must be arranged
c in the countercloskwise direction
c-------------------

      elseif(Iflow.eq.51) then

      open (4,file='rectangle.dat')

      read (4,*) NGL
      read (4,*) sizex
      read (4,*) sizey
      read (4,*) xcenter,ycenter
      read (4,*) theta
      read (4,*) Vx,Vy
      read (4,*) cr
      read (4,*)
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) 
      read (4,*) NE(1),RT(1)
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*) NE(4),RT(4)
      read (4,*) 
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c----------
c overide
c----------

      if(Itry.gt.1) then
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      end if

      V1x = sizex+xcenter
      V1y =-sizey+ycenter

      V2x = sizex+xcenter
      V2y = sizey+ycenter

      V3x =-sizex+xcenter
      V3y = sizey+ycenter

      V4x =-sizex+xcenter
      V4y =-sizey+ycenter

      theta = theta*pi
      cs = cos(theta)
      sn = sin(theta)

      P1x = V1x*cs+V1y*sn
      P1y =-V1x*sn+V1y*cs

      P2x = V2x*cs+V2y*sn
      P2y =-V2x*sn+V2y*cs

      P3x = V3x*cs+V3y*sn
      P3y =-V3x*sn+V3y*cs

      P4x = V4x*cs+V4y*sn
      P4y =-V4x*sn+V4y*cs

c---------------------------------
c place point vortex at the center
c of the rectangle
c---------------------------------

      Xpv = xcenter
      Ypv = ycenter

c-------------
c preparations
c-------------

      NSG   = 4
      Ic    = 0       ! collocation points counter
      sinit = 0.0D0   ! initialize arc length

c---
c side # 1  (right side)
c---

      Itp(1) = 1    ! straight elements
      Isym   = 1

      call elm_line
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,P1x,P1y
     +   ,P2x,P2y
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(1)+1
        XW(1,i) = Xe(i)
        YW(1,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Do i=1,NE(1)
        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx*ddx+ddy*ddy)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd = x0(Ic)-xpv    ! point vortex
        ypvd = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        if(Iwall.eq.1) then  ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        end if
c---

      End Do

      sinit = se(NE(1)+1)

c---
c side # 2  (top side)
c---

      Itp(2) = 1    ! straight segment
      Isym   = 1

      call elm_line 
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,P2x,P2y
     +   ,P3x,P3y
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(2)+1
        XW(2,i) = Xe(i)
        YW(2,i) = Ye(i)
      End Do

c---
c collocation points
c---

      Do i=1,NE(2)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv             ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then            ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        End If
c---

      End Do

      sinit = se(NE(2)+1)

c---
c  side # 3 (left side)
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1

      call elm_line
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,P3x,P3y
     +   ,P4x,P4y
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

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv           ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then          ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        End If
c---

      End Do

      sinit = se(NE(3)+1)

c---
c  side # 4
c---

      Itp(4) = 1    ! straight segment
      Isym   = 1

      call elm_line
     +
     +   (NE(4)
     +   ,RT(4)
     +   ,P4x,P4y
     +   ,P1x,P1y
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

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv        ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then      ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        End If
c---

      End Do

      Ncl = Ic          ! number of collocation points

c----------------
c TRIANGULAR BLOCK 
c
c Important:
c Segments arranged in the counterclockwise sense
c
c Normal vector points into the fluid
c-----------------

      Else If(Iflow.eq.52) then

      open (4,file='triangle.dat')

      read (4,*) NGL
      read (4,*) xfirst, yfirst
      read (4,*) xsecond,ysecond
      read (4,*) xthird, ythird
      read (4,*) Vx,Vy
      read (4,*) cr
      read (4,*) 
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) 
      read (4,*) NE(1),RT(1)
      read (4,*) NE(2),RT(2)
      read (4,*) NE(3),RT(3)
      read (4,*) 
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c----------
c over-ride
c----------

      If(Itry.gt.1) then
        cr      = cr_new                   ! graphics
        ycenter = ycenter_new              ! graphics
        yfirst  = yfirst  + ycenter_new    ! graphics
        ysecond = ysecond + ycenter_new    ! graphics
        ythird  = ythird  + ycenter_new    ! graphics
      End If

c---------------------------------
c place point vortex at the center
c of the triangle
c---------------------------------

      Xpv = (xfirst+xsecond+xthird)/3.0
      Ypv = (yfirst+ysecond+ythird)/3.0

c---
c preparations
c---

      NSG   = 3
      Ic    = 0    ! collocation point counter
      sinit = 0.   ! initialize arc length

c---
c side # 1
c---

      Itp(1) = 1    ! straight segment
      Isym   = 1

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

c---
c collocation points
c---

      Do i=1,NE(1)

        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv            ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then          ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        End If
c---

      End Do

      sinit = se(NE(1)+1)

c---
c side # 2
c---

      Itp(2) = 1    ! straight segment
      Isym   = 1

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

c---
c collocation points
c---

      Do i=1,NE(2)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv            ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then           ! image point vortex
         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)
        End If
c---

      End Do

      sinit = se(NE(2)+1)

c---
c side # 3
c---

      Itp(3) = 1    ! straight segment
      Isym   = 1

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

c---
c collocation points
c---

      Do i=1,NE(3)

        Ic = Ic + 1

        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = sm(i)

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv           ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then         ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        End If
c---

      End Do

      Ncl = Ic          ! number of collocation points

c------------------------------
c AIRFOIL
c
c Important:
c
c Points should be distributed
c in the counterclockwise direction
c around the airfoil
c------------------------------

      Else If(Iflow.eq.53) then

      open (4,file='airfoil.dat')

      read (4,*) NGL
      read (4,*) Vx,Vy
      read (4,*) cr
      read (4,*) 
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) 
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c----------
c over-ride
c----------

      If(Itry.gt.1) then
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      End If

c----------------
c points arranged
c in the counterclockwise sense
c
c point vortex placed at the
c contour centroid
c----------------

      Xpv = 0.0D0      ! point vortex
      Ypv = 0.0D0      ! point vortex

      read (4,*) 
      read (4,*) NE(1)    ! one segment consisting of linear elms
      Do i=1,NE(1)
        read (4,*) idle,xe(i),ye(i)
        xe(i) =          - 0.5D0*xe(i)  ! scale
        ye(i) =  ycenter + 0.5D0*ye(i)  ! scale
        Xpv = Xpv+xe(i)
        Ypv = Ypv+ye(i)
      End Do

      xe(NE(1)+1) = xe(1)   ! wrap around
      ye(NE(1)+1) = ye(1)

      Xpv = Xpv/NE(1)
      Ypv = Ypv/NE(1)

c--------
c prepare
c--------

      NSG = 1       ! one segments

      Itp(1) = 1    ! straight elements

      se(1) = 0.0D0   ! arc length at end points

      Do i=1,NE(1)
       xw(1,i) = Xe(i)
       yw(1,i) = Ye(i)
       se(i+1) = se(i)+Dsqrt((xe(i+1)-xe(i))**2+(ye(i+1)-ye(i))**2)
      End Do
      xw(1,NE(1)+1) = xw(1,1)
      yw(1,NE(1)+1) = yw(1,1)

      Do i=1,NE(1)                    ! middle points
        xm(i) = 0.5*(xw(1,i+1)+xw(1,i))
        ym(i) = 0.5*(yw(1,i+1)+yw(1,i))
      End Do

c-------------------
c collocation points
c-------------------

      Ncl = NE(1)

      Ic = 0    ! collocation point counter

      Do i=1,NE(1)

        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        s0(Ic) = 0.5D0*(se(i)+se(i+1))

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv              ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then            ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        End If
c---

      End Do

      Ncl = Ic          ! number of collocation points

c------------------------------
c RANKINE BODY
c
c Important:
c
c Points should be distributed
c in the counterclockwise direction
c around the airfoil
c------------------------------

      Else If(Iflow.eq.54) then

      open (4,file='rankine.dat')

      read (4,*) NGL
      read (4,*) Vx,Vy
      read (4,*) cr
      read (4,*)
      read (4,*) arank
      read (4,*) brank
      read (4,*)
      read (4,*) xcenter
      read (4,*) xcenter
      read (4,*) 
      read (4,*) Iwall
      read (4,*) wall
      read (4,*) 
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c----------
c over-ride
c----------

      If(Itry.gt.1) then
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      End If

c----------------
c points arranged
c in the counterclockwise sense
c
c point vortex placed at the
c contour centroid
c----------------

      Xpv = xcenter      ! point vortex
      Ypv = ycenter      ! point vortex

      read (4,*) 
      read (4,*) NE(1)    ! one segment consisting of linear elms

      AS = arank**2
      BS = brank**2

      dth = pi2/NE(1) 

      Do i=1,NE(1)+1
        pan = (i-1.0D0)*dth
        rad = Dsqrt(AS*BS/(AS*Dsin(pan)**2+BS*Dcos(pan)**2))
        te(i) = angle
        xe(i) = xcenter + rad*Dcos(pan)
        ye(i) = ycenter + rad*Dsin(pan)
      End Do

c--------
c prepare
c--------

      NSG = 1       ! one segments

      Itp(1) = 1    ! straight elements

      se(1) = 0.0D0   ! arc length at end points

      Do i=1,NE(1)
       xw(1,i) = Xe(i)
       yw(1,i) = Ye(i)
       se(i+1) = se(i)+Dsqrt((xe(i+1)-xe(i))**2+(ye(i+1)-ye(i))**2)
      End Do
      xw(1,NE(1)+1) = xw(1,1)
      yw(1,NE(1)+1) = yw(1,1)

      Do i=1,NE(1)                    ! middle points
        xm(i) = 0.5*(xw(1,i+1)+xw(1,i))
        ym(i) = 0.5*(yw(1,i+1)+yw(1,i))
      End Do

c-------------------
c collocation points
c-------------------

      Ncl = NE(1)

      Ic = 0    ! collocation point counter

      Do i=1,NE(1)

        Ic = Ic + 1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
	s0(Ic) = 0.5D0*(se(i)+se(i+1))

        ddx = xe(i+1)-xe(i)
        ddy = ye(i+1)-ye(i)
        ddl = Dsqrt(ddx**2+ddy**2)

        elml(Ic) = ddl
        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)
        dphidn0(Ic) = -Vx*vnx0(Ic)-Vy*vny0(Ic)

        xpvd  = x0(Ic)-xpv              ! point vortex
        ypvd  = y0(Ic)-ypv
        Rpvds = xpvd**2+ypvd**2
        dphidn0(Ic) = dphidn0(Ic) 
     +    -cr*(-ypvd*vnx0(Ic)+xpvd*vny0(Ic))/(pi2*Rpvds)

c---
        If(Iwall.eq.1) then            ! image point vortex

         xpvdi = xpvd
         ypvdi = y0(Ic)+ypv-2.0*wall
         Rpvdsi = xpvdi**2+ypvdi**2
         dphidn0(Ic) = dphidn0(Ic) 
     +    +cr*(-ypvdi*vnx0(Ic)+xpvdi*vny0(Ic))/(pi2*Rpvdsi)

        End If
c---

      End Do


c-----------
      End If       ! End of geometry module
c-----------

c-----
c Done
c-----

      Return
      End
