      subroutine body_ax_geo
     +
     +  (Ncl
     +  )

c==========================================
c FDLIB, CFDLAB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c==========================================

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
c  10 segments
c 128 elements per segment
c
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(129),Ye(129),Te(129),Se(129)
      Dimension Xm(128),Ym(128),Sm(28)

      Dimension    NE(10),   RT(10),  Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,129),yw(10,129),tw(10,129)

      Dimension X0(1280),Y0(1280),t0(1280),s0(1280)
      Dimension     dphidn0(1280)
      Dimension        arel(1280)
      Dimension        tnX0(1280),tnY0(1280)
      Dimension        vnX0(1280),vnY0(1280)

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw
      common/xxx03/actis,xcntr,ycntr

      common/xxx04/Vx,cr,Xlvr,Ylvr

      common/xxx05/X0,Y0,T0,S0,dphidn0
      common/xxx06/tnx0,tny0,vnx0,vny0,arel
      common/xxx07/xcenter,ycenter
      common/xxx08/xwmin,ywmin,xwmax,ywmax

      common/gr1/cr_new,ycenter_new,Itry   ! graphics

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

c--------
c  SPHERE
c--------

      If(Iflow.eq.50) then

      open (4,file='sphere.dat')

      read (4,*) NGL
      read (4,*) rad
      read (4,*) xcenter
      read (4,*) Vx       ! velocity of indident flow
      read (4,*) cr       ! line vortex ring strength
      read (4,*) 
      read (4,*) NE(1)    ! one segment consisting of arc elms
      read (4,*)
      read (4,*) xwmin,xwmax
      read (4,*) ywmin,ywmax

c----------
c over-ride
c----------

      If(Itry.gt.1) then               ! graphics
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      End If                           ! graphics

c--------------------------------
c place the lvr inside the center
c--------------------------------

      Xlvr = xcenter
      Ylvr = 0.5D0*rad

c--- 
c one semi-circular contour
c with evenly distributed arcs
c--- 

      ycenter = 0.0D0        ! sphere center is on the x axis

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

      NSG = 1           ! one circular segment
      actis(1) = rad
      xcntr(1) = xcenter
      ycntr(1) = ycenter

c---
c semicircular contour
c---

      Itp(1) = 2    ! circular arcs

      Do i=1,NE(1)+1 
       tw(1,i) = te(i)
       xw(1,i) = Xe(i)
       yw(1,i) = Ye(i)
      End Do

c--- 
c Collocation points
c--- 

      Ncl = NE(1)

      Do i=1,NE(1)        ! collocation points

        t0(i) = 0.5D0*(te(i)+te(i+1))
        x0(i) = xcenter+rad*Dcos(t0(i)) 
        y0(i) = ycenter+rad*Dsin(t0(i)) 
	s0(i) = 0.5D0*(se(i)+se(i+1))

        arel(i) = dth*rad*pi2*y0(i)
        tnx0(i) =-Dsin(t0(i))
        tny0(i) = Dcos(t0(i))
        vnx0(i) = tny0(i)
        vny0(i) =-tnx0(i)

        dphidn0(i) = -Vx*vnx0(i)

        Iopt = 1

        call lvr_fs
     +    (Iopt
     +    ,x0(i),y0(i)
     +    ,Xlvr,Ylvr
     +    ,ulvr,vlvr
     +    ,psi
     +    )
 
        dphidn0(i) = dphidn0(i)
     +              -cr*(ulvr*vnx0(i)+vlvr*vny0(i))

      End Do

c------------------------------------
c Flow past a triangular torus
c
c Important:
c
c Segments should be arranged in the
c counterclockwise sense
c
c Normal vector points into the fluid
c------------------------------------

      Else If(Iflow.eq.51) then

      open (4,file='torus_trgl.dat')

      read (4,*) NGL
      read (4,*) xfirst,yfirst     ! first vertex
      read (4,*) xsecond,ysecond   ! second vertex
      read (4,*) xthird,ythird     ! third vertex
      read (4,*) Vx
      read (4,*) cr                ! line vortex ring strength
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
        cr      = cr_new               ! graphics
        ycenter = ycenter_new          ! graphics
      End If

c------------------------------
c place the lvr at the centroid
c of the triangle
c------------------------------

      Xlvr = (xfirst+xsecond+xthird)/3.0D0
      Ylvr = (yfirst+ysecond+ythird)/3.0D0

c-------------
c preparations
c-------------

      NSG   = 3
      Ic    = 0       ! collocation point counter
      sinit = 0.0D0   ! initialize arc length

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

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        arel(Ic) = ddl*pi2*y0(Ic)

        dphidn0(Ic) = -Vx*vnx0(Ic)

        Iopt = 1

        call lvr_fs
     +
     +     (Iopt
     +     ,x0(Ic),y0(Ic)
     +     ,Xlvr,Ylvr
     +     ,ulvr,vlvr
     +     ,psi
     +     )
 
        dphidn0(Ic) = dphidn0(Ic)
     +              -cr*(ulvr*vnx0(Ic)+vlvr*vny0(Ic))

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

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        arel(Ic) = ddl*pi2*y0(Ic)

        dphidn0(Ic) = -Vx*vnx0(Ic)

        call lvr_fs
     +
     +    (Iopt
     +    ,x0(Ic),y0(Ic)
     +    ,Xlvr,Ylvr
     +    ,ulvr,vlvr
     +    ,psi
     +    )
 
        dphidn0(Ic) = dphidn0(Ic)
     +              -cr*(ulvr*vnx0(Ic)+vlvr*vny0(Ic))
 
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

        tnx0(Ic) = ddx/ddl
        tny0(Ic) = ddy/ddl
        vnx0(Ic) = tny0(Ic)
        vny0(Ic) =-tnx0(Ic)

        arel(Ic) = ddl*pi2*y0(Ic)

        dphidn0(Ic) = -Vx*vnx0(Ic)

        call lvr_fs
     +
     +   (Iopt
     +   ,x0(Ic),y0(Ic)
     +   ,Xlvr,Ylvr
     +   ,ulvr,vlvr
     +   ,psi
     +   )
 
        dphidn0(Ic) = dphidn0(Ic)
     +              -cr*(ulvr*vnx0(Ic)+vlvr*vny0(Ic))

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
