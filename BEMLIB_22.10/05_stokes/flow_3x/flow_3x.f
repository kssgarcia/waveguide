      program flow_3x

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------------
c Shear flow past a cavity, orifice, or protrusion
c on a plane wall,
c possibly in the presence of an upper parallel wall.
c
c  ^ y
c  |           |
c  |_____      |
c  |     |     |
c  -----------------> x
c   _____|     |
c  |           |
c  |           |
c
c The trace of the boundary in the xy plane
c consists of a number of segments
c that can be straight lines or circular arcs.
c
c Each segment is discretized into a number of elements
c with corresponding shapes.
c
c SYMBOLS:
c --------
c
c NSG: Number of segments
c
c Itp(i): Index for shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c NE(i): Number of elements on the ith segment 
c
c RT(i): Stretch ratio of elements on ith segment 
c
c X0, Y0: Coordinates of collocation points
c
c T0:     angle of a collocation point
c         subtended from a circular segment center
c
c Xop,Yop: coordinates of the element end-points on 
c          the base of a protuberance or opening
c          of a cavity
c
c opening: radius of the orifice or base
c
c fopx, fops, fopf: Fourier coefficents of the traction
c                   of the unperturbed flow
c                   over the opening
c
c NGL: Number of Gauss-Legendre points
c      for integration over each element
c
c shrt:	    shear rate of incident flow at the wall
c
c delta:    pressure gradient of incident flow
c
c fx, fs, ff:	polar cylindrical components of the first
c	        Fourier traction coefficients
c               (with capital letters in the references)
c
c----------------------------

      Implicit Double Precision (a-h,o-z) 

      Dimension Xe(129),Ye(129),Te(129),Se(129)
      Dimension Xm(128),Ym(128),Tm(128),Sm(128)

      Dimension Xop (129),Yop (129)
      Dimension fopx(129),fops(129),fopf(129)

      Dimension NE(10),RT(10),Itp(10),actis(10),xcntr(10),ycntr(10)

      Dimension  XW(10,129), YW(10,129),TW(10,129),SW(10,129)
      Dimension  fx(10,128), fs(10,128),ff(10,128)
      Dimension fty(10,128),fny(10,128),ftz(10,128)

      Dimension X0(500),Y0(500),T0(500),S0(0:500)
      Dimension tnX0(500),tnY0(500)

      Dimension AL(500,500),BL(500),SOL(500)    ! for the linear system

      Dimension xstr(900),ystr(900),zstr(900)   ! for streamlines

      Dimension ZZ(20),WW(20)

      Parameter (tol=0.00000001)

c--------------
c common blocks
c--------------

      common/VEL00/Iflow,NEop,NSG,NGL,NE,Itp
      common/VEL01/opening,Xop,Yop,fopx,fops,fopf
      common/VEL02/XW,YW,TW,fx,fs,ff
      common/VEL03/actis,xcntr,ycntr
      common/VEL04/visc,shrt,delta

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358 
      pi2 = 2.0*pi
      pi4 = 4.0*pi
      pi6 = 6.0*pi
      pi8 = 8.0*pi
      piH = 0.5*pi

      Null = 0
      None = 1

      zero = 0.0

c------------
c preferences
c------------

  93  Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) "  1 for spherical cavity"
      write (6,*) "  2 for spherical cavity in a channel"
      write (6,*) " 11 for cylindrical cavity"
      write (6,*) " 12 for cylindrical cavity in a channel"
      write (6,*) " 21 for orifice of zero thickness"
      write (6,*) " 22 for orifice of zero thickness in a channel"
      write (6,*) " 31 for orifice of finite thickness"
      write (6,*) " 32 for orifice of finite thickness in a channel"
      write (6,*) " 41 for spherical protrusion"
      write (6,*) " 42 for spherical protrusion in a channel"
      write (6,*) " 51 for cylindrical protrusion"
      write (6,*) " 52 for cylindrical protrusion in a channel"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " ---------"

      read (5,*) Iflow

      If(Iflow.eq.0) Go to 99

c---
c trap
c---

      if   (Iflow.ne.1.and.Iflow.ne.2
     + .and.Iflow.ne.11.and.Iflow.ne.12
     + .and.Iflow.ne.21.and.Iflow.ne.22
     + .and.Iflow.ne.31.and.Iflow.ne.32
     + .and.Iflow.ne.41.and.Iflow.ne.42
     + .and.Iflow.ne.51.and.Iflow.ne.52
     +     ) then
         write (6,*)
         write (6,*)  " This selection is not available;"
         write (6,*)  " Please try again"
         write (6,*)
         Go to 93
       end if
        
c---------------------------------------------------
c  element distribution
c
c  The input data will be read from individual files
c---------------------------------------------------

c=================
c spherical cavity
c=================

      if(Iflow.eq.1.or.Iflow.eq.2) then

       open (4,file='cvt_sph.dat')

        read (4,*)
        read (4,*)
        read (4,*) visc     ! viscosity
        read (4,*) shrt     ! shear rate
        read (4,*) delta    ! pressure gradient
        read (4,*) NGL      ! Gauss-Legendre points
        read (4,*) radius   ! cavity radius
        read (4,*) angle    ! cavity semi-angle
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2),truncw
        read (4,*) NEop,RTop
        read (4,*) 

        If(Iflow.eq.2) then    ! upper wall located at x = wall
         read (4,*) wall
         read (4,*) NE(3),RT(3),truncw_u
        End If

       close (4)

c---
c initialize and prepare
c---

      angle   = angle*pi
      opening = radius*sin(angle)

      Ic    = 0          ! counts collocation points
      sinit = 0.0        ! origin of arc length

c---
c cavity
c---

      Itp(1) = 2       ! circular segment
      Isym   = 0       ! non-symmetric element distribution

      actis(1) = radius
      angle1   = pi 
      angle2   = pi-angle

      xstart = radius*(cos(angle)-1.0)
      ystart = 0.0

      xcnt = xstart-radius*cos(angle1)   ! cavity center
      ycnt = ystart-radius*sin(angle1)

      xcntr(1) = xcnt
      ycntr(1) = ycnt

      call elm_arc        ! cavity
     +
     +  (NE(1)
     +  ,RT(1)
     +  ,Xcnt,Ycnt
     +  ,radius
     +  ,angle1,angle2
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,Te,se
     +  ,Xm,Ym,Tm,sm
     +  )

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
        s0  (Ic) = sm(i)
        tnx0(Ic) = sin(t0(Ic))
        tny0(Ic) =-cos(t0(Ic))
      End Do

      sinit = se(NE(1)+1)

      angle = angle/pi     ! unreduce for recording

c---
c  truncated wall
c---

      Itp(2) = 1           ! straight segments
      Isym   = 0           ! non-symmetric element distribution
      trunc  = truncw

      call elm_line
     +
     +  (NE(2)
     +  ,RT(2)
     +  ,zero,opening
     +  ,zero,trunc
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
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      NSG = 2               ! number of segments

c---
c  truncated upper wall
c---

      If(Iflow.eq.2) then

        Itp(3) = 1
        Isym   = 0
        sinit  = 0.0

        trunc  = truncw_u

        call elm_line
     +
     +    (NE(3)
     +    ,RT(3)
     +    ,wall,zero
     +    ,wall,trunc
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
          x0(Ic)   = xm(i)
          y0(Ic)   = ym(i)
          s0(Ic)   = sm(i)
          ddx      = xe(i+1)-xe(i)
          ddy      = ye(i+1)-ye(i)
          den      = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
        End Do

        NSG = 3

      End If

c=================
c cylindrical cavity 
c=================

      elseif(Iflow.eq.11.or.Iflow.eq.12) then

      open (4,file='cvt_cyl.dat')

        read (4,*) 
        read (4,*) 
        read (4,*) visc
        read (4,*) shrt
        read (4,*) delta
        read (4,*) NGL
        read (4,*) 
        read (4,*) depth
        read (4,*) opening
        read (4,*) NEop,RTop
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3),truncw
        read (4,*) 

        If(Iflow.eq.12) then
         read (4,*) wall
         read (4,*) NE(4),RT(4),truncw_u
        End If

      close (4)

c---
c initialize
c---

      depthm = -depth

      Ic    = 0          ! counts collocation points
      sinit = 0.0        ! origin of arc length

c---
c  bottom of the cavity
c---

      Itp(1) = 1           ! staright elements
      Isym   = 0           ! non-symmetric element distribution

      call elm_line            ! bottom of cavity
     +
     +    (NE(1)
     +    ,RT(1)
     +    ,depthm,zero
     +    ,depthm,opening
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
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      sinit = se(NE(1)+1)

c---
c  side of the cavity
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line            ! side of cavity
     +
     +  (NE(2)
     +  ,RT(2)
     +  ,depthm,opening
     +  ,zero  ,opening
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
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      sinit = se(NE(2)+1)

c---
c  truncated wall
c---

      trunc = truncw

      Itp(3) = 1
      Isym   = 0

      call elm_line           ! truncated wall
     +
     +  (NE(3)
     +  ,RT(3)
     +  ,zero,opening
     +  ,zero,trunc
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
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      NSG = 3               ! number of segments

c---
c  truncated upper wall
c---

      If(Iflow.eq.12) then

        trunc = truncw_u

        Itp(4) = 1
        Isym   = 0
        sinit  = 0.0

        call elm_line         ! truncated upper wall
     +
     +    (NE(4)
     +    ,RT(4)
     +    ,wall,zero
     +    ,wall,trunc
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
          x0(Ic)   = xm(i)
          y0(Ic)   = ym(i)
          s0(Ic)   = sm(i)
          ddx      = xe(i+1)-xe(i)
          ddy      = ye(i+1)-ye(i)
          den      = sqrt(ddx*ddx+ddy*ddy)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
        End Do

        NSG = 4

      End If

c=================
c orifice of zero thickness
c=================

      elseif(Iflow.eq.21.or.Iflow.eq.22) then

      open (4,file='orifice_0.dat')

        read (4,*)
        read (4,*) visc
        read (4,*) shrt
        read (4,*) delta
        read (4,*) NGL
        read (4,*)
        read (4,*) opening
        read (4,*) NEop,RTop
        read (4,*) 
        read (4,*) NE(1),RT(1),truncw
        read (4,*) 

        If(Iflow.eq.22) then
         read (4,*)
         read (4,*) wall
         read (4,*) NE(2),RT(2),truncw_u
        End If

      close (4)

c---
c initialize
c---

      Ic    = 0          ! counts collocation points
      sinit = 0.0        ! origin of arc length
      Itp(1) = 1       
      Isym   = 0     ! element distribution is not symmetric

      sinit  = 0.0

      call elm_line           ! truncated wall
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,zero,opening
     +   ,zero,truncw
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(1)+1
        XW(1,i) = Xe(i)
        YW(1,i) = Ye(i)
      End Do

      Ic = 0          ! counts collocation points

      Do i=1,NE(1)
        Ic       = Ic + 1
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      NSG = 1               ! number of segments

      If(Iflow.eq.22) then

        Itp(2) = 1
        Isym   = 0
        sinit  = 0.0

        call elm_line         ! truncated upper wall
     +
     +    (NE(2)
     +    ,RT(2)
     +    ,wall,zero
     +    ,wall,truncw_u
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
          den = sqrt(ddx**2+ddy**2)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
        End Do

        NSG = 2

      End If

c=================
c orifice of finite thickness
c=================

      elseif(Iflow.eq.31.or.Iflow.eq.32) then

      open (4,file='orifice.dat')

        read (4,*)
        read (4,*)
        read (4,*) visc
        read (4,*) shrt
        read (4,*) delta
        read (4,*) NGL
        read (4,*)
        read (4,*) depth
        read (4,*) opening
        read (4,*) NEop,RTop
        read (4,*) 
        read (4,*) NE(1),RT(1),truncw1
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3),truncw3
        read (4,*) 

        If(Iflow.eq.32) then
         read (4,*) wall
         read (4,*) NE(4),RT(4),truncw_u
        End If

      close (4)

c---
c initialize
c---

      depthm = -depth
      Ic     = 0          ! counts collocation points
      sinit  = 0.0

c---
c  lower side of the orifice
c---

      Itp(1) = 1
      Isym   = 0

      call elm_line            ! lower surface of orifice
     +
     +  (NE(1)
     +  ,RT(1)
     +  ,depthm,truncw1
     +  ,depthm,opening
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
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      sinit = se(NE(1)+1)

c---
c  side of the orifice
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line            ! side of orifice
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,depthm,opening
     +   ,zero  ,opening
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
        den = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      sinit = se(NE(2)+1)

c---
c  upper side of the orifice
c---

      Itp(3) = 1
      Isym   = 0

      call elm_line           ! truncated wall
     +
     +  (NE(3)
     +  ,RT(3)
     +  ,zero,opening
     +  ,zero,truncw3
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
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx**2+ddy**2)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      NSG = 3               ! number of segments

      If(Iflow.eq.33) then

        Itp(4) = 1
        Isym   = 0
        sinit  = 0.0

        call elm_line         ! truncated upper wall
     +
     +     (NE(4)
     +     ,RT(4)
     +     ,wall,zero
     +     ,wall,truncw_u
     +     ,sinit
     +     ,Isym
     +     ,Xe,Ye,se
     +     ,Xm,Ym,sm
     +     )

        Do i=1,NE(4)+1
          XW(4,i) = Xe(i)
          YW(4,i) = Ye(i)
        End Do

        Do i=1,NE(4)
          Ic = Ic + 1
          x0(Ic)   = xm(i)
          y0(Ic)   = ym(i)
          s0(Ic)   = sm(i)
          ddx      = xe(i+1)-xe(i)
          ddy      = ye(i+1)-ye(i)
          den      = sqrt(ddx*ddx+ddy*ddy)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
        End Do

        NSG = 4

      End If

c=====================
c spherical protrusion
c=====================

      elseif(Iflow.eq.41.or.Iflow.eq.42) then

      open (4,file='protr_sph.dat')

        read (4,*) 
        read (4,*) 
        read (4,*) visc
        read (4,*) shrt
        read (4,*) delta
        read (4,*) NGL
        read (4,*) radius
        read (4,*) angle
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2),truncw
        read (4,*) NEop,RTop

        If(Iflow.eq.42) then
         read (4,*) 
         read (4,*) wall
         read (4,*) NE(3),RT(3),truncw_u
        End If

      close (4)

c---
c initialize
c---
      Ic = 0          ! counts collocation points
      sinit  = 0.0

c---
c  protrusion
c---

      Itp(1) = 2
      Isym   = 0

      actis(1) = radius
      opening  = radius *sin(angle*pi)

      angle1 = 0.0      
      angle2 = angle*pi     ! in multiples of pi

      xstart = 1.0-cos(angle*pi)
      ystart = 0.0

      xcnt = xstart-radius*cos(angle1)
      ycnt = ystart-radius*sin(angle1)

      xcntr(1) = xcnt
      ycntr(1) = ycnt

      call elm_arc
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,Xcnt,Ycnt
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
        SW(1,i) = se(i)
      End Do

      Do i=1,NE(1)
        Ic = Ic+1
        x0(Ic) = xm(i)
        y0(Ic) = ym(i)
        t0(Ic) = tm(i)
        s0(Ic) = sm(i)
        tnx0(Ic) =-sin(t0(Ic))
        tny0(Ic) = cos(t0(Ic))
      End Do

      sinit = se(NE(1)+1)

c---
c truncated wall
c---

      trunc = truncw

      Itp(2) = 1
      Isym   = 0

      call elm_line         ! truncated wall
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,zero,opening
     +   ,zero,trunc
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(2)+1
        XW(2,i) = Xe(i)
        YW(2,i) = Ye(i)
        SW(2,i) = se(i)
      End Do

      Do i=1,NE(2)
        Ic = Ic+1
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      NSG = 2               ! number of segments

c---
c upper wall located at x = wall
c---

      if(Iflow.eq.42) then

        trunc = truncw_u
 
        Itp(3) = 1
        Isym   = 0
        sinit  = 0.0

        call elm_line         ! truncated upper wall
     +
     +    (NE(3)
     +    ,RT(3)
     +    ,wall,zero
     +    ,wall,trunc
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,se
     +    ,Xm,Ym,sm
     +    )

        Do i=1,NE(3)+1
          XW(3,i) = Xe(i)
          YW(3,i) = Ye(i)
          SW(3,i) = se(i)
        End Do

        Do i=1,NE(3)
          Ic = Ic + 1
          x0(Ic)   = xm(i)
          y0(Ic)   = ym(i)
          s0(Ic)   = sm(i)
          ddx      = xe(i+1)-xe(i)
          ddy      = ye(i+1)-ye(i)
          den      = sqrt(ddx*ddx+ddy*ddy)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
        End Do

        NSG = 3

      End If

c=======================
c cylindrical protrusion
c=======================

      elseif(Iflow.eq.51.or.Iflow.eq.52) then

      open (4,file='protr_cyl.dat')

        read (4,*) 
        read (4,*) 
        read (4,*) visc
        read (4,*) shrt
        read (4,*) delta
        read (4,*) NGL
        read (4,*) height
        read (4,*) opening
        read (4,*) 
        read (4,*) NE(1),RT(1)
        read (4,*) NE(2),RT(2)
        read (4,*) NE(3),RT(3),truncw
        read (4,*) NEop,RTop
        read (4,*) 

        if(Iflow.eq.52) then  ! upper wall
         read (4,*) wall
         read (4,*) NE(4),RT(4),truncw_u
        end if

      close (4)

c-----------
c initialize
c-----------

      Ic    = 0     ! count collocation points
      sinit = 0.0   ! initial arc length

c---
c top of the cylinder 
c---

      Itp(1) = 1
      Isym   = 0

      call elm_line 
     +
     +   (NE(1)
     +   ,RT(1)
     +   ,height,zero
     +   ,height,opening
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(1)+1
        XW(1,i) = Xe(i)
        YW(1,i) = Ye(i)
        SW(1,i) = se(i)
      End Do

      Do i=1,NE(1)
        Ic = Ic + 1
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      sinit = se(NE(1)+1)

c---
c side of the cylinder
c---

      Itp(2) = 1
      Isym   = 1

      call elm_line 
     +
     +   (NE(2)
     +   ,RT(2)
     +   ,height,opening
     +   ,zero  ,opening
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(2)+1
        XW(2,i) = Xe(i)
        YW(2,i) = Ye(i)
        SW(2,i) = se(i)
      End Do

      Do i=1,NE(2)
        Ic = Ic + 1
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      sinit = se(NE(2)+1)

c---
c truncated wall
c---

      trunc = truncw
      Itp(3) = 1
      Isym   = 0

      call elm_line           ! truncated wall
     +
     +   (NE(3)
     +   ,RT(3)
     +   ,zero,opening
     +   ,zero,trunc
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do i=1,NE(3)+1
        XW(3,i) = Xe(i)
        YW(3,i) = Ye(i)
        SW(3,i) = se(i)
      End Do

      Do i=1,NE(3)
        Ic = Ic + 1
        x0(Ic)   = xm(i)
        y0(Ic)   = ym(i)
        s0(Ic)   = sm(i)
        ddx      = xe(i+1)-xe(i)
        ddy      = ye(i+1)-ye(i)
        den      = sqrt(ddx*ddx+ddy*ddy)
        tnx0(Ic) = ddx/den
        tny0(Ic) = ddy/den
      End Do

      NSG = 3               ! number of segments

c---
c  truncated upper wall
c---

      if(Iflow.eq.52) then

        trunc = truncw_u
        Itp(4) = 1
        Isym  = 0
        sinit = 0.0

        call elm_line         ! truncated upper wall
     +
     +    (NE(4)
     +    ,RT(4)
     +    ,wall,zero
     +    ,wall,trunc
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
          x0(Ic)   = xm(i)
          y0(Ic)   = ym(i)
          s0(Ic)   = sm(i)
          ddx      = xe(i+1)-xe(i)
          ddy      = ye(i+1)-ye(i)
          den      = sqrt(ddx*ddx+ddy*ddy)
          tnx0(Ic) = ddx/den
          tny0(Ic) = ddy/den
        End Do

        NSG = 4

      End If

c---
      End If                ! End of geometry
c---

c===========================
c  End of element generation 
c===========================

      open (3,file="flow_3x.out1")
      open (9,file="flow_3x.out2")

c-------------------
c collocation points
c-------------------

      Ncl = Ic

      Ncl2 = 2*Ncl
      Ncl3 = 3*Ncl

      write (6,*) Ncl," Collocation points"

c--------------------------------------------
c print wall geometry to accompany streamlines
c--------------------------------------------

      write (3,102) NSG

      Do j=1,NSG

c      write (6,102) NE(j)+1
       write (3,102) NE(j)+1
       write (9,102) NE(j)+1

       Do i=1,NE(j)+1
c       write (6,102) i,YW(j,i),XW(j,i)
        write (9,102) i,YW(j,i),XW(j,i)
        write (3,102) i,YW(j,i),XW(j,i)
       End Do

      End Do

c---------------------------
c display the collocation points
c---------------------------

      write (6,*)
      write (6,*) Ncl," Collocation points:"
      write (6,*)

      Do i=1,Ncl
        write (6,102) i,Y0(i),X0(i)
      End Do

c--------
c prepare
c--------

      call gauss_leg (NGL,ZZ,WW)

c------------------------------------
c Element distribution along the 
c cavity opening, orifice opening,
c or protrusion base
c
c Evaluate the functions fx, fs, ff
c------------------------------------
  
      Isym = 0

      call elm_line
     +
     +   (NEop
     +   ,RTop
     +   ,zero,zero
     +   ,zero,opening
     +   ,zero
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

      Do I=1,NEop+1
        Xop(i)  = Xe(i)
        Yop(i)  = Ye(i)
        fopx(i) =           - delta*Yop(i)
        fops(i) = visc*shrt + delta*Xop(i)
        fopf(i) = visc*shrt + delta*Xop(i)
      End Do

c--------------------------
c optional printing session
c
c     write (6,*) " ------- "
c     write (6,*) " OPENING "
c     write (6,*) " ------- "
c
c     write (9,102) NEop+1
c     Do i = 1,NEop+1
c       write (6,102) i,Yop(i),Xop(i)
c       write (9,102) i,Xop(i),Yop(i)
c     End Do
c
c--------------------------

c----------------------------------------------------
c  GENERATE THE LINEAR SYSTEM
c
c  For flow past a cavity or an orifice,
c  solve equations (13) and (14) of Pozrikidis (1994)
c
c  For flow past a protrusion, solve equation (18)
c  of Pozrikidis J. Eng Math. 31, 29-42 (1997)
c----------------------------------------------------

c------------------------------
c Generate the influence matrix
c------------------------------

      Do I=1,Ncl    ! loop over collocation points

         NI = Ncl +I
        NII = Ncl2+I

        J = 0              ! counter

c       write (6,*) " Collocation point :",I,t0(i)

        Do k=1,NSG         ! loop over segments

        rad  = actis(k)
        xcnt = xcntr(k)
        ycnt = ycntr(k)

        Do L=1,NE(k)       ! loop over elements

         X1 = XW(K,L)
         Y1 = YW(K,L)
         T1 = TW(K,L)

         X2 = XW(K,L+1)
         Y2 = YW(K,L+1)
         T2 = TW(K,L+1)

         J   = J+1
         Ising = 0
         If(I.eq.J)  Ising = 1

         call flow_3x_slp 
     +
     +     (X0(I),Y0(I),t0(I)
     +     ,X1,Y1,T1
     +     ,X2,Y2,T2
     +     ,NGL
     +     ,Ising
     +     ,Itp(k)
     +     ,rad,xcnt,ycnt
     +     ,QRxx,QRxs,QRsx
     +     ,QRss,QRff,QIxf
     +     ,QIfx,QIsf,QIfs
     +     )

         NJ  = Ncl +J
         NJJ = Ncl2+J

         AL(  I,  J) =   QRxx
         AL(  I, NJ) =   QRxs
         AL(  I,NJJ) = - QIxf
         AL( NI,  J) =   QRsx
         AL( NI, NJ) =   QRss
         AL( NI,NJJ) = - QIsf
         AL(NII,  J) =   QIfx
         AL(NII, NJ) =   QIfs
         AL(NII,NJJ) =   QRff

         End Do

        End Do

c----------------------------
c compute the right-hand side
c----------------------------

        BL(  I) = 0.0
        BL( NI) = 0.0
        BL(NII) = 0.0

c--------------------------------------
c contribution from the opening or base
c--------------------------------------

        if(opening.gt.0.0001) then

        Ising = 0
        Itype = 1

        Do K=1,NEop     ! loop over opening elements

          X1 = Xop(K)
          Y1 = Yop(K)
          X2 = Xop(K+1)
          Y2 = Yop(K+1)

         call flow_3x_slp 
     +
     +       (X0(I),Y0(I),T0(I)
     +       ,X1,Y1,T1
     +       ,X2,Y2,T2
     +       ,NGL
     +       ,Ising
     +       ,Itype
     +       ,Rad,xcnt,ycnt
     +       ,QRxx,QRxs,QRsx
     +       ,QRss,QRff,QIxf
     +       ,QIfx,QIsf,QIfs
     +       )

          BL(  I) = BL(  I)
     +      +(QRxx*fopx(k)+QRxs*fops(k)-QIxf*fopf(k))
          BL( NI) = BL( NI)
     +      +(QRsx*fopx(k)+QRss*fops(k)-QIsf*fopf(k))
          BL(NII) = BL(NII)
     +      +(QIfx*fopx(k)+QIfs*fops(k)+QRff*fopf(k))

        End Do

        end if

       ! flow over a protrusion:

        if (Iflow.eq.41.or.Iflow.eq.42
     +  .or.Iflow.eq.51.or.Iflow.eq.52
     +      ) then       
          BL( NI) = BL( NI) 
     +         + pi8*(shrt*x0(i)+0.5*delta/visc*x0(i)**2)
          BL(NII) = BL(NII) 
     +         + pi8*(shrt*x0(i)+0.5*delta/visc*x0(i)**2)
        end if

c---------------
      End Do   ! over collocation points
c---------------

c------------------------
c solve the linear system
c------------------------

      write (6,*) 
      write (6,*) "Size of the linear system:",Ncl3
      write (6,*) 

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 1   ! pivoting enabled

      call gel
     +
     +   (Ncl3
     +   ,AL,BL,SOL
     +   ,Isym_g
     +   ,Iwlpvt
     +   ,Deter
     +   ,Istop
     +   )

c--------------------------
c display the linear system
c--------------------------
c
c      Do I=1,Ncl3
c        write (6,101) (AL(I,J),J=1,Ncl3),BL(I),SOL(I)
c        write (6,101) BL(i)
c      End Do

c------------------------
c distribute the solution
c------------------------

      K = 0        ! counter

      Do I=1,NSG

       Do J=1,NE(i)
        K = K+1
        fx(I,J) = SOL(K)
        fs(I,J) = SOL(K+Ncl)
        ff(I,J) = SOL(K+Ncl2)
       End Do

      End Do

c-------------------------------------------
c special treatment of flow over an orifice
c of zero thickness, because upper and lower
c surfaces are collapsed
c-------------------------------------------

      if(Iflow.eq.31.or.Iflow.eq.32) then

       Do J=1,NE(1)
        fx(1,J) = 0.5*fx(1,J)
        fs(1,J) = 0.5*fs(1,J)
        ff(1,J) = 0.5*ff(1,J)
       End Do

      end if

c----------------------------------
c compute the tangential components
c of the traction
c----------------------------------

      k = 0        ! counter

      Do i=1,NSG

       Do j=1,NE(i)
        k = k+1
        fty(i,j) = fx(i,j)*tnx0(k)+fs(i,j)*tny0(k)
        fny(i,j) =-fx(i,j)*tny0(k)+fs(i,j)*tnx0(k)
        ftz(i,j) = ff(i,j)
       End Do

      End Do

c--------------------------------------
c add shear stress of the incident flow
c over the wall
c--------------------------------------

c---
c set the segment number of upper wall
c---

      if(   Iflow.eq.01.or.Iflow.eq.02
     +  .or.Iflow.eq.41.or.Iflow.eq.42
     +  ) then
        Nwall = 2
      elseif( Iflow.eq.11.or.Iflow.eq.12
     +     .or.Iflow.eq.51.or.Iflow.eq.52
     +       ) then
        Nwall = 3
      end if

      const = visc*shrt

      Do j=1,NE(Nwall)
        fty(Nwall,j) = const + fty(Nwall,j) 
        ftz(Nwall,j) = const + ftz(Nwall,j)
      End Do

c--------------------
c record the solution
c--------------------

      write (6,*)
      write (6,*) "    X, Y, Arc Length, fx, fsigma, fphi,"
     +           ," ft_y, fn_y, ft_z "
      write (6,*)

      k = 0   ! collocation point counter

c----
      Do i=1,NSG
        write (6,109) NE(i),i
        write (3,109) NE(i),i
        Do j=1,NE(i)

         k = k+1

         sprint = s0(k)

         If(   Iflow.eq.01.or.Iflow.eq.02   ! spherical cavity
     +     .or.Iflow.eq.41.or.Iflow.eq.42   ! or protrusion
     +     )  sprint = sprint/pi

c        sprint = log10(abs(sprint-0.5))
c        yprint = log10(abs(fty(i,j)))
c        zprint = log10(abs(ftz(i,j)))
c        write (3,102) j,sprint,yprint,zprint

         write (6,102) j,X0(k),Y0(k),sprint
     +         ,fx(i,j),fs(i,j),ff(i,j)
     +         ,fty(i,j),fny(i,j),ftz(i,j)


         write (3,102) j,X0(k),Y0(k),sprint
     +         ,fx(i,j),fs(i,j),ff(i,j)
     +         ,fty(i,j),fny(i,j),ftz(i,j)

        End Do
      End Do

c----

c---------------------------------------
c Exact solution for flow over an orifice
c of zero thickness
c
c Displayed for comparison purposes
c---------------------------------------

       if(   Iflow.eq.01.or.Iflow.eq.02
     +   .or.Iflow.eq.31.or.Iflow.eq.32) then

       write (6,*)
       write (6,*) " Exact solution for simple shear flow"
       write (6,*) " over an orifice of zero thickness"
       write (6,*) " ---------------------------------"
       write (6,*)

       write (3,*) zero
       write (3,110) NE(2)

       cf = visc*shrt
       aas = opening**2

       Do J=1,NE(2)

        K  = Ncl-NE(2)+J
        RL = sqrt ((Y0(K)/opening)**2-1.0)
        AA = -atan(1.0/RL)+1.0/RL
        BB =  aas/(3.0*RL*Y0(K)**2)
        CC = cf*(AA+BB)/pi
        DD = cf*(AA-BB)/pi
        HH = CC+cf
        GG = DD+cf
        EE = 0.0

        sprint = S0(K)

        write (6,102) J,X0(K),Y0(K),Sprint,EE,CC,DD,HH,GG
        write (3,102) J,X0(K),Y0(K),Sprint,EE,CC,DD,HH,GG

       End Do

       write (3,102) Null

      end If

c----------------------------------------------
c Flow over a protuberance
c ------------------------
c
c Compute force and torque with respect to the z axis
c using equations (19) and (20) of Pozrikidis (1997)
c
c Torqo is torque with respect to the origin
c
c Torqc is torque with respect to the center 
c      of a spherical protrusion
c
c Will also compute the arc length
c and the surface area of the protrusion
c----------------------------------------------

      if(    Iflow.eq.41.or.Iflow.eq.42
     +   .or.Iflow.eq.51.or.Iflow.eq.52) then

       alngt = 0.0  ! arc length
       sarea = 0.0  ! surface area

       if(Iflow.eq.41.or.Iflow.eq.42) then
        Xctrq = Xcntr(1)
        Yctrq = Ycntr(1)
        imax = 1
       else
        Xctrq = 0.0
        Yctrq = 0.0
        imax = 2
       end if

       Force = 0.0
       Torqo = 0.0
       Torqc = 0.0

       Ic = 0

       Do i=1,imax
       Do j=1,NE(1)

        Ic=Ic+1
        contr = SW(i,j+1)-SW(i,j)
        alngt = alngt + contr
        fct   =         contr*Y0(Ic)
        sarea = sarea + fct*pi
        Force = Force + fct*(fs(i,j)+ff(i,j))

        Torqo = Torqo
     +   + fct * ( X0(Ic)*(fs(i,j)+ff(i,j)) 
     +            -Y0(Ic)*fx(i,j) )

        Torqc = Torqc 
     +   + fct * ( (x0(Ic)-Xctrq)*(fs(i,j)+ff(i,j))
     +            -(Y0(Ic)-Yctrq)* fx(i,j))

       End Do
       End Do

       Force = visc*pi*Force
       Torqo = visc*pi*Torqo
       Torqc = visc*pi*Torqc

c---
c normalize
c---

       fc = visc*pi
       alngt = alngt/fc
       sarea = sarea/fc

       Force = Force/fc
       Torqo = Torqo/fc
       Torqc = Torqc/fc

       write (6,900) sarea,alngt,Force,Torqo,Torqc

      End If

c----------------
c post-processing
c----------------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to draw streamlines"
      write (6,*) " 2 for an x-velocity profile"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c-----------------------
      if(menu.eq.1) then    ! will draw streamlines
c-----------------------

      write (6,*)
      write (6,*) " Enter maximum number of points"
      write (6,*) " before inquiring"
      write (6,*) " --------------- "
      read  (5,*) Mstr

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 2 for the second-order RK method"
      write (6,*) " 4 for the fourth-order RK method"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) IRK

      if(IRK.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to generate streamlines based on the"
      write (6,*) "   numerical solution"
      write (6,*) " 2 to generate streamlines using the"
      write (6,*) "   exact solution for flow over an orifice"
      write (6,*) "   of zero thickness"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) method

      if(method.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter the spatial step"
      write (6,*) " ----------------------"
      read  (5,*) Dl

      open (8,file="streamlines.dat")    ! file containing starting points

  22  Continue

      read (8,*) X00,Y00,Z00

      If(abs(X00-99).lt.tol) Go to 98      ! end of drawing
      If(abs(Y00-99).lt.tol) Go to 98
      If(abs(Z00-99).lt.tol) Go to 98

      Ycross = Y00    ! to be used for streamline stopping check

c---
      L = 1     ! local counter for inquiry
      K = 1     ! global counter

  20  Continue

      xstr(L) = X00
      ystr(L) = Y00
      zstr(L) = Z00

      S00 = sqrt(Y00*Y00+Z00*Z00)

      phi = acos(Y00/S00)

      If(method.eq.1) then

        call velocity 
     +
     +    (X00,S00,phi
     +    ,Ux1,Uy1,Uz1
     +    )

      else

        call davis 
     +
     +    (opening,shrt
     +    ,X00,S00,phi
     +    ,Ux1,Uy1,Uz1
     +    )

      end if

      step  = Dl/sqrt(Ux1*Ux1+Uy1*Uy1+Uz1*Uz1)  ! time step
      steph = 0.5D0*step

      Xsv = X00  ! save
      Ysv = Y00
      Zsv = Z00

c---
c second-order RK
c---

      if(IRK.eq.2) then

        X00 = Xsv + step * Ux1
        Y00 = Ysv + step * Uy1
        Z00 = Zsv + step * Uz1

        S00 = sqrt(Y00*Y00+Z00*Z00)
        phi = acos(Y00/S00)

        if(method.eq.1) then

           call velocity
     +
     +       (X00,S00,phi
     +       ,Ux2,Uy2,Uz2
     +       )
        else

           call davis 
     +
     +        (opening,shrt
     +        ,X00,S00,phi
     +        ,Ux2,Uy2,Uz2
     +        )

        end if

        X00 = Xsv + steph*(Ux1+Ux2)
        Y00 = Ysv + steph*(Uy1+Uy2)
        Z00 = Zsv + steph*(Uz1+Uz2)

      end if

c---
c fourth-order RK
c---

      If(IRK.eq.4) then

        X00 = Xsv + steph * Ux1
        Y00 = Ysv + steph * Uy1
        Z00 = Zsv + steph * Uz1

        S00 = sqrt(Y00*Y00+Z00*Z00)
        phi = acos(Y00/S00)

        if(method.eq.1) then

          call velocity 
     +
     +      (X00,S00,phi
     +      ,Ux2,Uy2,Uz2
     +      )

        else

           call davis 
     +
     +      (opening,shrt
     +      ,X00,S00,phi
     +      ,Ux2,Uy2,Uz2
     +      )

        end if

        X00 = Xsv + steph * Ux2
        Y00 = Ysv + steph * Uy2
        Z00 = Zsv + steph * Uz2

        S00 = sqrt(Y00*Y00+Z00*Z00)
        phi = acos(Y00/S00)

        if(method.eq.1) then

          call velocity 
     +
     +      (X00,S00,phi
     +      ,Ux3,Uy3,Uz3
     +      )

        else

           call davis 
     +
     +      (opening,shrt
     +      ,X00,S00,phi
     +      ,Ux3,Uy3,Uz3
     +      )

        end if

        X00 = Xsv + step * Ux3
        Y00 = Ysv + step * Uy3
        Z00 = Zsv + step * Uz3

        S00 = sqrt(Y00*Y00+Z00*Z00)
        phi = acos(Y00/S00)

        if(method.eq.1) then

           call velocity 
     +
     +       (X00,S00,phi
     +       ,Ux4,Uy4,Uz4
     +       )

        else

           call davis 
     +
     +       (opening,shrt
     +       ,X00,S00,phi
     +       ,Ux4,Uy4,Uz4
     +       )

        End If

        X00 = Xsv + step * (Ux1 +2.0*Ux2 +2.0*Ux3 +Ux4)/6.0
        Y00 = Ysv + step * (Uy1 +2.0*Uy2 +2.0*Uy3 +Uy4)/6.0
        Z00 = Zsv + step * (Uz1 +2.0*Uz2 +2.0*Uz3 +Uz4)/6.0
      End If

c------

      K = K+1
      L = L+1

      test = Ycross*Y00

      if(test.lt.0) then
        write (6,*)
        write (6,*) " Crossed the xz plane: I will stop"
        Go to 21
      end if

      If(K.lt.Mstr) Go to 20

      K = 1    ! reset local counter

      write (6,*) 
      write (6,*) " Continue this streamline ?"
      write (6,*) 
      write (6,*) " 0 for no, 1 for yes"
      write (6,*) "--------------------"

      read (5,*) Icon

      if(Icon.eq.1) Go to 20

  21  Continue

      xstr(L) = X00
      ystr(L) = Y00
      zstr(L) = Z00

c---------------------
c print the streamline
c---------------------

      write (6,*)
      write (6,*) " One streamline with ",L," points completed"
      write (6,*)

      write (9,104) L  ! number of points
      Do I=1,L
        write (9,104) I,ystr(I),xstr(I),zstr(I)
      End Do

c---
c If streamline lies off the xy plane,
c also print its three reflections
c---

      If(abs(Zstr(1)).gt.0.0000001) then

        write (9,104) L
        Do i=1,L
          YYY = -Ystr(i)
          write (9,104) i,YYY,Xstr(i),Zstr(i)
        End Do

        write (9,104) L
        Do i=1,L
          ZZZ = -Zstr(i)
          write (9,104) i,Ystr(i),Xstr(i),ZZZ
        End Do

        write (9,104) L
        Do i=1,L
          YYY = -Ystr(i)
          ZZZ = -Zstr(i)
          write (9,104) i,YYY,Xstr(i),ZZZ
        End Do

      End If

c------------

      Go to 22        ! Return for another streamline

c----------
      End If   ! end of drawing streamlines
c----------

      If(menu.eq.2) then

      open (9,file="flow_3x.out3")

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to use the numerical solution"
      write (6,*) " 2 to use the exact solution for flow "
      write (6,*) "   over an orifice of zero thickness"
      write (6,*) " 0 to quit"
      write (6,*) " --------- "
      read  (5,*) method

      If(method.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter Xmin Xmax"
      write (6,*) " ---------------"
      read  (5,*) Xmin,Xmax

      write (6,*)
      write (6,*) " How many points for the profile ?"
      write (6,*) " ---------------------------------"
      read  (5,*) Nprof

      Y00 = 0.0001   ! slightly off the axis to avoid singularities
      Z00 = 0.0

      Nprof1 = Nprof+1

c----------
c profiling
c----------

      DX0 = (Xmax-Xmin)/(Nprof+1.0-1.0)

      X00 = Xmin

      write (6,104) Nprof1

      Do i=1,Nprof1

        S00 = sqrt(Y00**2+Z00**2)
        arg = Y00/S00
        if(arg.gt. 1.0) arg =  0.999999999 
        if(arg.lt.-1.0) arg = -0.999999999 
        phi = acos(arg)

c---
        if(method.eq.1) then
c---
           call velocity 
     +       (X00,S00,phi
     +       ,Ux,Uy,Uz
     +       )
c---
        else if(method.eq.2) then
c---
           call davis 
     +       (opening,shrt
     +       ,X00,S00,phi
     +       ,Ux,Uy,Uz
     +       )

c---
        end if
c---

        X00 = X00 + DX0

        write (6,104) I,X00,Ux,Uy,Uz
        write (9,104) I,X00,Ux,Uy,Uz

      End Do

c----------
      End If     ! End of drawing a profile
c-----------------

  98  Continue

c------------------
c recording session
c for streamlines
c------------------

      write (9,104) Null
      write (9,*)
      write (9,200) Iflow
      write (9,201) visc
      write (9,202) shrt
      write (9,203) delta
      write (9,204) NGL
      write (9,205) radius
      write (9,206) angle
      write (9,207) truncw
      write (9,208) truncw_u
      write (9,*)

      Do i=1,NSG
        write (9,220) i,NE(i),RT(i)
      End Do

      write (9,221) NEop,RTop
      write (9,*)

      close (9)

  99  Continue

c------------------
c recording session
c for tractions
c------------------

      write (3,104) Null
      write (3,*)
      write (3,200) Iflow
      write (3,201) visc
      write (3,202) shrt
      write (3,203) delta
      write (3,204) NGL
      write (3,205) radius
      write (3,206) angle
      write (3,207) truncw
      write (3,208) truncw_u
      write (3,*)

      Do i = 1,NSG
       write (3,220) i,NE(i),RT(i)
      End Do

      write (3,221) NEop,RTop
      write (3,*)

      close (3)

c---
c Done
c---

 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
c109  Format (1x,i3,1x,"segment ",i2)
 109  Format (1x,i3,1x,i2)
 110  Format (1x,i3,1x,"Solution for an orifice of zero thickness")

 200  Format (" Iflow   = ",I2)
 201  Format (" visc    = ",f10.5)
 202  Format (" shrt    = ",f10.5)
 203  Format (" delta   = ",f10.5)
 204  Format (" NGL     = ",I3)
 205  Format (" radius  = ",f10.5)
 206  Format (" angle   = ",f10.5)
 207  Format (" truncw  = ",f10.5)
 208  Format (" truncw_u= ",f10.5)

 220  Format ("Segment",i2," NE = ",I3," RT =",f10.5)
 221  Format (" NEop = ",I3," RTop =",f10.5)

 900  Format (
     +        " sarea       =",f12.7,/,
     +        " Arc length  =",f12.7,/,
     +        " Drag        =",f12.7,/,
     +        " Torque_o    =",f12.7,/,
     +        " Torque_c    =",f12.7)

      Stop
      End
