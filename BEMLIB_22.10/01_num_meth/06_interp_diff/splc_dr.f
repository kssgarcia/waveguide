      program splc_dr

c=====================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=====================================

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c      Oxford University Press
c------------------------------------------------

c--------------------------------------
c Driver for cubic spline interpolation
c
c LEGEND:
c ------
c
c N: 	    Number of intervals
c xp(i):    X value of ith point
c fp(i):    f value of ith point
c a, b, c:  cubic spline coefficients
c
c kpl:  number of intervals for plotting the spline
c
c Capacity:
c ---------
c
c N = 512 intervals
c--------------------------------------

      Implicit Double Precision (a-h, o-z)

      Dimension xp(513),fp(513)
      Dimension  a(512), b(512),c(512)

      Dimension fpx(513),fpy(513)
      Dimension  ax(512), bx(512),cx(512)
      Dimension  ay(512), by(512),cy(512)

      Parameter (eps=0.00000000001)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

      Null = 0

c-----
c menu
c-----

 98   Continue

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to read data from file: splc.dat"
      write (6,*) " 2 for interpolating f(x) = 1/(1+25*x^2)"
      write (6,*) " 3 to read periodic data from file: splc_pr.dat"
      write (6,*) " 8 to reconstruct a closed line"
      write (6,*) "   based on data in file: splc_lc.dat"
      write (6,*) " 9 to reconstruct an open line"
      write (6,*) "   based on data in file: splc_lo.dat"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read (5,*) menu

      if(menu.eq.0) Go to 99

c-----------------------
      if(menu.eq.3.or.menu.eq.8) then
c-----------------------

       Itype = 9     ! periodic spline

c----------
       else
c----------

       write (6,*)
       write (6,*) " CHOOSE THE SPLINE TYPE"
       write (6,*)
       write (6,*) " Enter:"
       write (6,*)
       write (6,*) " 1 for clamped-end conditions "
       write (6,*) " 2 for a natural spline"
       write (6,*) " 0 to quit"
       write (6,*) "----------"
       read (5,*) Itype

       if(Itype.eq.0)  Go to 99

c----------
      end if
c----------

c-----------------------------
c clamped_end: specify the end-point slopes
c-----------------------------

      if(Itype.eq.1) then
        write (6,*)
        write (6,*) " Enter the slope at the first point"
        write (6,*) " ----------------------------------"
        read (5,*) slope1
        write (6,*)
        write (6,*) " Enter the slope at the last point"
        write (6,*) " ---------------------------------"
        read (5,*) slope2
      end if

c--------
c prepare
c--------

      open (7,file="splc.out")

c--------------------------------
c read or compute the data points
c--------------------------------

c-----------------------
      if(menu.eq.1) then
c-----------------------

        open (1,file="splc.dat")

          read (1,*) N1

          Do i=1,N1
           read (1,*) idle,xp(i),fp(i)
          End Do

        close (1)

        N = N1-1

        write (6,*)
        write (7,101) N1

        Do i=1,N1
          write(6,101) i,xp(i),fp(i)
          write(7,101) i,xp(i),fp(i)
        End Do

        aint = xp(1) +eps     ! limits of interpolation
        bint = xp(N1)-eps     ! limits of interpolation

c----------------------------
      else if(menu.eq.2) then
c----------------------------

        write (6,*)
        write (6,*) " Enter the end-points a, b"
        write (6,*) " -------------------------"
        read  (5,*) aint,bint

        write (6,*)
        write (6,*) " Enter the number of segments N"
        write (6,*) " ------------------------------"
        read  (5,*) N

        N1 = N+1
        h  = (bint-aint)/(N1-1.0D0)

        Do i=1,N1
         xp(i) = aint +(i-1.0D0)*h
         fp(i) = 1.0D0/(1.0D0+25.0D0*xp(i)**2)
        End Do

        write (7,101) N1
        write (6,*)

        Do i=1,N1
          write(7,101) i,xp(i),fp(i)
          write(6,101) i,xp(i),fp(i)
        End Do

c-----------------------
      elseif(menu.eq.3) then
c-----------------------

        write (6,*) "splc_dr: reading data from file splc_pr.dat"

        open (1,file="splc_pr.dat")

        read (1,*) N1

        Do i=1,N1
           read (1,*) idle,xp(i),fp(i)
        End Do

        close (1)

        N = N1-1

        write (6,*)
        write (6,*) "data:"
        write (6,*)
        write (7,101) N1

        Do i=1,N1
          write(6,101) i,xp(i),fp(i)
          write(7,101) i,xp(i),fp(i)
        End Do

        aint = xp(1) +eps     ! limits of interpolation
        bint = xp(N1)-eps     ! limits of interpolation

c----------------------------
      else if(menu.eq.8) then
c----------------------------

          write (6,*) " splc_dr: reading data from file: splc_lc.dat"

          open (1,file="splc_lc.dat")

          read (1,*) N

          Do i=1,N
           read (1,*) idle,fpx(i),fpy(i)
          End Do

          N1 = N+1

          fpx(N1) = fpx(1)     ! last point is the first point
          fpy(N1) = fpy(1)     ! last point is the first point

          close (1)

          write (6,*) " splc_dr: done reading"

c-----------------------------------------
      else if(menu.eq.9) then
c-----------------------------------------

          write (6,*) " splc_dr: reading data from file: splc_lo/dat"

          open (1,file="splc_lo.dat")

          read (1,*) N1

          Do i=1,N1
           read (1,*) idle,fpx(i),fpy(i)
          End Do

          N = N1-1

          close (1)

          write (6,*) " splc_dr: done reading"

c------------
      end if
c------------

c------------------------------------
c open or closed line:
c compute the polygonal arc length
c to be used as the interpolation variable
c------------------------------------

      if(menu.eq.8.or.menu.eq.9) then

        xp(1) = 0.0D0

        Do i=2,N1
         ia = i-1
         xp(i) = xp(ia)+ Dsqrt((fpx(i)-fpx(ia))**2
     +                        +(fpy(i)-fpy(ia))**2)
        End Do

        write (7,101) N1

        write (6,*)
        write (6,*) " splc_dr: x, fx, fy:"
        write (6,*)

        Do i=1,N1
          write(6,101) i,xp(i),fpx(i),fpy(i)
          write(7,101) i,xp(i),fpx(i),fpy(i)
        End Do

        aint = xp(1) +eps     ! limits of interpolation
        bint = xp(N1)-eps     ! limits of interpolation

c-----------
      end if
c-----------

c--------------------
c compute the splines
c--------------------

c===================
      if(menu.eq.1.or.menu.eq.2) then
c===================

c---
      if(Itype.eq.1) then   ! clamped
c---

      call splc_clm
     +
     +  (N
     +  ,xp,fp
     +  ,slope1,slope2
     +  ,a,b,c
     +  )

c---
      else if(Itype.eq.2) then   ! natural
c---

      call splc_nt
     +
     +  (N
     +  ,xp,fp
     +  ,a,b,c
     +  )

c---
       end if
c---

      write (6,*)
      write (6,*) " spline coefficients: a, b, c"
      write (6,*) " ----------------------------"
      write (6,*)

      Do i=1,N
        write (6,100) i,a(i),b(i),c(i)
      End Do

c===================
      end if
c===================

c===================
      if(menu.eq.3) then
c===================

      call splc_pr
     +
     +  (N
     +  ,xp,fp
     +  ,a,b,c
     +  )

      write (6,*)
      write (6,*) " spline coefficients: a, b, c"
      write (6,*) " ----------------------------"
      write (6,*)

      Do i=1,N
        write (6,100) i,a(i),b(i),c(i)
      End Do

c===================
      end if
c===================

c===================
      if(menu.eq.8.or.menu.eq.9) then     ! line
c===================

c---------
      if(Itype.eq.1) then   ! clamped-end
c---------

       call splc_clm
     +
     + (N
     + ,xp,fpx
     + ,slope1,slope2
     + ,ax,bx,cx
     + )

       call splc_clm
     +
     + (N
     + ,xp,fpy
     + ,slope1,slope2
     + ,ay,by,cy
     + )

c---------
      else if(Itype.eq.2) then   ! natural
c---------

       call splc_nt
     +
     + (N
     + ,xp,fpx
     + ,ax,bx,cx
     + )

       call splc_nt
     +
     + (N
     + ,xp,fpy
     + ,ay,by,cy
     + )

c---------
      else if(Itype.eq.9) then   ! periodic
c---------

       call splc_pr
     +
     + (N
     + ,xp,fpx
     + ,ax,bx,cx
     + )

       call splc_pr
     +
     + (N
     + ,xp,fpy
     + ,ay,by,cy
     + )

c-------
      end if
c-------

      write (6,*)
      write (6,*) " spline coefficients: a, b, c for x(s):"
      write (6,*) " --------------------------------------"
      write (6,*)

      Do i=1,N
        write (6,100) i,ax(i),bx(i),cx(i)
      End Do

      write (6,*)
      write (6,*) " spline coefficients: a, b, c for y(s):"
      write (6,*) " --------------------------------------"
      write (6,*)

      Do i=1,N
        write (6,100) i,ay(i),by(i),cy(i)
      End Do

c=============
      end if       ! close loop over menu
c=============

c--------------------------------------------------------------

c----------------
c post processing
c----------------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to interpolate at a point"
      write (6,*) " 2 to plot the spline"
      write (6,*) " 9 to return to the main menu"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Iopt

      if(Iopt.eq.0) Go to 99
      if(Iopt.eq.1) Go to 98

c-----------------------
      if(Iopt.eq.1) then    ! one point evaluation 
c-----------------------

  97  Continue

      write (6,*)
      write (6,*) " Please enter x_int; 99 to quit"
      write (6,*) "-------------------------------"
      read  (5,*) xint

      if(xint.eq.99) Go to 98

c---
c locate the host interpolation interval
c---

      Do i=1,N
        prod = (xint-xp(i)) * (xint-xp(i+1))
        if(prod.le.0.00) then
         l = i
         Go to 80
        end if
      End Do
 80   Continue

      xd = xint-xp(l)

c---
c evaluate the interpolating polynomial
c---

      if(menu.eq.1.or.menu.eq.3) then

         yint = ( (a(l)*xd + b(l) )*xd + c(l) )*xd + fp(l)

         write (6,102) xint,yint

      else if(menu.eq.2) then

         yint  = ( (a(l)*xd + b(l) )*xd + c(l) )*xd + fp(l)
         exact = 1.0D0/(1.0D0+25.0D0*xint**2)

         write (6,102) xint,yint,exact

      else if(menu.eq.8.or.menu.eq.9) then

         xx = ( (ax(l)*xd+ bx(l) )*xd+ cx(l) )*xd + fpx(l)
         yy = ( (ay(l)*xd+ by(l) )*xd+ cy(l) )*xd + fpy(l)
         write (6,102) xint,xx,yy

      end if

      Go to 97

c----------------------------
      else if(Iopt.eq.2) then  ! generate data for plotting
c----------------------------

      write (6,*)
      write (6,*) " Enter the number of plotting points"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read (5,*) kpl

      if(kpl.eq.0) Go to 99

      kpl1 = kpl+1

      write (7,101) kpl1

c---
c plotting interval etc
c---

      xint = aint
      dxint = (bint-aint)/(kpl-1.0D0+1.0D0)

c---
c loop over plotting points
c---

      Do j=1,kpl1

c---
c find the host interval
c---

        Do i=1,N
          prod = (xint-xp(i)) * (xint-xp(i+1))
          if(prod.le.0.00) then
            l = i
            Go to 81
          end if
        End Do

 81     Continue

        xd = xint-xp(l)

c---
c compute the spline
c---

      if(menu.eq.1.or.menu.eq.3) then 

        yint = ( (a(l)*xd+ b(l) )*xd+ c(l) )*xd + fp(l)

        write (7,101) j,xint,yint
        write (6,101) j,xint,yint

      else if(menu.eq.2) then ! data

        yint = ( (a(l)*xd+ b(l) )*xd+ c(l) )*xd + fp(l)

        exact = 1.0D0/(1.0D0+25.0D0*xint**2)

        write (7,101) j,xint,yint,exact
        write (6,101) j,xint,yint,exact

      else if(menu.eq.8.or.menu.eq.9) then   ! a line

        xx  = ( (ax(l)*xd+ bx(l) )*xd+ cx(l) )*xd + fpx(l)
        yy  = ( (ay(l)*xd+ by(l) )*xd+ cy(l) )*xd + fpy(l)

        write (7,101) j,xint,xx,yy
        write (6,101) j,xint,xx,yy

      end if

       xint = xint+dxint

      End Do

c-----------
      end if
c-----------

      Go to 98  ! repeat

c-----
c done
c-----

 99   continue

      write (7,101) Null
      close (7)

 100  Format (1x,i3,3(1x,f10.5),2x,3(1x,f10.5))
 101  Format (1x,I3,10(1x,f10.5))
 102  Format (1x,10(1x,f10.5))

      stop
      end
