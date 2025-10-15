      program sgf_2d_1p_dr

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c Driver for the singly-periodic Green's function
c of two-dimensional Stokes flow: sgf_2d_1p
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c------
c flags
c------

      Idesing  = 0   ! will not subtract off the Stokeslet
      Iopt     = 2   ! will compute G, P, and T

c-------
c launch
c-------

 98   Continue

      write (6,*)
      write (6,*) " Please enter: "
      write (6,*)
      write (6,*) " 1 for one evaluation"
      write (6,*) " 2 to test integral identities on a circle"
      write (6,*) " 3 to test integral identities on a periodic line"
      write (6,*) " 0 to quit "
      write (6,*) " ----------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter the period"
      write (6,*) " ----------------"
      read  (5,*) period

      write (6,*) 
      write (6,*) " Enter the coordinates of a point force: x0, y0"
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x0,y0

      if(x0.eq.99) Go to 99
      if(y0.eq.99) Go to 99

c---------------
c one evaluation
c---------------

      if(menu.eq.1) then

      write (6,*) 
      write (6,*) " Please enter:"
      write (6,*) " the coordinates of the velocity point: x,y"
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99
  
      call sgf_2d_1p
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,period
     +   ,Idesing
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,Px,Py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

c---------------------------------

      write (6,*)  "-----------------------------"
      write (6,*)  " Velocity Green function: Gij"
      write (6,*)
      write (6,100) Gxx,Gxy
      write (6,100) Gyx,Gyy
      write (6,*)
      write (6,*)  " Pressure Green function: Pj"
      write (6,*)
      write (6,100) Px,Py
      write (6,*)
      write (6,*)  " Stress Green function: Tijk"
      write (6,*)
      write (6,100) Txxx,Txxy,Tyxx,Tyxy
      write (6,100) Txyx,Txyy,Tyyx,Tyyy
      write (6,*)
      write (6,*)  "----------------------------"
      write (6,*)

c----------------------------------------
c Test of integral identities on a circle
c see Pozrikidis (1992)
c----------------------------------------

      Else If(menu.eq.2) then

      write (6,*) 
      write (6,*) " Enter: "
      write (6,*) 
      write (6,*) " the circle center coordinates (xcnt, ycnt)"
      write (6,*) " the radius of the test circle"
      write (6,*) " the number of integration points "
      write (6,*) " --------------------------------"
      read (5,*) xcnt,ycnt,rad,mint

      dth = pi2/mint

      sgfx = 0.0D0   ! test on the continuity equation
      sgfy = 0.0D0

      sm11 = 0.0D0   ! test on force balance
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)
        x  = xcnt + rad*cs
        y  = ycnt + rad*sn

        vnx = cs      ! unit normal vector
        vny = sn

        call sgf_2d_1p
     +
     +    (Iopt
     +    ,X,Y
     +    ,X0,Y0
     +    ,period
     +    ,Idesing
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

        sgfx = sgfx + vnx*Gxx + vny*Gyx
        sgfy = sgfy + vnx*Gxy + vny*Gyy 

        sm11 = sm11 + Txxx*vnx + txxy*vny
        sm12 = sm12 + Tyxx*vnx + tyxy*vny
        sm21 = sm21 + Txyx*vnx + txyy*vny
        sm22 = sm22 + Tyyx*vnx + tyyy*vny

      End Do

      cf = dth*rad/pi4

      sm11 = sm11 * cf
      sm12 = sm12 * cf
      sm21 = sm21 * cf
      sm22 = sm22 * cf

      write (6,*)
      write (6,*)  " should be 0:"
      write (6,*)
      write (6,100) sgfx,sgfy
      write (6,*)
      write (6,*)  " should be 0 or -1:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22

c---
c Second test
c---

      sm11 = 0.0D0
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)
        x = xcnt + rad*cs
        y = ycnt + rad*sn

        vx = cs
        vy = sn

        call sgf_2d_1p
     +
     +     (Iopt
     +     ,X0,Y0
     +     ,X,Y
     +     ,period
     +     ,Idesing
     +     ,Gxx,Gxy
     +     ,Gyx,Gyy
     +     ,px,py
     +     ,Txxx,Txxy,Tyxx,Tyxy
     +     ,Txyx,Txyy,Tyyx,Tyyy
     +     )

        sm11 = sm11 + txxx*vx + txyx*vy
        sm12 = sm12 + tyxx*vx + tyyx*vy
        sm21 = sm21 + txxy*vx + txyy*vy
        sm22 = sm22 + tyxy*vx + tyyy*vy

      End Do

      sm11 = sm11 * cf
      sm12 = sm12 * cf
      sm21 = sm21 * cf
      sm22 = sm22 * cf


      write (6,*)
      write (6,*)  " should be 0 or 1:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22

c--------------------------------------------
c Test integral identities on a periodic line
c--------------------------------------------

      Else If(menu.eq.3) then

      write (6,*)
      write (6,*) " Test on a sinusoidal line "
      write (6,*)
      write (6,*) " Enter the y coordinate of the line"
      write (6,*) "       the amplitude of the sinusoid"
      write (6,*) " -----------------------------------"
      read  (5,*) yline,amp

      write (6,*)
      write (6,*) " Enter number of integration points "
      write (6,*) " -----------------------------------"
      read (5,*) mint

      wn = pi2/period    ! wave number

      sgfx = 0.0D0
      sgfy = 0.0D0

      sm11 = 0.0D0
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      dx = period/mint

      Do i=1,mint

        x = (i-1.0D0)*dx
        arg = wn*x
        y = yline+amp*sin(arg)

        vx =  amp*wn*cos(arg)   ! normal vector downward
        vy = -1.0

        call sgf_2d_1p
     +
     +    (Iopt
     +    ,X,Y
     +    ,X0,Y0
     +    ,period
     +    ,Idesing
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

        sgfx = sgfx + vx*Gxx + vy*Gyx
        sgfy = sgfy + vx*Gxy + vy*Gyy

        sm11 = sm11 + txxx*vx + txxy*vy
        sm12 = sm12 + tyxx*vx + tyxy*vy
        sm21 = sm21 + txyx*vx + txyy*vy
        sm22 = sm22 + tyyx*vx + tyyy*vy

      End Do

      cf  = dx/pi4

      sm11 = sm11 * cf
      sm12 = sm12 * cf
      sm21 = sm21 * cf
      sm22 = sm22 * cf

      write (6,*)
      write (6,*)  " should be 0:"
      write (6,*)
      write (6,100) sgfx,sgfy
      write (6,*)
      write (6,*)  " should be 0 or 0.5:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22

c---
c Second test
c---

      sm11 = 0.0D0
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      Do i=1,mint

        x = (i-1.0D0)*dx
        arg = wn*x
        y  = yline+amp*sin(arg)

        vx =  amp*wn*cos(arg)    ! normal vector downward
        vy = -1.0D0

        call sgf_2d_1p
     +
     +     (Iopt
     +     ,X0,Y0
     +     ,X,Y
     +     ,period
     +     ,Idesing
     +     ,Gxx,Gxy
     +     ,Gyx,Gyy
     +     ,px,py
     +     ,Txxx,Txxy,Tyxx,Tyxy
     +     ,Txyx,Txyy,Tyyx,Tyyy
     +     )

        sm11 = sm11 + txxx*vx + txyx*vy
        sm12 = sm12 + tyxx*vx + tyyx*vy
        sm21 = sm21 + txxy*vx + txyy*vy
        sm22 = sm22 + tyxy*vx + tyyy*vy

      End Do

      sm11 = sm11 * cf
      sm12 = sm12 * cf
      sm21 = sm21 * cf
      sm22 = sm22 * cf

      write (6,*)
      write (6,*)  " should be 0 or -0.5"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22
      write (6,100)

c--------
      End If
c--------

      Go to 98

  99  Continue

c-----
c Done
c-----

 100  Format (5(1x,f15.10))
 101  Format (1x,i3,4(1x,f15.10))

      Stop
      End
