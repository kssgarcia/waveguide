      program sgf_2d_1p_w_dr

c========================================
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------
c Driver for the 2D periodic Green's function
c of Stokes flow in a semi-infinite domain bounded
c by a plane wall located at y = wall
c
c SYMBOLS:
c -------
c
c  RL: period
c
c  Iselect  = 1  compute only the Green's function
c                for the velocity
c          ne 1  compute the Green's function
c                for the velocity, pressure, stress
c------------------------------------------

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

      Ising   = 0   ! will not desingularize
      Iselect = 2   ! evaluate all Green's functions

c---

      write (6,*) 
      write (6,*) " Enter " 
      write (6,*) 
      write (6,*) " 1 for one evaluation"
      write (6,*) " 2 to test of integral identities on a circle"
      write (6,*) " 3 to test of integral identities on a line  "
      write (6,*) " 0 to quit "
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

      write (6,*) 
      write (6,*) " wall is located at y = wall: please enter: wall"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) wall
c     wall = -0.10D0
      if(wall.eq.99) Go to 99

      write (6,*) 
      write (6,*) " please enter the period"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) RL
c     RL = 1.3D0
      if(RL.eq.99) Go to 99

 98   Continue          ! target for repeat

      write (6,*) 
      write (6,*) " Please enter the coordinates"
      write (6,*) " of the point force: x0, y0"
      write (6,*) " 99 for either one to quit"
      write (6,*) " -------------------------"
      read  (5,*) x0,y0
c     x0 = 0.50D0
c     y0 = 0.56D0

      if(x0.eq.99) Go to 99
      if(y0.eq.99) Go to 99

c-----------------------
      if(menu.eq.1) then       ! one evaluation
c-----------------------

      write (6,*) 
      write (6,*) " Please enter: x, y"
      write (6,*) " 99 for either one to quit"
      write (6,*) " -------------------------"
      read  (5,*) x,y

      if(x.eq.99) Go to 99
      if(y.eq.99) Go to 99
  
      call sgf_2d_1p_w 
     +
     +  (Iselect
     +  ,X,Y
     +  ,X0,Y0
     +  ,wall
     +  ,RL
     +  ,Ising
     +  ,Gxx,Gxy
     +  ,Gyx,Gyy
     +  ,Px,Py
     +  ,Txxx,Txxy,Tyxx,Tyxy
     +  ,Txyx,Txyy,Tyyx,Tyyy
     +  )

      write (6,*)  " -----------------------"
      write (6,*)  " Velocity Green function (G):"
      write (6,*)
      write (6,100) Gxx,Gxy
      write (6,100) Gyx,Gyy
      write (6,*)
      write (6,*)  " Pressure (P):"
      write (6,*)
      write (6,100) Px,Py
      write (6,*)
      write (6,*)  " Stress Green function (T):"
      write (6,*)
      write (6,100) Txxx,Txxy,Tyxx,Tyxy
      write (6,100) Txyx,Txyy,Tyyx,Tyyy
      write (6,*)  " -----------------------"
      write (6,*)

c-----------
      else if(menu.eq.2) then        ! test on a circle
c-----------

      write (6,*)
      write (6,*) " Enter "
      write (6,*)
      write (6,*) " (a) the circle center coordinates (xcnt,ycnt)"
      write (6,*) " (b) the circle radius "
      write (6,*) " (c) the number of integration points "
      write (6,*) " ---------------------------------"
      read (5,*) xcnt,ycnt,rad,mint

      dth = pi2/mint

      sgfx = 0.0D0   ! test on continuity equation
      sgfy = 0.0D0

      sm11 = 0.0D0   ! test of force balance
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)
        x = xcnt + rad*cs
        y = ycnt + rad*sn

        vx = cs              ! normal vector outward
        vy = sn

        call sgf_2d_1p_w 
     +
     +   (Iselect
     +   ,X,Y
     +   ,X0,Y0
     +   ,wall
     +   ,RL
     +   ,Ising
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,Px,Py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

        sgfx = sgfx + vx*Gxx + vy*Gyx
        sgfy = sgfy + vx*Gxy + vy*Gyy 

        sm11 = sm11 + txxx*vx + txxy*vy
        sm12 = sm12 + tyxx*vx + tyxy*vy
        sm21 = sm21 + txyx*vx + txyy*vy
        sm22 = sm22 + tyyx*vx + tyyy*vy

      end do

      cf = dth*rad/pi4

      sm11 = sm11 * cf
      sm12 = sm12 * cf
      sm21 = sm21 * cf
      sm22 = sm22 * cf

      write (6,*)
      write (6,*)  " Should be 0:"
      write (6,*)
      write (6,100) sgfx,sgfy
      write (6,*)
      write (6,*)  " Should be 0 or -1:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22

c---
c second test
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

        call sgf_2d_1p_w 
     +
     +    (Iselect
     +    ,X0,Y0
     +    ,X,Y
     +    ,wall
     +    ,RL
     +    ,Ising
     +    ,Gxx,Gxy
     +    ,Gyx,Gyy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

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
      write (6,100)

c----------------------------
      else if(menu.eq.3) then     ! test of identities on a line
c----------------------------

      write (6,*)
      write (6,*) " rest on a sinusoidal line "
      write (6,*)
      write (6,*) " Enter the y coordinate of the line"
      write (6,*) "       the amplitude of the sinusoid"
      write (6,*) " -----------------------------------"
      read  (5,*) yline,amp

      write (6,*)
      write (6,*) " Enter the number of integration points"
      write (6,*) " --------------------------------------"
      read  (5,*) mint

      wn = pi2/RL    ! wave number

      sgfx = 0.0D0
      sgfy = 0.0D0

      sm11 = 0.0D0
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      dx = RL/mint

      Do i=1,mint

        x   = (i-1.0D0)*dx
        arg = wn*x
        y   = yline+amp*Dsin(arg)

        vx =  amp*wn*Dcos(arg)   ! normal vector toward the wall
        vy = -1.0D0

        call sgf_2d_1p_w 
     +
     +      (Iselect
     +      ,X,Y
     +      ,X0,Y0
     +      ,wall
     +      ,RL
     +      ,Ising
     +      ,Gxx,Gxy
     +      ,Gyx,Gyy
     +      ,Px,Py
     +      ,Txxx,Txxy,Tyxx,Tyxy
     +      ,Txyx,Txyy,Tyyx,Tyyy
     +      )

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
      write (6,*)  " should be 0 or -1:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22

c---
c second test
c---

      sm11 = 0.0D0
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      Do i=1,mint

        x   = (i-1.0D0)*dx
        arg = wn*x
        y   = yline+amp*Dsin(arg)

        vx =  amp*wn*Dcos(arg)    ! normal vector toward the wall
        vy = -1.0D0

        call sgf_2d_1p_w 
     +
     +   (Iselect
     +   ,X0,Y0
     +   ,X,Y
     +   ,wall
     +   ,RL
     +   ,Ising
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,Px,Py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

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
      write (6,100)

c--------
      end if
c--------

      Go to 98

c-----
c done
c-----

  99  Continue

 100  Format (5(1x,f15.10))
 101  Format (1x,i3,4(1x,f15.10))

      stop
      end
