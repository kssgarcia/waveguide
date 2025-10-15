      program lgf_2d_1p_dr

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=========================================

c---------------------------------------
c Driver for the singly periodic
c Green's function of Laplace's equation
c in two dimensions
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      Iopt = 2  ! need G and gradient

c-------
c launch
c-------

 98   Continue

      write (6,*)
      write (6,*) " Enter: "
      write (6,*)
      write (6,*) " 1 for one evaluation"
      write (6,*) " 2 to test integral identities on a circle"
      write (6,*) " 3 to test integral identities on a periodic line"
      write (6,*) " 0 to quit "
      write (6,*) " ----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Please enter the period"
      write (6,*) " -----------------------"
      read  (5,*) RL

      write (6,*) 
      write (6,*) " Please enter the coordinates"
      write (6,*) " of one singularity: x0, y0"
      write (6,*) 
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x0,y0

      if(x0.eq.99) Go to 99
      if(y0.eq.99) Go to 99

c---------------
c One evaluation
c---------------

      if(menu.eq.1) then

      write (6,*) 
      write (6,*) " Enter the coordinates of the field point"
      write (6,*) " ----------------------------------------"
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99
  
      call lgf_2d_1p
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,RL
     +   ,G
     +   ,Gx,Gy
     +   )

c------------------------------
c uncomment to subtract off the
c free-space Green's function
c if desired
c---
c
c     Drs = (X-X0)**2+(Y-Y0)**2
c     G   = G + log(Drs)/pi4
c
c------------------------------

      write (6,*)  " ---------------------------"
      write (6,*)  " Green function and gradient"
      write (6,100) G
      write (6,100) Gx,Gy
      write (6,*)  " ---------------------------"

c----------------------------------------
c test of integral identities on a circle
c----------------------------------------

      elseif(menu.eq.2) then

      write (6,*) 
      write (6,*) " Enter the center of the test circle (xc, yc)"
      write (6,*) " --------------------------------------------"
      read (5,*) xcnt,ycnt

      write (6,*) 
      write (6,*) " Enter the radius of the test circle"
      write (6,*) " -----------------------------------"
      read (5,*) rad

      write (6,*) 
      write (6,*) " Enter the number of integration points"
      write (6,*) " --------------------------------------"
      read (5,*) mint

      dth = pi2/mint

      sm1 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)

        x = xcnt + rad*cs
        y = ycnt + rad*sn

        vnx = cs    ! normal vector pointing outward
        vny = sn

        call lgf_2d_1p
     +
     +     (Iopt
     +     ,X,Y
     +     ,X0,Y0
     +     ,RL
     +     ,G
     +     ,Gx,Gy
     +     )

        sm1 = sm1 + Gx*vnx + Gy*vny

      End Do

      cf = dth*rad

      sm1 = sm1 * cf

      write (6,*)
      write (6,*) " Should be 0 or -1:"
      write (6,*)
      write (6,100) sm1
      write (6,*)

c-----------------------------------------------
c Test of integral identities on a periodic line
c-----------------------------------------------

      elseif(menu.eq.3) then

       write (6,*)
       write (6,*) " The line will be a sinusoidal wave"
       write (6,*)
       write (6,*) " Enter the mean y location of the line"
       write (6,*) " -------------------------------------"
       read  (5,*) yline

       write (6,*)
       write (6,*) " Enter the wave amplitude"
       write (6,*) " ------------------------"
       read  (5,*) amp

       write (6,*)
       write (6,*) " Enter the number of integration points"
       write (6,*) " --------------------------------------"
       read  (5,*) mint

       wn = pi2/RL   ! wave number

       dx = RL/mint

       sm1 = 0.0D0

       Do i=1,mint

        x   = (i-1.0D0)*dx
        arg = wn*x
        y   = yline+amp*Dsin(arg)

        vnxdldx = -amp*wn*cos(arg)   ! normal vector pointing upward
        vnydldx =  1.0D0

        call lgf_2d_1p
     +
     +    (Iopt
     +    ,X,Y
     +    ,X0,Y0
     +    ,RL
     +    ,G
     +    ,Gx,Gy
     +    )

        sm1 = sm1 + Gx*vnxdldx + Gy*vnydldx

      End Do

      sm1 = sm1 * dx

      write (6,*)
      write (6,*) " Should be 0.5 or -0.5:"
      write (6,*)
      write (6,100) sm1
      write (6,*)

c-----------
      End If
c----------

      Go to 98

c-----
c Done
c-----

  99  Continue

 100  Format (5(1x,f15.10))
 101  Format (1x,i3,4(1x,f15.10))

      Stop
      End
