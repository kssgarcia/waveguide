      program ode

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c            C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c       Oxford University Press
c------------------------------------------------

c------------------------------
c Solve a system of first-order 
c ODEs by several methods
c------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(10),f(10),xsave(10)
      Dimension xtmp1(10),xtmp2(10),xtmp3(10),xtmp4(10)
      Dimension ftmp1(10),ftmp2(10),ftmp3(10)

c--------------
c common blocks
c--------------

      common/menu1/a
      common/menu2/rk,r,b
      common/menu3/rnu
      common/menu4/rlam,pinf,epsilon

      common/pii/pih,pi

c----------
c constants
c----------

      null = 0

      pi = 3.14159 265358 D0
      pih = 0.5D0*pi

c------
c input
c------

 98   Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for the linear ode dx/dt = a x "
      write (6,*) " 2 for the Lorenz system"
      write (6,*) " 3 for the Renardy system"
      write (6,*) " 4 for the Reyleigh-Plesset equation"
      write (6,*) " 5 for problem 10.2.3"
      write (6,*) " 6 for trajectories of two spheres"
      write (6,*) "   in simple shear flow"
      write (6,*) " 0 to quit "
      write (6,*) " ----------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

c-----
c trap
c-----

      If(    menu.ne.1
     +  .and.menu.ne.2
     +  .and.menu.ne.3
     +  .and.menu.ne.4
     +  .and.menu.ne.5
     +  .and.menu.ne.6
     +   ) then

       write (6,*)
       write (6,*) " Sorry this selection is invalid;"
       write (6,*)
       write (6,*) " Please try again"
       write (6,*)

       Go to 98

      End If

c----------
c Inquiries
c----------

      write (6,*) 
      write (6,*) " Select the integration method"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) "  0 to quit "
      write (6,*) "  1 for Euler's method"
      write (6,*) "  2 for RK2"
      write (6,*) "  4 for RK4"
      write (6,*) " 10 for the mid-point method"
      write (6,*) " 34 for RKF34"
      write (6,*) " ------------"
      read  (5,*) method

      If(method.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter the time step"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Dt

      If(Dt.lt.0.00000000001) Go to 99

      write (6,*) 
      write (6,*) " How many steps before pausing?"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) npause

      If(npause.eq.0) Go to 99

c-------------------------
      If(method.eq.10) then   ! need to start up
c-------------------------

       write (6,*) 
       write (6,*) " Enter the number of subdivisions for start-up"
       write (6,*) " ---------------------------------------------"
       read  (5,*) nst

       Dts = Dt/nst
       Dt2 = 2.0D0*Dt

c-----------
      End If
c-----------

  96  Continue   ! target for repeat

c-----------------------
      If(menu.eq.1) then
c-----------------------

       n = 1    ! one equation

       write (6,*)
       write (6,*) " Please enter the value of a"
       write (6,*) " ---------------------------"
       read  (5,*) a

       write (6,*)
       write (6,*) " Please enter the initial value x(1)"
       write (6,*) " -----------------------------------"
       read  (5,*) x(1)

c-----------------------------
      Else If (menu.eq.2) then
c-----------------------------

       n = 3     ! three equations

       write (6,*)
       write (6,*) " Please enter the values of k, r, b"
       write (6,*) " ----------------------------------"
       read  (5,*) rk,r,b

       write (6,*)
       write (6,*) " Enter the initial values x(1), x(2), x(3)"
       write (6,*) " -----------------------------------------"
       read  (5,*) x(1),x(2),x(3)

c----------------------------
      Else If(menu.eq.3) then
c----------------------------

       n = 1    ! one equation

       write (6,*)
       write (6,*) " Please enter the value of nu"
       write (6,*) " ----------------------------"
       read  (5,*) rnu

       write (6,*)
       write (6,*) " Please enter the initial value x(1)"
       write (6,*) " -----------------------------------"
       read  (5,*) x(1)

c-----------------------------
      Else If (menu.eq.4) then
c-----------------------------

       n = 2     ! two equations

       write (6,*)
       write (6,*) " Please enter lambda, p_inf, epsilon"
       write (6,*) " -----------------------------------"
c      read  (5,*) rlam,pinf,epsilon

       rlam = 1.25D0
       pinf = 100.0D0
       epsilon = 0.10D0

       write (6,*)
       write (6,*) " Enter the initial values x(1), x(2)"
       write (6,*) " -----------------------------------"
c      read  (5,*) x(1),x(2)

       x(1) = 1.0D0
       x(2) = 0.0D0

c-----------------------------
      Else If (menu.eq.5) then
c-----------------------------

       n = 2        ! two equations

       x(1) = 1.0D0

       write (6,*)
       write (6,*) " Please enter x'(0)"
       write (6,*) " ------------------"
       read  (5,*) x(2)

c-----------------------------
      Else If (menu.eq.6) then
c-----------------------------

       n = 3        ! three equations

       write (6,*)
       write (6,*) " Please enter x,y,z at t=0"
       write (6,*) " -------------------------"
       read  (5,*) x(1),x(2),x(3)

       x(1)=2.0D0*x(1)
       x(2)=2.0D0*x(2)
       x(3)=2.0D0*x(3)

c-----------
      End If 
c-----------

c--------
c prepare
c--------

      open (1,file="PLOTDAT")

      time = 0.0D0

      kstep  = 1   ! global step counter
      ipause = 0   ! batch step counter

      If(menu.eq.6) then
       write (1,100) kstep,time,(0.5D0*x(i),i=1,n)
       write (6,100) kstep,time,(0.5D0*x(i),i=1,n)
      Else
       write (1,100) kstep,time,(x(i),i=1,n)
       write (6,100) kstep,time,(x(i),i=1,n)
      End If

c----------
c launching
c----------

      If(method.eq.34) call velocity (menu,x,time,f)

      Dth = 0.50D0*Dt

 97   Continue

c----------------------
c EULER EXPLICIT METHOD
c----------------------

      If(method.eq.1) then

      call velocity (menu,x,time,f)

      Do i=1,n
       x(i) = x(i)+Dt*f(i)
      End Do

      kstep = kstep+1
      time  = time + Dt

      If(menu.eq.6) then
       write (1,100) kstep,time,(0.5D0*x(i),i=1,n)
       write (6,100) kstep,time,(0.5D0*x(i),i=1,n)
      Else
       write (1,100) kstep,time,(x(i),i=1,n)
       write (6,100) kstep,time,(x(i),i=1,n)
      End If

c-----------
c RK2 METHOD
c-----------

      Else If(method.eq.2) then

      Do i=1,n
       xsave(i) = x(i)
      End Do

      call velocity
     +
     +    (menu
     +    ,x
     +    ,time
     +    ,f
     +    )

      Do i=1,n
       x(i) = xsave(i)+Dt*f(i)
      End Do

      call velocity
     +
     +    (menu
     +    ,x
     +    ,time+Dt
     +    ,ftmp1
     +    )

      Do i=1,n
       velav = 0.5D0*(f(i)+ftmp1(i))
       x(i)  = xsave(i)+Dt*velav
      End Do

      kstep = kstep+1

      time = time + Dt

      If(menu.eq.6) then
       write (1,100) kstep,time,(0.5D0*x(i),i=1,n)
       write (6,100) kstep,time,(0.5D0*x(i),i=1,n)
      Else
       write (1,100) kstep,time,(x(i),i=1,n)
       write (6,100) kstep,time,(x(i),i=1,n)
      End If

c-----------
c RK4 METHOD
c-----------

      Else If(method.eq.4) then

      Do i=1,n
       xsave(i) = x(i)
      End Do

      call velocity
     +
     +    (menu
     +    ,x
     +    ,time
     +    ,f
     +    )

      Do i=1,n
       x(i) = xsave(i)+Dth*f(i)
      End Do

      call velocity
     +
     +    (menu
     +    ,x
     +    ,time+Dth
     +    ,ftmp1
     +    )

      Do i=1,n
       x(i) = xsave(i)+Dth*ftmp1(i)
      End Do

      call velocity
     +
     +    (menu
     +    ,x
     +    ,time+Dth
     +    ,ftmp2
     +    )

      Do i=1,n
       x(i) = xsave(i)+Dt*ftmp2(i)
      End Do

      call velocity
     +
     +    (menu
     +    ,x
     +    ,time+Dt
     +    ,ftmp3
     +    )

      Do i=1,n
       velav = (f(i)+2.0D0*ftmp1(i)+2.0D0*ftmp2(i)+ftmp3(i))/6.0D0
       x(i)  = xsave(i)+Dt*velav
      End Do

      kstep = kstep+1

      time = time + Dt

      If(menu.eq.6) then
       write (1,100) kstep,time,(0.5D0*x(i),i=1,n)
       write (6,100) kstep,time,(0.5D0*x(i),i=1,n)
      Else
       write (1,100) kstep,time,(x(i),i=1,n)
       write (6,100) kstep,time,(x(i),i=1,n)
      End If

c------------------
c MID-POINT METHOD
c-----------------

      Else If(method.eq.10) then

c------------------------
      If(kstep.eq.1) then     ! start up
c------------------------

       Do i=1,n
        xsave(i) = x(i)
       End Do

       call velocity (menu,x,time,f)

       Do j = 1,nst
        Do i = 1,n
         x(i) = x(i)+Dts*f(i)
        End Do
        time = time + Dts
       End Do

c------------------------
      Else      ! further steps
c------------------------

       call velocity (menu,x,time,f)

       Do i = 1,n
        tmp  = x(i)
        x(i) = xsave(i)+Dt2*f(i)
        xsave(i) = tmp
       End Do

       time = time + Dt

c-----------
      End If
c-----------
 
      kstep  = kstep+1

      If(menu.eq.6) then
       write (1,100) kstep,time,(0.5D0*x(i),i=1,n)
       write (6,100) kstep,time,(0.5D0*x(i),i=1,n)
      Else
       write (1,100) kstep,time,(x(i),i=1,n)
       write (6,100) kstep,time,(x(i),i=1,n)
      End If

c-------------
c RKF34 METHOD
c-------------

      Else If(method.eq.34) then

      Do i = 1,n
       xtmp1(i) = x(i)+0.25D0*Dt*f(i)
      End Do

      tmp = time+0.25D0*Dt

      call velocity (menu,xtmp1,tmp,ftmp1)

      Do i = 1,n
       xtmp2(i) = x(i)+Dt*(-189.0D0*f(i)
     +                     +729.0D0*ftmp1(i))/800.0D0
      End Do

      tmp = time + 27.0D0/40.0D0 * Dt

      call velocity (menu,xtmp2,tmp,ftmp2)

      Do i = 1,n
       xtmp3(i) = x(i)+Dt*(214.0D0*f(i)
     +                     +27.0D0*ftmp1(i)
     +                    +650.0D0*ftmp2(i))/891.0D0
      End Do

      tmp = time + Dt

      call velocity (menu,xtmp3,tmp,ftmp3)

      Do i = 1,n
       xtmp4(i) = x(i)+Dt*(533.0D0*f(i)
     +                   +1600.0D0*ftmp2(i)
     +                    - 27.0D0*ftmp3(i))/2106.0D0
      End Do

      Do i = 1,n
       x(i) = xtmp3(i)
       f(i) = ftmp3(i)
      End Do

      kstep  = kstep+1
      time   = time + Dt

      If(menu.eq.6) then
       write (1,100) kstep,time,(0.5D0*x(i),i=1,n)
       write (6,100) kstep,time,(0.5D0*x(i),i=1,n)
      Else
       write (1,100) kstep,time,(x(i),i=1,n)
       write (6,100) kstep,time,(x(i),i=1,n)
      End If

c-----------
      End If
c-----------

      Ipause = Ipause + 1

      If(ipause.ne.npause) Go to 97

      ipause = 0

      write (6,*)
      write (6,*) " Continue the integration ?"
      write (6,*) " 0 for NO, 1 for YES"
      write (6,*) " -------------------"
      read  (5,*) Icon

      If(Icon.eq.1) Go to 97

      Go to 98
c     Go to 96

c-----
c Done
c-----

 99   Continue

      write (1,100) null
      close (1)

 100  Format (1x,i3,1x,f8.5,10(1x,f9.6))

      Stop
      End
