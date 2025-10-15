      program nol1
 
c============================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c------------------------------------------------
c This program accompanies the book:
c               C. Pozrikidis
c Numerical Computation in Science and Engineering
c          Oxford University Press
c------------------------------------------------

c------------------------------------------------
c Computation of a zero of one nonlinear equation
c by several methods
c
c The objective function is defined in file: fun1
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)
    
      common/menu/a,b,OH,OOT,RT,P,T,al,eps_con,rl,rlp
      common/pii/pi
 
c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      Null = 0

c--------------------
c select the function
c--------------------

  98  Continue

      write (6,*)
      write (6,*) ' SOLVE A NONLINEAR ALGEBRAIC EQUATION'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) '  1 for spheroidal equation'
      write (6,*) '  2 for the column problem'
      write (6,*) '  3 for the CSTR problem'
      write (6,*) '  4 for the equation of state problem'
      write (6,*) '  5 for f(x) = exp(x)+cosx-2.52'
      write (6,*) '  6 for f(x) = 3*x**2+...'
      write (6,*) '  7 for f(x) = 9-2*x**2+5*x**3-x**4'
      write (6,*) '  8 for f(x) = lnx + a* x - b'
      write (6,*) '  9 for Stokes flow in a corner'
      write (6,*) ' 10 for f(x) = lnx + 3-3.1 x**2'
      write (6,*) ' 11 for f(x) = x*lnx'
      write (6,*) 
      write (6,*) '101 for antisymmetric Stokes flow in a corner'
      write (6,*) '102 for symmetric Stokes flow in a corner'
      write (6,*) 
      write (6,*) '110 for the Marshall & Trowbridge equation'
      write (6,*) 
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read (5,*) menu

      if(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) ' Choose the method'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for Newton'
      write (6,*) ' 2 for Muller'
      write (6,*) ' 3 for secant'
      write (6,*) ' 4 for Newton with third-order convergence'
      write (6,*) ' 0 to quit   '
      write (6,*) ' ---------'
      read  (5,*) method

c     method=1

      if(method.eq.0) Go to 99

c------
c input
c------

      if(menu.eq.1) then

        RT  = 2.0D0
        RT  = sqrt(RT)
        OOT = 1.0D0/3.0D0
        OH  = 0.50D0

        write (6,*)
        write (6,*) ' Enter the sphericity constant'
        write (6,*) ' -----------------------------'
        read  (5,*) A

      else if(menu.eq.4) then

        write (6,*)
        write (6,*) ' Enter P and T'
        write (6,*) ' -------------'
        read  (5,*) P,T

      else if(menu.eq.6) then

        write (6,*)
        write (6,*) ' Enter epsilon for continuation'
        write (6,*) ' ------------------------------'
        read  (5,*) eps_con

      else if(menu.eq.8) then

        write (6,*)
        write (6,*) ' Please enter a and b'
        write (6,*) ' --------------------'
        read (5,*) a,b

      else if(menu.eq.9) then

        al=0.25D0*pi

        write (6,*)
        write (6,*) ' Please enter lambda'
        write (6,*) ' --------------------'
        read (5,*) rl

      end if

c------------------
c a solution branch
c------------------

      if(menu.eq.101.or.menu.eq.102) then

        eps    = 0.000001  
        Niter  = 20 

        if(menu.eq.101) then
          N      = 128 
          x      = 1.5 
          alpha1 = 180.0
          alpha2 = 73.0
        elseif(menu.eq.102) then
          N      = 128
          x      = 1.5 
          alpha1 = 180.0
          alpha2 = 78.0
        end if

         N1      = N+1
         alpha1  = alpha1/180.0*pi
         alpha2  = alpha2/180.0*pi
         dal     = (alpha2-alpha1)/(N1-1.0)

         open (2,file="nol1.out")
         write (2,100) N1

         italk = 0

         Do i=1,N1
          al = alpha1 + (i-1.0)*dal
          if(method.eq.1) call newton1_2 (menu,Niter,x,eps,italk)
          if(method.eq.2) call muller    (menu,Niter,x,Iflag,italk)
          if(method.eq.3) call secant    (menu,Niter,x,Iflag,italk)
          if(method.eq.4) call newton1_3 (menu,Niter,x,eps,italk)
          aa = al/pi
          write (6,100) i,aa,x
          write (2,101)   aa,x
        End Do

        write (2,100) Null
        close (2)

c------------------
c a solution branch
c------------------

      elseif(menu.eq.110) then

        Niter = 10 
        eps   = 0.00001  
        N     = 2*128
        rlp1  = 0.00001
        rlp2  = 100.0

        N1    = N+1
        dlp   = (rlp2-rlp1)/(N1-1.0)

        x     = 0.001
        italk = 0

        open (2,file="nol1.out")
        write (2,100) N1

        rlp = rlp1

         Do i=1,N1
c         rlp = rlp1+(i-1.0)*dlp
          rlp=1.1*rlp
          if(method.eq.1) call newton1_2 (menu,Niter,x,eps,italk)
          if(method.eq.2) call muller    (menu,Niter,x,Iflag,italk)
          if(method.eq.3) call secant    (menu,Niter,x,Iflag,italk)
          if(method.eq.4) call newton1_3 (menu,Niter,x,eps,italk)
          aa = al/pi
          write (6,100) i,rlp,x
          write (2,101)   rlp,x
        End Do

        write (2,100) Null
        close (2)

c-------------
c one solution
c-------------

      else

      write (6,*)
      write (6,*) ' Enter the maximum number of iterations'
      write (6,*) ' --------------------------------------'
c     read  (5,*) Niter
      Niter = 10

c---
      if(method.eq.1.or.method.eq.4) then
        write (6,*) 
        write (6,*) ' Enter epsilon for computation of df/dx'
        write (6,*) ' --------------------------------------'
c       read  (5,*)  eps
        eps=0.0001
      end if

        write (6,*) 
        write (6,*) ' Enter the initial guess'
        write (6,*) ' -----------------------'

        read (5,*) x
c       x=0.3
        italk = 1;

        if(method.eq.1) call newton1_2 (menu,Niter,x,eps,italk)
        if(method.eq.2) call muller    (menu,Niter,x,Iflag,italk)
        if(method.eq.3) call secant    (menu,Niter,x,Iflag,italk)
        if(method.eq.4) call newton1_3 (menu,Niter,x,eps,italk)

c-----------
      end if
c-----------

      Go to 98  ! repeat

c-----
c finishing
c-----

 99   Continue

      write (6,*)
      write (6,*) " Thank you for your business"

c-----
c Done
c-----

 100  Format (2X,I6,2X,F20.10,2X,F15.10)
 101  Format (1X,F20.10,2X,F20.10)

      stop
      end
