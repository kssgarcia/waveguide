      program nol

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------------
c This program accompanies the book:
c
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c------------------------------------------------
c solve a system of nonlinear algebraic equations
c using Newton's or Broyden's method
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(20),Jac(20,20),Bra(20,20)
    
      common/pii/pi
      common/menu12/al
      common/menu5/a,b,c
 
c----------
c constants
c----------

      pi = 3.14159 265358 D0

      Null = 0

c--------
c prepare
c--------

      open (2,file="nol.out")

c--------
c prepare
c--------

      Ispeak = 0   ! silent mode
      Ispeak = 1   ! verbose mode

c-------
c launch
c-------

 98   Continue

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 1 for anti-symm Stokes flow in a corner '
      write (6,*) ' 2 for      symm Stokes flow in a corner '
      write (6,*) ' 3 for      a test case with n=2'
      write (6,*) ' 4 for      a test case with n=3        '
      write (6,*) ' 5 for      three touching circles'
      write (6,*) ' 6 for      the tubular reactor problem'
      write (6,*) '            equations (4.1.21)         '
      write (6,*) ' 7 for      equations (4.5.10) (n=2)    '
      write (6,*) ' 8 for      some equations (n=2)    '
      write (6,*) ' 0 to quit'
      write (6,*) '----------' 
      read  (5,*) menu

      If(menu.eq.0) Go to 99

c------------------
c choose the method
c------------------

      write (6,*)
      write (6,*) ' Choose the method'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for Newton iterations'
      write (6,*) ' 2 for Broyden iterations' 
      write (6,*) ' 0 to quit'
      write (6,*) '----------' 
      read (5,*) method

      if(method.eq.0) Go to 99

c---------------------
c How many equations ?
c---------------------

      if(menu.eq.1) n = 2
      if(menu.eq.2) n = 2
      if(menu.eq.3) n = 2
      if(menu.eq.4) n = 3
      if(menu.eq.5) n = 2
      if(menu.eq.6) n = 2
      if(menu.eq.7) n = 2
      if(menu.eq.8) n = 2

      write (6,*) 
      write (6,*) ' nonlinear: system size: ',n
      write (6,*) 

c------------------------
c definition of equations
c and commentary
c------------------------

      if(menu.eq.1) then

        eps   = 0.00001 !epsilon for Newton iterations
        Niter = 100     !maximum number of iterations
        x(1)  = 3.3     !initial guess for xi and eta (antisymmetric)
        x(2)  = 1.75 
        m     = 128    !number of intervals
        ai    = 73.0   ! alpha(initial) in degrees
        af    = 0.01    ! alpha(final) in degrees

      elseif(menu.eq.2) then

        eps   = 0.00001 !epsilon for Newton iterations
        Niter = 100     !maximum number of iterations
        x(1)  = 7.5    !initial guess for xi and eta (antisymmetric)
        x(2)  = 0.5
        m     = 128    !number of intervals
        ai    = 78.0   ! alpha(initial) in degrees
        af    = 0.01    ! alpha(final) in degrees

      elseif(menu.eq.3) then

        write (6,*)
        write (6,*) ' A solution is: 1, 2'
        write (6,*) ' -------------------'

      elseif(menu.eq.4) then

        write (6,*)
        write (6,*) ' A solution is : 0.4981, -0.1996, -0.5288'
        write (6,*) ' ----------------------------------------'

      elseif(menu.eq.5) then

        write (6,*)
        write (6,*) ' Please enter a, b'
        write (6,*) ' -----------------'
        read  (5,*) a,b
        c = 0.10*a

      elseif(menu.eq.7) then

        write (6,*)
        write (6,*) ' A solution is: 2, 2'
        write (6,*) ' -------------------'

      elseif(menu.eq.8) then

        write (6,*)
        write (6,*) ' Choose initial guess: 3.0, 3.0'
        write (6,*) ' ------------------------------'

      End If

      if(menu.eq.1.or.menu.eq.2) Go to 996

c--------------
c initial guess
c--------------

      write (6,*)
      write (6,*) ' Enter the components '
      write (6,*) ' of the starting vector:'
      write (6,*) ' -----------------------'

      read  (5,*) (x(i),i=1,n)

      write (6,*)
      write (6,*) ' Please enter maximum number of iterations'
      write (6,*) ' -----------------------------------------'
      read (5,*) Niter

c-----------------
c Newton's method:
c----------------

      If(method.eq.1) then

        write (6,*)
        write (6,*) ' Newton iterations:'
        write (6,*) ' Jacobian will be computed by finite differences'
        write (6,*) ' Please enter epsilon '
        write (6,*) ' ---------------------'
        read (5,*) eps

      End If

c--------------------------------
c Broyden method:
c
c Initial approximation to the inverse
c of the Jacobian
c
c Here we choose the initial guess
c to be the identity matrix
c--------------------------------

      if(method.eq.2) then

       Do i=1,n
         Do j=1,n
          Bra(i,j) = 0.0D0
         End Do
         Bra(i,i) = 1.0D0
       End Do

      end if

c---------------------------------
c Stokes flow in a corner
c
c compute a whole solution branch
c--------------------------------

  996 Continue

      if(menu.eq.1.or.menu.eq.2)then

        m1 = m+1
        write (2,100) m1

        ai = ai/180.0*pi
        af = af/180.0*pi
        dal = (af-ai)/(m1-1.0)

c---
c loop over alpha
c---

        Do i=1,m1

          al = ai+(i-1.0D0)*dal

          if(method.eq.1) then

            call newton 
     +
     +        (n
     +        ,menu
     +        ,Niter
     +        ,eps
     +        ,Ispeak
     +        ,Jac
     +        ,x
     +        ,Iflag
     +        )

          elseif(method.eq.2) then

           call broyden
     +
     +        (n
     +        ,menu
     +        ,Niter
     +        ,Ispeak
     +        ,x
     +        ,Bra
     +        ,Iflag
     +        )

          End If

          aa  = al/pi
          al2 = 2.0D0*al
          RLr = 1.0D0+x(1)/al2
          RLi =       x(2)/al2

          if(abs(Rli).gt.0.00001) then
           rho = pi/RLi
           sig = pi*(RLr-1.0)/RLi
          else
           rho = 0.0
           sig = 0.0
          endif

c         write (2,100) i,aa,RLr,RLi
          write (2,101)   aa,x(1),x(2),RLr,RLi,rho,sig
c         write (2,100) i,aa,x(1),x(2),RLr,RLi,rho,sig
          write (6,100) i,aa,x(1),x(2),RLr,RLi,rho,sig
c         write (6,100) i,aa,RLr,RLi

        End Do

      Go to 99

      end if

c---------------------
c Compute one solution
c---------------------

      if(method.eq.1) then

        call newton 
     +
     +      (n
     +      ,menu
     +      ,Niter
     +      ,eps
     +      ,Ispeak
     +      ,Jac
     +      ,x
     +      ,Iflag
     +      )

      elseif(method.eq.2) then

        call broyden 
     +
     +      (n,menu
     +      ,Niter
     +      ,Ispeak
     +      ,x
     +      ,Bra
     +      ,Iflag
     +      )

      end if

      Go to  98

c-----
c done
c-----

  99  Continue

      write (2,100) Null
      close (2)

 100  Format (2X,I3,9(1X,F9.6))
 101  Format (9(1X,F9.6))
 
      stop
      end
