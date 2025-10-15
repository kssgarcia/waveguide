      program newton2_dr

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
c Solve two nonlinear algebraic equations
c using Newton's method
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(2)
    
      common/pii/pi
 
c----------
c constants
c----------

      pi = 3.14159 265358 D0

      Null = 0

c--------
c prepare
c--------

      open (2,file="newton2.out")

c--------
c prepare
c--------

      Ispeak = 0   ! silent mode
      Ispeak = 1   ! verbose mode

      eps = 0.001D0;

c-------
c launch
c-------

 98   Continue

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 1 for some equations'
      write (6,*) ' 2 for some other equations'
      write (6,*) ' 0 to quit'
      write (6,*) '----------' 
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c--------------
c initial guess
c--------------

     
      if(menu.eq.1) then
       write (6,*) ' A solution is: 1, 2'
      end if

      write (6,*)
      write (6,*) ' Enter the components '
      write (6,*) ' of the starting vector:'
      write (6,*) ' -----------------------'

      read  (5,*) x(1),x(2)


      write (6,*)
      write (6,*) ' Please enter maximum number of iterations'
      write (6,*) ' -----------------------------------------'
      read (5,*) Niter

c---------------------
c Compute one solution
c---------------------

        call newton2 
     +
     +    (menu
     +    ,Niter
     +    ,eps
     +    ,Ispeak
     +    ,x
     +    ,Iflag
     +    )

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
