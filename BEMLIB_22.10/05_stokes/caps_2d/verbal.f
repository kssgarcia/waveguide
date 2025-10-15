      subroutine verbal
     +
     + (Iflow
     + ,rho1,rho2
     + ,vs1,vs2
     + ,tinit
     + ,cinit
     + ,Ds
     + ,betas
     + ,elst
     + ,gac
     + ,shrt
     + ,ant_c1,ant_c2
     + ,a11,a12
     + ,a21,a22
     + ,Iglut
     + ,Max1,Max2
     + ,roexp,mexp,mexp_foam
     + ,wall
     + ,U1,U2
     + ,pg
     + ,RL
     + ,h,Ngfww,IQPD
     + ,thet0
     + ,IS_slp,IS_dlp
     + ,NGL
     + ,thmax
     + ,spmin,spmax
     + ,Iread
     + ,Ishape
     + ,aob,req
     + ,awob,bwob
     + ,xdc,udc
     + ,NSG
     + ,Isym
     + ,Isolve
     + ,eps,Nter,Idfl,JGS
     + ,Iadj,Eadj,Norm
     + ,Nprint
     + ,IRK
     + ,Dt
     + ,Istop
     + )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c-----------
c initialize
c-----------

      Istop = 0

c-----------
c preferences
c-----------

      write (6,*)
      write (6,*) '   MENU OF FLOWS'
      write (6,*)
      write (6,*) ' Enter '
      write (6,*)
      write (6,*) '  1  for a solitary capsule'
      write (6,*) '     in infinite straining flow '
      write (6,*) '  2  for a solitary capsule'
      write (6,*) '     in infinite simple shear flow '
      write (6,*) '  5  for a doubly periodic array of capsules'
      write (6,*) '     in simple shear flow'
      write (6,*) '  6  for a doubly periodic expanding foam'
      write (6,*) ' 10  for a solitary capsule above a wall'
      write (6,*) ' 11  for a periodic file of capsules'
      write (6,*) '     above a plane wall'
      write (6,*) ' 20  for periodic channel flow'
      write (6,*) '  0  to quit'
      write (6,*) ' -----------'
      read  (5,*)   Iflow

      if(Iflow.eq.0) then
        Istop = 1
        Go to 99
      end if

c---
c physical properties
c---

      if(Iflow.gt.5) then

       write (6,*)
       write (6,*) ' Enter the density of the ambient fluid'
       write (6,*) ' --------------------------------------'
       read  (5,*)  rho1

       write (6,*)
       write (6,*) ' Enter the density of the capsule fluid'
       write (6,*) ' --------------------------------------'
       read  (5,*)  rho2

      end if

      write (6,*)
      write (6,*) ' Enter the viscosity of the ambient fluid'
      write (6,*) ' ----------------------------------------'
      read  (5,*)   vs1

      write (6,*)
      write (6,*) ' Enter the viscosity of the capsule fluid'
      write (6,*) ' ----------------------------------------'
      read  (5,*)   vs2

      write (6,*)
      write (6,*) ' Enter the initial surface tension'
      write (6,*) ' ---------------------------------'
      read  (5,*)  tinit

      write (6,*)
      write (6,*) ' Enter the initial surfactant concentration'
      write (6,*) ' ------------------------------------------'
      read  (5,*)  cinit

      write (6,*)
      write (6,*) ' Enter the surfactant diffusivity'
      write (6,*) ' --------------------------------'
      read  (5,*)  Ds

      write (6,*)
      write (6,*) ' Enter beta: sensitivity of surface tension'
      write (6,*) '             to surfactant concentration'
      write (6,*) ' -----------------------------------------'
      read  (5,*)  betas

      write (6,*)
      write (6,*) ' Enter the interface modulus of elasticity'
      write (6,*) ' -----------------------------------------'
      read  (5,*)  elst

      write (6,*)
      write (6,*) ' Enter the acceleration of gravity '
      write (6,*) ' ----------------------------------'
      read  (5,*)  gac

c---

      if   (Iflow.eq.1.or.Iflow.eq.2
     +  .or.Iflow.eq.5.or.Iflow.eq.10
     +  .or.Iflow.eq.11
     +  ) then

        write (6,*)
        write (6,*) ' Enter the shear rate or rate of extension '
        write (6,*) ' ------------------------------------------'
        read (5,*)   shrt

      end if

c---

      if(Iflow.eq.1) then

        write (6,*)
        write (6,*) ' Enter the coefficient c_1'
        write (6,*) '       of the straining flow'
        write (6,*) ' -----------------------------'
        read  (5,*)   ant_c1

        write (6,*)
        write (6,*) ' Enter the coefficient c_2'
        write (6,*) '       of the straining flow'
        write (6,*) ' -----------------------------'
        read  (5,*)   ant_c2

      end if

c---------------------
c doubly-periodic flow
c---------------------

      if(Iflow.eq.5.or.Iflow.eq.6) then 

       write (6,*)
       write (6,*) ' Enter the x and y coordinates'
       write (6,*) ' of the first lattice vector'
       write (6,*) ' -----------------------------'
       read  (5,*)   a11,a12

       write (6,*)
       write (6,*) ' Enter the x and y coordinates '
       write (6,*) ' of the second the lattice vector'
       write (6,*) ' --------------------------------'
       read  (5,*)   a21,a22

       write (6,*)
       write (6,*) ' Enter 0 to compute the green function'
       write (6,*) '       1 to interpolate'
       write (6,*) ' -------------------------------------'
       read  (5,*)   Iglut

       if(Iglut.eq.0) then

       write (6,*)
       write (6,*) ' Enter the truncation limits Max1 and MAx2'
       write (6,*) '       for the Ewald summation of the'
       write (6,*) '       doubly-periodic Green function'
       write (6,*) '       in real and reciprocal space'
       write (6,*) ' ---------------------------------------'
       read  (5,*)   Max1,Max2

       end if

      end if

c---------------
c expanding bubble
c---------------

      if(vs2.lt.0.000001) then

      write (6,*)
      write (6,*) ' Please enter the rate of expansion'
      write (6,*) ' ----------------------------------'
      read  (5,*)  roexp

       if(Iflow.ne.6) then
        write (6,*)
        write (6,*) ' Please enter the method of solution '
        write (6,*) ' ------------------------------------'
        read  (5,*)  mexp
       Else
        write (6,*)
        write (6,*) ' Please enter the method of solution '
        write (6,*) ' ------------------------------------'
        read  (5,*)  mexp_foam
       end if

      end if

c-------------------
c semi-infinite flow
c-------------------

      if(Iflow.eq.10.or.Iflow.eq.11) then

        write (6,*) ' Enter the y location of the wall'
        write (6,*) ' -------------------------------'
        read  (5,*)  wall

      end if

c-------------
c channel flow
c--------------

      if(Iflow.eq.20) then

        write (6,*) ' Enter the channel half-width h'
        write (6,*) ' ------------------------------'
        read  (5,*)   h

        write (6,*) ' Choose the conditions of the flow'
        write (6,*)
        write (6,*) ' Enter  0 for fixed pressure drop'
        write (6,*) '        1 for foxed flow rate '
        write (6,*) ' -------------------------------'
        read  (5,*)   IQPD

        write (6,*) ' Enter the truncation limit for summing'
        write (6,*) ' the Green function for channel flow'
        write (6,*) ' --------------------------------------'
        read  (5,*)   Ngfww

        write (6,*) ' Enter lower and upper wall velocities'
        write (6,*) ' -------------------------------------'
        read (5,*)  U1,U2

        write (6,*) ' Enter the negative of the presure gradient'
        write (6,*) ' ------------------------------------------'
        read (5,*)   pg

        write (6,*)
        write (6,*) ' Enter the channel inclination'
        write (6,*) '       angle in multiples of pi'
        write (6,*) ' ------------------------------'
        read (5,*)  thet0

      end if

c---------------------
c singly-periodic flow
c---------------------

      if(menu.eq.11.or.menu.eq.20  ! a single file
     +  ) then

       write (6,*) ' Enter the period: L '
       write (6,*) ' --------------------'
       read  (5,*)   RL

      end if

c--------------------------------
c further input
c related to the numerical method
c--------------------------------

      write (6,*)
      write (6,*) ' Enter 1 to desingularize'
      write (6,*) '       the single-layer potential'
      write (6,*) '       0 for no'
      write (6,*) ' --------------'
      read  (5,*)   IS_slp

      write (6,*)
      write (6,*) ' Enter 1 to desingularize'
      write (6,*) '       the double-layer potential'
      write (6,*) '       0 for no'
      write (6,*) ' --------------'
      read  (5,*)   IS_dlp

      write (6,*)
      write (6,*) ' Enter the number of integration quadature'
      write (6,*) ' points over each arc'
      write (6,*)
      write (6,*) ' Choose from: 1,2,3,4,5,6,12,20; default is 20'
      write (6,*) ' ---------------------------------------------'
      read  (5,*)   NGL

      write (6,*)
      write (6,*) ' Enter the maximum angle subtended by an arc'
      write (6,*) '                          in multiples of pi'
      write (6,*) ' -------------------------------------------'
      read  (5,*)   thmax

      write (6,*)
      write (6,*) ' Enter the minimum and maximum '
      write (6,*) '              point separation'
      write (6,*) ' -----------------------------'
      read  (5,*)   spmin,spmax

c---
c Initial conditions
c---

      write (6,*)
      write (6,*) '  Enter 0 to generate the initial shape'
      write (6,*) '          of the interface'
      write (6,*) '        1 to read it from file: caps_2d.inp'
      write (6,*) ' -------------------------------------------'
      read  (5,*)   Iread

      write (6,*)
      write (6,*) ' Choose the initial or unstressed shape'
      write (6,*)
      write (6,*) ' Enter 1 for a circular or elliptical shape'
      write (6,*) '       2 for a wobbly shape'
      write (6,*) ' ---------------------------------------'
      read  (5,*)   Ishape

      if(Ishape.eq.1) then

        write (6,*)
        write (6,*) ' Circular or elliptical capsule'
        write (6,*)
        write (6,*) ' Enter the axis ratio'
        write (6,*) ' --------------------'
        read  (5,*)  aob

        write (6,*)
        write (6,*) ' Enter the equivalent capsule radius'
        write (6,*) ' -----------------------------------'
        read  (5,*)  req

      else

        write (6,*)
        write (6,*) ' Wobbly capsule'
        write (6,*)
        write (6,*) ' Enter the wobbly coefficients a and b'
        write (6,*) ' -------------------------------------'
        read  (5,*)  awob,bwob

      end if

      write (6,*)
      write (6,*) ' Enter the x-y position of the capsule center'
      write (6,*) ' --------------------------------------------'
      read  (5,*) xdc,ydc

      if(Iread.eq.0) then

       write (6,*)
       write (6,*) ' Enter the number of points'
       write (6,*) '       around the interface'
       write (6,*) ' --------------------------'
       read  (5,*) NSG

      end if

      write (6,*)
      write (6,*) ' Is shape symmetric ?'
      write (6,*)
      write (6,*) '  Enter 0 for no '
      write (6,*) '        2 for dual symmetry wrt x and y axes'
      write (6,*) ' -------------------------------------------'
      read  (5,*)   Isym

c---
c integral equation
c---

      write (6,*)
      write (6,*) ' Integral Equation'
      write (6,*)
      write (6,*) ' Enter 1 for a direct solution'
      write (6,*) '       2 for iterative solution'
      write (6,*) ' ------------------------------'
      read  (5,*)   Isolve

      If(Isolve.eq.2) then

       write (6,*)
       write (6,*) ' Enter: (a) Iteration tolerance '
       write (6,*) '        (b) Maximum Number of Iterations '
       write (6,*) '        (c) 0, 1 or 2 for Deflation '
       write (6,*) '        (d) 1 for Jacobi, 2 for Gauss Siedel '
       write (6,*) ' -------------------------------------------'
       read  (5,*)  eps,Nter,Idfl,JGS

      end if

c---
c related to time stepping
c---

      write (6,*)
      write (6,*) ' Adjust time step according to'
      write (6,*) '        maximum curvature ?'
      write (6,*)
      write (6,*) ' Enter 0 for NO; 1 or 2 for yes'
      write (6,*) ' ------------------------------'
      read  (5,*)  Iadj

      if(Iadj.ne.0) then
        write (6,*)
        write (6,*) ' Please enter the adjustment constant'
        write (6,*) ' ------------------------------------'
        read  (5,*)  Eadj
      end if

      write (6,*)
      write (6,*) ' Enter 1 to normalize the capsule area '
      write (6,*) '         after each step'
      write (6,*) '       0 for No'
      write (6,*) ' -----------------------------------'
      read (5,*)  Norm

      write (6,*)
      write (6,*) ' After how many steps of printout ? '
      write (6,*) ' -----------------------------------'
      read  (5,*)  Nprint

      write (6,*)
      write (6,*) ' Time integration by the Runge-Kutta method'
      write (6,*)
      write (6,*) ' Enter 1 for the first-order method'
      write (6,*) '       2 for the second-order method'
      write (6,*) ' -----------------------------------'
      read  (5,*)   IRK

      write (6,*)
      write (6,*) ' Please enter the time step'
      write (6,*) ' --------------------------'
      read  (5,*)   Dt

      write (6,*)
      write (6,*) ' Point advancement policy'
      write (6,*)
      write (6,*) ' Enter 0 for the marker points to move with'
      write (6,*) '         the fluid velocity'
      write (6,*) '       1 for the marker points to move with'
      write (6,*) '         the velocity component normal to '
      write (6,*) '         the interface '
      write (6,*) ' ------------------------------------------'
      read  (5,*)   Move

c-----
c done
c-----

   99 Continue

      return
      end
