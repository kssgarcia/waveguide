      subroutine verbal
     +
     + (Iflow
     + ,cinit
     + ,tinit
     + ,Diff
     + ,betas
     + ,gac
     + ,th0
     + ,IS_slp,IS_dlp
     + ,thmax,SPmin,SPmax
     + ,Iread,Ishape
     + ,aob,req,awob,bwob
     + ,xdc,ydc
     + ,Isolve
     + ,eps,Nter,Idfl,JGS
     + ,Norm,Nprint,IRK,Dt,Move
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

      Dimension  rhod(49),  vsd(49)
      Dimension  elst(49)
      Dimension cinit(49),tinit(49),Diff(49),betas(49)

      Dimension xdc(49),ydc(49)
      Dimension aob(49),req(49),awob(49),bwob(49)

      Dimension  NSG(49),NSG1(49),NSG2(49)

c---
c various common blocks
c---

      common/ancR2/vs1,vsd
      common/ancR4/rho1,rhod,elst
      common/ancR5/RL
      common/ancR6/shrt,U1,U2,pg

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI1/Ndrops,NSG,NSG1,NSG2
      common/ancI2/NGFww,IQPD
      common/ancI3/NGL

c---
c doubly-periodic flow
c---

      common/aaaa/a11,a12,a21,a22
      common/bbbb/b11,b12,b21,b22
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

      common/sgf2d2p/Iglut

c------------
c preferences
c------------

      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) '  0  to quit'
      write (6,*) '  1  for uniaxial elongational flow'
      write (6,*) '  2  for simple shear flow'
      write (6,*) '  5  for a doubly periodic suspension'
      write (6,*) '     in simple shear flow'
      write (6,*) ' 10  for semi-infinite simple shear flow'
      write (6,*) '     above a plane wall'
      write (6,*) ' 11  for a periodic suspension'
      write (6,*) '     in semi-infinite simple shear flow'
      write (6,*) '     above a plane wall'
      write (6,*) ' 20  for periodic Couette channel flow '
      write (6,*) ' 21  for periodic Poiseuille channel flow '
      write (6,*) ' 22  for periodic gravity-driven channel flow'
      write (6,*) ' --------------------------------------------'
      read  (5,*)  Iflow

      If(Iflow.eq.0) stop

      write (6,*)
      write (6,*) ' Enter the number of drops'
      write (6,*) ' -------------------------'
      read  (5,*)  Ndrops

c--------------------
c physical properties
c--------------------

      write (6,*)
      write (6,*) ' Enter the density of the ambient fluid'
      write (6,*) ' --------------------------------------'
      read  (5,*)  rho1

      write (6,*)
      write (6,*) ' Enter the drop densities, one by one'
      write (6,*) ' ------------------------------------'
      read  (5,*)  (rhod(i),i=1,Ndrops)

      write (6,*)
      write (6,*) ' Enter the viscosity of the ambient fluid'
      write (6,*) ' ----------------------------------------'
      read  (5,*)  vs1

      write (6,*)
      write (6,*) ' Enter the drop viscosities, one by one'
      write (6,*) ' --------------------------------------'
      read  (5,*)  (vsd(i),i=1,Ndrops)

      write (6,*)
      write (6,*) ' Enter the initial surface tensions, one by one'
      write (6,*) ' ----------------------------------------------'
      read  (5,*) (tinit(i),i=1,Ndrops)

      write (6,*)
      write (6,*) ' Enter the initial surfactant concentrations'
      write (6,*) ' -------------------------------------------'
      read  (5,*) (cinit(i),i=1,Ndrops)

      write (6,*)
      write (6,*) ' Enter the surfactant diffusivity'
      write (6,*) ' over each drop interface'
      write (6,*) ' --------------------------------'
      read  (5,*) (Diff(i),i=1,Ndrops)

      write (6,*)
      write (6,*) ' Enter the surfactant sensitivity parameter: beta'
      write (6,*) '       for each drop'
      write (6,*) ' ------------------------------------------------'
      read  (5,*) (betas(i),i=1,Ndrops)

      write (6,*)
      write (6,*) ' Enter the interface elasticity modulii'
      write (6,*) ' --------------------------------------'
      read  (5,*) (elst(i),i=1,Ndrops)

      write (6,*)
      write (6,*) ' Enter the acceleration of gravity '
      write (6,*) ' ----------------------------------'
      read  (5,*) gac

c----------------------------------
c elongational or simple shear flow
c----------------------------------

      If   (Iflow.eq.1.or.Iflow.eq.2
     +  .or.Iflow.eq.5.or.Iflow.eq.10
     +  .or.Iflow.eq.11
     +  ) then

         If(Iflow.eq.1) then
          write (6,*)
          write (6,*) ' Enter the rate of elongation '
          write (6,*) ' -----------------------------'
         Else
          write (6,*)
          write (6,*) ' Enter the shear rate '
          write (6,*) ' ---------------------'
         End If

        read  (5,*)   shrt

      End If


c---------------------
c doubly-periodic flow
c---------------------

      If(Iflow.eq.5) then

        write (6,*) 
        write (6,*) ' Enter the x and y coordinates '
        write (6,*) ' of the first lattice vector'
        write (6,*) ' ---------------------------'
        read  (5,*)   a11,a12

        write (6,*)
        write (6,*) ' Enter the x and y coordinates'
        write (6,*) ' of the second lattice vector'
        write (6,*) ' -----------------------------'
        read  (5,*)   a21,a22

        write (6,*)
        write (6,*) ' Enter 0 to compute the green function'
        write (6,*) '       1 to interpolate'
        write (6,*) ' -------------------------------------'
        read  (5,*)   Iglut

        If(Iglut.eq.0) then

        write (6,*)
        write (6,*)
        write (6,*) ' Enter truncation limits for Ewald summation'
        write (6,*) ' of the doubly-periodic Green function'
        write (6,*) ' -------------------------------------------'
        read  (5,*)   Max1,Max2

        End If

      End If

c-------------------
c semi-infinite flow
c-------------------

      If(Iflow.eq.10.or.Iflow.eq.11) then

        write (6,*) 
        write (6,*) ' Enter the y location of the wall'
        write (6,*) ' -------------------------------'
        read  (5,*)  wall

      End If

c------------
c channel flow
c-------------

      If(Iflow.eq.20) then

        write (6,*) 
        write (6,*) ' Enter lower and upper wall velocities'
        write (6,*) ' -------------------------------------'
        read  (5,*)  U1,U2

        write (6,*) 
        write (6,*) ' Enter the Presure Gradient pg'
        write (6,*) ' -----------------------------'
        read  (5,*)   pg

        write (6,*)
        write (6,*) ' Enter the channel inclination angle'
        write (6,*) ' in multiples of pi '
        write (6,*) ' -------------------'
        read  (5,*) th0

      End If

c---------------------
c singly-periodic flow
c---------------------

      If(menu.eq.11.or.menu.eq.20) then

        write (6,*)
        write (6,*) ' Enter the period: L '
        write (6,*) ' --------------------'
        read  (5,*)   RL

      End If


c-------------
c channel flow
c-------------

      If(Iflow.eq.20) then

        write (6,*)
        write (6,*) ' Enter the channel half-width H'
        write (6,*) ' ------------------------------'
        read  (5,*)  h

        write (6,*)
        write (6,*) ' Enter the N for summation of the gf '
        write (6,*) ' ------------------------------------'
        read  (5,*) NGFww

        write (6,*)
        write (6,*) ' Choose the Green function '
        write (6,*)
        write (6,*) ' Enter 0 for zero pressure drop'
        write (6,*) '       1 for zero flow rate '
        write (6,*) ' ------------------------------'
        read  (5,*) IQPD

      End If


c--------------
c further input
c--------------

      write (6,*)
      write (6,*) ' Enter 1 to desingularize'
      write (6,*) '       the single-layer potential'
      write (6,*) '       0 for no'
      write (6,*) ' ------------------------------------------'
      read  (5,*)   IS_slp

      write (6,*)
      write (6,*) ' Enter 1 to desingularize'
      write (6,*) '       the double-layer potential'
      write (6,*) '       0 for no'
      write (6,*) ' ------------------------------------------'
      read  (5,*)   IS_dlp

      write (6,*)
      write (6,*) ' Enter the number of integration quadature'
      write (6,*) ' points over each arc'
      write (6,*)
      write (6,*) ' Choose from: 2,6,12,20; default is 20'
      write (6,*) ' ------------------------------------------'
      read  (5,*)   NGL

      write (6,*)
      write (6,*) ' Enter the maximum angle subtended by an arc'
      write (6,*) '                          in multiples of pi'
      write (6,*) ' -------------------------------------------'
      read  (5,*)   thmax

      write (6,*)
      write (6,*) ' Enter the minimum and maximum '
      write (6,*) '       point separation'
      write (6,*) ' ------------------------------'
      read  (5,*)   SPmin,SPmax
c-------------------
c Initial conditions
c-------------------

      write (6,*)
      write (6,*) '  Enter 0 to generate the initial shape'
      write (6,*) '          of the interface'
      write (6,*) '        1 to read it from file: em_2d.inp'
      write (6,*) ' -------------------------------------------'
      read  (5,*)   Iread

      write (6,*) ' Choose the initial shapes'
      write (6,*)
      write (6,*) ' Enter 1 for circular or elliptical'
      write (6,*) '       2 for wobbly'
      write (6,*) ' ---------------------------------------'
      read  (5,*)   Ishape

      If(Ishape.eq.1) then

        write (6,*) ' Circular or elliptical drops'
        write (6,*)
        write (6,*) ' Enter the drop axis ratio, one by one'
        write (6,*) ' -------------------------------------'
        read  (5,*)  (aob(i),i=1,Ndrops)

        write (6,*)
        write (6,*) ' Enter the equivalent drop radii'
        write (6,*) ' one by one'
        write (6,*) ' --------------------------------'
        read  (5,*)  (req(i),i=1,Ndrops)

      Else

        write (6,*)
        write (6,*) ' Wobbly drops'
        write (6,*)
        write (6,*) ' Enter the coefficients a and b'
        write (6,*) ' one by one'
        write (6,*) ' ------------------------------'
        read  (5,*)  (awob(i),i=1,Ndrops)
        read  (5,*)  (bwob(i),i=1,Ndrops)

      End If

      write (6,*)
      write (6,*) ' Enter the x and y position of the drop centers'
      write (6,*) ' one by one'
      write (6,*) ' ----------------------------------------------'

      Do i=1,Ndrops
       read (5,*) xdc(i),ydc(i)
      End Do

      If(Iread.eq.0) then

       write (6,*)
       write (6,*) ' Enter the number of points'
       write (6,*) ' around each interface'
       write (6,*) '              Maximum is 64'
       write (6,*) ' --------------------------'

       Do i=1,Ndrops
        read (5,*) NSG(i)
       End Do

      End If

c---
c integral equation
c---

      write (6,*)
      write (6,*) 'Boundary Integral Equation'
      write (6,*)
      write (6,*) ' Enter 1 for direct solution'
      write (6,*) '       2 for iterative solution'
      write (6,*) ' ------------------------------'
      read  (5,*)   Isolve

      If(Isolve.eq.2) then

       write (6,*)
       write (6,*) ' Enter Iteration tolerance '
       write (6,*) '       Maximum Number of Iterations '
       write (6,*) '       0, 1, or 2 for Deflation '
       write (6,*) '       1 for Jacobi, 2 for Gauss Siedel '
       write (6,*) ' ---------------------------------------'
       read  (5,*)   eps,Nter,Idfl,JGS

      End If

c---
c related to time stepping
c---

      write (6,*)
      write (6,*) ' Enter 1 to normalize the drop'
      write (6,*) '         area after each step'
      write (6,*) '       0 for No'
      write (6,*) ' -----------------------------'
      read (5,*)  Norm

      write (6,*)
      write (6,*) ' After how many steps of Printout ? '
      write (6,*) ' -----------------------------------'
      read  (5,*)  Nprint

      write (6,*)
      write (6,*) ' Time integration by the Runge-Kutta method'
      write (6,*)
      write (6,*) ' Enter 1 for first-order method'
      write (6,*) '       2 for second-order method'
      write (6,*) ' --------------------------------'
      read  (5,*)  IRK

      write (6,*)
      write (6,*) ' Please enter the time step'
      write (6,*) ' --------------------------'
      read  (5,*)   Dt

      write (6,*)
      write (6,*) ' Point advancement policy:'
      write (6,*)
      write (6,*) ' Enter 0 for the marker points to move with'
      write (6,*) '         the fluid velocity'
      write (6,*) '       1 for the normal velocity'
      write (6,*) ' ------------------------------------------'
      read  (5,*)   Move

c-----
c Done
c-----

      Return
      End
