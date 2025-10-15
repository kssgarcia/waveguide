      subroutine verbal
     +
     +  (Iflow
     +
     +  ,Ioctaicos
     +  ,ndiv
     +
     +  ,req
     +  ,boa,coa
     +  ,cx,cy,cz
     +  ,phi1,phi2,phi3
     +
     +  ,mint
     +  ,NGL
     +
     +  ,wall
     +
     +  ,a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +
     +  ,Max1,Max2
     +
     +  ,shrt
     +
     +  ,vs1,vs2
     +
     +  ,tinit
     +  ,cinit
     +
     +  ,Ds
     +
     +  ,Isurf
     +  ,betas
     +  ,psis
     +  ,Ismeth
     +
     +  ,Nter
     +  ,tol
     +  ,Idfl
     +
     +  ,Iread
     +
     +  ,Norm
     +  ,Isym_xy
     +
     +  ,IRK
     +  ,Dt
     +
     +  ,Nprint_xy
     +  ,Nprint_xyz
     +
     +  ,Move
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c----------------------------
c Input inquiries for drop_3d
c----------------------------

      Implicit double precision (a-h,o-z)

      write (6,*)
      write (6,*) " Choose the flow"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for flow in infinite space"
      write (6,*) " 2 for flow bounded by a plane wall"
      write (6,*) " 3 for triply periodic flow"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Iflow

      if(Iflow.eq.0) stop

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to discretize an octahedron"
      write (6,*) " 2 to discretize an icosahedron"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ioctaicos

      if(Ioctaicos.eq.0) stop

 88   Continue

      write (6,*)
      write (6,*) " Enter the level of triangulation"
      write (6,*)
      write (6,*) " Choose from 0, 1, 2, 3"
      write (6,*)
      write (6,*) " Enter 99 to quit"
      write (6,*) " ----------------"
      read  (5,*) Ndiv 

      if(Ndiv.eq.99) stop

      if(Ndiv.gt.3) then
        write (6,*) 'Too high; try again'
        Go to 88
      End If

      write (6,*)
      write (6,*) " Enter the equivalent drop radius"
      write (6,*) " --------------------------------"
      read  (5,*) req

      write (6,*)
      write (6,*) " The drop has an ellipsoidal initial shape"
      write (6,*) " with x,y,z semi-axis: a,b,c"
      write (6,*)
      write (6,*) " Enter the axes ratios b/a and c/a"
      write (6,*) " ---------------------------------"
      read  (5,*) boa,coa

      write (6,*)
      write (6,*) " Enter coordinates of the ellipsoid center"
      write (6,*) " -----------------------------------------"
      read  (5,*) cx,cy,cz

      write (6,*)
      write (6,*) " Enter the three rotation angles"
      write (6,*) "       about the x,y,z axes, in multiples of pi"
      write (6,*) " -----------------------------------------------"
      read  (5,*) phi1,phi2,phi3

      write (6,*)
      write (6,*) " Non-singular integration over each triangle"
      write (6,*)
      write (6,*) " Will use the m-point rule."
      write (6,*) 
      write (6,*) " Please enter m"
      write (6,*) 
      write (6,*) "   choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) "   enter 0 to quit"
      write (6,*) " -----------------------------------"
      read  (5,*) mint

      if(mint.eq.0) Stop

      write (6,*)
      write (6,*) " Singular integration over each triangle"
      write (6,*)
      write (6,*) " Will use the m-point Gauss-Legendre quadrature"
      write (6,*) 
      write (6,*) " Please enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) NGL

      if(NGL.eq.0) Stop

c----
c wall bounded flow
c----

      if(Iflow.eq.2) then

       write (6,*)
       write (6,*) " Wall is located at y = wall"
       write (6,*)
       write (6,*) " Please enter: wall"
       write (6,*) " ------------------"
       read  (5,*) wall

      end if

c----
c triply-periodic flow
c----

      if(Iflow.eq.3) then

       write (6,*) " Triply-periodic flow:"
       write (6,*)
       write (6,*) ' Enter the x, y, and z coordinates'
       write (6,*) '       of the first lattice vector'
       write (6,*) ' ---------------------------------'
       read  (5,*)  a11,a12,a13

       write (6,*) ' Repeat for the second lattice vector'
       write (6,*) ' ------------------------------------'
       read  (5,*)  a21,a22,a23

       write (6,*) ' Repeat for the third lattice vector'
       write (6,*) ' -----------------------------------'
       read  (5,*)  a31,a32,a33

       write (6,*)
       write (6,*) ' Enter the truncation limits Max1 and MAx2'
       write (6,*) '       for the Ewald summation of the'
       write (6,*) '       doubly-periodic Green function'
       write (6,*) '       in real and reciprocal space'
       write (6,*) ' ---------------------------------------'
       read  (5,*) Max1,Max2

      end if

      write (6,*)
      write (6,*) " Simple shear flow"
      write (6,*)
      write (6,*) " Enter the shear rate" 
      write (6,*) " --------------------"
      read  (5,*) shrt

      write (6,*)
      write (6,*) " Enter the viscosity of the ambient fluid"
      write (6,*) " ----------------------------------------"
      read  (5,*) vs1

      write (6,*)
      write (6,*) " Enter the viscosity of the drop"
      write (6,*) " -------------------------------"
      read  (5,*) vs2

      write (6,*)
      write (6,*) " Enter the initial surface tension"
      write (6,*) " ---------------------------------"
      read  (5,*) tinit

      write (6,*)
      write (6,*) " Enter the initial surfactant concentration"
      write (6,*) " ------------------------------------------"
      read  (5,*) cinit

      write (6,*)
      write (6,*) " Enter the surfactant surface diffusivity"
      write (6,*) " ----------------------------------------"
      read  (5,*) Ds

      write (6,*)
      write (6,*) " Choose the surface equation of state"
      write (6,*)
      write (6,*) " Enter 1 for Henry"
      write (6,*) "       2 for Langmuir"
      write (6,*) " --------------------"
      read  (5,*) Isurf

      write (6,*)
      write (6,*) " Enter the maximum number of iterations "
      write (6,*) "       in solving the integral equation"
      write (6,*) " --------------------------------------"
      read  (5,*) Nter

      write (6,*)
      write (6,*) " Enter the error tolerance "
      write (6,*) "       in solving the integral equation"
      write (6,*) " --------------------------------------"
      read  (5,*) tol

      write (6,*)
      write (6,*) " Enable deflation ? "
      write (6,*)
      write (6,*) " Enter 0 for no"
      write (6,*) "       1 for deflation of one eigenvalue "
      write (6,*) " ----------------------------------------"
      read  (5,*) Idfl

      write (6,*)
      write (6,*) " Enter 0 to generate the initial shape"
      write (6,*) "       1 to read from file: drop_3d.inp "
      write (6,*) "         (restart)"
      write (6,*) " ---------------------------------------"
      read  (5,*) Iread

      write (6,*)
      write (6,*) " Normalize volume after each step ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      read  (5,*) Norm

      write (6,*)
      write (6,*) " Exploit symmetry wr to the xy plane ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      read  (5,*) Isym_xy

      write (6,*)
      write (6,*) " Time integration"
      write (6,*)
      write (6,*) " Enter 1 for the Euler explicit method"
      write (6,*) "       2 for the RK2 method          "
      write (6,*) " -------------------------------------"
      read  (5,*) IRK

      write (6,*)
      write (6,*) " Enter the time step Dt"
      write (6,*) " ----------------------"
      read  (5,*) Dt

      write (6,*)
      write (6,*) " Will record a profile in the xy plane"
      write (6,*) "      after N steps; please enter N"
      write (6,*) " -------------------------------------"
      read  (5,*) Nprint_xy

      write (6,*)
      write (6,*) " Will record a 3D profile"
      write (6,*) "      after N steps; please enter N"
      write (6,*) " ----------------------------------"
      read  (5,*) Nprint_xyz

      write (6,*)
      write (6,*) ' Select point advancement policy '
      write (6,*)
      write (6,*) ' Enter 0 to for the total velocity '
      write (6,*) '       1 to for the normal velocity '
      write (6,*) ' -----------------------------------------'
      read  (5,*) Move

c-----
c Done
c-----

      return
      end
