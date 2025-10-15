      subroutine verbal
     +
     + (Iflow
     + ,ndiv
     + ,boa,coa
     + ,req
     + ,cxp,cyp,czp
     + ,phi1,phi2,phi3
     + ,mint
     + ,NGL
     + ,visc
     + ,Uinf
     + ,shrt
     + ,Iprec
     + ,Ireg
     + )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

      Implicit Double Precision (a-h,o-z)

      common/var/wall

      common/sgfr/a11,a12,a13,a21,a22,a23,a31,a32,a33
     +           ,b11,b12,b13,b21,b22,b23,b31,b32,b33
     +           ,ew,cell_vlm

      common/sgfi/Max1,Max2,method

c------------
c preferences
c------------


 93   Continue

      write (6,*)
      write (6,*) " Enter 1 to triangulate an octrahedron"
      write (6,*) "       2 to triangulate an icosahedron"
      write (6,*) "       0 to quit"
      write (6,*) " ---------------"
      read  (5,*) Ioctaicos
      write (6,*)

      if(Ioctaicos.eq.0) stop

      write (6,*)
      write (6,*) " Enter the level of triangulation: 0, 1, 2, 3"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) ndiv
      write (6,*)

      if(ndiv.eq.99) stop

      if(ndiv.gt.3) then
       write (6,*) 'ndiv too high; please try again'
       Go to 93
      end if

      write (6,*)
      write (6,*) " The particle will be an ellipsoid"
      write (6,*) " with semi-axes: a, b, c"
      write (6,*)
      write (6,*) " Please enter the axes ratios b/a and c/a"
      write (6,*) " ----------------------------------------"
      read  (5,*) boa,coa

      write (6,*)
      write (6,*) " Please enter the equivalent particle radius"
      write (6,*) " defined with respect to the particle volume"
      write (6,*) " -------------------------------------------"
      read  (5,*) req
      write (6,*)

      write (6,*) " Enter coordinates of the particle center"
      write (6,*) " ----------------------------------------"
      read  (5,*) cxp,cyp,czp

      write (6,*)
      write (6,*) " The particle will be rotated about"
      write (6,*) " the x, y, and z axes"
      write (6,*)
      write (6,*) " Enter the three rotation angles"
      write (6,*) "              in multiples of pi"
      write (6,*) " -------------------------------"
      read  (5,*) phi1,phi2,phi3

      write (6,*)
      write (6,*) " Non-singular integration over each triangle"
      write (6,*) " will be carried out by the m-point triangle rule."
      write (6,*)
      write (6,*) " Please enter m"
      write (6,*)
      write (6,*) " Choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) " Enter 0 to quit"
      write (6,*) " ----------------------------------"
      read  (5,*) mint

      If(mint.eq.0) stop

      write (6,*)
      write (6,*) " Singular integrals over each triangle"
      write (6,*) " will be computed by the polar integration rule"
      write (6,*) " using the m-point double Gauss-Legendre quadrature"
      write (6,*)
      write (6,*) " Please enter m"
      write (6,*)
      write (6,*) " Choose from 2, 6, 12, 20"
      write (6,*) " Enter 0 to quit"
      write (6,*) " ------------------------"
      read  (5,*) NGL

      If(NGL.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter the fluid viscosity"
      write (6,*) " -------------------------"
      read  (5,*) visc

c--------
      If(Iflow.eq.2) then
c--------

        write (6,*)
        write (6,*) " The wall is located at x = wall"
        write (6,*)
        write (6,*) " Please enter: wall "
        write (6,*) " -------------------"
        read  (5,*) wall

c--------
      Else If(Iflow.eq.3) then ! triply periodic flow
c--------

        write (6,*)
        write (6,*) " Enter the coordinates"
        write (6,*) " of the first lattice base vector"
        write (6,*) " --------------------------------"
        read  (5,*) a11,a12,a13

        write (6,*)
        write (6,*) " Enter the coordinates"
        write (6,*) " of the second lattice base vector"
        write (6,*) " ---------------------------------"
        read  (5,*) a21,a22,a23

        write (6,*)
        write (6,*) " Enter the coordinates"
        write (6,*) " of the third lattice base vector"
        write (6,*) " --------------------------------"
        read  (5,*) a31,a32,a33

        write (6,*)
        write (6,*) " Enter Max1, Max2 for summation"
        write (6,*) "       of the Green funcion "
        write (6,*) "       in real and reciprocal space  "
        write (6,*) " ------------------------------------"
        read  (5,*) Max1,Max2

c---------
      Else If(Iflow.eq.4.or.Iflow.eq.5) then  ! doubly periodic flow
c---------

        write (6,*)
        write (6,*) " Please enter the coordinates"
        write (6,*) " of the first lattice base vector"
        write (6,*) " --------------------------------"
        read  (5,*) a11,a12

        write (6,*)
        write (6,*) " Please enter the coordinates"
        write (6,*) " of the second lattice base vector"
        write (6,*) " --------------------------------"
        read  (5,*) a21,a22

        write (6,*)
        write (6,*) " Choose the Method of computing"
        write (6,*) " the Green funcion "
        write (6,*)
        write (6,*) " Enter 1 for Fourier series"
        write (6,*) "       2 for the Pozrikidis method"
        write (6,*) " ---------------------------------"
        read  (5,*) method

c---
        If(method.eq.1) then

         write (6,*)
         write (6,*) " Enter the truncation limit for summing"
         write (6,*) "       the Green funcion "
         write (6,*) "       in wave number space  "
         write (6,*) " ------------------------------"
         read  (5,*) Max2

        Else If(method.eq.2) then

         write (6,*)
         write (6,*) " Enter truncation limits for summing"
         write (6,*) "       the Green funcion "
         write (6,*) "       in real and wavenumber space "
         write (6,*) " -----------------------------------"
         read  (5,*) Max1,Max2

         write (6,*)
         write (6,*) " Enter epsilon for numerical differentiation"
         write (6,*) " -------------------------------------------"
         read  (5,*) eps

        End If
c---

        If(Iflow.eq.4) then
         write (6,*)
         write (6,*) " Enter the velocity of the incident flow"
         write (6,*) " ---------------------------------------"
         read  (5,*) Uinf
        End If

        If(Iflow.eq.5) then
         write (6,*)
         write (6,*) " Enter the shear rate far above the lattice"
         write (6,*) " ------------------------------------------"
         read  (5,*) shrt
        End If

c-----------
      End If
c-----------

      write (6,*)
      write (6,*) " Precondition the linear system?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      read  (5,*) Iprec

      write (6,*)
      write (6,*) " Regularize the linear system?"
      write (6,*)
      write (6,*) " Enter 0 for no"
      write (6,*) "       1 to set the last unknown"
      write (6,*) "       2 to set the projection"
      write (6,*) "       3 to deflate the integral equation"
      write (6,*) " ----------------------------------------"
      read  (5,*) Ireg

c-----
c Done
c-----

  99  Return

      End
