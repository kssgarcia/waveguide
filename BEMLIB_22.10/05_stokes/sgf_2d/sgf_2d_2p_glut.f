      program sgf_2d_2p_glut

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------------------
c Generate a look up table for the nonsingular
c part of the doubly-periodic
c Green's function of two-dimensional Stokes flow
c with vanishing flow rate across the face of a cell.
c
c GLUT: Green's Look Up Table
c
c Will tabulate, and interpolate for, the difference G-S
c where G is the Green's function and S is the Stokeslet,
c and then evaluate the Green's function by interpolation
c to confirm the accuracy of the computation.
c
c Will take a11 = 1.0, a12 = 0.0, 
c                      a22 = 1.0,
c
c and will interpolate with respect to: a21,
c                                       x-x0,
c                                       y-y0
c
c LEGEND:
c ------
c
c Iselect = 1:   will compute G only
c Iselect = 2:   will compute G and T
c
c-----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c-------------------------------------
c Look Up Tables
c for interpolation
c
c GLUT stands for Green Look Up Table
c-------------------------------------

      Dimension glut_xx(-64:64,0:64,-64:64)
      Dimension glut_xy(-64:64,0:64,-64:64)
      Dimension glut_yy(-64:64,0:64,-64:64)

c---
c common statement to pass
c to the interpolation routine
c---

      common/glut_r/glut_xx,glut_xy,glut_yy
      common/glut_i/Na21,Nxxx,Nyyy

c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      piq  = 0.25D0*pi
      pih  = 0.50D0*pi
      pi2  = 2.00D0*pi
      pi4  = 4.00D0*pi
      pi6  = 6.00D0*pi
      pi8  = 8.00D0*pi
      srpi = sqrt(pi)
      
c---
c initialize
c---

      Iselect = 1   ! will compute G only

c---
c input
c---

      write (6,*)
      write (6,*) " Please enter Max1, Max2 for summation"
      write (6,*) "        in real and reciprocal space."
      write (6,*)
      write (6,*) " Must be less than 9"
      write (6,*) " -------------------"
      read  (5,*) Max1,Max2

      write (6,*) 
      write (6,*) "Will produce a look up table with respect to:"
      write (6,*) 
      write (6,*) "  1) a21  in the range [-0.5, 0.5]"
      write (6,*) "  2) x-x0 in the range [-0.5, 0.5]"
      write (6,*) "  3) y-y0 in the range [ 0.0, 0.5]"
      write (6,*) 
      write (6,*) " Please enter half the number of divisions"
      write (6,*) "        for each tabulation variable"
      write (6,*) " -----------------------------------------"
      read  (5,*) Na21,Nxxx,Nyyy

c--------
c prepare
c--------

      x0 = 0.0D0        ! one point force at the origin
      y0 = 0.0D0

      a11 = 1.0D0       ! fixed base vectors
      a12 = 0.0D0 
      a22 = 1.0D0  

c---
c tabulation steps
c---

      Da21 = 0.50D0/Na21
      Dxxx = 0.50D0/Nxxx
      Dyyy = 0.50D0/Nyyy

c--------------------
c loop over a21, x, y
c--------------------

c--------------------
      Do k=-Na21,Na21         ! loop for a21
c--------------------

       a21 = k*Da21

c---
c prepare
c---

       call sgf_2d_2p_ewald
     +
     +  (a11,a12,a21,a22
     +  ,b11,b12,b21,b22
     +  ,ew,tau
     +  )

       call sgf_2d_2p_qqq
     +
     +   (b11,b12
     +   ,b21,b22
     +   ,Max2,ew
     +   )

c-----------------
       Do j=0,Nyyy               ! loop over y
c-----------------

         y = y0+j*Dyyy

c-----------------------
         Do i=-Nxxx,Nxxx         ! loop over x
c-----------------------

         x = x0+i*Dxxx

         If(i.eq.0.and.j.eq.0) then
           x = x + 0.0000001
           y = y + 0.0000001
         End If

         call sgf_2d_2p 
     +
     +     (Iselect
     +     ,x,y
     +     ,x0,y0
     +     ,a11,a12,a21,a22
     +     ,b11,b12,b21,b22
     +     ,ew,tau
     +     ,Max1,Max2
     +     ,Gxx,Gxy
     +     ,Gyx,Gyy
     +     ,px,py
     +     ,txxx,txxy,tyxx,tyxy
     +     ,txyx,txyy,tyyx,tyyy
     +     )

c---
c subtract off the Stokeslet
c before tabulating
c---

          Dx = x-x0
          Dy = y-y0

          rs   = Dx**2+Dy**2
          sing = 0.5D0*log(rs)

          Glut_xx(i,j,k) = Gxx + sing - Dx*Dx/rs
          Glut_xy(i,j,k) = Gxy        - Dx*Dy/rs
          Glut_yy(i,j,k) = Gyy + sing - Dy*Dy/rs
         
         End Do
       End Do
      End Do

c--------------------------------
c End of look-up table generation
c--------------------------------

c--------------------------------------------
c Record the data in file: sgf_2d_2p_glut.out
c--------------------------------------------
 
      write (6,*) 
      write (6,*) " Record the interpolation tables in file:"
      write (6,*) "        sgf_2d_2p_glut.out ?"
      write (6,*) 
      write (6,*) " Enter 0 for NO, 1 for YES"
      write (6,*) "--------------------------"
      read  (5,*) Irecord

      If(Irecord.eq.1) then

       open (3,file="sgf_2d_2p_glut.out")
 
       write (3,*) Na21,Nxxx,Nyyy

       Do k=-Na21,Na21
         Do j=0,Nyyy
          Do i=-Nxxx,Nxxx

            write (3,104) Glut_xx(i,j,k),Glut_xy(i,j,k)
     +                   ,Glut_yy(i,j,k)

          End Do
         End Do
       End Do
 
       close (3)

      End If

c---------------------------------------

c----
c Evaluate the Green's function at a point
c by trilinear interpolation
c---

      write (6,*)
      write (6,*) ' Will interpolate'
      write (6,*)

 95   Continue

      write (6,*)
      write (6,*) ' Enter (x,y) of the field point'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) x,y

      If(x.eq.0) Go to 99
      If(y.eq.0) Go to 99

      write (6,*)
      write (6,*) ' Enter the component a21'
      write (6,*) ' -----------------------'
      read  (5,*) a21

      call sgf_2d_2p_int
     +
     +   (x,y
     +   ,x0,y0
     +   ,a21
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   )

      write (6,*) ' --------------------------------'
      write (6,*) ' Green function for the velocity'
      write (6,*)
      write (6,*) ' Interpolated values:'
      write (6,*)
      write (6,100) Gxx,Gxy
      write (6,100) Gyx,Gyy

c---
c Direct evaluation for verification
c---

      call sgf_2d_2p_ewald
     +
     +   (a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,tau
     +   )

      call sgf_2d_2p_qqq
     +
     +   (b11,b12
     +   ,b21,b22
     +   ,Max2,ew
     +   )

      call sgf_2d_2p 
     +
     +   (Iselect
     +   ,x,y
     +   ,x0,y0
     +   ,a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,tau
     +   ,Max1,Max2
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,txxx,txxy,tyxx,tyxy
     +   ,txyx,txyy,tyyx,tyyy
     +   )

      write (6,*)
      write (6,*) ' Directly computed values:'
      write (6,*)
      write (6,100) Gxx,Gxy
      write (6,100) Gyx,Gyy
      write (6,*)
      write (6,*) ' ---------------------------'

      Go to 95

c-----
c Done
c-----

 99   Continue

 100  Format (8(1x,f15.10))
 104  Format (8(1x,f15.10))

      Stop
      End
