      program sgf_2d_1p_ww_glut

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------------------
c Generate a look up table for the nonsingular
c part of the singly-periodic
c Green's function of two-dimensional Stokes flow
c in a two-dimensional channel with parallel-sided walls
c located at y= +- h
c
c glutww_ij: Green's Look Up Table
c
c Will tabulate, and interpolate for, the difference:
c
c  G-S1-S2-S3
c
c where:
c
c G  is the Green's function,
c S1 is the upper-wall Green's function
c S2 is the lower-wall Green's function
c S3 is the Stokeslet
c
c and then evaluate the Green's function by interpolation
c to confirm the accuracy of the computation.
c
c Interpolation is done with respect to: 
c
c  y0, y, x-x0
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

      Dimension glutww_xx(0:64,-64:64,-64:64)
      Dimension glutww_xy(0:64,-64:64,-64:64)
      Dimension glutww_yx(0:64,-64:64,-64:64)
      Dimension glutww_yy(0:64,-64:64,-64:64)

c---
c common statements to pass
c to the interpolation routine
c---

      common/glutww_r/glutww_xx,glutww_xy,glutww_yx,glutww_yy
      common/glutww_i/Nwwy0,Nwwy,Nwwx

c-----------
c initialize
c-----------

      Iselect = 1   ! will compute G only

      NGF = 10

c---
c input
c---

c     write (6,*)
c     write (6,*) " Please enter the period L"
c     write (6,*) " -------------------------"
c     read  (5,*) RL

      RL = 2.0

c     write (6,*)
c     write (6,*) " Please enter the channel semi-height h"
c     write (6,*) " --------------------------------------"
c     read  (5,*) h

      h = 1.0

c     write (6,*) " Please enter IQPD"
c     write (6,*) " -----------------"
c     read  (5,*) IQPD

      IQPD = 0

      write (6,*) 
      write (6,*) "Will produce a look up table with respect to:"
      write (6,*) 
      write (6,*) "  1) y0   in the range [-h, h]"
      write (6,*) "  2) y    in the range [-h, h]"
      write (6,*) "  3) x-x0 in the range [ 0, L/2]"
      write (6,*) 
      write (6,*) " Please enter the number of divisions"
      write (6,*) "        for each tabulation variable"
      write (6,*) 
      write (6,*) "        Maximum: 64"
      write (6,*) " -----------------------------------------"
      read  (5,*) Nwwy0,Nwwy,Nwwx

c-----------------
c set the point force at the origin of the x axis
c and prepare
c-----------------

      x0 = 0.0D0

      RLH = 0.5D0*RL

c---
c tabulation steps
c---

      Dy0 =   h/Nwwy0
      Dy  =   h/Nwwy
      Dx  = RLH/Nwwx

c----------------------
c loop over y0, y, x-x0
c----------------------

       write (6,*) " Will count up to ",Nwwy0

c--------------------
      Do k=-Nwwy0,Nwwy0         ! loop for y0
c--------------------

       write (6,*) k

       y0 = k*Dy0

       If(k.eq.-Nwwy0) y0 = -h+0.000001D0
       If(k.eq. Nwwy0) y0 =  h-0.000001D0

c-----------------
       Do j=-Nwwy,Nwwy             ! loop over y
c-----------------

       y = j*Dy

       If(j.eq.-Nwwy) y = -h+0.000001D0
       If(j.eq. Nwwy) y =  h-0.000001D0

c-----------------------
         Do i=0,Nwwx         ! loop over x
c-----------------------

         x = x0+i*Dx

         If(i.eq.0.and.j.eq.k) then
           x = x + 0.000001D0
           y = y + 0.000001D0
         End If

       call sgf_2d_1p_ww
     +
     +       (Iselect
     +       ,IQPD
     +       ,x,y
     +       ,x0,y0
     +       ,RL
     +       ,NGF
     +       ,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

c---
c subtract off the single-wall Green's functions
c before tabulating
c---

      Iss = 0

      call sgf_2d_1p_w
     +
     +       (Iselect
     +       ,x,y
     +       ,x0,y0
     +       ,-h
     +       ,RL
     +       ,Iss
     +       ,Gxx1,Gxy1
     +       ,Gyx1,Gyy1
     +       ,px1,py1
     +       ,Txxx1,Txxy1,Tyxx1,Tyxy1
     +       ,Txyx1,Txyy1,Tyyx1,Tyyy1
     +       )

c     Gxx1 = 0.0D0
c     Gxy1 = 0.0D0
c     Gyx1 = 0.0D0
c     Gyy1 = 0.0D0

      call sgf_2d_1p_w
     +
     +       (Iselect
     +       ,x,y
     +       ,x0,y0
     +       ,h
     +       ,RL
     +       ,Iss
     +       ,Gxx2,Gxy2
     +       ,Gyx2,Gyy2
     +       ,px2,py2
     +       ,Txxx2,Txxy2,Tyxx2,Tyxy2
     +       ,Txyx2,Txyy2,Tyyx2,Tyyy2
     +       )

c     Gxx2 = 0.0D0
c     Gxy2 = 0.0D0
c     Gyx2 = 0.0D0
c     Gyy2 = 0.0D0

c---
c add back one stokeslet
c to take it off from one of the
c single-wall Green's functions
c---

      xx   = x-x0
      yy   = y-y0
      rs   = xx**2+yy**2
      rlog = 0.5D0*log(rs)
      
      Gxx3 = -rlog + xx*xx/rs
      Gxy3 =         xx*yy/rs
      Gyx3 =         xx*yy/rs
      Gyy3 = -rlog + yy*yy/rs

c     Gxx3 = 0.0D0
c     Gxy3 = 0.0D0
c     Gyx3 = 0.0D0
c     Gyy3 = 0.0D0

      Glutww_xx(i,j,k) = Gxx - Gxx1 - Gxx2 + Gxx3
      Glutww_xy(i,j,k) = Gxy - Gxy1 - Gxy2 + Gxy3
      Glutww_yx(i,j,k) = Gyx - Gyx1 - Gyx2 + Gyx3
      Glutww_yy(i,j,k) = Gyy - Gyy1 - Gyy2 + Gyy3
         
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
      write (6,*) "        sgf_2d_1p_ww_glut.out ?"
      write (6,*) 
      write (6,*) " Enter 0 for NO, 1 for YES"
      write (6,*) "--------------------------"
      read  (5,*) Irecord

c     Irecord = 0

      If(Irecord.eq.1) then

       open (3,file="sgf_2d_1p_ww.table")
 
       write (3,*) Nwwy0,Nwwy,Nwwx

       Do k=-Nwwy0,Nwwy0
         Do j=-Nwwy,Nwwy
          Do i=0,Nwwx

            write (3,104) Glutww_xx(i,j,k)
     +                   ,Glutww_xy(i,j,k)
     +                   ,Glutww_yx(i,j,k)
     +                   ,Glutww_yy(i,j,k)

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

      write (6,*)
      write (6,*) ' Please enter y0'
      write (6,*) '      99 to quit'
      write (6,*) ' ---------------'
      read  (5,*) y0

      If(y0.eq.99) Go to 99

 95   Continue

      write (6,*)
      write (6,*) ' Enter (x,y) of the field point'
      write (6,*) '          99 for either to quit'
      write (6,*) ' ------------------------------'
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99

c---
c Direct evaluation for verification
c---

       call sgf_2d_1p_ww
     +
     +       (Iselect
     +       ,IQPD
     +       ,x,y
     +       ,x0,y0
     +       ,RL
     +       ,NGF
     +       ,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

      write (6,*) ' --------------------------------'
      write (6,*) ' Green function for the velocity'
      write (6,*)
      write (6,*)
      write (6,*) ' Directly computed values:'
      write (6,*)
      write (6,100) Gxx,Gxy
      write (6,100) Gyx,Gyy

c---
c interpolation
c---

       call sgf_2d_1p_ww_int
     +
     +       (Iselect
     +       ,x,y
     +       ,x0,y0
     +       ,RL
     +       ,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

      write (6,*)
      write (6,*) ' Interpolated values:'
      write (6,*)
      write (6,100) Gxx,Gxy
      write (6,100) Gyx,Gyy
      write (6,*)
      write (6,*) ' --------------------------------'


      Go to 95

c-----
c Done
c-----

 99   Continue

 100  Format (8(1x,f15.10))
 104  Format (8(1x,f15.10))

      Stop
      End
