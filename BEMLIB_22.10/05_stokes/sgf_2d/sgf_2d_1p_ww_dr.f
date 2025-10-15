      program sgf_2d_1p_ww_dr

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c---------------------------------------------------------
c Driver for the 2D periodic
c Green's function of Stokes flow
c between two parallel flat plates
c
c The first  plate is located at y = -h
c The second plate is located at y =  h
c
c SYMBOLS:
c -------
c
c RL:  period
c N:   summation limit for the Green's function
c
c IQPD = 0  :  Pressure Drop (PD) = 0, Flow rate Q is finite
c IQPD = 1  :  Pressure Drop (PD) is finite, and Q = 0
c
c Iselect = 1 to compute G
c Iselect = 2 to compete G and T
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension u1(900),u2(900),u3(900)

c---
c Green Look Up Tables (GLUT) for interpolation
c---

      Dimension glutww_xx(0:64,-64:64,-64:64)
      Dimension glutww_xy(0:64,-64:64,-64:64)
      Dimension glutww_yx(0:64,-64:64,-64:64)
      Dimension glutww_yy(0:64,-64:64,-64:64)

      common/glutww_r/glutww_xx,glutww_xy,glutww_yx,glutww_yy
      common/glutww_i/Nwwy0,Nwwy,Nwwx

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c----------------
c read parameters
c----------------

      open (9,file="sgf_2d_1p_ww.dat")

       read (9,*) h
       read (9,*) RL
       read (9,*) NGF
       read (9,*) IQPD

      close (9)

c-------
c launch
c-------

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for one value "
      write (6,*) " 2 to test integral identities on G-T "
      write (6,*) " 3 for a y-profile and flow rate of G"
      write (6,*) " 4 for a y-profile and flow rate of T"
      write (6,*) " 5 for a x-profile of G"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(Menu.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) "  1 for direct evaluation"
      write (6,*) "  2 for interpolation"
      write (6,*) "    (Look up tables will be read from file:"
      write (6,*) "    sgf_2d_1p_ww.glut)"
      write (6,*) "  0 to quit"
      write (6,*) " ----------"
      read  (5,*) method

      if(method.eq.0) Go to 99

c--------
c prepare
c--------

      Iselect = 2     ! compute G, p, and T

      open (1,file="sgf_2d_1p_ww.out")

c------------------------------
c read the interpolation tables
c from file: sgf_2d_1p_ww_glut.dat
c------------------------------

      if(method.eq.2) then

       open (4,file="sgf_2d_1p_ww.glut")

       write (6,*)
       write (6,*)  " Reading the look-up tables"
       write (6,*)  " This may take a while"
       write (6,*)  " Please twiddle your thumbs in anticipation"
       write (6,*)

       read (4,*) Nwwy0,Nwwy,Nwwx

       Do k=-Nwwy0,Nwwy0
         Do j=-Nwwy,Nwwy
          Do i=0,Nwwx
            read (4,*) Glutww_xx(i,j,k),Glutww_xy(i,j,k)
     +                ,Glutww_yx(i,j,k)
     +                ,Glutww_yy(i,j,k)
          End Do
        End Do
       End Do

       close (4)

      end if

c-------------
  98  Continue     ! target for repeat
c-------------

      write (6,*) 
      write (6,*) " Enter the point-force coordinates: x0, y0"
      write (6,*) 
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x0,y0

      if(x0.eq.99) Go to 99
      if(y0.eq.99) Go to 99

c---------------
c one evaluation
c---------------

      if(menu.eq.1) then

       write (6,*) 
       write (6,*) " Enter field-point coordinates: x, y"
       write (6,*) 
       write (6,*) " 99 for either to quit"
       write (6,*) " ---------------------"
       read  (5,*) x,y

       If(x.eq.99) Go to 99
       If(y.eq.99) Go to 99

       call sgf_2d_1p_ww 
     +
     +   (Iselect
     +   ,IQPD
     +   ,x,y
     +   ,x0,y0
     +   ,RL,NGF,h
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

       write (6,*) ' --------------------------------'
       write (6,*) ' Direct evaluation'
       write (6,*)
       write (6,*) ' Green function for the velocity:'
       write (6,*)
       write (6,100) Gxx,Gxy
       write (6,100) Gyx,Gyy
       write (6,*)
       write (6,*) ' Green function for the pressure:'
       write (6,*)
       write (6,100) px,py
       write (6,*)
       write (6,*)
       write (6,*) ' Green function for the stress:'
       write (6,*)
       write (6,100) Txxx,Txxy
       write (6,100) Tyxx,Tyxy
       write (6,*)
       write (6,100) Txyx,Txyy
       write (6,100) Tyyx,Tyyy
       write (6,*) ' --------------------------------'

       if(method.eq.2) then

       call sgf_2d_1p_ww_int
     +
     +   (Iselect
     +   ,x,y
     +   ,x0,y0
     +   ,RL,h
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

       write (6,*) ' --------------------------------'
       write (6,*) ' Interpolation from lookup tables:'
       write (6,*)
       write (6,*) ' Green function for the velocity:'
       write (6,*)
       write (6,100) Gxx,Gxy
       write (6,100) Gyx,Gyy
       write (6,*)
       write (6,*) ' --------------------------------'

       end if

      end if

c-------------------------
c Test integral identities
c-------------------------

      if(Menu.eq.2) then

  96  Continue

      write (6,*)
      write (6,*) ' Will test integral identities for a test circle'
      write (6,*)
      write (6,*) ' Enter the coordinates of the center'
      write (6,*) '       of the test circle'
      write (6,*)
      write (6,*) '       99 for either to quit'
      write (6,*) ' ---------------------------'
      read  (5,*) xcnt,ycnt

      if(xcnt.eq.99) Go to 99
      if(ycnt.eq.99) Go to 99

      write (6,*)
      write (6,*) ' Enter the radius of the test circle'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) rad

      if(rad.eq.0) Go to 99

      write (6,*)
      write (6,*) ' Will integrate using the trapezoidal rule'
      write (6,*)
      write (6,*) ' Enter the number of intervals'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) mint

      if(mint.eq.0) Go to 99

c---
c initialize
c---

      sgfx = 0.0D0
      sgfy = 0.0D0

      dth = pi2/mint

      sm11 = 0.0D0
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = cos(th)
        sn = sin(th)
        x  = xcnt + rad*cs
        y  = ycnt + rad*sn
        vnx = cs             ! normal vector
        vny = sn

        if(method.eq.1) then

        call sgf_2d_1p_ww
     +
     +       (Iselect
     +       ,IQPD
     +       ,x,y
     +       ,x0,y0
     +       ,RL,NGF,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )
        else

       call sgf_2d_1p_ww_int
     +
     +       (Iselect
     +       ,x,y
     +       ,x0,y0
     +       ,RL,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

        end if

        sgfx = sgfx + vnx*Gxx + vny*Gyx
        sgfy = sgfy + vnx*Gxy + vny*Gyy
  
        sm11 = sm11 + Txxx*vnx + Txxy*vny
        sm12 = sm12 + Tyxx*vnx + Tyxy*vny
        sm21 = sm21 + Txyx*vnx + Txyy*vny
        sm22 = sm22 + Tyyx*vnx + Tyyy*vny

      End Do

      cf = dth*rad/pi4

      sgfx = sgfx * cf
      sgfy = sgfy * cf
      sm11 = sm11 * cf
      sm12 = sm12 * cf
      sm21 = sm21 * cf
      sm22 = sm22 * cf

      write (6,*)
      write (6,*) " Should be zero:"
      write (6,*)
      write (6,100) sgfx,sgfy
      write (6,*)
      write (6,*) " Should be either zero"
      write (6,*) " or negative of the identity matrix:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22
      write (6,*)

      Go to 96

      end if

c--------------------------------------
c y-profile of x-velocity and flow rate
c--------------------------------------

      if(menu.eq.3) then

      write (6,*)
      write (6,*) " Enter the x position of the profile"
      write (6,*) " -----------------------------------"
      read  (5,*) x

      write (6,*) " Enter the number of profile intervals"
      write (6,*) " -------------------------------------"
      read  (5,*) mint

      dy = 2.0D0*h/mint

      Do i=1,mint+1

        y = -h+(i-1.0D0)*dy

        if(method.eq.1) then

        call sgf_2d_1p_ww 
     +
     +   (Iselect
     +   ,IQPD
     +   ,x,y
     +   ,x0,y0
     +   ,RL,NGF,h
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

        else

        call sgf_2d_1p_ww_int
     +
     +   (Iselect
     +   ,x,y
     +   ,x0,y0
     +   ,RL,h
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

        end if

        u1(i) = Gxx
        u2(i) = Gxy

        write (6,101) i,y,u1(i),u2(i)

      End Do

c---
c compute the flow rate
c---

      Flow1 = 0.0D0
      Flow2 = 0.0D0

      Do i=1,mint
        Flow1 = Flow1 + u1(i)+u1(i+1)
        Flow2 = Flow2 + u2(i)+u2(i+1)
      End Do

      Flow1 = 0.5D0*Flow1*dy
      Flow2 = 0.5D0*Flow2*dy

      write (6,*)
      write (6,*) " x and y flow rate"
      write (6,*)
      write (6,102) Flow1
      write (6,103) Flow2

      end if

c-------------------------------------------
c y-profile of x-velocity and flow rate of T
c-------------------------------------------

      if(menu.eq.4) then

      write (6,*)
      write (6,*) " Enter the position (x, y)"
      write (6,*) " -------------------------"
      read  (5,*) x,y

      write (6,*) " Enter the x0 position"
      write (6,*) " ---------------------"
      read  (5,*) x0

      write (6,*) " Enter the number of integration intervals"
      write (6,*) " -----------------------------------------"
      read  (5,*) mint

      dy = (2.0*h)/mint

      Do i=1,mint+1

        y0 = -h+(i-1.0)*dy

        if(method.eq.1) then

        call sgf_2d_1p_ww 
     +
     +       (Iselect
     +       ,IQPD
     +       ,x,y
     +       ,x0,y0
     +       ,RL,NGF,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

        else

        call sgf_2d_1p_ww_int
     +
     +       (Iselect
     +       ,x,y
     +       ,x0,y0
     +       ,RL,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

        end if

        u1(i) = Txxx
        u2(i) = Txxy
        u3(i) = Tyxy

        write (6,101) i,y0,u1(i),u2(i),u3(i)

      End Do

      Flow1 = 0.0
      Flow2 = 0.0
      Flow3 = 0.0

      Do i=1,mint
        Flow1 = Flow1 +u1(i)+u1(i+1)
        Flow2 = Flow2 +u2(i)+u2(i+1)
        Flow3 = Flow3 +u3(i)+u3(i+1)
      End Do

      Flow1 = 0.50*Flow1*dy
      Flow2 = 0.50*Flow2*dy
      Flow3 = 0.50*Flow3*dy

      write (6,102) Flow1
      write (6,103) Flow2
      write (6,104) Flow3

      end if

c------------------
c x-profile of G_ij
c------------------

      if(menu.eq.5) then

      write (6,*)
      write (6,*) " Enter the y position of the profile"
      write (6,*) " -----------------------------------"
      read  (5,*) y

      write (6,*)
      write (6,*) " Enter the x_start and x_end"
      write (6,*) " ---------------------------"
      read  (5,*) x_start,x_end

      write (6,*) " Enter the number of points"
      write (6,*) " --------------------------"
      read  (5,*) mint

      dx = (x_end-x_start)/mint

      Do i=1,mint+1

        x = x_start+(i-1.0D0)*dx

        if(method.eq.1) then

        call sgf_2d_1p_ww 
     +
     +       (Iselect
     +       ,IQPD
     +       ,x,y
     +       ,x0,y0
     +       ,RL,NGF,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

        else

        call sgf_2d_1p_ww_int
     +
     +       (Iselect
     +       ,x,y
     +       ,x0,y0
     +       ,RL,h
     +       ,Gxx,Gxy
     +       ,Gyx,Gyy
     +       ,px,py
     +       ,Txxx,Txxy,Tyxx,Tyxy
     +       ,Txyx,Txyy,Tyyx,Tyyy
     +       )

        end if

        write (6,101) i,x,Gxx,Gxy,Gyx,Gyy
        write (1,101) i,x,Gxx,Gxy,Gyx,Gyy

      End Do

      end if

c---------------
c return to menu
c---------------

      Go to 98

  99  Continue

      close (1)

c-----
c Done
c-----

 100  format (10(1x,f12.6))
 101  format (1x,i3,10(1x,f12.08))
 102  format (" Flow rate 1 =",f20.10)
 103  format (" Flow rate 2 =",f20.10)
 104  format (" Flow rate 3 =",f20.10)
 105  format (" Flow rate 4 =",f20.10)
 106  format (" Flow rate 5 =",f20.10)
 108  format (8(1x,f15.10))

      Stop
      End
