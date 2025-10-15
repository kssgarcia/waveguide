      program sgf_2d_2p_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------------------
c Driver for the doubly-periodic Green's function
c of two-dimensional Stokes flow
c
c The Green's function yields zero
c flow rate across the face of each flow cell
c
c Important note:
c --------------
c
c In the case of direct evaluation,
c the subroutines
c
c ewald_2d, qqq_2d, vvv_2d,
c
c must be called
c before sgf_2d_2p, as explained in the User Guide
c
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c Lookup tables for interpolation
c---

      Dimension glut_xx(-64:64,0:64,-64:64)
      Dimension glut_xy(-64:64,0:64,-64:64)
      Dimension glut_yy(-64:64,0:64,-64:64)

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

      srpi = Dsqrt(pi)
      
c-----------
c input data
c-----------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) "  1 for one evaluation"
      write (6,*) "  2 for a test of integral identities on a circle"
      write (6,*) "  0 to quit"
      write (6,*) " ----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) "  1 for direct evaluation"
      write (6,*) "  2 for interpolation"
      write (6,*)
      write (6,*) "     Interpolation implemented for G_ij only"
      write (6,*) "     and for base vectors: a11=1 a12=0 a22 = 1;"
      write (6,*) "     rescale for different cases"
      write (6,*) "     Look up tables will be read from file:"
      write (6,*) "     sgf_2d_2p_glut.inp"
      write (6,*) " ----------------------"
c     read  (5,*) method
      method = 1

      if(method.eq.1) Iopt = 2   ! will compute G and T
      if(method.eq.2) Iopt = 1   ! will compute G only

c------------------------------
c read the interpolation tables
c from file: sgf_2d_2p_glut.dat
c------------------------------

      if(method.eq.2) then

       open (4,file="sgf_2d_2p_glut.inp")

       write (6,*)
       write (6,*)  " Reading the look-up tables"
       write (6,*)  " This may take a while"
       write (6,*)  " Please twiddle your thumbs in anticipation"
       write (6,*)

       read (4,*) Na21,Nxxx,Nyyy

       Do k=-Na21,Na21
         Do j=0,Nyyy
          Do i=-Nxxx,Nxxx
            read (4,104) Glut_xx(i,j,k),Glut_xy(i,j,k)
     +                  ,Glut_yy(i,j,k)
          End Do
        End Do
       End Do

       close (4)

      End If

c-------------
  98  Continue     ! target for repeat
c-------------

      write (6,*)
      write (6,*) ' Enter the coordinates of a point force'
      write (6,*)
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
c     read  (5,*) x0,y0
      x0 = 0.0D0
      y0 = 0.0D0

      if(x0.eq.99) Go to 99
      if(y0.eq.99) Go to 99

c----------------------------
c input for direct evaluation
c----------------------------

      if(method.eq.1) then

       write (6,*)
       write (6,*) " Enter the coordinates of the"
       write (6,*) "       first lattice vector"
       write (6,*) "  99 to quit"
       write (6,*) " -----------"
c      read  (5,*) a11,a12
       a11 = 1.0D0
       a12 = 0.0D0

       if(a11.eq.99) Go to 99
       if(a12.eq.99) Go to 99

       write (6,*)
       write (6,*) " Enter the coordinates of the"
       write (6,*) "        second lattice vector"
       write (6,*) "        99 to quit"
       write (6,*) " ----------------------------"
c      read  (5,*) a21,a22
       a21 = 0.0D0
       a22 = 1.0D0

       if(a21.eq.99) Go to 99
       if(a22.eq.99) Go to 99

c--------------------------------------------
c ewald_2d will produce the reciprocal vectors
c and the optimal value of xi
c--------------------------------------------

       call sgf_2d_2p_ewald
     +
     +    (a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,tau
     +    )

       ewsave = ew

       write (6,*) 
       write (6,*) " Reciprocal lattice vectors:"
       write (6,*) " --------------------------"
       write (6,100) b11,b12
       write (6,100) b21,b22

       write (6,*) 
       write (6,*) " Enter the ewald parameter xi"
       write (6,*) 
       write (6,*) " 99 for the optimal value ",ewsave
       write (6,*) " 0 to quit"
       write (6,*) " ---------"
       read  (5,*) ew

       if(ew.eq.0 ) Go to 99
       if(ew.eq.99) ew = ewsave

       write (6,*)
       write (6,*) " Enter the summation limits max1, max2"
       write (6,*) " in real and reciprocal space"
       write (6,*)
       write (6,*) " Must be less than 9"
       write (6,*) " -------------------"
c      read  (5,*) max1,max2
       max1 = 6;
       max2 = 6;

c--------------------------------------------
c  qqq_2d will produce an array used to compute
c  the sum of the velocity Green's function 
c  in reciprocal space
c--------------------------------------------

       call sgf_2d_2p_qqq
     +
     +   (b11,b12
     +   ,b21,b22
     +   ,max2,ew
     +   )

c-----------------------------------------
c  vvv_2d will produce an array used to compute
c  the sum of the stress Green's function
c  in reciprocal space
c-----------------------------------------

       if(Iopt.eq.2) then

        call sgf_2d_2p_vvv
     +
     +    (b11,b12
     +    ,b21,b22
     +    ,max2,ew
     +    )

       end if

      end if

c-------------------------
c input for interpolation
c-------------------------

      If(method.eq.2) then

       a11 = 1.0D0
       a12 = 0.0D0
       a22 = 1.0D0

       write (6,*) 
       write (6,*) " Will take: a11 = 1.0, a12 = 0.0"
       write (6,*) "                       a22 = 1.0"
       write (6,*) 
       write (6,*) " Please enter: a21"
       write (6,*) " 99 to quit"
       write (6,*) " ----------"
       read  (5,*) a21

       If(a21.eq.99) Go to 99

       write (6,*)
       write (6,*) " Coordinates of the lattice vectors:"
       write (6,*) " -----------------------------------"
       write (6,100) a11,a12
       write (6,100) a21,a22

       write (6,*)
       write (6,*) " One point force located at the origin"
       write (6,*)

      End If

c---------------
c One evaluation
c---------------

      If(menu.eq.1) then
 
       write (6,*)
       write (6,*) " Enter field-point coordinates x, y"
       write (6,*)
       write (6,*) " 99 for either to quit"
       write (6,*) " ---------------------"
       read  (5,*) x, y

       If(x.eq.99) Go to 99
       If(y.eq.99) Go to 99

       If(method.eq.1) then

       call sgf_2d_2p 
     +
     +       (Iopt
     +       ,x,y
     +       ,x0,y0
     +       ,a11,a12,a21,a22
     +       ,b11,b12,b21,b22
     +       ,ew,tau
     +       ,max1,max2
     +       ,gxx,gxy
     +       ,gyx,gyy
     +       ,px,py
     +       ,txxx,txxy,tyxx,tyxy
     +       ,txyx,txyy,tyyx,tyyy
     +       )

       Else

       call sgf_2d_2p_int
     +
     +      (x,y
     +      ,x0,y0
     +      ,a21
     +      ,Gxx,Gxy
     +      ,Gyx,Gyy
     +      )

       End If

       write (6,*) ' --------------------------------'
       write (6,*) ' Green function for the velocity:'
       write (6,*)
       write (6,100) Gxx,Gxy
       write (6,100) Gyx,Gyy
       write (6,*)
       write (6,*) ' Green function for the pressure:'
       write (6,*)
       write (6,100) px,py
       write (6,*)
       write (6,*) ' Green function for the stress:'
       write (6,*)
       write (6,100) Txxx,Txxy
       write (6,100) Tyxx,Tyxy
       write (6,100)
       write (6,100) Txyx,Txyy
       write (6,100) Tyyx,Tyyy
       write (6,*) ' --------------------------------'

      End If

c----------------------------
c Test of integral identities
c----------------------------

      If(menu.eq.2) then
 
 96   Continue

      write (6,*)
      write (6,*) ' Enter the coordinates of the center'
      write (6,*) ' of the test circle'
      write (6,*) ' 99 for either to quit'
      write (6,*) ' ---------------------------'
      read  (5,*) xcnt,ycnt

      If(xcnt.eq.99) Go to 99
      If(ycnt.eq.99) Go to 99

      write (6,*)
      write (6,*) ' Enter the radius of the test circle'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) rad

      If(rad.eq.0) Go to 99

      write (6,*)
      write (6,*) ' Will integrate using the trapezoidal rule'
      write (6,*) ' Enter the number of intervals'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) mint

      If(mint.eq.0) Go to 99

c---
c initialize
c---

      sgfx = 0.0
      sgfy = 0.0

      If(method.eq.1) then
       sm11 = 0.0D0
       sm12 = 0.0D0
       sm21 = 0.0D0
       sm22 = 0.0D0
      End If

      dth = pi2/mint

c---
c trapezoidal rule
c---

      Do i=1,mint

        th  = (i-1.0D0)*dth
        cs  = Dcos(th)
        sn  = Dsin(th)
        x   = xcnt + rad*cs
        y   = ycnt + rad*sn
        vnx = cs              ! normal vector
        vny = sn              ! normal vector

        if(method.eq.1) then

         call sgf_2d_2p
     +
     +     (Iopt
     +     ,x,y
     +     ,x0,y0
     +     ,a11,a12,a21,a22
     +     ,b11,b12,b21,b22
     +     ,ew,tau
     +     ,max1,max2
     +     ,gxx,gxy
     +     ,gyx,gyy
     +     ,px,py
     +     ,txxx,txxy,tyxx,tyxy
     +     ,txyx,txyy,tyyx,tyyy
     +      )

        Else

         call sgf_2d_2p_int
     +
     +      (x,y
     +      ,x0,y0
     +      ,a21
     +      ,Gxx,Gxy
     +      ,Gyx,Gyy
     +      )

        End If

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
      write (6,*) " Should be zero or negative of the identity matrix:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22

      Go to 96

c-----------
      End If
c-----------

      Go to 98     ! return to repeat

c-----
c Done
c-----

 99   Continue

 100  Format (8(1x,f15.10))
 104  Format (8(1x,f15.10))

      Stop
      End
