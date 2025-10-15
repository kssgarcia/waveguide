      program sgf_2d_fs_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------
c Driver for Green's function of
c Stokes flow in free space
c
c Iopt = 1 computes G
c Iopt = 2 computes G and T
c------------------------------
 
      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      Iopt = 2

c-------
c launch
c-------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' 1 for one evaluation'
      write (6,*) ' 2 to test integral identities'
      write (6,*) '------------------------------'
      read  (5,*) menu

      If(menu.eq.0) Go to 99

  98  Continue

      write (6,*)
      write (6,*) ' Enter the coordinates of the point force'
      write (6,*) ' 99 to quit'
      write (6,*) ' ----------'
      read  (5,*) x0,y0

      If(x0.eq.99) Go to 99
      If(y0.eq.99) Go to 99

c---------------
c One evaluation
c---------------

      If(menu.eq.1) then

      write (6,*)
      write (6,*) ' Enter the coordinates of the velocity point'
      write (6,*) ' 99 to quit'
      write (6,*) ' ----------'
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99

      call sgf_2d_fs
     +
     +  (Iopt
     +  ,x,y
     +  ,x0,y0
     +  ,Gxx,Gxy
     +  ,Gyx,Gyy
     +  ,px,py
     +  ,Txxx,Txxy,Tyxx,Tyxy
     +  ,Txyx,Txyy,Tyyx,Tyyy
     +  )

      write (6,*) ' --------------------------------'
      write (6,*) ' Green function and stress tensor'
      write (6,*)
      write (6,*)  " velocity Green function"
      write (6,*)
      write (6,100) Gxx,Gxy
      write (6,100) Gyx,Gyy
      write (6,*)
      write (6,*)  " stress tensor"
      write (6,*)
      write (6,100) Txxx,Txxy
      write (6,100) Tyxx,Tyxy
      write (6,*)
      write (6,100) Txyx,Txyy
      write (6,100) Tyyx,Tyyy
      write (6,*) ' --------------------------------'

c-------------------------------
c Test of integral properties
c by integrating around a circle
c-------------------------------

      Else If(menu.eq.2) then

      write (6,*)
      write (6,*) " Enter the coordinates of the center"
      write (6,*) " of the test circle"
      write (6,*) " ------------------"
      read (5,*) xcnt,ycnt

      write (6,*)
      write (6,*) " Enter the radius of the test circle"
      write (6,*) " -----------------------------------"
      read (5,*) rad

      write (6,*)
      write (6,*) " Enter the number of integration points "
      write (6,*) " ---------------------------------------"
      read (5,*) mint

      sgfx = 0.0D0
      sgfy = 0.0D0

      sm11 = 0.0D0
      sm12 = 0.0D0
      sm21 = 0.0D0
      sm22 = 0.0D0

      dth = pi2/mint

      Do i=1,mint

       th = (i-1.0D0)*dth
       cs = Dcos(th)
       sn = Dsin(th)
       x  = xcnt + rad*cs
       y  = ycnt + rad*sn

       call sgf_2d_fs
     +
     +   (Iopt
     +   ,x,y
     +   ,x0,y0
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

       vnx = cs
       vny = sn

       sgfx = sgfx + vnx*Gxx + vny*Gyx
       sgfy = sgfy + vnx*Gxy + vny*Gyy

       sm11 = sm11 + Txxx*vnx + Txxy*vny
       sm12 = sm12 + Tyxx*vnx + Tyxy*vny
       sm21 = sm21 + Txyx*vnx + Txyy*vny
       sm22 = sm22 + Tyyx*vnx + Tyyy*vny

      End Do

      cf  = dth*rad/pi4

      sm11 = sm11 * cf
      sm12 = sm12 * cf
      sm21 = sm21 * cf
      sm22 = sm22 * cf


      write (6,*)
      write (6,*) " Should be 0:"
      write (6,*)
      write (6,100) sgfx,sgfy
      write (6,*)
      write (6,*) " Should be 0 or -1:"
      write (6,*)
      write (6,100) sm11,sm12
      write (6,100) sm21,sm22
      write (6,100)

c--------
      End If
c--------

      Go to 98

c-----
c Done
c-----

  99  Continue
 
 100  Format (3(2x,f15.10))
 110  Format (1x,i3,1x,f15.10,1x,f20.5)

      Stop
      End
