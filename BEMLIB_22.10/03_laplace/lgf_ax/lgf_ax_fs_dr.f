      program lgf_ax_fs_dr

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=========================================

c-----------------------------------------
c Driver for asixymmetric Green's function
c of Laplace's equation in free space
c-----------------------------------------
 
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

      Iopt = 2         ! need G and gradient

c-------
c launch
c-------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 0 to quit'
      write (6,*) ' 1 for one evaluation'
      write (6,*) ' 2 to test the integral identity'
      write (6,*) '--------------------------------'
      read  (5,*) menu

      If(menu.eq.0) Go to 99

  98  Continue

      write (6,*)
      write (6,*) ' Enter (x0, s0) of the singular point'
      write (6,*)
      write (6,*) ' 99 to quit'
      write (6,*) ' ----------'
      read  (5,*) x0,s0

      If(x0.eq.99) Go to 99
      If(sigma0.eq.99) Go to 99

c----------------------------
      If(menu.eq.1) then
c----------------------------

c---------------
c One evaluation
c---------------

      write (6,*)
      write (6,*) ' Enter (x,s) of the field point'
      write (6,*)
      write (6,*) ' 99 to quit'
      write (6,*) ' ----------'
      read  (5,*) x,s

      If(x.eq.99) Go to 99
      If(sigma.eq.99) Go to 99

      call lgf_ax_fs
     +
     +   (Iopt
     +   ,x,s
     +   ,x0,s0
     +   ,G
     +   ,Gx,Gs
     +   )

      write (6,*) ' ---------------------------'
      write (6,*) ' Green function and gradient'
      write (6,*) ' ---------------------------'
      write (6,*)
      write (6,100) G
      write (6,*)
      write (6,100) Gx,Gs
      write (6,*)

c----------------------------
      Else If(menu.eq.2) then
c----------------------------

c-------------------------------
c Test of integral properties
c by integrating around a circle
c
c sm1 should be 0 or -1.0
c-------------------------------


      write (6,*)
      write (6,*) " Enter the center of the test circle (x, sig)"
      write (6,*) ' --------------------------------------------'
      read  (5,*) xcnt,scnt

      write (6,*) " Enter the radius of the test circle"
      write (6,*) " -----------------------------------"
      read  (5,*) rad

      write (6,*) " Enter the number of integration points "
      write (6,*) " ---------------------------------------"
      read  (5,*) mint

      dth = pi2/mint

      sm1 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)
        x = xcnt + rad*cs
        s = scnt + rad*sn

        call lgf_ax_fs
     +
     +    (Iopt
     +    ,x,s
     +    ,x0,s0
     +    ,G
     +    ,Gx,Gs
     +    )

       vnx = cs           ! unit normal vector
       vns = sn

       sm1 = sm1 + (Gx*vnx + Gs*vns) * s

      End Do

      cf = dth*rad

      sm1 = sm1 * cf

      write (6,*)
      write (6,*) " Should be 0 or -1: "
      write (6,*)
      write (6,100) sm1
      write (6,*)

c-----------
      End If
c-----------

      Go to 98

c-----
c Done
c-----

  99  Continue
 
 100  Format (3(2x,f15.10))
 110  Format (1x,i3,1x,f15.10,1x,f20.5)

      Stop
      End
