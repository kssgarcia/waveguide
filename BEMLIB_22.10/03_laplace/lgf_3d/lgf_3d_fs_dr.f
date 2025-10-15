      program lgf_3d_fs_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------
c Driver for free space Green's function 
c of Laplace's equation
c---------------------------
 
      Implicit Double Precision (a-h,o-z)

      Dimension xx(50),yy(50),zz(50),ww(50)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi
      pi8 = 8.0D0*pi

c-----------
c initialize
c-----------

      Iopt = 2

c----------
c launching
c----------

      write (6,*)
      write (6,*) ' Enter 0 to quit'
      write (6,*) '       1 for one evaluation'
      write (6,*) '       2 for a test of integral identities'
      write (6,*) '         on a sphere'
      write (6,*) ' -----------------------------------------'
      read  (5,*) menu

      If(menu.eq.0) Go to 99

  98  Continue

      write (6,*)
      write (6,*) ' Enter the (x,y,z) coordinates'
      write (6,*) '          of the singular point'
      write (6,*) '                     99 to quit'
      write (6,*) ' ------------------------------'
      write (6,*)
      read  (5,*) x0,y0,z0

      If(x0.eq.99) Go to 99
      If(y0.eq.99) Go to 99
      If(z0.eq.99) Go to 99

c---------------
c One evaluation
c---------------

      If(menu.eq.1) then

       write (6,*)
       write (6,*) ' Enter the (x,y,z) coordinates'
       write (6,*) '            of the field point'
       write (6,*) '                    99 to quit'
       write (6,*) ' -----------------------------'
       read  (5,*) x,y,z

       If(x.eq.99) Go to 99
       If(y.eq.99) Go to 99
       If(z.eq.99) Go to 99

       call lgf_3d_fs
     +
     +   (Iopt
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

       write (6,*) ' ----------------------------'
       write (6,*) ' Green function and gradient '
       write (6,100)
       write (6,100) G
       write (6,100) 
       write (6,100) Gx,Gy,Gz
       write (6,*) ' ----------------------------'

      End If

c----------------------------
c Test of integral properties
c----------------------------

      If(menu.eq.2) then

       write (6,*)
       write (6,*) ' Enter the coodrinates of the center'
       write (6,*) '                  of the test sphere'
       write (6,*) '                          99 to quit'
       write (6,*) ' -----------------------------------'
       read  (5,*) xcnt,ycnt,zcnt

       If(xcnt.eq.99) Go to 99
       If(ycnt.eq.99) Go to 99
       If(zcnt.eq.99) Go to 99

       write (6,*)
       write (6,*) ' Enter the radius of the test sphere'
       write (6,*) '                           0 to quit'
       write (6,*) ' -----------------------------------'
       read  (5,*) rad

       If(rad.eq.0) Go to 99

       write (6,*)
       write (6,*) ' Choose the number of quadrature points '
       write (6,*)
       write (6,*) ' Select from 6, 18'
       write (6,*) ' Enter 0 to quit'
       write (6,*) ' -----------------'
       read  (5,*) Nint

       If(mint.eq.0) Go to 99

       call gauss_sph
     +
     +   (mint
     +   ,xx,yy,zz,ww
     +   )

c---
c integrate over the sphere
c---

       sm1 = 0.0D0

       Do i=1,mint

        vnx = xx(i)
        vny = yy(i)
        vnz = zz(i)

        x = xcnt+rad*vnx
        y = ycnt+rad*vny
        z = zcnt+rad*vnz

        call lgf_3d_fs
     +
     +    (Iopt
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,G
     +    ,Gx,Gy,Gz
     +    )

        sm1 = sm1 + (vnx*Gx + vny*Gy + vnz*Gz)*ww(i)

       End Do

       cf = pi4*rad**2
 
       sm1 = sm1 * cf

       write (6,*) 
       write (6,*)  " Should be -1 or 0"
       write (6,*) 
       write (6,100) sm1
 
c------------
       End If
c------------

      Go to 98

c-----
c Done
c-----

  99  Continue
 
 100  Format (3(2x,f15.10))

      Stop
      End
