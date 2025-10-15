      program sgf_2d_sph_blyth_dr

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c
c With some minor cosmetic
c edits by Marc G. Blyth
c----------------------------------------

c----------------------------------------
c Driver for Stokes Green's function in a
c domain bounded internally by a sphere.
c
c SYMBOLS:
c -------
c
c (xc,yc,zc) Center of the sphere
c a          Radius of the sphere
c
c----------------------------------------
 
      Implicit Double Precision (a-h,o-z)

      Dimension xx(50),yy(50),zz(50),ww(50)

c---------
c constants
c----------

      pi = 3.14159 265358 D0

c------------
c preferences
c------------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for one evaluation'
      write (6,*) ' 2 for a test of the continuity equation'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'

      read  (5,*) menu
      If(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) ' Enter the coordinates of the sphere center'
      write (6,*) ' ------------------------------------------'
      read  (5,*) xc,yc,zc
      
      write (6,*)
      write (6,*) ' Enter the sphere radius'
      write (6,*) ' -----------------------'
      read  (5,*) a

c-----------------
c return to repeat
c-----------------

 97   Continue

      write (6,*)
      write (6,*) ' Enter the coordinates of the point force'
      write (6,*) ' ----------------------------------------'
      read  (5,*) x0,y0,z0
      
c---
c make sure point
c is not inside the sphere
c---

      rad = sqrt(x0**2+y0**2+z0**2)

      If(rad.le.a) then
       write (6,*) ' Point is in the sphere; please try again'
       Go to 97
      End If

c---------------
c one evaluation
c---------------

  98  Continue

      If(menu.eq.1) then

          write (6,*)
          write (6,*) ' Enter: '
          write (6,*)
          write (6,*) ' 1 for a point on the sphere'
          write (6,*) ' 2 for a point off the sphere'
          write (6,*) ' 0 to quit'
          write (6,*) ' ---------'
          read  (5,*) ion

        If(ion.eq.0) Go to 99

c-------------
        If(ion.eq.1) then

          write (6,*)
          write (6,*) ' Enter the angles: theta,  phi'
          write (6,*) ' in multiples of pi'
          write (6,*) ' ------------------'
          read  (5,*) theta,phi

          theta = theta*pi
          phi   = phi  *pi

          x = xc + a*cos(theta)*cos(phi)
          y = yc + a*sin(theta)*cos(phi)
          z = zc + a*sin(phi)

c----
c Added by MGB 14/03/2011
c
c Problem if vector x-x0 is collinear with x-xc and x
c lies on the sphere's surface. Specifically, h7 and h9
c in {sgf_3d_sph.f} involve divisions by zero. The problem is purely
c numerical: if x is shifted slightly
c we get zero Green's function as required. 
c
c Remedy: Set G to zero manually at this problem point.
c
c First determine if vectors (x-xc) and (x0-xc) are
c collinear
c----

          xh1 = x-x0
          xh2 = y-y0
          xh3 = z-z0

          xb1 = x-xc
          xb2 = y-yc
          xb3 = z-zc

          curl1 = xh2*xb3 - xb2*xh3
          curl2 = xh3*xb1 - xb3*xh1
          curl3 = xh1*xb3 - xb1*xh3

          check = dsqrt(curl1**2 + curl2**2 + curl3**2)

c          print*,curl1,curl2,curl3

          if (check.lt.1.d-6) then
c---
c If so, adjust manually
c---
             gxx = 0.d0
             gxy = 0.d0
             gxz = 0.d0
             gyx = 0.d0
             gyy = 0.d0
             gyz = 0.d0
             gzx = 0.d0
             gzy = 0.d0
             gzz = 0.d0
             goto 451
          end if

c------------------
c End of MGB's edits
c------------------

          write (6,*)
          write (6,*) ' x, y, z coordinates of the point are:'
          write (6,100) x,y,z
          write (6,*)

        Else

          write (6,*) 
          write (6,*) ' Enter the (x,y,z) of the velocity point'
          write (6,*) ' --------------------------------------'
          read  (5,*) x,y,z

        End If
c-------------

        call sgf_3d_sph
     +
     +    (x,y,z
     +    ,x0,y0,z0
     +    ,xc,yc,zc
     +    ,a
     +    ,gxx,gxy,gxz
     +    ,gyx,gyy,gyz
     +    ,gzx,gzy,gzz
     +    )

 451    write (6,*) ' --------------'
        write (6,*) ' Green function'
        write (6,*) ' --------------'

        write (6,100) gxx,gxy,gxz
        write (6,100) gyx,gyy,gyz
        write (6,100) gzx,gzy,gzz

      End If

c--------------------------------
c test of the continuity equation
c--------------------------------

      If(menu.eq.2) then

        write (6,*)
        write (6,*) ' Enter the center of the test sphere'
        write (6,*)
        write (6,*) ' 99 for one coordinate to quit'
        write (6,*) ' -----------------------------'
        read  (5,*) xcnt,ycnt,zcnt


        If(xcnt.eq.99) Go to 99
        If(ycnt.eq.99) Go to 99
        If(zcnt.eq.99) Go to 99

        write (6,*)
        write (6,*) ' Enter the radius of the test sphere'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) rad

        If(rad.eq.0) Go to 99

        write (6,*)
        write (6,*) ' Choose the integration formula'
        write (6,*) ' select from 6, 18'
        write (6,*) ' -----------------'
        read  (5,*) mint

        call gauss_sph (mint,xx,yy,zz,ww)

        sum1 = 0.0D0
        sum2 = 0.0D0
        sum3 = 0.0D0

        Do i=1,mint

         x = xcnt+rad*xx(i)
         y = ycnt+rad*yy(i)
         z = zcnt+rad*zz(i)

         call sgf_3d_sph
     +
     +        (x,y,z
     +        ,x0,y0,z0
     +        ,xc,yc,zc,a
     +        ,gxx,gxy,gxz
     +        ,gyx,gyy,gyz
     +        ,gzx,gzy,gzz
     +        )

         vnx = xx(i)
         vny = yy(i)
         vnz = zz(i)

         sum1 = sum1 + (gxx*vnx+gyx*vny+gzx*vnz)*ww(i)
         sum2 = sum2 + (gxy*vnx+gyy*vny+gzy*vnz)*ww(i)
         sum3 = sum3 + (gxz*vnx+gyz*vny+gzz*vnz)*ww(i)

c        sum1 = sum1 + (gxx*vnx+gxy*vny+gxz*vnz)*ww(i)
c        sum2 = sum2 + (gyx*vnx+gyy*vny+gyz*vny)*ww(i)
c        sum3 = sum3 + (gzx*vnx+gzy*vny+gzz*vny)*ww(i)

        End Do
        write (6,100) sum1,sum2,sum3

      End If

c---
      Go to 98
c---
 
  99  Continue

c-----
c done
c-----

 100  Format (3(2x,f15.10))
 110  format(1x,A,20f12.6)

      Stop
      End
