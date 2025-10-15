      program sgf_2d_sph_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c----------------------------------------
c Driver for Stokes Green's function
c outside a apshere
c
c SYMBOLS:
c -------
c
c (xc,yc,zc) Center of the sphere
c a          Radius of the sphere
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
      write (6,*) 'Enter:'
      write (6,*)
      write (6,*) '1 for one evaluation'
      write (6,*) '2 for a test of the continuity equation'
      write (6,*) '0 to quit'
      write (6,*) '---------'

      read  (5,*) menu
      If(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) 'Enter the coordinates of the center'
      write (6,*) 'of the sphere: xc, yc,zc'
      write (6,*) '------------------------'
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
      write (6,*) 'Enter the coordinates of the point force'
      write (6,*) 'x0, y0, z0'
      write (6,*) '----------'
      read  (5,*) x0,y0,z0

c---
c make sure the point force
c is not in the sphere
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

          write (6,*) 'Enter:'
          write (6,*)
          write (6,*) '1 for a point on the sphere'
          write (6,*) '2 for a point off the sphere'
          write (6,*) '0 to quit'
          write (6,*) '---------'
          read  (5,*) ion

        If(ion.eq.0) Go to 99

c-------------
        If(ion.eq.1) then

          write (6,*)
          write (6,*) ' Enter the angles: theta,  phi'
          write (6,*) ' in multiples of pi'
          write (6,*) ' -------------------------------'
          read  (5,*) theta,phi

          theta = theta*pi
          phi   = phi  *pi

          x = xc + a*cos(theta)*cos(phi)
          y = yc + a*sin(theta)*cos(phi)
          z = zc + a*sin(phi)

          write (6,*)
          write (6,*) 'The x, y, z coordinates of the point force are:'
          write (6,100) x,y,z
          write (6,*)

        Else

          write (6,*) 
          write (6,*) ' Enter (x,y,z) of the velocity point'
          write (6,*) ' -----------------------------------'
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

      write (6,*) ' --------------'
      write (6,*) ' Green function'
      write (6,*) ' --------------'

      write (6,100) gxx,gxy,gxz
      write (6,100) gyx,gyy,gyz
      write (6,100) gzx,gzy,gzz

      End If

c--------------------------------
c test of the continuity equation
c--------------------------------

      if(menu.eq.2) then

        write (6,*)
        write (6,*) ' Enter the center of the test sphere'
        write (6,*) ' 99 for any coordinate to quit'
        write (6,*) ' -----------------------------'
        read  (5,*) xcnt,ycnt,zcnt

        if(xcnt.eq.99) Go to 99
        if(ycnt.eq.99) Go to 99
        if(zcnt.eq.99) Go to 99

        write (6,*)
        write (6,*) ' Enter the radius of the test sphere'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) rad

        if(rad.eq.0) Go to 99

        write (6,*)
        write (6,*) ' Choose the integration formula'
        write (6,*) ' Select from 6, 18'
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
     +      (x,y,z
     +      ,x0,y0,z0
     +      ,xc,yc,zc,a
     +      ,gxx,gxy,gxz
     +      ,gyx,gyy,gyz
     +      ,gzx,gzy,gzz
     +      )

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
 
  99  continue

c-----
c Done
c-----

 100  Format (3(2x,f15.10))

      stop
      end
