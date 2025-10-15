      subroutine lgf_3d_2p_int
     +
     +    (x,y,z
     +    ,x0,y0,z0
     +    ,a11,a12
     +    ,a21,a22
     +    ,g
     +    ,gx,gy,gz
     +    )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c all rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------
c Trilinear interpolation of the doubly-periodic
c Green's function  of the three-dimensional Laplace
c equation from look-up tables. 
c 
c
c Symbols:
c --------
c
c x,y,z    : coordinates of the evaluation or field point
c x0,y0,z0 : coordinates of a point force
c 
c glut     : green's function look-up table
c glut_x   : x-component of the gradient of the green's function
c glut_y   : y-component of the gradient of the green's function
c glut_z   : z-component of the gradient of the green's function
c
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension glut  (-64:64,-64:64,0:64)
      Dimension glut_x(-64:64,-64:64,0:64)
      Dimension glut_y(-64:64,-64:64,0:64)
      Dimension glut_z(-64:64,-64:64,0:64)

      common/glut_r/glut,glut_x,glut_y,glut_z
      common/glut_i/nxxx,nyyy,nzzz

c----------
c constants
c----------

      pi  = 3.14159 265358
      pi4 = 4.0*pi

c--------
c prepare
c--------

      xx = x-x0
      yy = y-y0
      zz = z-z0

c--------------------------------------------
c bring the interpolation variables xx and yy
c into to the master cell:
c
c -0.50<xx<0.50      -0.50<yy<0.50
c--------------------------------------------

      imove = nint(yy)
      xx = xx - imove*a21   ! move along the second base vector
      yy = yy - imove*a22

      imove = nint(xx)
      xx = xx - imove*a11   ! move along the first base vector
      yy = yy - imove*a12

c---
c checks
c
c de-activate after testing
c if desired
c---

      If(zz.gt.0.5.or.zz.lt.0.0) then
        write (6,*)
        write (6,*) " lgf_3d_2p_int:"
        write (6,*) zz," zz exceeded the designated range"
        stop
      End If

      If(xx.gt.0.5 .or.xx.lt.-0.5) then
        write (6,*)
        write (6,*) " lgf_3d_2p_int:"
        write (6,*) xx," xx exceeded the designated range"
        stop
      End If

      If(yy.gt.0.5.or.yy.lt.-0.5) then
        write (6,*)
        write (6,*) " message from lgf_3d_2p_int"
        write (6,*)
        write (6,*) yy," yy exceeded the designated range"
        stop
      End If

c---
c compute the cell indices
c---

      xtn = 2.0*nxxx*xx   ! remember: -0.5<xx<0.5
      ytn = 2.0*nyyy*yy   ! remember: -0.5<yy<0.5
      ztn = 2.0*nzzz*zz   ! remember:  0.0<zz<0.5

      ict = int(xtn)
      jct = int(ytn)
      kct = int(ztn)

      if(xtn.gt.0) then
       ict1 = ict + 1
      else
       ict1 = ict - 1
      end if

      if(ytn.gt.0) then
       jct1 = jct + 1
      else
       jct1 = jct - 1
      end if

      kct1 = kct + 1

c---
c compute the interpolation coefficients
c---

      pct = abs(xtn-ict)
      qct = abs(ytn-jct)
      rct =    (ztn-kct)

      c00 = (1.0-rct)*(1.0-qct)
      c10 = (1.0-rct)*qct
      c01 = rct*(1.0-qct)
      c11 = rct*qct

c--
c trilinear interpolation
c---

      g   =  c00*((1.0-pct)*glut  (ict, jct, kct )
     +                +pct *glut  (ict1,jct, kct ))
     +     + c10*((1.0-pct)*glut  (ict, jct1,kct )
     +                +pct *glut  (ict1,jct1,kct ))
     +     + c01*((1.0-pct)*glut  (ict, jct, kct1)
     +                +pct *glut  (ict1,jct, kct1))
     +     + c11*((1.0-pct)*glut  (ict, jct1,kct1)
     +                +pct *glut  (ict1,jct1,kct1))

      gx  =  c00*((1.0-pct)*glut_x(ict, jct, kct )
     +                +pct *glut_x(ict1,jct, kct ))
     +     + c10*((1.0-pct)*glut_x(ict, jct1,kct )
     +                +pct *glut_x(ict1,jct1,kct ))
     +     + c01*((1.0-pct)*glut_x(ict, jct, kct1)
     +                +pct *glut_x(ict1,jct, kct1))
     +     + c11*((1.0-pct)*glut_x(ict, jct1,kct1)
     +                +pct *glut_x(ict1,jct1,kct1))

      gy  =  c00*((1.0-pct)*glut_y(ict, jct, kct )
     +                +pct *glut_y(ict1,jct, kct ))
     +     + c10*((1.0-pct)*glut_y(ict, jct1,kct )
     +                +pct *glut_y(ict1,jct1,kct ))
     +     + c01*((1.0-pct)*glut_y(ict, jct, kct1)
     +                +pct *glut_y(ict1,jct, kct1))
     +     + c11*((1.0-pct)*glut_y(ict, jct1,kct1)
     +                +pct *glut_y(ict1,jct1,kct1))

      gz  =  c00*((1.0-pct)*glut_z(ict, jct, kct )
     +                +pct *glut_z(ict1,jct, kct ))
     +     + c10*((1.0-pct)*glut_z(ict, jct1,kct )
     +                +pct *glut_z(ict1,jct1,kct ))
     +     + c01*((1.0-pct)*glut_z(ict, jct, kct1)
     +                +pct *glut_z(ict1,jct, kct1))
     +     + c11*((1.0-pct)*glut_z(ict, jct1,kct1)
     +                +pct *glut_z(ict1,jct1,kct1))

c-----
c done
c-----

 100  Format (10(1x,f15.10))

      return
      end 
