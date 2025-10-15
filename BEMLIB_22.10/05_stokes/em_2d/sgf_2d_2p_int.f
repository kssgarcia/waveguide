      subroutine sgf_2d_2p_int
     +
     +   (x,y
     +   ,x0,y0
     +   ,a21
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------
c trilinear interpolation of Green's function 
c from look-up tables.
c
c Interpolation is done with respect to: 
c
c a21, x-x0, y-y0
c                               
c for a11 = 1.0, a12 = 0.0, a22 = 1.0
c
c SYMBOLS:
c -------
c
c (x,y): coordinates of the field point
c (x0,y0): coordinates of a point force
c
c NOTE:
c ----
c
c The f77 intrinsic function nint(a)
c produces the nearest integer.
c              -------
c
c Thus, int(9.7) = 9 but nint(9.7) = 10
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension glut_xx(-64:64,0:64,-64:64)
      Dimension glut_xy(-64:64,0:64,-64:64)
      Dimension glut_yy(-64:64,0:64,-64:64)

      common/glut_r/glut_xx,glut_xy,glut_yy
      common/glut_i/Na21,Nxxx,Nyyy

c---
c set the fixed components of the base vectors
c---

      a11 = 1.0D0
      a12 = 0.0D0
      a22 = 1.0D0

c---
c compute the interpolation variables: xx and yy
c---

      xx = x-x0
      yy = y-y0

c---------------------------------------
c bring a21 into the range [-0.50, 0.50]
c---------------------------------------

      a21 = a21 - nint(a21)*a11

c---------------------------------------
c bring the interpolation variables xx and yy
c into the master cell:
c
c -0.50<xx<0.50      0.0<yy,0.50
c---------------------------------------

      Imove = nint(yy)
      xx = xx - Imove*a21   ! move along the second base vector
      yy = yy - Imove*a22

      Imove = nint(xx)
      xx = xx - Imove*a11   ! move along the first base vector
      yy = yy - Imove*a12

c---------------------------------------
c If yy < 0, reflect the pair (xx,yy) 
c            with respect to the origin
c            to bring it within the master cell
c---------------------------------------

      If(yy.lt.0) then
        xx = -xx
        yy = -yy
      End If

c---
c checks
c
c De-activate after testing
c if desired
c---

      if(a21.gt.0.5.or.a21.lt.-0.5) then
        write (6,*)
        write (6,*) " Message from sgf_2d_sp_int"
        write (6,*)
        write (6,*) a21," a21 exceeded the designated range"
        stop
      end if

      if(xx.gt.0.5.or.xx.lt.-0.5) then
        write (6,*)
        write (6,*) " Message from sgf_2d_sp_int"
        write (6,*)
        write (6,*) xx," xx exceeded the designated range"
        stop
      end if

      if(yy.gt.0.5.or.yy.lt.0.0) then
        write (6,*)
        write (6,*) " Message from sgf_2d_sp_int"
        write (6,*)
        write (6,*) yy," yy exceeded the designated range"
        stop
      end if

c---
c compute the cell indices
c---

      atN = 2.0D0*Na21*a21  ! remember: -0.5<a21<0.5
      xtN = 2.0D0*Nxxx*xx   ! remember: -0.5< xx<0.5
      ytN = 2.0D0*Nyyy*yy   ! remember:  0.0< yy<0.5

      kct = int(atN)
      ict = int(xtN)
      jct = int(ytN)

      if(atn.gt.0) then
       kct1 = kct + 1
      else
       kct1 = kct - 1
      end if

      if(xtn.gt.0) then
       ict1 = ict + 1
      else
       ict1 = ict - 1
      end if

      jct1 = jct + 1

c---
c compute the interpolation coefficients
c---

      rct = abs(atN-kct)
      pct = abs(xtN-ict)
      qct =     ytN-jct

      c00 = (1.0D0-rct)*(1.0D0-qct)
      c10 = (1.0D0-rct)*qct
      c01 = rct*(1.0D0-qct)
      c11 = rct*qct

c--
c trilinear interpolation
c---

      Gxx =  c00*((1.0D0-pct)*Glut_xx(ict, jct, kct )
     +                  +pct *Glut_xx(ict1,jct, kct ))
     +     + c10*((1.0D0-pct)*Glut_xx(ict, jct1,kct )
     +                  +pct *Glut_xx(ict1,jct1,kct ))
     +     + c01*((1.0D0-pct)*Glut_xx(ict, jct, kct1)
     +                  +pct *Glut_xx(ict1,jct, kct1))
     +     + c11*((1.0D0-pct)*Glut_xx(ict, jct1,kct1)
     +                  +pct *Glut_xx(ict1,jct1,kct1))

      Gxy =  c00*((1.0D0-pct)*Glut_xy(ict, jct, kct )
     +                  +pct *Glut_xy(ict1,jct, kct ))
     +     + c10*((1.0D0-pct)*Glut_xy(ict, jct1,kct )
     +                  +pct *Glut_xy(ict1,jct1,kct ))
     +     + c01*((1.0D0-pct)*Glut_xy(ict, jct, kct1)
     +                  +pct *Glut_xy(ict1,jct, kct1))
     +     + c11*((1.0D0-pct)*Glut_xy(ict, jct1,kct1)
     +                  +pct *Glut_xy(ict1,jct1,kct1))

      Gyy =  c00*((1.0D0-pct)*Glut_yy(ict, jct, kct )
     +                  +pct *Glut_yy(ict1,jct, kct ))
     +     + c10*((1.0D0-pct)*Glut_yy(ict, jct1,kct )
     +                  +pct *Glut_yy(ict1,jct1,kct ))
     +     + c01*((1.0D0-pct)*Glut_yy(ict, jct, kct1)
     +                  +pct *Glut_yy(ict1,jct, kct1))
     +     + c11*((1.0D0-pct)*Glut_yy(ict, jct1,kct1)
     +                  +pct *Glut_yy(ict1,jct1,kct1))

c---
c Add the free space Green's function
c---

      rs  = xx**2 + yy**2
      tmp = 0.5D0*log(rs)

      Gxx = Gxx - tmp + xx*xx/rs
      Gxy = Gxy       + xx*yy/rs
      Gyy = Gyy - tmp + yy*yy/rs
      Gyx = Gxy

c-----
c done
c-----

 100  Format (10(1x,f15.10))

      return
      end
