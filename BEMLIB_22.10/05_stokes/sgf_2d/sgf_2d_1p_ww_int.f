      subroutine sgf_2d_1p_ww_int
     +
     +  (Iselect
     +  ,x,y
     +  ,x0,y0
     +  ,RL
     +  ,h
     +  ,Gxx,Gxy
     +  ,Gyx,Gyy
     +  ,px,py
     +  ,Txxx,Txxy,Tyxx,Tyxy
     +  ,Txyx,Txyy,Tyyx,Tyyy
     +  )

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c Trilinear interpolation of Green's function 
c from look-up tables.
c
c Interpolation is done with respect to: 
c
c y0, y, x-x0
c                               
c SYMBOLS:
c -------
c
c (x,  y) : coordinates of the field point
c (x0, y0): coordinates of a point force
c
c NOTE:
c ----
c
c The f77 intrinsic function nint(a)
c produces the nearest integer.
c              -------
c
c Thus, int(9.7) = 9 but nint(9.7) = 10
c
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension glutww_xx(0:64,-64:64,-64:64)
      Dimension glutww_xy(0:64,-64:64,-64:64)
      Dimension glutww_yx(0:64,-64:64,-64:64)
      Dimension glutww_yy(0:64,-64:64,-64:64)

      common/glutww_r/glutww_xx,glutww_xy,glutww_yx,glutww_yy
      common/glutww_i/Nwwy0,Nwwy,Nwwx

c--------
c prepare
c--------

      RLH = 0.5D0*RL

c------------------------------
c bring the interpolation point
c to the right of x0
c------------------------------

      xx = x-x0

      Do i=1,100
       If(xx.gt.0) Go to 1
       xx = xx+RL
      End Do

  1   Continue

c---------------------------------------
c compute the interpolation coefficients
c---------------------------------------

      xqr  = mod(xx+RLH,RL) - RLH
      sgn  = sign(1.0D0,xqr)
      xqr  = sgn*xqr
      xtN  = Nwwx*xqr/RLH
      ict  = int(xtN)
      ict1 = ict+1

      yqr  = y
      ytN  = Nwwy*yqr/h
      jct  = int(ytN)
      incr = nint(sign(1.0D0,ytN))
      jct1 = jct+incr

      zqr  = y0
      ztN  = Nwwy0*zqr/h
      kct  = int(ztN)
      incr = nint(sign(1.0D0,ztN))
      kct1 = kct+incr

      pct  = xtN - ict
      pctc = 1.0D0-pct
      qct  = abs(ytN-jct)
      rct  = abs(ztN-kct)

      c00 = (1.0D0-rct)*(1.0D0-qct)
      c10 = (1.0D0-rct)*qct
      c01 = rct*(1.0D0-qct)
      c11 = rct*qct

c     write (6,*) ict,jct,kct

c------------------------
c trilinear interpolation
c------------------------

      Gxx =  c00*(pctc*Glutww_xx(ict, jct, kct )
     +           +pct *Glutww_xx(ict1,jct, kct ))
     +     + c10*(pctc*Glutww_xx(ict, jct1,kct )
     +           +pct *Glutww_xx(ict1,jct1,kct ))
     +     + c01*(pctc*Glutww_xx(ict, jct, kct1)
     +           +pct *Glutww_xx(ict1,jct, kct1))
     +     + c11*(pctc*Glutww_xx(ict, jct1,kct1)
     +           +pct *Glutww_xx(ict1,jct1,kct1))

      Gxy =  c00*(pctc*Glutww_xy(ict, jct, kct )
     +           +pct *Glutww_xy(ict1,jct, kct ))
     +     + c10*(pctc*Glutww_xy(ict, jct1,kct )
     +           +pct *Glutww_xy(ict1,jct1,kct ))
     +     + c01*(pctc*Glutww_xy(ict, jct, kct1)
     +           +pct *Glutww_xy(ict1,jct, kct1))
     +     + c11*(pctc*Glutww_xy(ict, jct1,kct1)
     +           +pct *Glutww_xy(ict1,jct1,kct1))

      Gyx =  c00*(pctc*Glutww_yx(ict, jct, kct )
     +           +pct *Glutww_yx(ict1,jct, kct ))
     +     + c10*(pctc*Glutww_yx(ict, jct1,kct )
     +           +pct *Glutww_yx(ict1,jct1,kct ))
     +     + c01*(pctc*Glutww_yx(ict, jct, kct1)
     +           +pct *Glutww_yx(ict1,jct, kct1))
     +     + c11*(pctc*Glutww_yx(ict, jct1,kct1)
     +           +pct *Glutww_yx(ict1,jct1,kct1))

      Gyy =  c00*(pctc*Glutww_yy(ict, jct, kct )
     +           +pct *Glutww_yy(ict1,jct, kct ))
     +     + c10*(pctc*Glutww_yy(ict, jct1,kct )
     +           +pct *Glutww_yy(ict1,jct1,kct ))
     +     + c01*(pctc*Glutww_yy(ict, jct, kct1)
     +           +pct *Glutww_yy(ict1,jct, kct1))
     +     + c11*(pctc*Glutww_yy(ict, jct1,kct1)
     +           +pct *Glutww_yy(ict1,jct1,kct1))

      Gxy = sgn * Gxy
      Gyx = sgn * Gyx

c---
c Add the complementary components
c---

      Iss = 0

      call sgf_2d_1p_w
     +
     +       (Iselect
     +       ,x,y
     +       ,x0,y0
     +       ,-h          ! location of the lower wall
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
     +       ,h        ! location of the upper wall
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

      yy  = yqr-zqr
      rs  = xqr**2 + yy**2
      rlog = 0.5D0*log(rs)

c     write (6,*) xqr,yy

      Gxx3 = - rlog + xqr**2/rs
      Gxy3 =        + xqr*yy/rs * sgn
      Gyx3 =        + xqr*yy/rs * sgn
      Gyy3 = - rlog + yy**2 /rs

c     Gxx3 = 0.0D0
c     Gxy3 = 0.0D0
c     Gyx3 = 0.0D0
c     Gyy3 = 0.0D0

      Gxx = Gxx + Gxx1 + Gxx2 - Gxx3
      Gxy = Gxy + Gxy1 + Gxy2 - Gxy3
      Gyx = Gyx + Gyx1 + Gyx2 - Gyx3
      Gyy = Gyy + Gyy1 + Gyy2 - Gyy3

c-----
c Done
c-----

 100  Format (10(1x,f15.10))

      Return
      End
