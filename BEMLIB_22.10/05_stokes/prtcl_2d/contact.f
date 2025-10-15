      subroutine contact
     +
     +   (a11,a21,xc1,yc1,th1
     +   ,a12,a22,xc2,yc2,th2
     +   ,Tmin,xmin,ymin
     +   ,Iloc
     +   )

c-------------------------------------
c Find the point on ellipse numbered 2
c closest to another ellipse numbered 1
c
c SYMBOLS
c ------
c
c first ellipse:
c
c a11:       first  axis
c a21:       second axis
c xc1, yc1: coordinates of the center
c th1:       inclination
c
c second ellipse:
c
c a12:        first  axis 
c a22:        second axis
c xc2, yc2:  coordinates of the center
c th2:        inclination of the second ellipse
c
c xmin, ymin: coordinates of the minimum
c
c If Iloc=0, the point of minimum separation is in the exterior
c            of the first ellipse
c
c If Iloc=1, the point of minimum separation is in the interior
c            of the first ellipse
c
c---------------------------------------
 
      Implicit Double Precision (a-h,o-z)

      Parameter (Ndiv=48,Nter=20,eps=0.0001)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

c-----------
c Initialize
c-----------

      Iloc = 0

c--------
c prepare
c--------

      a11s = a11*a11
      a21s = a21*a21

      a12s = a12*a12
      a22s = a22*a22

      cs1 = dcos(th1)
      sn1 = dsin(th1)

      cs2 = dcos(th2)
      sn2 = dsin(th2)

c---------------------------
c search for a rough minimum
c of the shape function
c by comparison
c---------------------------

      dth = pi2/Ndiv

      Dmin = 100000000.0 
      Tmin = 0.0

      Do i=1,Ndiv

       th = (i-1.0D0)*dth

       cs = cos(th)
       sn = sin(th)

       xx = xc2 + a12*cs*cs2 - a22*sn*sn2
       yy = yc2 + a12*cs*sn2 + a22*sn*cs2

       fnc = ( (xx-xc1)*cs1+(yy-yc1)*sn1 )**2/a11s
     +     + (-(xx-xc1)*sn1+(yy-yc1)*cs1 )**2/a21s
     +     - 1.0D0

c      write (6,100) i,th,fnc

       If(fnc.lt.Dmin) then
        Dmin = fnc
        Tmin = th
       End If

      End Do

      Tmin_rough = Tmin

c     write (6,*)
c     write (6,*)  " Rough minimum:",Tmin,Dmin
c     write (6,*)

c     Go to 2

c---------------------------
c Compute the true minimum
c by Newton's method
c---------------------------

      Iter = 1

  1   Continue

      call fnc_min
     +
     +   (a11,a21,xc1,yc1,th1,a11s,a21s,cs1,sn1
     +   ,a12,a22,xc2,yc2,th2,a12s,a22s,cs2,sn2
     +   ,Tmin,fnc0
     +   )

c     write (6,100) Iter,Tmin,fnc

      Tmin_save = Tmin

      Tmin = Tmin+eps

      call fnc_min
     +
     +   (a11,a21,xc1,yc1,th1,a11s,a21s,cs1,sn1
     +   ,a12,a22,xc2,yc2,th2,a12s,a22s,cs2,sn2
     +   ,Tmin,fnc1
     +   )

      corr = -eps*fnc0/(fnc1-fnc0)

      Tmin = Tmin_save + corr

      If(Iter.eq.Nter) then
        write (6,*)
        write (6,*) " contact: number of Newton iterations: ",Nter
        write (6,*) "          exceeded"
        write (6,*)
        write (6,100) Iter,Tmin,fnc
        write (6,*)
        Tmin = Tmin_rough
        Go to 2
      End If

      If(abs(corr).gt.0.0000001) then
        Iter = Iter+1
        Go to 1
      End If

  2   Continue

c---------------------------------------------
c compute the coordinates of the minimum point
c---------------------------------------------

       cs = cos(Tmin)
       sn = sin(Tmin)

       xmin = xc2 + a12*cs*cs2 - a22*sn*sn2
       ymin = yc2 + a12*cs*sn2 + a22*sn*cs2

       fnc = ( (xmin-xc1)*cs1+(ymin-yc1)*sn1 )**2/a11s
     +     + (-(xmin-xc1)*sn1+(ymin-yc1)*cs1 )**2/a21s
     +     - 1.0D0

       If(fnc.lt.0) Iloc = 1

c-----
c Done
c-----

 100  Format (1x,i4,2(1x,f10.5))
 101  Format (2(1x,f10.5))

      Return
      End

c===============================================

      subroutine fnc_min 
     +
     +   (a11,a21,xc1,yc1,th1,a11s,a21s,cs1,sn1
     +   ,a12,a22,xc2,yc2,th2,a12s,a22s,cs2,sn2
     +   ,t,fnc
     +   )

c---------------------------------
c Evaluate the first derivative of
c the location function
c---------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c begin
c---
      cs = cos(t)
      sn = sin(t)

      xx = xc2 + a12*cs*cs2 - a22*sn*sn2
      yy = yc2 + a12*cs*sn2 + a22*sn*cs2

      dxdt = -a12*sn*cs2 - a22*cs*sn2
      dydt = -a12*sn*sn2 + a22*cs*cs2 

      fnc = ( (xx-xc1)*cs1+(yy-yc1)*sn1 )/a11s
     +      * (dxdt*cs1+dydt*sn1)
     +    + (-(xx-xc1)*sn1+(yy-yc1)*cs1 )/a21s
     +      * (-dxdt*sn1+dydt*cs1)

c-----
c Done
c-----

      Return
      End
