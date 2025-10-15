      subroutine ker_kei_0
     +
     +   (Iopt
     +   ,X
     +   ,KER0, KEI0
     +   ,KER0D,KEI0D  ! derivatives
     +   )

c--------------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c--------------------------------------------

c--------------------------------------------
c Computation of Kelvin functions ker0, kei0,
c and their derivatives ker0d, ker0i
c
c Set Iopt = 0 to compute ker, kei
c     Iopt = 1 to compute ker, kei, kerd, keid
c
c Note that ker and kei are defined in terms of
c ber and bei
c--------------------------------------------

      Implicit Double Precision (a-h, o-z)

      Double Precision KER0,KEI0,KER0D,KEI0D

      Parameter (eps=0.00001)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      piq = 0.25D0*pi
      pih = 0.50D0*pi

c--------
      if(X.le.8.00D0) then
c--------

      XH = 0.5D0*X
      Y  = X/8.0D0

      Y2  =   Y*Y
      Y4  =  Y2*Y2
      Y6  =  Y4*Y2
      Y8  =  Y6*Y2
      Y10 =  Y8*Y2
      Y12 = Y10*Y2
      Y14 = Y12*Y2
      Y16 = Y14*Y2
      Y18 = Y16*Y2
      Y20 = Y18*Y2
      Y22 = Y20*Y2
      Y24 = Y22*Y2
      Y26 = Y24*Y2
      Y28 = Y26*Y2

      BER0 = 1.0 - 64.0       * Y4  + 113.77777774 * Y8
     +          - 32.36345652 * Y12 +   2.64191397 * Y16
     +          -  0.08349609 * Y20 +   0.00122552 * Y24
     +          -  0.00000901 * Y28

      BEI0 =      16.0        * Y2  - 113.77777774 * Y6
     +          + 72.81777742 * Y10 -  10.56765779 * Y14
     +          +  0.52185615 * Y18 -   0.01103667 * Y22
     +          +  0.00011346 * Y26

      KER0 = - DLOG(XH) * BER0 + piq*BEI0 - 0.57721566
     +          - 59.05819744 * Y4  + 171.36272133 * Y8
     +          - 60.60977451 * Y12 +   5.65539121 * Y16
     +          -  0.19636347 * Y20 +   0.00309699 * Y24
     +          -  0.00002458 * Y28

      KEI0 = - DLOG(XH) * BEI0 - piq*BER0 + 6.76454936 * Y2
     +          - 142.91827687 * Y6  + 124.23569650 * Y10
     +          -  21.30060904 * Y14 +   1.17509064 * Y18
     +          -   0.02695875 * Y22 +   0.00029532 * Y26

c---- compute the derivatives:

      If(Iopt.ne.0) then

      BER0D = X*( - 4.0        * Y2  + 14.22222222 * Y6
     +           - 6.06814810 * Y10 +  0.66047849 * Y14
     +           - 0.02609253 * Y18 +  0.00045957 * Y22
     +           - 0.00000394 * Y26 )

      BEI0D = X*(   0.5              - 10.66666666 * Y4
     +          + 11.37777772 * Y8  -  2.31167514 * Y12
     +          +  0.14677204 * Y16 -  0.00379386 * Y20
     +          +  0.00004609 * Y24 )

      KER0D = - DLOG(XH) * BER0D - BER0/X + piq*BEI0D
     +      + X* ( -  3.69113734 * Y2  + 21.42034017 * Y6
     +             - 11.36433272 * Y10 +  1.41384780 * Y14
     +             -  0.06136358 * Y18 +  0.00116137 * Y22
     +             -  0.00001075 * Y26 )

      KEI0D = - DLOG(XH) * BEI0D - BEI0/X - piq*BER0D
     +      + X* (    0.21139217       - 13.39858846 * Y4
     +             + 19.41182758 * Y8  -  4.65950823 * Y12
     +             +  0.33049424 * Y16 -  0.00926707 * Y20
     +             +  0.00011997 * Y24 )

      End If

c---------
      else
c---------

      call ker_kei_0_large (X,KER0,KEI0)
      
      eps2 = 2.0*eps

      ! compute the derivatives
      ! by finite differences:

      if(Iopt.ne.0) then

       X1 = X + eps
       call ker_kei_0_large (X1,WER1,WEI1)
       X2 = X - eps
       call ker_kei_0_large (X2,WER2,WEI2)

       KER0D = (WER1-WER2)/eps2
       KEI0D = (WEI1-WEI2)/eps2

       end if

c-------
      end if
c-------

c=====
c done
c=====

 100  Format (5(1x,f12.8))

      return
      end
