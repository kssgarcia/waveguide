      subroutine ber_bei_0 
     +
     +  (Iopt
     +  ,X
     +  ,ber_0,bei_0
     +  ,ber_0_p,bei_0_p
     +  )

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c========================================
c Evaluate the Kelvin functions:
c
c  ber_0, bei_0
c
c  and their first derivatives:
c
c  ber_0_p, bei_0_p
c
c  If Iopt.ne.1 derivatives are not computed
c
c A&S: Abramowitz and Stegun
c========================================

      Implicit Double Precision (a-h,o-z)

c-----------------------
      if(X.le.8.00) then
c-----------------------

c---
c Use formulas (9.11.1) and (9.11.2) of A&S, p. 384
c---

      XH  = 0.5D0*X
      Y   = X/8.0D0

      Y2  = Y**2
      Y4  = Y2**2
      Y6  = Y4 *Y2
      Y8  = Y6 *Y2
      Y10 = Y8 *Y2
      Y12 = Y10*Y2
      Y14 = Y12*Y2
      Y16 = Y14*Y2
      Y18 = Y16*Y2
      Y20 = Y18*Y2
      Y22 = Y20*Y2
      Y24 = Y22*Y2
      Y26 = Y24*Y2
      Y28 = Y26*Y2

      ber_0 =  1.0
     +
     +   - 64.0        * Y4  + 113.77777774 * Y8
     +   - 32.36345652 * Y12 +   2.64191397 * Y16
     +   -  0.08349609 * Y20 +   0.00122552 * Y24
     +   -  0.00000901 * Y28

      bei_0 =
     +
     +  16.0        * Y2  - 113.77777774 * Y6
     +  + 72.81777742 * Y10 -  10.56765779 * Y14
     +  +  0.52185615 * Y18 -   0.01103667 * Y22
     +  +  0.00011346 * Y26

c---
c  Use formulas (9.11.5) and (9.11.6) of A&S, p. 384
c---

      if(iopt.ne.0) then

      ber_0_p = X*( - 4.0        * Y2  + 14.22222222 * Y6
     +             - 6.06814810 * Y10 +  0.66047849 * Y14
     +             - 0.02609253 * Y18 +  0.00045957 * Y22
     +             - 0.00000394 * Y26 )

      bei_0_p = X*(   0.5              - 10.66666666 * Y4
     +             + 11.37777772 * Y8  -  2.31167514 * Y12
     +             +  0.14677204 * Y16 -  0.00379386 * Y20
     +             +  0.00004609 * Y24 )
      end if

c---------
      else
c---------

      write (6,*)
      write (6,*) " ber_bei_0: functions not yet implemented"
      write (6,*)

      stop

c---------
      End If
c---------

c-----
c done
c-----

 100  Format (5(1x,f12.8))

      return
      end
