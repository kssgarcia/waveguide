      function erfuni(y)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c-----------------------
c inverse error function
c-----------------------

      Implicit Double Precision (a-h,o-z)
 
c-------------------------
c constants and parameters
c-------------------------

      maxit = 20
      eps = 0.0000001 D0
      pi = 3.14159265358 D0

      const = Dsqrt(pi)/2.0D0

c-----------------
c very small value
c-----------------

      if(abs(y).le.eps) then
        erfi = Const*y
        iter = 1
        return
      end if

c----------------------------
c Newton iterations otherwise
c----------------------------

      if(abs(y).lt.1.0) then

        erfi  = const*abs(y)
        y0    = erfun(0.9D0*erfi)
        Derfi = 0.1D0*erfi

        Do iter=1,maxit

         y1  = erfun(erfi)
         Dy1 = Dabs(y)-y1
         If (Dabs(Dy1).lt.eps) Go to 99
         Dy0   = y1-y0
         Derfi = Derfi*Dy1/Dy0
         y0    = y1
         erfi  = erfi + Derfi
         if(abs(Derfi/erfi).lt.eps) Go to 99

        End Do

      end if

c-----------
c unfeasible
c-----------

      Iter = 0     ! Did not converge
      erfi = 0.0D0

      Return

  99  Continue

      if(y.lt.0.0) erfuni = -erfuni

c-----
c Done
c-----

      Return
      End
