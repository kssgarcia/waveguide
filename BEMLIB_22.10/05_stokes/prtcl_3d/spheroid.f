      program spheroid

c----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c----------------------------------------
c Computation of the force resistance tensor
c for the translation of a prolate spheroid
c
c Forulas on page 268 of Pozrikidis (1997)
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c---------
c constants
c---------

      pi = 3.14157 265358 D0
      oot = 1.0D0/3.0D0

c----------
c launching
c----------

  97  Continue

      write (6,*)
      write (6,*) " The particle is a prolate spheroid"
      write (6,*) " with axes: a, b where b<a"
      write (6,*)
      write (6,*) " Please enter the axes ratio b/a"
      write (6,*) " (less than 1)"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) boa

      if(boa.eq.0) Go to 99

      if(boa.gt.1.0) then
       write (6,*)
       write (6,*) "Should be b/a < 1; please try again"
       write (6,*)
       Go to 97
      end if

      write (6,*)
      write (6,*) " Please enter the equivalent radius"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) req
      write (6,*)

      if(req.eq.0) Go to 99

c----------
c intialize
c----------

      R11 = 0.0
      R12 = 0.0
      R13 = 0.0

      R21 = 0.0
      R22 = 0.0
      R23 = 0.0

      R31 = 0.0
      R32 = 0.0
      R33 = 0.0

c----------------------
c compute the semi-axes
c----------------------

      scale = req/(boa*boa)**oot

      a = scale
      b = scale*boa

c-------
c launch
c-------

      c  = Dsqrt(a*a-b*b) 
      e  = c/a
      es = e*e
      tmp = log((1.0+e)/(1.0-e))
 
      R11 = 8.0/3.0 *e*a*es/(-2.0*e+(1.0+es)*tmp)
      R22 = 8.0/3.0 *e*a*2.0*es/(2.0*e-(1.0-3.0*es)*tmp)
      R33 = R22

c-----
c print
c------

      write (6,*) " -------------------------------"
      write (6,*)  " Force resistance for the translation"
      write (6,*)  " of a prolate spheroid"
      write (6,*)  " Spheroid is symmetric around the x axis"
      write (6,*)  " Force = - 6 pi* viscosity * R * U"
      write (6,*)  " The matrix R is given by:"
      write (6,*)
      write (6,100) R11,R12,R13
      write (6,100) R21,R22,R23
      write (6,100) R31,R32,R33
      write (6,*) " -------------------------------"

      Go to 97

c-----
c Done
c-----

 99   Continue

 100  Format (1x,4(2x,f12.8))

      Stop
      End
