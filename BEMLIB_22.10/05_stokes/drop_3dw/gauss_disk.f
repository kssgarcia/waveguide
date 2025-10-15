      subroutine gauss_disk (n,xi,eta,w)

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c------------------------------------------------
c This program accompanies the book:
c
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c----------------------------------------
c Abscissas and weights for 
c Gaussian Integration over a circular disk
c of unit radius
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xi(30),eta(30),w(30)

c---
c constants
c---

      pi = 3.14159 265358D0

c------
c query
c------

      if(   n.ne.  4
     +  .and.n.ne. 5
     +  .and.n.ne.21
     +  ) then
         write (6,*) ' gauss_disk: number of gaussian points'
         write (6,*) '             not available'
         write (6,*) '             Will take 21 points'
         n = 7
      end if      

c--------------------
      if(n.eq.4) then
c--------------------

       xi(1) = 0.5D0
      eta(1) = 0.5D0
        w(1) = 0.25D0

       xi(2) =-0.5D0
      eta(2) = 0.5D0
        w(2) = 0.25D0

       xi(3) =-0.5D0
      eta(3) =-0.5D0
       w(3) = 0.25D0

       xi(4) = 0.5D0
      eta(4) =-0.5D0
        w(4) = 0.25D0

c-------------------------
      else if(n.eq.5) then
c-------------------------

       xi(1) = 1.0D0
      eta(1) = 0.0D0
        w(1) = 0.125D0

       xi(2) = 0.0D0
      eta(2) = 1.0D0
        w(2) = 0.125D0

       xi(3) = -1.0D0
      eta(3) =  0.0D0
        w(3) =  0.125D0

       xi(4) =  0.0D0
      eta(4) = -1.0D0
        w(4) =  0.125D0

       xi(5) = 0.0D0
      eta(5) = 0.0D0
        w(5) = 0.5D0

c--------------------------
      else if(n.eq.21) then
c--------------------------

      six  = 6.0D0
      rsix = Dsqrt(six)
      a    = sqrt((6.0D0-rsix)/10.0D0)
      b    = sqrt((6.0D0+rsix)/10.0D0)
      c    = (16.0D0+rsix)/360.0D0
      d    = (16.0D0-rsix)/360.0D0
      e    = 0.2D0*pi

      Do i=1,10
        cs = Dcos(i*e)
        sn = Dsin(i*e)
         xi(i) = a*cs
        eta(i) = a*sn
          w(i) = c
        j = i+10
         xi(j) = b*cs
        eta(j) = b*sn
         w (j) = d
      End Do

       xi(21) = 0.0D0
      eta(21) = 0.0D0
        w(21) = 1.0D0/9.0D0

c-----------
      end if
c-----------

c---
c scale the weights
c---

      Do i=1,n
       w(i) = w(i)*pi
      End Do

c-----
c done
c-----

      return
      end
