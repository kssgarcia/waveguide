      subroutine sgf_3d_2p_aux 
     +
     +  (a,b,ew,eps
     +  ,AB
     +  ,ABD
     +  ,ABDD
     +  )

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c-------------------------------------------------
c Evaluation of equation (4.6) and its derivatives
c The derivatives are computed by central differences
c
c Definitions:
c -----------
c
c a = delta
c b = k/ew
c-------------------------------------------------

      Implicit double precision (a-h,o-z)

c--------
c prepare
c--------

      eps2 = 2.0D0*eps   ! for numerical differentiation

      ew2 = ew*ew
      ew3 = ew2*ew

      b2 = b*b
      b3 = b2*b

      d2 = 1.0D0+b2
      d  = sqrt(d2)
      d5 = d*d*d*d*d

c---
c begin evaluations
c---

      a2 = a*a
      c  = a*b
      ad = a*d

      AB = (1.0D0+c)*exp(-c)/b3
     +   + (2.0D0-b2+a*(2.0D0-b2)*d+a2*d2)*exp(-ad)/d5 

c-------

      a  = a+eps
      a2 = a*a
      c  = a*b
      ad = a*d

      AB2= (1.0D0+c)*exp(-c)/b3
     +   + (2.0D0-b2+a*(2.0-b2)*d+a2*d2)*exp(-ad)/d5 

c-------

      a  = a-2.0D0*eps
      a2 = a*a
      c  = a*b
      ad = a*d

      AB1= (1.0+c)*exp(-c)/b3
     +   + (2.0-b2+a*(2.0-b2)*d+a2*d2)*exp(-ad)/d5 

c---------

      a = a+eps      ! rectify

c---
c apply finite differences
c---

      ABD  = ew *(AB2-AB1)/(2.0D0*eps)
      ABDD = ew2*(AB2-2.0D0*AB+AB1)/eps**2

c---
c finish up
c---

      fc = 0.5D0/ew3

      AB   = fc*AB
      ABD  = fc*ABD 
      ABDD = fc*ABDD

c-----
c done
c-----

  100 Format (9(1x,f20.10))

      return
      end
