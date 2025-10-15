      function exp_int(x)

c===========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c------------------------------------------------
c This program accompanies the book:
c
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c----------------------------------
c Polynomial and rational approximations
c for the exponential integral
c
c Abromowitz and Stegun, p. 231
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      if(x.le.1.0) then

         a0 = -0.57721 566
         a1 =  0.99999 193
         a2 = -0.24991 055
         a3 =  0.05519 968
         a4 = -0.00976 004
         a5 =  0.00107 857

         x2 = x*x
         x3 = x*x2
         x4 = x*x3
         x5 = x*x4

         expint = -log(x)+a0 + a1*x + a2*x2 + a3*x3 +
     +             a4*x4 + a5*x5

c--------
c lower-order accuracy
c
c        expint = -log(x)
c        sign = 1.0
c        fact = 1.0
c        power = 1.0
c        Do i = 1,10
c          fact = fact*i
c          power = power*x
c          expint = expint + sign*power/(i*fact)
c          sign = -sign
c        End Do
c--------

      else

c--------
c lower-order accuracy
c
c        a1 = 2.334733 
c        a2 = 0.250621 
c        b1 = 3.330657 
c        b2 = 1.681534 
c        x2 = x*x
c        expint = (x2+a1*x+a2)/(x2+b1*x+b2)/(x*exp(x))
c--------

         a1 = 8.57332 87401
         a2 =18.05901 69730
         a3 = 8.63476 08925
         a4 =  .26777 37343

         b1 = 9.57332 23454
         b2 =25.63295 61486
         b3 =21.09965 30827
         b4 = 3.95849 69228

         x2 = x*x
         x3 = x*x2
         x4 = x*x3

         expint = (x4+a1*x3+a2*x2+a3*x+a4)
     +           /(x4+b1*x3+b2*x2+b3*x+b4)/(x*exp(x))

c-----------
      end if
c-----------

      return
      end
