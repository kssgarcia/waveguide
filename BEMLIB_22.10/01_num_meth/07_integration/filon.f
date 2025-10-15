      subroutine filon 
     +
     +  (menu
     +  ,a,b,k
     +  ,n
     +  ,rint_cs,rint_sn
     +  )

c===========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c------------------------------------------------
c This program accompanies the book:
c                C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c---------------------------------------------------
c Computes the integral of an oscillatory integrand,
c as described in the text
c
c SYMBOLS:
c -------
c
c n: number of intervals
c menu: index for the function to be integrated
c a,b:  end-points of interval of integration
c
c rint_cs: integral of cosine
c rint_sn: integral of sine
c
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k
      Dimension x(901),f(901)

c--------
c prepare
c--------

      n1 = n+1
      h  = (b-a)/n    ! step size

c---
c integrand evaluation
c---

      Do i=1,n1
        x(i) = a + (i-1)*h
        If(menu.eq.1) then
           f(i) = 1.0
        Else If(menu.eq.2) then
           f(i) = x(i)
        Else If(menu.eq.3) then
           f(i) = x(i)**2
        Else If(menu.eq.4) then
           f(i) = x(i)**3
        Else If(menu.eq.5) then
           f(i) = exp(-x(i))
        End If
      End Do

c---
c  formulae (7.6.3)
c---

      c1 = 0.5D0*(f(n1)*cos(k*x(n1))+f(1)*cos(k*x(1)))
      s1 = 0.5D0*(f(n1)*sin(k*x(n1))+f(1)*sin(k*x(1)))

      Do i=3,n-1,2
        c1 = c1 + f(i)*cos(k*x(i))
        s1 = s1 + f(i)*sin(k*x(i))
      End Do

      c2 = 0.0D0
      s2 = 0.0D0

      Do i=2,n,2
        c2 = c2+f(i)*cos(k*x(i))
        s2 = s2+f(i)*sin(k*x(i))
      End Do

c---
c compute alpha, beta, gamma
c---

      u = k*h

      ws = f(n1)*sin(k*x(n1))-f(1)*sin(k*x(1))
      wc = f(n1)*cos(k*x(n1))-f(1)*cos(k*x(1))

      alpha  = (u**2+0.5*u*sin(2*u)-2.0*sin(u)**2)/u**3
      beta  = 2.0*(u*(1.0+cos(u)**2)-sin(2*u))/u**3
      gamma = 4.0*(sin(u)-u*cos(u))/u**3

c----------
c integrals
c----------

      rint_cs = h*( alpha*ws + beta*c1 + gamma*c2)
      rint_sn = h*(-alpha*wc + beta*s1 + gamma*s2)

c-----
c Done
c-----

      Return
      End
