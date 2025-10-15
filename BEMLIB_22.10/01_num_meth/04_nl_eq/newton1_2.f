      subroutine newton1_2 (menu,Niter,x,eps,italk)

c=======================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement
c=======================================

c------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c        Oxford University Press
c------------------------------------------------

c-----------------------------------
c Solution of one nonlinear equation
c by the quadratic Newton method
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (tol=0.0000001)

c---
c launch
c---

      Do i=1,Niter
        call newton1_fun (menu,x,f)
        x = x + eps          ! derivative by finite differences
        call newton1_fun (menu,x,f1)
        x = x - eps          ! reset
        Df = (f1-f)/eps
        Dx = -f/Df
        x  = x+Dx
        if(italk.eq.1) write (6,100) i,x,f
        if(abs(Dx).le.tol) Go to 99
      End Do

c---
c done
c---
      write (6,*)
      write (6,*) "Did not converge in ",Niter," iterations"

  99  continue

 100  Format (2X,I3,2X,F15.10,2X,F15.10)

      return
      end
