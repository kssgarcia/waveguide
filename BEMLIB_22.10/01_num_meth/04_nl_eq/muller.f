      subroutine muller (menu,Niter,xs
     +
     + ,iflag
     + ,italk
     +)

c======================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c======================================

c------------------------------------------------
c This program accompanies the book:
c
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c----------------------------------------------
c  Locate a root of a single function of one variable by
c  iteration using a parabola fit between three consecutive points.
c
c  Function evaluation is done in the user-defined function
c  'fun'
c  which must be modified for each problem (see the driver)
c
c  Algorithm (4.6.6)
c
c  SYMBOLS
c  -------
c
c  xs .... initial guess for root
c  root .. returned value for approximated root
c  tol ... convergence criterion
c  Niter ... max allowable iterations
c  iflag . flag for convergence of iteration
c          returns 0 if convergence fails)
c
c  x ...... root estimate at next step
c  xk ..... x(k), root estimate at current step
c  xk1 .... x(k-1), root estimate from previous step
c  xk2 .... x(k-2), root estimate from next previous step
c  del ... interval for evaluation of initial secant
c
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (del=0.01,tol=0.00000001)

c-----------
c initialize
c-----------

      Iflag  = 0
      Icount = 1   ! counter

c-----------
c first step
c-----------

      xk2 = xs-2.0D0*del
      xk1 = xs-del

      call muller_fun (menu,xk2,fk2)
      call muller_fun (menu,xk1,fk1)

      xk  = xs

c--------
c iterate
c--------

      Do while(icount.le.niter)

       call muller_fun (menu,xk,fk)

       if(italk.eq.1) write (6,100) icount,xk,fk

       If(abs(fk).lt.tol) then
         Iflag = 1
         xs  = xk
         Go to 98
       End If

       q = (xk-xk1)/(xk1-xk2)
       a = q*fk-q*(q+1.0D0)*fk1+q**2*fk2
       b = (1.0D0+2.0*q)*fk-(1.0D0+q)**2*fk1+q**2*fk2
       c = (1.0D0+q)*fk
       d = b**2-4.0D0*a*c

c---------------------
       if(d.gt.0) then
c---------------------

         d1 = b+Dsqrt(d)
         d2 = b-Dsqrt(d)

         if(abs(d1).gt.abs(d2)) then
            fact = 2.0D0*c/d1
         else
            fact = 2.0D0*c/d2
         end if

c----------
       else
c----------

         f    = Dsqrt(abs(d)) 
         fact = 2.0D0*c*b/(b**2+f**2)

c------------
       end if
c------------

       x   = xk-(xk-xk1)*fact
       xk2 = xk1
       fk2 = fk1
       xk1 = xk
       fk1 = fk
       xk  = x
       icount = icount+1

      End Do

c---
c done iterating
c---

      xs = x    ! did not converge

 98   Continue

c-----
c done
c-----

 100  Format (2X,I3,2X,F15.10,2X,F15.10)

      Return
      End
