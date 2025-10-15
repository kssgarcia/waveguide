      subroutine secant 
     +    (menu
     +    ,Niter
     +    ,xs
     +    ,Iflag
     +    ,italk
     +    )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c--------------------------------------------
c  Approximate a root of a single function 
c  of one variable by
c  a newton type iteration in which
c  the derivative is evaluated
c  at each step by a secant of the function.  
c
c  Iflag returns 0 if convergence fails
c--------------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c        C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c  ___   external variables   __________________
c
c  xs .... initial guess for root
c  tol ... convergence criterion
c  Niter ... max allowable iterations
c  Iflag . flag for convergence of iteration
c
c  ___   internal variables   __________________
c
c  xk1 ... x(k-1), root estimate from previous step
c  xk .... x(k), root estimate at current step
c  x ..... x(k+1), root estimate at next step
c  del ... interval for evaluation of initial secant
c  ______________________________________________

      Implicit Double Precision (a-h,o-z)

      Parameter (del=0.01,tol=0.0000001)

c--------
c prepare
c--------

      Iflag  = 0
      Icount = 1

      xk1 = xs-del
      call secant_fun (menu,xk1,fk1)
      xk  = xs

c--------
c iterate
c--------

       Do while(Icount.le.niter)

        call secant_fun (menu,xk,fk)

        if(italk.eq.1) write (6,100) Icount,xk,fk

        if(abs(fk).lt.tol) then
          Iflag = 1
          xs  = xk
          Go to 98
        end If

        slope = (fk-fk1)/(xk-xk1)  

        If(abs(slope).lt.tol) Go to 98

        x   = xk-fk/slope
        xk1 = xk
        fk1 = fk
        xk  = x
        Icount = Icount+1

      End Do

c-----
c done
c-----

      xs = x   ! did not converge

 98   Continue

 100  Format (2X,I3,2X,F15.10,2X,F15.10)

      Return
      End
