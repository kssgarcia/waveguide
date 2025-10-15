      subroutine newton1_3 
     +
     +    (menu
     +    ,Niter
     +    ,x,eps
     +    ,italk
     +    )

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c========================================

c------------------------------------------------
c This program accompanies the book:
c
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c--------------------------------------------
c Modified Newton method with cubic convergence
c for one equation
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (tol=0.0000001)

c---
c launch
c---

      Do i=1,Niter

        call newton1_fun (menu,x,f)      ! derivative
        x1 = x + eps
        call newton1_fun (menu,x1,f1)

        fD = (f1-f)/eps
        Dx = - f/fD
        xt = x + Dx

        call newton1_fun (menu,xt,ft)
        Dx = - (f+ft)/fD
        x  = x + Dx

        if(italk.eq.1) write (6,100) i,x,f

        if(abs(Dx).le.tol) Go to 99

      end do

      write (6,*)
      write (6,*) "newton1_3: did not converge in ",Niter," iterations"

   99 Continue

c---
c Done
c---

 100  Format (2X,I3,2X,F15.10,2X,F15.10)

      Return
      End
