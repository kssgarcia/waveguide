      subroutine newton2
     +
     +   (menu
     +   ,Niter
     +   ,eps
     +   ,Ispeak
     +   ,x
     +   ,Iflag
     +   )

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------------
c This program accompanies the book:
c
c         C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c    Oxford University Press
c------------------------------------------------

c---------------------------------------------
c  Newton's method for two nonlinear equations
c
c  SYMBOLS:
c  --------
c
c  eps:	      small number for computing the Jacobian
c   	      by numerical differentiation
c  Dx: 	      correction vector
c  tol:	      accuracy
c  Iflag:     will set equal to 1 if something is wrong
c
c  Ispeak = 0: silent mode
c           1: verbose mode
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(2),f(2),f1(2),dx(2),rhs(2)
      Double precision Jac(2,2)

      Parameter (tol=0.0000001,relax=1.00)

c-------
c message
c--------

      if(Ispeak.eq.1) then

       write (6,*) 
       write (6,*) " Performing the Newton iterations"
       write (6,*) " --------------------------------"

       write (6,*) 
       write (6,*) ' Iteration No, x(i), f(i)'
       write (6,*) 

      end if

c-----------
c initialize
c-----------

      Iflag = 0
      n = 2

c---------------------
c start the iterations
c---------------------

      Do 1 Iter=1,Niter

       call newton2_fun (menu,x,f)

       if(Ispeak.eq.1) then 
         write (6,100) Iter,(x(i),i=1,n),(f(i),i=1,n)
       end if

c---------------------
c Compute the Jacobian
c by numerical differentiation
c---------------------

      Do j=1,2
        x(j) = x(j)+eps      ! perturb
        call newton2_fun (menu,x,f1)
        x(j) = x(j)-eps      ! reset
        Do i=1,2
          Jac(i,j) = (f1(i)-f(i))/eps
        End Do
      End Do

c-------------------
c print the Jacobian
c-------------------
c
c     Do i=1,2
c      write (6,101) (Jac(i,j),j=1,n)
c     End Do

c---
c
c solve the equation: Jac . Dx = - f
c for the correction vector Dx
c by Cramer's rule
c---

      b1  = -f(1)
      b2  = -f(2)
      Det = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)
      dx(1) = (b1*Jac(2,2)-Jac(1,2)*b2)/Det
      dx(2) = (b2*Jac(1,1)-Jac(2,1)*b1)/Det

c--------
c correct
c--------

      x(1) = x(1) + relax*dx(1)
      x(2) = x(2) + relax*dx(2)

c-------
c escape
c-------

      if(abs(dx(1)).gt.tol) Go to 1
      if(abs(dx(2)).gt.tol) Go to 1

      Go to 99

  1   Continue     ! target for iterations

c-------------
c unsuccessful
c-------------

      write (6,*) "newton2: iterations did not converge"
      write (6,*) "             in ",Niter," iterations"
      Iflag = 1
      return

c-----------
c successful
c-----------

 99   Continue

c--------
c confirm
c--------

      Iter = Iter+1

      call newton2_fun (menu,x,f)

      if(Ispeak.eq.1) then
        write (6,100) Iter,(x(i),i=1,n),(f(i),i=1,n)
      end if

c-----
c done
c-----

 100  Format (1X,I3,8(1X,F8.5),/,4X,8(1X,F8.5))
 101  Format (18(1X,F8.3))

      return
      end
