      subroutine newton 
     +
     +   (n
     +   ,menu
     +   ,Niter
     +   ,eps
     +   ,Ispeak
     +   ,Jac
     +   ,x
     +   ,Iflag
     +   )

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c---------------------------------------------
c  Newton's method for a system of n equations
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

      Dimension x(20),f(20),f1(20),dx(20),rhs(20)
      Double precision Jac(20,20)

      Parameter (tol=0.00001,relax=1.00)

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

c---------------------
c start the iterations
c---------------------

      Do 1 Iter=1,Niter

       call nol_fun (menu,x,f)

       if(Ispeak.eq.1) then 
         write (6,100) Iter,(x(i),i=1,n),(f(i),i=1,n)
       endif

c---------------------
c compute the Jacobian
c by numerical differentiation
c---------------------

c     if(Iter==1) then

      Do j=1,n

        x(j) = x(j)+eps      ! perturb

        call nol_fun (menu,x,f1)

        x(j) = x(j)-eps      ! reset

        Do i=1,n
          Jac(i,j) = (f1(i)-f(i))/eps
        End Do

      End Do

c     end if

c-------------------
c print the Jacobian
c-------------------
c
c     Do i=1,n
c      write (6,101) (Jac(i,j),j=1,n)
c     End Do

c---
c
c Solve equation: Jac . Dx = - f
c
c for the correction vector Dx
c
c Cramer's rule for 2 or 3 equations
c
c Gauss elimination for a higher number
c of equations
c---

c-------------------
      if(n.eq.2) then   ! Cramer's rule
c-------------------

        b1  = -f(1)
        b2  = -f(2)
        Det = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)

        dx(1) = (b1*Jac(2,2)-Jac(1,2)*b2)/Det
        dx(2) = (b2*Jac(1,1)-Jac(2,1)*b1)/Det

c-------------------
      else if(n.eq.3) then  ! Cramer's rule
c-------------------

        A11 = Jac(1,1)
        A12 = Jac(1,2)
        A13 = Jac(1,3)
        A21 = Jac(2,1)
        A22 = Jac(2,2)
        A23 = Jac(2,3)
        A31 = Jac(3,1)
        A32 = Jac(3,2)
        A33 = Jac(3,3)

        b1 = -f(1)
        b2 = -f(2)
        b3 = -f(3)

        call cramer33 
     +
     +    (A11,A12,A13
     +    ,A21,A22,A23
     +    ,A31,A32,A33
     +    ,b1,b2,b3
     +    ,dx(1),dx(2),dx(3)
     +    )

c---------
      else     ! Gauss elimination
c---------

        Do i=1,n
          rhs(i) = -f(i)
        End Do

        Isym_gel = 0  ! system is non-symetric
        Iwlpvt   = 1  ! pivoting enabled

        call gel
     +
     +    (n,Jac,rhs,dx
     +    ,Isym_gel,Iwlpvt
c    +    ,l,u
     +    ,det
     +    ,Istop
     +    )

        if(Istop.eq.1) Go to 99

c---------
      end if
c---------

c--------
c correct
c--------

      Do i=1,n
        x(i) = x(i)+relax*dx(i)
      End Do

c-------
c escape
c-------

      Do i=1,n
       If(abs(dx(i)).gt.tol) Go to 1
      End Do

      Go to 99

  1   Continue           ! Target for iterations

c-------------
c unsuccessful
c-------------

      write (6,*) "Newton: solution did not converge"
      write (6,*) "        in ",Niter," iterations"
      Iflag = 1
      Return

c-----------
c Successful
c-----------

 99   Continue

c--------
c confirm
c--------

      Iter = Iter+1

      call nol_fun (menu,x,f)

      if(Ispeak.eq.1)  then
        write (6,100) Iter,(x(i),i=1,n),(f(i),i=1,n)
      endif

c-----
c Done
c-----

 100  Format (1X,I3,8(1X,F8.5),/,4X,8(1X,F8.5))
 101  Format (18(1X,F8.3))

      Return
      End
