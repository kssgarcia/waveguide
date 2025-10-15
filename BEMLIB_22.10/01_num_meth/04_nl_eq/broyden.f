      subroutine broyden 
     +
     +   (n,menu
     +   ,Niter
     +   ,Ispeak
     +   ,x
     +   ,a
     +   ,Iflag
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c    Oxford University Press, 1998
c------------------------------------------------

c-------------------------------------------------
c Roots of a system of n nonlinear equations by the
c iterative method of Broyden (Math. Comp. 19, 1965). 
c
c An initial guess of the inverse jacobian is supplied 
c along with the initial guess of the solution vector.  
c
c Method is discussed in Section 4.6 of Pozrikidis (1998)
c
c The program returns Iflag=0 if convergence failed.
c
c This version sequentially reduces the parameter lambda
c searching for a reduction in the vector function for
c a fixed iteration of the inverse jacobian.
c
c SYMBOLS
c -------
c	
c  a .... inverse jacobian matrix
c  n .... size of a (no. of non-linear equations to be solved)
c  x .... solution of the non-linear system
c  f .... function vector
c  Niter .. maximum number of iterations allowed
c  tol .. convergence tolerance
c  Iflag  flag for convergence (1 = no convergence)
c       
c  x1 ... temporary storage of x vector
c  dx ... updating increment to solution vector
c  lambda . relaxation parameter
c  u .... vector for updating inverse jacobian
c  v .... vector for updating inverse jacobian
c
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision a(20,20),x(20),x1(20),f(20),f1(20)
      Double Precision dx(20),df(20),u(20),v(20)

      Parameter (tol=0.000000001)

c-----------
c initialize
c-----------

      Iflag  = 1

c---------------------------------------------------
c evaluate the vector function at the initial guess
c---------------------------------------------------

      call nol_fun (menu,x,f1)

c-----
c print
c------

      If(Ispeak.eq.1) then

       write (6,*)
       write (6,*) ' Iteration No, x(i), f(i)'
       write (6,*)

       write (6,100) Icount,(x(i),i=1,n),(f1(i),i=1,n)

      End If

c---------------------------------------------
c  calculate the update increment 
c  to the solution vector
c  according to eq. (4.6.19) with omega = 1.0
c---------------------------------------------

      Do i=1,n
        dx(i) = 0.0D0
        Do j=1,n
         dx(i) = dx(i)+a(i,j)*f1(j)
        End Do
        x(i) = x(i) - dx(i)
      End Do

      call nol_fun (menu,x,f)

      Icount = 1

c------
c print
c------

      If(Ispeak.eq.1) then
        write (6,100) Icount,(x(i),i=1,n),(f(i),i=1,n)
      End If

c-----------------
c first correction
c-----------------

      Do i=1,n
       df(i) = f(i)-f1(i)
      End do

c---------------------------------
c  iterate for the solution vector 
c  while updating both the solution
c  vector and the inverse of the jacobian
c---------------------------------

      Do while (Iflag.eq.1
     +         .and.Icount.le.Niter)    !  outer loop for solution

       Icount = Icount+1  
 
c---
c update the inverse jacobian 
c according to equation (4.6.18)
c
c par:	denominator
c---

      par = 0.0D0

      Do i=1,n
        Do j=1,n
         par = par+dx(i)*a(i,j)*df(j)
        End Do
      End Do

c--
c Guard against round-off error
c---

      If(abs(par).lt.tol**2) return

c---
c Numerator
c---

      Do i=1,n
        u(i) = 0.0D0
        v(i) = 0.0D0
        Do j=1,n
         u(i) = u(i)+a(i,j)*df(j)
         v(i) = v(i)+dx(i)*a(i,j)
        End Do
        u(i) = u(i)-dx(i)
      End Do

      Do i=1,n
        Do j=1,n
          a(i,j) = a(i,j)-(u(i)*v(j))/par
        End Do
      End Do

c-----------------------------------------
c rename the solution and function vectors
c-----------------------------------------

      Do i=1,n
       x1(i) = x(i)
       f1(i) = f(i)
      End Do

c---
c Compute updated solution vector 
c so that the vector function is
c monotonically decreasing.
c---

      alambda = 1.0D0

      alen1 = 0.0D0
 
      Do i=1,n
       alen1 = alen1+f1(i)**2
      End Do

      jflag  = 0
      jcount = 0
 
c---
      Do while(jflag.eq.0.and.jcount.le.Niter)
c---
 
         jcount=jcount+1  
 
         Do i=1,n
           dx(i) = 0.0 
           Do j   = 1,n
            dx(i) = dx(i)-alambda*a(i,j)*f1(j)
           End Do
           x(i) = x1(i)+dx(i)
         End Do
 
c---------------------
c ckeck the magnitude
c of the vector function 
c-----------------------
 
         call nol_fun (menu,x,f)

         alen = 0.0D0

         Do i=1,n
           alen = alen+f(i)**2
         End Do

         If(alen.lt.alen1) then
          jflag=1
         Else
          alambda = 0.9*alambda
         End If

c---
      End Do              !  end of loop for function evaluation
c---

c------
c print
c-----
   
      If(Ispeak.eq.1) then
        write (6,100) icount,(x(i),i=1,n),(f(i),i=1,n)
      End If

c---
c  compute df vector for next iteration
c---

      Do i=1,n
        df(i) = f(i)-f1(i)
      End Do

c----------------------
c check for convergence
c----------------------

      error = 0.0D0

      Do i=1,n
        error = error +f(i)**2
      End Do

      error = Dsqrt(error)/n

      If(error.lt.tol) Iflag = 0

c----------------------
      End Do               !  end of outer loop for solution vector
c----------------------

c-----
c Done
c-----

 100  Format (1X,I3,8(1X,F9.5),/,4X,8(1X,F9.5))

      Return
      End
