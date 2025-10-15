      subroutine power 
     +
     +    (n,a,x
     +    ,Inv,max,eps
     +    ,ev,icount
     +    ,Istop
     +    )

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
c                   C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c             Oxford University Press
c------------------------------------------------

c-------------------------------------------------------
c Computation of the leading eigenvalue of a real matrix
c by the power method according to algorithm (5.5.4).
c
c Additional eigenvalues are found by the
c repeated application of the power method
c in combination with matrix deflation. 
c
c
c SYMBOLS:
c --------
c
c  a .... square matrix
c  n .... size (rows/columns) of matrix a
c  x .... initial guess for eigenvector 
c  max .. maximum number of iterations allowed
c  eps .. accuracy
c
c  Istop: flag for convergence (1 indicates no convergence)
c
c  inv = 1 ... for regular iterations
c  inv = 2 ... for inverse iterations
c
c  u .... approximated eigenvector
c  save . old eigenvector for assessing convergence
c  ev.... approximated eigenvalue
c  evs... sequence of approximated eigenvalues
c
c--------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(512,512),x(512),u(512),save(512)

      Dimension evs(1000)

c-----------
c initialize
c-----------

      Icount = 0  ! counter
      Istop  = 0

c-------------------------
c normalize initial vector
c-------------------------

      par = 0.0D0

      Do i=1,n
        par = par+x(i)**2
      End Do

      par = Dsqrt(par)

      Do i=1,n
        x(i) = x(i)/par
      End Do

c----------------------
c begin the projections
c----------------------

      write (6,*) 
      write (6,*) " Power method"
      write (6,*) 
      write (6,*) " Iteration, eigenvalue, extrapolated"
      write (6,*)

      Do 97 iter =1,max            !  loop over iterations

      Icount = Icount+1

c----
c transform the vector by the matrix a
c---

c----------------------
      If(Inv.eq.1) then               ! regular iteration
c----------------------

       Do j=1,n
         u(j) = 0.0D0
         Do k=1,n
          u(j) = u(j)+a(j,k)*x(k)
         End Do
       End Do

c---------------------------
      Else If(Inv.eq.2) then     ! inverse iteration
c---------------------------

        Isym   = 0 ! set to 0 to allow for pivoting
        Iwlpvt = 1 ! pivoting enabled

        call gel_power
     +
     +    (n,a,x,u
     +    ,Isym
     +    ,Iwlpvt
c    +    ,l,u
     +    ,det
     +    ,Istop
     +    )

c-----------
      End If
c-----------

c-------------------------------
c compute approximate eigenvalue
c-------------------------------

      ev = 0.0D0

      Do j=1,n
       ev = ev + u(j)*x(j)
      End Do

      evs(Icount) = ev
 
c-------------------------
c normalize the new vector
c------------------------

      par = 0.0D0

      Do i=1,n
        par = par+u(i)**2
      End Do

      fc = ev*Dsqrt(par)/(abs(ev))

      Do i=1,n
       u(i) = u(i)/fc
      End Do

c----------------------
c check for convergence
c----------------------

      If(iter.ge.2) then

          sum=0.0
          Do j=1,n
           sum = sum+(u(j)-save(j))**2
          End Do
          sum = sqrt(sum)

          If(sum.lt.eps) then      ! converged
           Do j=1,n
            x(j) = u(j)
           End Do
           Return
          End If

      End If

c---------------
c update vectors
c---------------

      Do j=1,n
        save(j) = u(j)
           x(j) = u(j)
      End Do

c---------------------
c Aitken extrapolation
c---------------------

       If(icount.gt.2) then
         e0 = evs(icount-2)
         e1 = evs(icount-1)
         e2 = evs(icount)
         ea = (e0*e2-e1**2)/(e2-2.0*e1+e0)
         write (6,100) icount,evs(icount),ea
       Else
         write (6,100) icount,evs(icount)
       End If
              
 97   Continue                !  end of iteration loop

c---------------------------------------------
c If it did not converge after max iterations,
c                         set: Istop=1
c---------------------------------------------

      Istop = 1

c-----
c Done
c-----

 100  Format (1x,i3,3(1x,f20.10))

      Return
      End
