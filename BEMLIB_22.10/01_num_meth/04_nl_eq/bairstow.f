      subroutine bairstow 
     +
     +   (a
     +   ,npoly
     +   ,x
     +   ,Iflag
     +   )

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------------
c This program accompanies the book:
c
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c-------------------------------------------------
c  Bairstow's method for finding all roots
c  of an n'th degree 
c  polynomial with real coefficients.  
c
c  Equations are derived in Section 4.7
c  of Pozrikidis (1998)
c
c  A truly crash-free code may require user input
c  for the initial guesses for r and s
c  for some of the reduced-order
c  polynomials
c
c  Failure to converge is indicated by Iflag
c  returned as 0
c
c  SYMBOLS
c  -------
c       
c  npoly  order of the polynomial
c
c  a .... vector of polynomial coefficients
c         a(1) is coefficient of x**n
c
c  x .... roots of polynomial, possibly complex:
c  x(i,1) real part of x(i)
c  x(i,2) imaginary part of x(i)
c
c  Iflag: Flag for success (0 indicates failure)
c	
c  b .... vector of order n-2 polynomial coefficients
c  c .... vector of partial derivatives
c
c  max... maximum number of Newton iterations
c  relax: relaxation factor in Newton iterations
c  eps:   accuracy
c
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(20),b(20),c(20),x(20,2)

      Parameter (eps=0.0000000001,max=200,relax=1.0)

c----------------------------------------
c start outer loop for synthetic devision
c
c order of descendant polynomial
c is reduced by 2 after each pass
c----------------------------------------

      n = npoly             !   order of original polynomial

      Do while(n.gt.2)     !  loop over reduced order polynomials

      Iflag  = 0
      jcount = 0   ! iteratin counter

      r = 1.0D0    ! initial guess
      s = 1.0D0    ! initial guess
  
c---
c inner loop for Newton iteration
c on synthetic division
c---

      b(1) = a(1)
      c(1) = 0.0D0
      c(2) = b(1)
  
      Do while (Iflag.eq.0)           ! loop until convergence    

       jcount = jcount+1    ! count iterations

       b(2) = a(2)+r*b(1)

       Do i=3,n+1
        b(i) = a(i)+r*b(i-1)+s*b(i-2)
       End Do

       c(3) = b(2)+r*c(2)

       Do i=4,n+1
         c(i) = b(i-1)+r*c(i-1)+s*c(i-2)
       End Do

c---
c check for convergence and increment (r,s)
c---

       if(    abs(b(n)  ).lt.eps
     +   .and.abs(b(n+1)).lt.eps) then
         Iflag = 1
       else
         dt =  c(n)**2-c(n-1)*c(n+1)
         dr = -( b(n)*c(n)-c(n-1)*b(n+1) )/dt
         ds = -( c(n)*b(n+1)-c(n+1)*b(n) )/dt
         r  = r + relax*dr
         s  = s + relax*ds
       end if 

       if(jcount.eq.max) return  ! failed to converge
   
      End Do   !  end of inner loop for solution of (r,s)

c---
c given (r,s), find the next pair 
c of roots of the polynomial
c---

      p = 1.0D0

      call bairstow_root (p,-r,-s,x,n)
 
c---
c reduce the size of the working polynomial
c and update the coefficient vector
c---
 
      n = n-2
 
      Do i=1,n+1
       a(i) = b(i)
      End Do

      End Do       !  end of outer loop over polynomial order 

c--------------------------
c compute the final root(s)
c--------------------------

      if(n.eq.1) then    ! odd-degree polynomial   

        x(1,1) = -a(2)/a(1)
        x(1,2) = 0.0D0

      else if(n.eq.2) then   ! even-degree polynomial

        call bairstow_root(a(1),a(2),a(3),x,n)

      end if

c------------------------------------
c set convergence flag to "converged"
c------------------------------------

      Iflag=1

c-----
c done
c-----

      return
      end
