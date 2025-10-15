      subroutine power_sym 
     +
     +    (n,a,x,max
     +    ,l1,l2
     +    ,Icount
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
c C. Pozrikidis
c
c ``Numerical Computation in Science and Engineering''
c
c Oxford University Press, 1998
c------------------------------------------------

c------------------------------------------------
c  Compute the two largest eigenvalues
c  of a real symmetric matrix
c
c  If the subroutine returns without converging, 
c  either l1d or l2d may have became smaller
c  than allowed. 
c
c  SYMBOLS:
c  -------
c
c  a .... real symmetric square matrix
c  n .... size (rows/columns) of matrix a
c  x .... initial guess for eigen vector
c  max .. maximum number of iterations allowed
c  l1 ... dominant eigen value
c  l2 ... 2nd largest eigen value
c  Istop: flag for convergence (1 indicates no-convergence)
c
c  x1 ... j+1 approximation to eigen vector
c  x2 ... j+2 approximation to eigen vector
c  mu_01  
c  mu_11  Rayleigh-Schwartz quotients
c  mu_12  
c
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),x(10),x1(10),x2(10)

      Double Precision l1,l2,mu_01,mu_11,mu_12
      Double Precision l1n,l1d,l2n,l2d

      Parameter (eps=0.00000001)

c-----------
c initialize
c-----------

      Istop  = 0
      Icount = 0    ! counter of projections

c---
c normalize the input vector
c---

      sum = 0.0D0
      Do i=1,n
        sum = sum+x(i)*x(i)
      End Do

      Do i=1,n
       x(i) = x(i)/Dsqrt(sum)
      End Do

c----------------------------------------------
c begin iteration loop for eigenvector updating
c----------------------------------------------

      write (6,*)
      write (6,*) " Estimated eigenvalues:"
      write (6,*)

      Do 98 i=1,max          !  begin outer loop

        icount = icount + 1

c---
c transform the approximate eigenvector twice
c---

       Do j=1,n
         x1(j)=0.0D0
         Do k=1,n
           x1(j)=x1(j)+a(j,k)*x(k)
         End Do
       End Do
	
       Do j=1,n
         x2(j)=0.0D0
         Do k=1,n
          x2(j)=x2(j)+a(j,k)*x1(k)
         End Do
       End Do

c---
c compute the inner products for the R-S quotients
c---
 
       rs_00 = 0.0D0
       rs_01 = 0.0D0
       rs_11 = 0.0D0
       rs_12 = 0.0D0
	 
       Do j=1,n
         rs_00=rs_00 +  x(j)* x(j)
         rs_01=rs_01 +  x(j)*x1(j)
         rs_11=rs_11 + x1(j)*x1(j)
         rs_12=rs_12 + x1(j)*x2(j)
       End Do
	
c---
c compute R-S quotients
c using eqs (5.5.11)-(5.5.13)
c---
 
       mu_01 = rs_01/rs_00
       mu_11 = rs_11/rs_01
       mu_12 = rs_12/rs_11
 
c---
c compute eigenvalue estimates
c---
 
       save1 = l1
       save2 = l2

       l1n = (mu_12-mu_11)**2
       l1d = mu_12-2.0D0*mu_11+mu_01

       l1  = mu_12-l1n/l1d

       l2n = mu_12-mu_11
       l2d = mu_11-mu_01

       l2  = l1*l2n/l2d

       write (6,100) icount,l1,l2
 
c-----
c stopping check
c-----

       err1 = abs(save1 - l1)
       err2 = abs(save2 - l2)

       If(err1.lt.eps
c    +   .and.err2.lt.eps
     +   ) Return    ! success
 
c----------------------------------------------
c normalize vector x2 and check for convergence
c----------------------------------------------
 
        sum = 0.0
        Do j=1,n
         sum = sum+x2(j)*x2(j)
        End Do

        Do j=1,n
          x2(j) = x2(j)/Dsqrt(sum)
        End Do
	  
        sum=0.0
        Do j=1,n
         sum = sum+(x2(j)-x(j))*(x2(j)-x(j))
        End Do

c       If(sum.lt.eps) then
c         Return
c       End If

c-------------
c update vector x for next iterations
c-------------

       Do j=1,n
        x(j) = x2(j)
       End Do
	
  98  Continue              !  end of outer loop

c---------------------------------------
c if no convergence, then return istop=1
c---------------------------------------

      Istop = 1

 100  Format (1x,i3,2(1x,f20.10))

      Return
      End
