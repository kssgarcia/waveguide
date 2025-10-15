      subroutine hermite_card (n,m,xp,x,q)

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
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c 1998
c------------------------------------------------

c-----------------------------------------
c  Evaluates cardinal polynomials for higher-order 
c  Hermite interpolation
c  with (m-1) matched derivatives.  
c
c  Up to 4th order derivatives in this implementation
c
c  Derivatives of the auxiliary polynomials r
c  are computed by finite differences with a user specified
c  interval (eps),using 2nd order five point centered differences. 
c
c  SYMBOLS:
c  --------
c
c  x .... x coordinate for interpolation
c  xp ... x coordinate of prescribed data
c  n .... number of prescribed points (n+1)
c  m .... number of conditions at each prescribed point
c  r .... auxiliary polynomial
c  rp ... auxiliary polynomial evaluated at prescribed point
c  q .... cardinal polynomial 
c  eps .. interval for numerical computation of derivatives
c
c  fact . factorial
c  prod . product of n+1 terms in auxiliary polynomials
c  t .... temporary r
c  s .... integer counter for derivative #
c
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xp(50),r(50,0:4),rp(50,0:4,0:4),q(50,0:4),t(-2:2)
      Dimension m(10)
      Integer s

      Parameter (eps=0.00001)

c--------
c prepare
c--------

      n1 = n+1

      epst = 2.0D0*eps
      eps2 = eps**2
      eps3 = eps**3
      eps4 = eps**4

c-----------------------------------------------------------
c  Algorithm (6.7.11)
c
c  Evaluate auxiliary polynomial at interpolation point and
c  prescribed data (and adjacent points for determination of 
c  derivatives).
c-----------------------------------------------------------

      Do i=1,n1
        Do k=0,m(i)-1

          p = x
          r(i,k) = rik(p,xp,m,n,i,k) 
    
c---
c  compute derivatives of auxiliary polynomials
c  numerically
c---
  
         Do l=-2,2
           p = xp(i)+l*eps
           t(l) = rik(p,xp,m,n,i,k)
         End Do 
	   
         rp(i,k,0)=t(0)

         Do s=1,m(i)-1
          If(s.eq.1) then
           rp(i,k,1)=(t(1)-t(-1))/epst   
          Else If(s.eq.2) then
           rp(i,k,2)=(t(-1)-2.0*t(0)+t(1))/eps2  
          Else If(s.eq.3) then
           rp(i,k,3)=(-t(-2)+2.0*t(-1)-2.0*t(1)+t(2))/(2.0*eps3)
          Else If(s.eq.4) then
           rp(i,k,4)=(t(-2)-4.0*t(-1)+6.0*t(0)-4.0*t(1)+t(2))/eps4
          End If
         End Do
	  
        End Do  !  loop over k
      End Do   !  loop over i

c---
c evaluate the cardinal polynomial q recursively 
c using algorithm (6.7.12)
c---
	
      Do i=1,n1
       Do k=m(i)-1,0,-1

        q(i,k) = r(i,k)

        If(k.ne.m(i)-1) then
          sum=0.0D0
          Do s=k+1,m(i)-1
           q(i,k)=q(i,k)-rp(i,k,s)*q(i,s)
          End Do
        End If

       End Do
      End Do

c-----
c Done
c-----

      Return
      End

c============================================
	
      function rik(x,xp,m,n,i,k)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---
c evaluate r_i,k(x) defined in algorithm (6.7.11)
c---

      Implicit Double Precision (a-h,o-z)
      Dimension m(10),xp(50)
	
c---
c compute factorial of k
c---
	
      If(k.eq.0) then
        fact = 1.0D0
      Else
        fact =1.0D0
        Do l=1,k
         fact =fact*real(l)
        End Do
      End If

c---
c evaluate the auxillary polynomial
c---

      prod = 1.0D0
	
      Do j=1,n+1
        If(j.ne.i) then
          prod = prod*((x-xp(j))/(xp(i)-xp(j)))**m(j)
        End If
      End Do

      rik = ((x-xp(i))**k)*prod/fact

c-----
c Done
c-----

      End
