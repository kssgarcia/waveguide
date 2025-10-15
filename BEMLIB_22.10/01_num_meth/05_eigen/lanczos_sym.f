      subroutine lanczos_s 
     +
     +  (n
     +  ,a
     +  ,p
     +  ,alpha,beta
     +  ,Istop
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------
c Transform a real symmetric matrix into a similar 
c tridiagonal matrix by the method of Lanczos.
c
c Algorithm 5.7.1
c
c vector p(1,i) is imported from the driver.
c-------------------------------------------

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

c  _________________   external variables   ________________
c
c  a .... square matrix
c  n .... size (rows/columns) of matrix a
c  p .... set of orthogonal vectors
c  b .... tridiagonal matrix
c  l .... evaluation point for characteristic polynomial
c  f .... scalar lanczos function evaluated at l
c  Istop flag for successful completion (0 = completion)
c       
c  _________________   internal variables   ________________
c
c  q .... orthogonal matrix
c  alpha  coefficient for orthogonalization
c  beta . coefficient for orthogonalization
c  _________________________________________________________

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),p(10,10)
      Dimension alpha(10),beta(10)

      Parameter (eps=0.00000001)
	
c-------
c prepare
c--------

      Iflag = 0

c-----------------------------------
c normalize the initial vector p(1,i)
c-----------------------------------

      sum = 0.0D0

      Do i=1,n
        sum=sum+p(1,i)**2
      End Do

      Do i=1,n
        p(1,i)=p(1,i)/Dsqrt(sum)
      End Do

c------------------------
c compute the first alpha
c------------------------

      alpha(1)=0.0D0

      Do i=1,n
       Do j=1,n
         alpha(1)=alpha(1)+p(1,i)*a(i,j)*p(1,j)
       End Do
      End Do

c---------------
c compute p(2,i)
c---------------

      Do i=1,n
        p(2,i)=0.0D0
        Do j=1,n
          p(2,i)=p(2,i)+a(i,j)*p(1,j)
        End Do
        p(2,i)=p(2,i)-alpha(1)*p(1,i)
      End Do

c-----------------------
c compute the first beta
c-----------------------

      sum=0.0D0

      Do i=1,n
         sum=sum+p(2,i)**2
      End Do
      beta(1)=Dsqrt(sum)

      If(abs(beta(1)).lt.eps) then
       Istop = 1
       Return
      End If

c----------------------------
c normalize the second vector
c----------------------------

      Do i=1,n
        p(2,i)=p(2,i)/beta(1)
      End Do

c-----------------------------
c Loop over krylov projections
c-----------------------------

      Do 95 i=2,n         !  outer loop over krylov iterations

c---
c compute alpha
c---

      alpha(i)=0.0D0

      Do j=1,n
        Do k=1,n
         alpha(i)=alpha(i)+p(i,j)*a(j,k)*p(i,k)
        End Do
      End Do

c---
c  compute the components of the p-vector
c---

      Do j=1,n
        p(i+1,j)=0.0D0
        Do k=1,n
          p(i+1,j)=p(i+1,j)+a(j,k)*p(i,k)
        End Do
        p(i+1,j)=p(i+1,j)-alpha(i)*p(i,j)-beta(i-1)*p(i-1,j)
      End Do

c---
c compute beta
c---

      sum=0.0D0

      Do j=1,n
        sum=sum+p(i+1,j)**2
      End Do

      beta(i) = Dsqrt(sum)

      If(abs(beta(i)).lt.eps) then
        Istop = 1
        Return
      End If

c---
c normalize the current p vector
c----

      Do j=1,n
        p(i+1,j)=p(i+1,j)/beta(i)
      End Do

  95  Continue            !  end of outer loop

c-----
c Done
c-----

      Return
      End
