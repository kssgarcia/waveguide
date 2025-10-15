      subroutine lanczos 
     +
     +  (n,a,p,s
     +  ,alpha,beta,gamma
     +  ,Istop
     +  )

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
c        C. Pozrikidis
c Numerical Computation in Science and Engineering
c    Oxford University Press,1998
c------------------------------------------------

c---------------------------------
c Tridiagonalize an arbitrary non-symmetric matrix
c by the method of Lanczos
c
c The vectors p(1,i) and s(1,i) are imported from driver.
c
c Algorithm 5.8.1
c
c  SYMBOLS:
c  -------
c
c  a .... arbitrary (non-symmetric) matrix
c  n .... size (rows/columns) of matrix a
c  p .... set of normalized vectors orthogonal to s
c  s .... set of vectors orthogonal to p
c  Istop. flag for successful completion (0 = completion)
c
c  alpha:  scale factor for vectors p, s
c  beta:   scale factor for vectors p, s
c  gamma:  scale factor for vectors p, s
c
c  Greek variables will be elements of the tridiagonal matrix
c
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),p(10,10),s(10,10)
      Dimension alpha(10),beta(10),gamma(10)

      Parameter (eps=0.0000000001)

c-----
c flag
c-----

      Istop = 0

c-----------------------------
c normalize the first p-vector
c-----------------------------

      sum = 0.0D0

      Do i=1,n
        sum = sum+p(1,i)**2
      End Do

      sum = Dsqrt(sum)

      Do i=1,n
        p(1,i) = p(1,i)/sum
      End Do

c---------------------------------------------
c scale the first s-vector so that p(1).s(1)=1
c---------------------------------------------

      sum = 0.0D0

      Do i=1,n
        sum = sum+p(1,i)*s(1,i)
      End Do

      Do i=1,n
        s(1,i) = s(1,i)/sum
      End Do

c------------------------
c compute the first alpha
c------------------------

      alpha(1) = 0.0D0

      Do i=1,n
        Do j=1,n
         alpha(1) = alpha(1)+s(1,i)*a(i,j)*p(1,j)
        End Do
      End Do

c---------------------------------------
c compute the second vector p (unscaled)
c---------------------------------------

      Do i=1,n
        p(2,i) = -alpha(1)*p(1,i)
        Do j=1,n
          p(2,i) = p(2,i)+a(i,j)*p(1,j)
        End Do
      End Do

c----------------------
c compute the first beta
c-----------------------

      sum = 0.0D0

      Do i=1,n
        sum = sum+p(2,i)**2
      End Do

      beta(1) = Dsqrt(sum)

c-----
c trap
c-----

      If(abs(beta(1)).lt.eps) then
        Istop = 1
        write (6,*) " lanczos: Returning from position 1"
        Return
      End If

c---------------------------------------
c compute the second vector s (uncsaled)
c---------------------------------------

      Do i=1,n
        s(2,i) = -alpha(1)*s(1,i)
        Do j=1,n
         s(2,i) = s(2,i)+a(j,i)*s(1,j)
        End Do
      End Do

c-----------------------
c compute the first gamma
c------------------------

      sum = 0.0D0

      Do i=1,n
        sum = sum+p(2,i)*s(2,i)
      End Do

      gamma(1) = sum/beta(1)

c---
c trap
c---

      If(abs(gamma(1)).lt.eps) then
        Istop = 1
        Return
      End If

c--
c scale the second vectors p and s
c---

      Do i=1,n
        p(2,i) = p(2,i)/beta(1)
        s(2,i) = s(2,i)/gamma(1)
      End Do

c---
c  begin outer loop for vectors x and y,
c  and parameters alpha, beta, gamma
c---

      Do 95 i=2,n              !  begin loop over krylov iteration

c---
c compute alpha
c---

      alpha(i) = 0.0D0

      Do j=1,n
        Do k=1,n
         alpha(i) = alpha(i)+s(i,j)*a(j,k)*p(i,k)
        End Do
      End Do

c---
c compute unscaled vector p
c---

      Do j=1,n
        p(i+1,j) = -alpha(i)*p(i,j)-gamma(i-1)*p(i-1,j)
        Do k=1,n
         p(i+1,j) = p(i+1,j)+a(j,k)*p(i,k)
        End Do
      End Do

c---
c  compute beta
c---

      sum = 0.0D0

      Do j=1,n
        sum = sum+p(i+1,j)**2
      End Do

      beta(i) = sqrt(sum)

c---
c trap
c---

      If(abs(beta(i)).lt.eps) then
       Istop = 1
       Return
      End If
 
c---
c  compute unscaled vector s
c---

      Do j=1,n
        s(i+1,j) = -alpha(i)*s(i,j)-beta(i-1)*s(i-1,j)
        Do k=1,n
          s(i+1,j) = s(i+1,j)+a(k,j)*s(i,k)
        End Do
      End Do

c---
c compute gamma
c---

      sum = 0.0D0

      Do j=1,n
        sum = sum+p(i+1,j)*s(i+1,j)
      End Do

      gamma(i)=sum/beta(i)
 
c---
c trap
c---

      If(abs(gamma(i)).lt.eps) then
       Istop = 1
       Return
      End If

c---
c scale vectors p and s
c---

      Do j=1,n
        p(i+1,j) = p(i+1,j)/beta(i)
        s(i+1,j) = s(i+1,j)/gamma(i)
      End Do

c-------------
  95  Continue                   !  end of outer loop
c-------------
       
c-----
c Done
c-----

      Return
      End
