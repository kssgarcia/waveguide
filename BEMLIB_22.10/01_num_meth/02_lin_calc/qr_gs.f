       subroutine qr_gs (a,n,index,q,r)

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
c          C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press, 1998
c------------------------------------------------

c------------------------------------------------
c QR Factorization 
c by Gram-Schmidt orthogonalization
c
c If index = 1 will do the orthogonal Gram-Schmidt
c If index = 2 will do the orthoNORMAL Gram-Schmidt
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),q(128,128),r(128,128)
      Dimension u(128,128),v(128,128),alpha(128,128)

c     Dimension ver(128,128)

c------------------------
c copy into the v matrix
c------------------------

      Do i=1,n
        Do j=1,n
         v(j,i) = a(j,i)
        End Do
      End Do

c------------------------
c orthogonal Gram-Schmidt
c------------------------

      If(index.eq.1) then

c---
c initialize the coefficient matrix alpha
c----

      Do i=1,N
        Do j=1,N
          alpha(i,j) = 0
        End Do
        alpha(i,i) = 1.0D0
      End Do

c---
c first Vector
c---

      Do k=1,N
        u(k,1) = v(k,1)
      End Do

c---
c further Vectors
c---

      Do j=2,N

        Do k=1,N
          u(k,j) = v(k,j)
        End Do

        Do i=1,j-1

          t1 = 0.0D0
          t2 = 0.0D0

          Do l=1,N
            t1 = t1 + u(l,i)*v(l,j)
            t2 = t2 + u(l,i)*u(l,i)
          End Do

          alpha(i,j) = t1/t2

          Do k=1,N
            u(k,j) = u(k,j)-alpha(i,j)*u(k,i)
          End Do

        End Do

      End Do

c-----------------------
c Printing and verifying
c-----------------------
c
c     write (6,*)
c     write (6,*) " Matrix A"
c     write (6,*)
c     Do i=1,N
c       write (6,100) (v(i,j),j=1,N)
c     End Do
c
c     write (6,*)
c     write (6,*) " Matrix Qpr"
c     write (6,*)
c     Do i=1,N
c       write (6,100) (u(i,j),j=1,N)
c     End Do
c
c     write (6,*)
c     write (6,*) " Matrix Rpr"
c     write (6,*)
c     Do i=1,N
c       write (6,100) (alpha(i,j),j=1,N)
c     End Do
c
c---------
c  Verify that Qpr Rpr = A
c---
c     Do i=1,N
c       Do j=1,N
c           ver(i,j) = 0.
c           Do k=1,N
c           ver(i,j) = ver(i,j)+u(i,k)*alpha(k,j)
c           End Do
c        End Do
c     End Do
c---------
c Printing
c---------
c     write (6,*)
c     write (6,*) " Verify A = Qpr Rpr"
c     write (6,*)
c     Do i=1,N
c       write (6,100) (ver(i,j),j=1,N)
c     End Do
c---------

c-----------------
c QR decomposition
c----------------

      Do i=1,N

        dn = 0.0D0

        Do j=1,N
          dn = dn+u(j,i)**2
        End Do

        dn = sqrt(dn)

        Do j=1,N
          q(j,i) = u(j,i)/dn
          r(i,j) = alpha(i,j)*dn
        End Do

      End Do

c---------
c Printing
c---------
c----------------
c     write (6,*)
c     write (6,*) " Matrix Q"
c     write (6,*) " --------"
c     Do  i=1,N
c       write (6,100) (q(i,j),j=1,N)
c     End Do
c----------------

      End If     ! End of orthogonal Gram-Scmidt

c-------------------------
c orthonormal Gram-Schmidt
c-------------------------

      If(index.eq.2) then

c---
c First Vector
c---

       den = 0.0D0

       Do  k=1,n
        den = den + v(k,1)**2
       End Do

       den = Dsqrt(den)

       Do k=1,n
        q(k,1) = v(k,1)/den
       End Do

c---
c Further Vectors
c---

      Do j=2,N

        Do k=1,M
          q(k,j) = v(k,j)
        End Do

        Do i=1,j-1

          tmp = 0.0D0

          Do l=1,M
            tmp = tmp + q(l,i)*v(l,j)
          End Do

          Do k=1,M
            q(k,j) = q(k,j)-tmp*q(k,i)
          End Do

        End Do

        den = 0.0D0

        Do k=1,M
          den = den + q(k,j)**2
        End Do

        den = sqrt(den)

        Do k=1,M
         q(k,j) = q(k,j)/den
        End Do

      End Do

c---------
c Printing
c---------

c----------------
c     write (6,*)
c     write (6,*) " Matrix Q again"
c     write (6,*) " --------------"
c     Do  i=1,N
c       write (6,100) (q(i,j),j=1,N)
c     End Do
c----------------

      End If     ! End of orthonormal Gram-Scmidt

c-----
c Done
c-----

 100  Format (100(1x,f10.5))

      Return
      End
