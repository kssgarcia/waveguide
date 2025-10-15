       program gram_schmidt

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c----------------------------------------
c Gram-Schmidt orthogonalization of a set of N 
c M-dimensional vectors, where N is less than M
c QR Decomposition of A
c----------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension u(100,100),v(100,100),w(100,100),a(100,100)
      Dimension b(100,100)

c---
c read the vectors
c---

      open (3,file="matrix_gs.dat")
       read (3,*) M,N
       Do i=1,N
        read (3,*) (v(j,i),j=1,M)
       End Do
      close (3)

c---
c initialize the coefficient matrix alpha
c----

      Do i=1,M
        Do j=1,N
          a(i,j) = 0.0D0
        End Do
        a(i,i) = 1.0D0
      End Do

c------------------------
c Orthogonal Gram-Schmidt
c------------------------

c---
c First vector
c---

      Do k=1,M
        u(k,1) = v(k,1)
      End Do

c---
c Further vectors
c---

      Do j=2,N

       Do k=1,M
         u(k,j) = v(k,j)
       End do

       Do i=1,j-1

         tmp1 = 0.0D0
         tmp2 = 0.0D0

         Do l=1,N
           tmp1 = tmp1 + u(l,i)*v(l,j)
           tmp2 = tmp2 + u(l,i)*u(l,i)
         End Do

         a(i,j) = tmp1/tmp2

         Do k=1,M
          u(k,j) = u(k,j)-a(i,j)*u(k,i)
         End Do

       End Do

      End Do

c-----------------------
c Printing and verifying
c-----------------------

      write (6,*)
      write (6,*) " Columns of v"
      write (6,*) " ------------"
      write (6,*)

      Do i=1,N
        write (6,100) (v(i,j),j=1,M)
      End Do

      write (6,*)
      write (6,*) " Columns of u (Matrix Q_pr)"
      write (6,*) " --------------------------"
      write (6,*)

      Do i=1,N
        write (6,100) (u(i,j),j=1,M)
      End Do

      write (6,*)
      write (6,*) " Matrix of alphas, R_pr"
      write (6,*) " ----------------------"
      write (6,*)

      Do i=1,N
        write (6,100) (a(i,j),j=1,M)
      End Do

      write (6,*)
      write (6,*) " Verify Qpr R_pr = V"
      write (6,*) " -------------------"
      write (6,*)

      Do i=1,N
       Do j=1,M
         b(i,j) = 0.0D0
         Do k=1,N
          b(i,j) = b(i,j)+u(i,k)*a(k,j)
         End Do
       End Do
      End Do

      Do i=1,N
        write (6,100) (b(i,j),j=1,M)
      End do

c-----------------
c QR decomposition
c-----------------

      Do i=1,N

        dn = 0.0D0

        Do j=1,N
          dn = dn+u(j,i)*u(j,i)
        End Do

        dn = Dsqrt(dn)

        Do j=1,N
          u(j,i) = u(j,i)/dn
          a(i,j) = a(i,j)*dn
        End Do

      End Do

      write (6,*)
      write (6,*) " Matrix Q"
      write (6,*) " --------"
      write (6,*)

      Do i=1,N
        write (6,100) (u(i,j),j=1,M)
      End Do

      write (6,*)
      write (6,*) " Matrix R"
      write (6,*) " --------"
      write (6,*)

      Do i=1,N
        write (6,100) (a(i,j),j=1,M)
      End Do

      write (6,*)
      write (6,*) " Verify Q R = V"
      write (6,*) " --------------"
      write (6,*)

      Do i=1,N
       Do j=1,N
        b(i,j) = 0.0
        Do k=1,N
         b(i,j) = b(i,j)+u(i,k)*a(k,j)
        End Do
       End Do
      End Do

      Do i=1,N
        write (6,100) (b(i,j),j=1,M)
      End do

      write (6,*)
      write (6,*) " Verify Q QT = I"
      write (6,*) " ---------------"
      write (6,*)

      Do i=1,N
       Do j=1,N
         b(i,j) = 0.
         Do k=1,N
          b(i,j) = b(i,j)+u(i,k)*u(j,k)
         End Do
       End Do
      End Do

      Do i=1,N
        write (6,100) (b(i,j),j=1,M)
      End Do

c-------------------------
c Orthonormal Gram-Schmidt
c-------------------------

c---
c First vector
c---

       den = 0.0D0

       Do k=1,M
        den = den + v(k,1)*v(k,1)
       End Do

       den = sqrt(den)

       Do k=1,M
        w(k,1) = v(k,1)/den
       End Do

c---
c Further vectors
c---

      Do j=2,N

        Do k=1,M
         w(k,j) = v(k,j)
        End Do

        Do i=1,j-1
           tmp = 0.0D0
           Do l=1,M
            tmp = tmp+w(l,i)*v(l,j)
           End Do
           Do k=1,M
            w(k,j) = w(k,j)-tmp*w(k,i)
           End Do
        End Do

        den = 0.0D0

        Do k=1,M
          den = den + w(k,j)*w(k,j)
        End Do

        den = sqrt(den)

        Do k=1,M
          w(k,j) = w(k,j)/den
        End Do

      End Do

c---------
c printing
c---------

      write (6,*)
      write (6,*) " Matrix of w"
      write (6,*) " -----------"

      Do i=1,N
        write (6,100) (w(i,j),j=1,M)
      End Do

c-----
c Done
c-----

 100  Format (100(1x,f10.5))

      Stop
      End
