        program qr_dec

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c               C. Pozrikidis
c ''Numerical Computation in Science and Engineering''
c       Oxford University Press, 1998
c------------------------------------------------

c----------------------------------
c QR decomposition of a real matrix
c by three methods
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),q(128,128),r(128,128),v(128,128)

  97  Continue

      write (6,*)
      write (6,*) " Choose the method"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 for Gram-Schmidt"
      write (6,*) " 2 for reflection (Householder)"
      write (6,*) " 3 for rotation   (Givens)"
      write (6,*) " --------------------------"
      Read  (5,*) Method

      If(Method.eq.0) Go to 99

      If(Method.eq.1) then
        write (6,*)
        write (6,*) "   Enter:"
        write (6,*)
        write (6,*) "   0 to quit"
        write (6,*) "   1 for orthognal Gram-Schmidt"
        write (6,*) "   2 for orthonormal Gram-Schmidt"
        write (6,*) "   ------------------------------"
        read (5,*) index
      End If

c----------------
c read the matrix
c----------------

      open (unit=8,file='matrix.dat')

       read (8,*) n
       Do i=1,n
         read(8,*) (a(i,j),j=1,n)
       End Do

      close(8)

c---------------------
c Do the decomposition
c---------------------

      If(Method.eq.1) call qr_gs     (a,n,index,q,r)
      If(Method.eq.2) call qr_reflex (a,n,q,r)
      If(Method.eq.3) call qr_rot    (a,n,q,r)

c------
c print
c------

      write (6,*)
      write (6,*) " Original Matrix A:"
      write (6,*) " -----------------"

      Do i=1,n
        write (6,100) (a(i,j),j=1,n)
      End Do 

      write (6,*)
      write (6,*) " Orthogonal Matrix Q:"
      write (6,*) " -------------------"

      Do i=1,n
        write(6,100) (q(i,j),j=1,n)
      End Do 

      write (6,*)
      write (6,*) " Right-Triangular Matrix R:"
      write (6,*) " -------------------------"

      Do i=1,n
        write(6,100) (r(i,j),j=1,n)
      End Do 

c-------
c Verify
c-------

      Do i=1,n
        Do j=1,n
          v(i,j) = 0.0D0
          Do m=1,n
            v(i,j) = v(i,j)+q(i,m)*r(m,j)
           End Do
         End Do
      End Do

      write (6,*)
      write (6,*) " Product QR (should be equal to A)"
      write (6,*) " ---------------------------------"
      write (6,*)

      Do i=1,n
         write(6,100) (v(i,j),j=1,n)
      End Do 

      Go to 97

c-----
c Done
c-----

 99   Continue

 100  Format(20(3x,f10.5))

      Stop
      End
