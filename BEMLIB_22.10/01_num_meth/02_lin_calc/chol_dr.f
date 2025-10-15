      program chol_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c      Oxford University Press
c------------------------------------------------

c------------------------------------------------
c Cholesky decomposition of a symmetric and
c positive-definite matrix by column or row
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)
      Double Precision l(10,10)
      Dimension a(10,10),v(10,10)

c-------
c launch
c-------

  98  Continue       ! return to repeat"

      write (6,*)
      write (6,*) "Choose the method"
      write (6,*)
      write (6,*) "Enter 0 to quit"
      write (6,*) "      1 for decomposition by column"
      write (6,*) "      2 for decomposition by row"
      write (6,*) "--------------------------------"
      read  (5,*) method

      If(method.eq.0) Go to 99

c----------------
c read the matrix
c----------------

      open (unit=2,file='matrix_sym.dat')

      read (2,*) n
      Do i=1,n
        read (2,*) (a(i,j),j=1,n)
      End Do

      close(2)

c----------------------------
c Carry out the decomposition
c----------------------------
    
      If(method.eq.1) call chol_c (l,a,n)
      If(method.eq.2) call chol_r (l,a,n)

c---------
c printing
c---------

      write (6,*)
      write (6,*) " Matrix A"
      write (6,*)

      Do i=1,n
       write(6,100) (a(i,j),j=1,n)
      End Do 

      write (6,*)
      write (6,*) " Matrix L"
      write (6,*)

      Do i=1,n
        write(6,100) (l(i,j),j=1,n)
      End Do 

c------------------------
c Verify: L*L(transp) = A
c------------------------

      write (6,*)
      write (6,*) " Matrix L*LT"
      write (6,*)

      Do i=1,n
         Do j=1,n
           v(i,j) = 0.0D0
           Do m=1,n
             v(i,j) = v(i,j)+l(i,m)*l(j,m)
           End Do
         End Do
       End Do

      Do i=1,n
        write(6,100) (v(i,j),j=1,n)
      End Do 

      Go to 98

c-----
c Done
c-----

 99   Continue

100   Format(20(3x,f10.7))

      Stop
      End
