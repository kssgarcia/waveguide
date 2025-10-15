      program gel_inv_dr

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------------
c This program accompanies the book:
c
c               C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c------------------------------------------------

c------------------------------------------------
c Compute the inverse of a square matrix A
c by gauss elimination with row pivoting
c
c SYMBOLS:
c -------
c
c Iwlpvt: I will pivot; 1 for yes
c Isym  : Exploit symmetry ? 1 for yes
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   a(500,500),u(500,500)
      Double Precision       l(500,500)
      Dimension             ai(500,500)
      Dimension         verify(500,500)

c-----------
c preferences
c------------

 98   Continue

      write (6,*)
      write (6,*) " Is the matrix symmetric? (0 for No, 1 for Yes)"
      write (6,*) " Enter 2 to quit"
      write (6,*) " ----------------"
      read  (5,*) Isym

      if(Isym.eq.2) Go to 99

c-----------------------

      if(Isym.eq.0) then  ! symmetric system

       open (8, file="mat.dat")  ! input file

       write (6,*)
       write (6,*) " Enable row pivoting? (0 for No, 1 for Yes)"
       write (6,*) " ------------------------------------------"
       read  (5,*) Iwlpvt

      else     !  pivoting is disabled for a symmetric system

       open (8, file="mat_sym.dat")  ! input file

       Iwlpvt = 0 

       write (6,*) " gel_inv_dr: pivoting disabled"

      end if

c------------------------

c----------------
c read the matrix
c----------------

      read (8,*) n

      Do i=1,n
        read  (8,*) (a(i,j),j=1,n)
      End Do

      close (8)

c------------------
c find the solution
c------------------

      call gel_inv 
     +
     +   (n
     +   ,a
     +   ,ai               ! inverse
     +   ,Isym,Iwlpvt
     +   ,l,u
     +   ,det
     +   ,Istop
     +   )

c-----------------
c printing session
c-----------------

      write (6,*)
      write (6,*) " matrix A:"
      write (6,*) " ---------"
      write (6,*)

      Do i=1,n
       write (6,101) (a(i,j),j=1,n)
      End Do

      write (6,*)
      write (6,*) " matrix L:"
      write (6,*) " ---------"
      write (6,*)

      Do i=1,n
       write (6,101) (l(i,j),j=1,n)
      End Do

      write (6,*)
      write (6,*) " matrix U:"
      write (6,*) " ---------"
      write (6,*)

      Do i=1,n
       write (6,101) (u(i,j),j=1,n)
      End Do

      write (6,*) 
      write (6,*) " inverse of A:"
      write (6,*) " ------------"
      write (6,*) 

      Do i=1,n
       write (6,101) (ai(i,j),j=1,n)
      End Do

      write (6,*) 
      write (6,*) " Determinant:"
      write (6,*) " ------------"
      write (6,102) det

c-------
c verify
c-------

      Do i=1,n
       Do j=1,n
         verify(i,j)=0.0D0
         Do k=1,n
          verify(i,j) = verify(i,j)+a(i,k)*ai(k,j)
         End Do
       End Do
      End Do

      write (6,*) 
      write (6,*) " A times AI"
      write (6,*) " ----------"
      write (6,*) 

      Do i=1,n
       write (6,101) (verify(i,j),j=1,n)
      End Do

c-------
c repeat
c-------

      Go to 98

c-----
c done
c-----

 99   Continue

      write (6,*) "Thank you for using our services"

 100  Format (1x,i3,1x,f10.6)
 101  Format (20(1x,f8.4))
 102  Format (f20.5)

      stop
      end
