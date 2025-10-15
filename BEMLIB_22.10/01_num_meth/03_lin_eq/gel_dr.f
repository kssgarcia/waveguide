      program gel_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press, 1998
c------------------------------------------------

c------------------------------------------------
c Solves the linear system Ax = b 
c by Gauss Elimination
c with row pivoting
c
c SYMBOLS:
c -------
c
c Iwlpvt: I will pivot;      1 for yes
c Isym  : Exploit symmetry ? 1 for yes
c
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(500,500),b(500),u(500,500)
      Double Precision l(500,500)
      Dimension x(500)

c------------
c preferences
c------------

 98   Continue

      write (6,*)
      write (6,*) " Is system symmetric ?"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 for no"
      write (6,*) " 1 for yes "
      write (6,*) " 2 to quit"
      write (6,*) " ----------"
      read  (5,*) Isym

      If(Isym.eq.2) Go to 99

c---------------------
c non-symmetric system
c---------------------

      If(Isym.eq.0) then 

       open (8, file="mat_vec.dat")  ! input file

       write (6,*)
       write (6,*) " Is row pivoting enabled ?"
       write (6,*)
       write (6,*) " 0 for no "
       write (6,*) " 1 for yes"
       write (6,*) " ---------"
       read  (5,*) Iwlpvt

c-----------------
c symmetric system
c-----------------

      Else     ! pivoting is disabled for a symmetric system

       open (8, file="mat_s_vec.dat")  ! input file

       Iwlpvt = 0 

       write (6,*) " gel_dr: pivoting disabled"

      End If

c--------------------------------
c read matrix and right-hand side 
c--------------------------------

      read (8,*) n

      Do i=1,n
        read  (8,*) (a(i,j),j=1,n)
      End Do

      read (8,*) 
      read (8,*) (b(j),j=1,n)

      close (8)

c------------------
c find the solution
c------------------

      call gel
     +
     +  (n
     +  ,a,b,x
     +  ,Isym
     +  ,Iwlpvt
     +  ,l,u
     +  ,det
     +  ,Istop
     +  )

c-----------------
c printing session
c-----------------

      write (6,*)
      write (6,*) " matrix A:"
      write (6,*) " --------"
      write (6,*)

      Do i=1,n
       write (6,101) (a(i,j),j=1,n)
      End Do

      write (6,*)
      write (6,*) " Right-hand side:"
      write (6,*) " ---------------"
      write (6,*)

      write (6,101) (b(j),j=1,n)

      write (6,*)
      write (6,*) " matrix L:"
      write (6,*) " ---------"
      write (6,*)

      Do i=1,n
       write (6,101) (l(i,j),j=1,n)
      End Do

      write (6,*)
      write (6,*) " matrix U:"
      write (6,*) " --------"
      write (6,*)

      Do i=1,n
       write (6,101) (u(i,j),j=1,n)
      End Do

      write (6,*) 
      write (6,*) " Solution vector:"
      write (6,*) " ---------------"
      write (6,*) 

      Do i=1,n
        write (6,100) i,x(i)
      End Do

      write (6,*) 
      write (6,*) " Determinant:"
      write (6,*) " -----------"
      write (6,102) Det
   
c---
      write (6,*) 
      write (6,*) " Residuals:"
      write (6,*) " ----------"
      write (6,*) 

      Do i=1,n
        res = b(i)
        Do j=1,n
         res = res - a(i,j)*x(j)
        End Do
        write (6,100) i,res
      End Do

c-------
c repeat
c-------

      Go to 98

c-----
c Done
c-----

 99   Continue

      write (6,*)  "gel: thank you"

 100  Format (1x,i3,1x,f10.6)
 101  Format (20(1x,f8.4))
 102  Format (f20.5)

      Stop
      End
