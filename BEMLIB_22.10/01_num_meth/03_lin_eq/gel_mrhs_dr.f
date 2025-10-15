      program gel_mrhs_dr

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

c-------------------------------------
c Solves the linear system Ax = b 
c with multiple right-hand sides
c by Gauss Elimination
c
c computes the inverse of the matrix A
c
c SYMBOLS:
c -------
c
c iwlpvt: I will pivot; 1 for yes
c isym  : Exploit symmetry ? 1 for yes
c
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(500,500),rhs(500,500),u(500,500)
      Double Precision l(500,500)
      Dimension x(500,500),ver(500,500)

c---
c preferences
c---

 98   Continue

      write (6,*)
      write (6,*) " Is system symmetric ?"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 for no  "
      write (6,*) " 1 for yes "
      write (6,*) " 2 to quit"
      write (6,*) " ---------"
      read  (5,*) Isym

      If(Isym.eq.2) Go to 99

c-----------------------

      If(Isym.eq.0) then  ! symmetric system

       open (8, file="mat_vec_many.dat")  ! input file

       write (6,*)
       write (6,*) " Is row pivoting enabled ?"
       write (6,*) " 1 for yes "
       write (6,*) " -------------------------"
       read  (5,*) iwlpvt

      Else     ! pivoting is disabled for a symmetric system

       open (8, file="mat_s_vec_many.dat")  ! input file

       Iwlpvt = 0

       write (6,*) " gel_mrhs_dr: pivoting is disabled"

      End If

c---------------------------------
c read matrix and right-hand sides
c---------------------------------

      read (8,*) n
      Do i=1,n
        read  (8,*) (a(i,j),j=1,n)
      End Do

      read (8,*) nrhs

      Do j = 1,nrhs
        read (8,*) (rhs(i,j),i=1,n)
      End Do

      close(8)

c------------------
c find the solution
c------------------

      call gel_mrhs 
     +
     +    (n,a
     +    ,nrhs,rhs
     +    ,x
     +    ,isym,iwlpvt
     +    ,l,u
     +    ,det
     +    ,Istop
     +    )

c-----------------
c printing session
c-----------------

      write (6,*)
      write (6,*) " matrix A"
      write (6,*) " --------"
      write (6,*)

      Do i=1,n
       write (6,101) (a(i,j),j=1,n)
      End Do

      write (6,*)
      write (6,*) " Right-hand sides"
      write (6,*) " ---------------"
      write (6,*)

      Do j=1,nrhs
       write (6,101) (rhs(i,j),i=1,n)
      End Do

      write (6,*)
      write (6,*) " matrix L"
      write (6,*) " --------"
      write (6,*)

      Do i=1,n
       write (6,101) (l(i,j),j=1,n)
      End Do

      write (6,*)
      write (6,*) " matrix U"
      write (6,*) " --------"
      write (6,*)

      Do i = 1,n
       write (6,101) (u(i,j),j=1,n)
      End Do

      write (6,*) 
      write (6,*) " Solution vectors"
      write (6,*) " ----------------"
      write (6,*) 

      Do i=1,n
        write (6,100) i,(x(i,j),j=1,nrhs)
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

      Do ll=1,nrhs
       Do i=1,n

        res = rhs(i,ll)

        Do j=1,n
         res = res - a(i,j)*x(j,ll)
        End Do
        write (6,100) i,res

       End Do
       write (6,*)
      End Do

c----------------
c matrix inverse
c----------------

c---
c put the identity into the
c right-hand sides
c---

      nrhs = n
      Do j=1,n
       Do i=1,n
        rhs(i,j) = 0.0D0
       End Do
       rhs(j,j) = 1.0D0
      End Do

c---
c compute the inverse
c of matrix A
c---

      call gel_mrhs 
     +
     +    (n,a
     +    ,nrhs,rhs
     +    ,x
     +    ,isym,iwlpvt
     +    ,l,u
     +    ,det
     +    ,Istop
     +    )

      write (6,*) 
      write (6,*) " Matrix inverse"
      write (6,*) " --------------"

      Do i=1,n
        write (6,101) (x(i,j),j=1,nrhs)
      End Do

c-------------------
c verify the inverse
c-------------------

      Do i=1,n
       Do j=1,n
        ver(i,j) = 0.0D0
        Do k = 1,n
         ver(i,j) = ver(i,j)+a(i,k)*x(k,j)
        End Do
       End Do
      End Do

      write (6,*) 
      write (6,*) " Product A * A(inv)"
      write (6,*) " ------------------"

      Do i=1,n
        write (6,101) (ver(i,j),j=1,n)
      End Do

c---
c repeat
c---

      Go to 98

c-----
c Done
c-----

 99   Continue

      write (6,*)
      write (6,*) "Thank you for running me"

 100  Format (1x,i3,100(1x,f10.6))
 101  Format (20(1x,f8.4))
 102  Format (f20.7)

      Stop
      End
