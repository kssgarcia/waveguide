      program det_hess

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c-----------------------------------
c Determinant of a Hessenberg matrix
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10)

  97  Continue    ! return to repeat

      write (6,*)
      write (6,*) " Enter"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " 1 for lower Hessenberg"
      write (6,*) " 2 for upper Hessenberg"
      write (6,*) "-----------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

c---
c read the matrix
c---

      If(menu.eq.1) open (3,file="matrix_hl.dat")
      If(menu.eq.2) open (3,file="matrix_hu.dat")

      read (3,*) n
      Do i = 1,n
        read (3,*) (a(i,j),j=1,n)
      End Do

      close (3)

c---
c call hessenberg
c---

      If(menu.eq.1) call det_hess_l (n,a,det)
      If(menu.eq.2) call det_hess_u (n,a,det)

      write (6,100) det

      Go to 97

c-----
c Done
c-----

   99 Continue

  100 Format("Determinant = ", f10.5)

      Stop
      End 
