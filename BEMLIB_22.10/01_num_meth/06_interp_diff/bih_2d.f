      program bih_2d

c====================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c====================================

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c     Oxford University Press
c------------------------------------------------

c----------------------------------------
c Compute the biharmonic of a function of
c two variables by finite differences
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c------
c input
c------

c     write (6,*) 
c     write (6,*) " Choose the method"
c     write (6,*) 
c     write (6,*) " Enter: "
c     write (6,*) 
c     write (6,*) " 0 to quit"
c     write (6,*) " 1 for the five-point formula"
c     write (6,*) " 2 for the diagonal five-point formula"
c     write (6,*) " 3 for the nine-point formula"
c     write (6,*) " ----------------------------"
c     read  (5,*) method

      method = 1

      If(method.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter x and y:"
      write (6,*) "---------------"
      read  (5,*) x, y

c----
c return to repeat
c----

 98   Continue

      write (6,*) 
      write (6,*) " Please enter Dx (= Dy)"
      write (6,*) " 0 to quit"
      write (6,*) " ----------"
      read  (5,*) Dx

      If(Dx.eq.0) Go to 99

      Dy = Dx

c----------------------
c compute the laplacian
c at 5 points
c----------------------

      call lapl_2d
     +
     +  (method
     +  ,x,y
     +  ,Dx
     +  ,Dy
     +  ,Rlap0
     +  )

      write (6,101) Rlap0

      call lapl_2d
     +
     +  (method
     +  ,x-Dx,y
     +  ,Dx
     +  ,Dy
     +  ,Rlap1
     +  )

      call lapl_2d
     +
     +  (method
     +  ,x+Dx,y
     +  ,Dx
     +  ,Dy
     +  ,Rlap2
     +  )

      call lapl_2d
     +
     +  (method
     +  ,x,y-Dy
     +  ,Dx
     +  ,Dy
     +  ,Rlap3
     +  )

      call lapl_2d
     +
     +  (method
     +  ,x,y+Dy
     +  ,Dx
     +  ,Dy
     +  ,Rlap4
     +  )

c---
c compute the laplacian of the laplacian
c---

      bih = (Rlap1+Rlap2+Rlap3+Rlap4-4.0D0*Rlap0)/Dx**2

      write (6,100) bih

c-----------------
c return to repeat
c-----------------

      Go to 98

c-----
c done
c-----

 99   Continue

      write (6,*)
      write (6,*) " thank you for running me"
      write (6,*)

 100  Format (" Biharmonic = ",f22.12)
 101  Format (" Laplacian  = ",f22.12)

      Stop
      End
