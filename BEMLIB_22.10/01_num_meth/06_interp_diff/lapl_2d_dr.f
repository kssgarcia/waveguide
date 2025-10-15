      program lapl_2d_dr

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
c           C. Pozrikidis
c
c ``Numerical Computation in Science and Engineering''
c
c       Oxford University Press, 1998
c------------------------------------------------

c----------------------------------------
c Compute the laplacian of a function of
c two variables 
c
c Section 6.12
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension A(3,3),F(3,3)

c------
c input
c------

 98   Continue

      write (6,*) 
      write (6,*) " Choose the method"
      write (6,*) 
      write (6,*) " Enter: "
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 for the five-point formula"
      write (6,*) " 2 for the diagonal five-point formula"
      write (6,*) " 3 for the nine-point formula"
      write (6,*) " ----------------------------"
      read  (5,*) method

      If(method.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter x and y:"
      write (6,*) "---------------"
      read  (5,*) x, y

      If(method.eq.1) then

         write (6,*) 
         write (6,*) " Please enter Dx and Dy"
         write (6,*) " 0 to quit"
         write (6,*) " ---------"
         read  (5,*) Dx, Dy

         If(Dx.eq.0) Go to 99
         If(Dy.eq.0) Go to 99

      Else

         write (6,*) 
         write (6,*) " Please enter Dx (= Dy)"
         write (6,*) " 0 to quit"
         write (6,*) " ----------"
         read  (5,*) Dx
         If(Dx.eq.0) Go to 99
         Dy = Dx

      End If

c----------------------
c compute the laplacian
c----------------------

      call lapl_2d
     +
     +  (method
     +  ,x,y
     +  ,Dx
     +  ,Dy
     +  ,Rlap
     +  )

      write (6,100) Rlap

c---
c return to repeat
c---

      Go to 98

c-----
c Done
c-----

 99   Continue

      write (6,*)
      write (6,*) " thank you for running me"
      write (6,*)

 100  Format (" Laplacian = ",f22.12)

      Stop
      End
