      subroutine lapl_2d
     +
     +  (method
     +  ,x,y
     +  ,Dx
     +  ,Dy
     +  ,Rlap
     +  )

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press
c------------------------------------------------

c----------------------------------------
c Compute the laplacian of a function of
c two variables by finite differences
c
c Section 6.12
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension A(3,3),F(3,3)

c-------------------------------------------
c function definition  must be placed on top
c-------------------------------------------

      func(x,y) = (x*y*Dcos(x)-(1.5D0+2.0D0*y+y**2)*Dsin(x)) * Dexp(-y)
c     write (6,*) " function: f(x,y) = ( x y cosx ... ) exp(-y)"

c     func(x,y) = (x*y*Dcos(x)+(1.5D0-2.0D0*y+y**2)*Dsin(x)) * Dexp(y)
c     write (6,*) " function: f(x,y) = ( x y cosx ... ) exp(y)"

c     func(x,y) = exp(x)*sin(y)

c     func(x,y) = exp(x)*cos(y)

c     func(x,y) = x**2+y**2
c     write (6,*) " function: f(x,y) = x^2+y^2"

c     func(x,y) = cos(x)*cos(y)
c     write (6,*) " function: f(x,y) = cosx cosy"

c------------------------------------
c generate the differentiation matrix
c------------------------------------

      if(method.eq.1) then  ! five-point formula

         beta = (Dx/Dy)**2

         A(1,1) =  0.0D0
         A(1,2) =  beta
         A(1,3) =  0.0D0
         A(2,1) =  1.0
         A(2,2) = -2.0D0*(1.0D0+beta)
         A(2,3) =  1.0D0
         A(3,1) =  0.0D0
         A(3,2) =  beta
         A(3,3) =  0.0D0
         cf     =  Dx*Dx

      else if(method.eq.2) then

         A(1,1) =  1.0
         A(1,2) =  0.0
         A(1,3) =  1.0
         A(2,1) =  0.0
         A(2,2) = -4.0
         A(2,3) =  0.0
         A(3,1) =  1.0
         A(3,2) =  0.0
         A(3,3) =  1.0
         cf     =  2.0*Dx*Dx

      else if(method.eq.3) then

         A(1,1) =   1.0D0
         A(1,2) =   4.0D0
         A(1,3) =   1.0D0
         A(2,1) =   4.0D0
         A(2,2) = -20.0D0
         A(2,3) =   4.0D0
         A(3,1) =   1.0D0
         A(3,2) =   4.0D0
         A(3,3) =   1.0D0
         cf     =   6.0D0*Dx*Dx

      end if

c--------------------
c function evaluation
c--------------------

      f(3,1) = func(x-Dx,y-Dy)
      f(2,1) = func(x-Dx,y)
      f(1,1) = func(x-Dx,y+Dy)
c---
      f(3,2) = func(x,y-Dy)
      f(2,2) = func(x,y)
      f(1,2) = func(x,y+Dy)
c---
      f(3,3) = func(x+Dx,y-Dy)
      f(2,3) = func(x+Dx,y)
      f(1,3) = func(x+Dx,y+Dy)

c----------------------
c compute the laplacian
c----------------------

      sum = 0.0D0

      Do i=1,3
        Do j=1,3
        sum = sum + A(i,j)*f(i,j) 
        End Do
      End Do

      Rlap = sum/cf

c-----
c Done
c-----

      return
      end
