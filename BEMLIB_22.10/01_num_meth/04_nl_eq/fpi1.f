      program fixed_point_iter

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c    C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c----------------------------------------
c One-point iterations for one equation
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(0:2000)

c----------
c constants
c----------

      null = 0
      three = 3.0

      open (3,file="fpi1.out")

  95  Continue

c---
c read input
c---

      write (6,*) 
      write (6,*) "          MENU"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " 1 for g(x) = log(3) + 2*log(x)"
      write (6,*) " 2 for g(x) = exp(x/2)/sqrt(3)"
      write (6,*) " 3 for g(x) = exp(-x)/3"
      write (6,*) " 4 for g(x) = exp(-x)-2x"
      write (6,*) " 5 for g(x) = l x(1-x) (Logistic mapping)"
      write (6,*) " 6 for g(x) = 1 - 0.2 * x**3"
      write (6,*) " -----------------------------------------"
      Read  (5,*) menu

      If(menu.eq.0)  Go to 99

      If(menu.eq.5) then
        write (6,*) 
        write (6,*) " Enter lamda"
        write (6,*) " -----------"
        read  (5,*) rl
      End If

      write (6,*)
      write (6,*) " Perform Thiele-Aitken extrapolation ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      Read  (5,*) Iextr

      write (6,*)
      write (6,*) " Enter initial guess"
      write (6,*) " -------------------"
      Read  (5,*) x(0)

c---
c run
c---
      i = 0   ! counter
      k = 0   ! counter

  98  Continue

      xx = x(k)
      i = i+1
      k = k+1

c---
      If(Menu.eq.1) then
        x(k) = log(three) + 2.0*Dlog(xx)

      Else If(Menu.eq.2) then
        x(k) = exp(0.50*xx) / sqrt(three)

      Else If(Menu.eq.3) then
        x(k) = exp(-xx) / 3.0

      Else If(Menu.eq.4) then
        x(k) = exp(-xx)-2.0*xx

      Else If(Menu.eq.5) then
        x(k) = rl*xx*(1.0-xx)

      Else If(Menu.eq.6) then
        x(k) = 1.0-0.2*xx**3

      End If
        
c---
      If(Iextr.eq.1) then

       If(k.ge.2) then 
        xh = x(k) -(x(k)-x(k-1))**2/(x(k)-2.0*x(k-1)+x(k-2))
        write (6,100) k,x(k),xh
       Else
        write (6,100) k,x(k)
       End If

      Else

        write (6,100) k,x(k)
        write (3,100) k,x(k)

      End If
c---
      If(i.lt.10) Go to 98

      i = 0
      write (6,*)
      write (6,*) " Continue the iterations ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      Read  (5,*) Icon

      If(Icon.eq.1) Go to 98

c---
      Go to 95

c---
c Done
c---

  99  Continue

      write (3,100) Null
      Close (3)

 100  Format (1x,i3,2(1x,f15.10))

      Stop
      End
