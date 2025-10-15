      program fp2_iter

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------
c  Fixed-point iterations for two equations
c
c  SYMBOLS:
c  -------
c
c  Iask: number of iterations before pausing
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (Iask=20)

c------
c input
c------

 97   Continue

      write (6,*) 
      write (6,*) "   SYSTEM MENU"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 for system (4.1.4) of the text"
      write (6,*) " 2 for system (4.1.21) of the text"
      write (6,*) " ---------------------------------"
      Read  (5,*) menu

      If(menu.eq.0) Go to 99

c----------------
c Initial guesses
c----------------

      If(menu.eq.1) then

       write (6,*) 
       write (6,*) " Enter the initial guesses for x1, x2"
       write (6,*) " ------------------------------------"
       read  (5,*) x1,x2

      Else If(menu.eq.2) then

       write (6,*) 
       write (6,*) " Enter the initial guesses for x1, x2"
       write (6,*) 
       write (6,*) " 0.5 and 0.02 recommended"
       write (6,*) " ------------------------"
       read  (5,*) x1,x2

      End If

c---------------------------------
c perform the iterations
c
c will pause every Iask iterations
c---------------------------------

      k = 0   ! global iteration counter
      i = 0   ! local counter

  98  Continue

      k = k+1
      i = i+1

c-----------------------
      If(Menu.eq.1) then
c-----------------------

        x1n = 2.0*exp(-x1-x2)
        x2n = 0.1*exp(-x2)/x1

        x1 = x1n
        x2 = x2n

c-----------------------
      Else If(Menu.eq.2) then
c-----------------------

c---
c a sleuth of bad choices
c---
c       x1n = 2.0*(0.312*(1.0-x1-x2)**2/(0.50*x1+x2)+x2)
c       x1n = 2.0*(0.312*(1.0-x1-x2)**2/(0.50*x1-x2)+x2)
c       x2n = 0.480*(1.0-x1-x2)*(0.5*x1-x2)/x2-0.5*x1
c       x2n = 1.0-x1-x2*(0.5*x1+x2)/(0.480*(0.5*x1-x2))
c       tmp = (0.5*x1-x2)*(0.5*x1+x2)/0.312
c       x1n = 1.0-x2-sqrt(tmp)
c       x2n = 0.480*(1.0-x1-x2)*(0.5*x1-x2)/(0.5*x1+x2)

c---
c something nice
c
c solve first equation for x1 inder the square
c
c solve second equation for x2 on the left-hand side
c and then modify a bit
c---

        tmp = (0.5*x1-x2)*(0.5*x1+x2)/0.312

        x1n = 1.0-x2-sqrt(tmp)
        x2n = (0.480*(1.0-x1-x2)*(0.5*x1-x2)/(0.5*x1+x2)+x2)/2.0

        x1 = x1n
        x2 = x2n

c-----------
      End If
c-----------

      write (6,100) k,x1,x2

      If(i.lt.Iask) Go to 98

      i = 0

      write (6,*) 
      write (6,*) " Continue the iterations ?"
      write (6,*) 
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      read  (5,*) Icon

      If(Icon.eq.1) Go to 98

c--------------------------------------------
c Confirm solution by computing functions f_i
c--------------------------------------------

c-----------------------
      If(Menu.eq.1) then
c-----------------------

        f1 = x1*exp(x1+x2)-2.0
        f2 = x1*x2 - 0.10*exp(-x2)

        write (6,*)
        write (6,*) "Functions f1 and f2 at the root:"
        write (6,*)
        write (6,101) f1,f2

c----------------------------
      Else If(Menu.eq.2) then
c----------------------------

        f1 = (0.5*x1-x2)*(0.50*x1+x2)-0.312*(1.0-x1-x2)**2
        f2 = x2*(0.5*x1+x2)-0.480*(1.0-x1-x2)*(0.5*x1-x2)

        write (6,*)
        write (6,*) "Functions f1 and f2 at the root:"
        write (6,*)
        write (6,101) f1,f2

c-----------
      End If
c-----------

c-----------------
c return to repeat
c-----------------

      Go to 97

c-----
c Done
c-----

  99  Continue

 100  Format (1x,i3,2(1x,f15.10))
 101  Format (2(2x,f15.10))

      Stop
      End
