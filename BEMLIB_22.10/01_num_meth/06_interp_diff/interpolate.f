      program interpolate

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c---------------------------------------------

c----------------------------------------
c Mono-variate polynomial interpolation
c by several methods 
c----------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension x(50),f(50)
      Dimension a(50,50),c(0:05,0:10,10)

c--------------
c read the data
c--------------

      open(3,file="interpolate.dat")

         read (3,*) n1
         n = n1-1
         Do i=1,n1
           read (3,*) k,x(i),f(i)
         End Do

      close (3)

 98   Continue

c-----------------
c display the data
c-----------------

      write (6,*)
      write (6,*) " Data points:"
      write (6,*) " ------------"
      write (6,*)

      Do i=1,n1
        write (6,100) i,x(i),f(i)
      End Do

c---
c proceed to interpolate
c---

      write (6,*) 
      write (6,*) " Choose the interpolation method"
      write (6,*) 
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1  for Lagrange "
      write (6,*) "       10 for inverse Lagrange "
      write (6,*) "       2  for Neville "
      write (6,*) "       3  for Aitken  "
      write (6,*) "       4  for Newton"
      write (6,*) "       5  for Fornberg"
      write (6,*) " -------------------------------"
      read  (5,*) method

      If(method.eq.0) Go to 99

      If(method.ne.10) then
        write (6,*)
        write (6,*) " Enter the value of x for interpolation"
        write (6,*) " --------------------------------------"
        read  (5,*) xint
      Else
        xint = 0.0
      End If

c-------------------------
      If(method.eq.1) then
c-------------------------

         call Lagrange (xint,n,x,f,fint)

         write (6,*)
	 write (6,105) fint
         write (6,*)

c-------------------------
      Else If(method.eq.10) then
c-------------------------

         Do i=1,n1
           save = x(i)
           x(i) = f(i)
           f(i) = save
         End Do

         call Lagrange (xint,n,x,f,fint)

         write (6,*)
         write (6,105) fint
         write (6,*)

c-------------------------
      Else If(method.eq.2) then
c-------------------------

         call Neville (xint,n,x,f,fint,a)

         write (6,*)
         write (6,105) fint
         write (6,*)
         write (6,*) "Neville Table"
         write (6,*) "------------"
         write (6,*)

         Do i=1,n1
           write (6,102) (a(i,j),j=1,n+2-i)
         End Do

c-------------------------
      Else If(method.eq.3) then
c-------------------------

        call Aitken (xint,n,x,f,fint,a)

         write (6,*)
         write (6,105) fint
         write (6,*)

         write (6,*) "Aitken Table"
         write (6,*) "------------"
         write (6,*)

         Do i=1,n1
          write (6,102) (a(i,j),j=1,i)
         End Do

c-------------------------
      Else If(method.eq.4) then
c-------------------------

         call Newton (xint,n,x,f,fint,a)

         write (6,*)
         write (6,105) fint
         write (6,*)
         write (6,*) "Newton Table"
         write (6,*) "------------"
         write (6,*)

         Do i=1,n1
          write (6,102) (a(i,j),j=1,n+2-i)
         End Do

c-------------------------
      Else If(method.eq.5) then
c-------------------------

      write (6,*)
      write (6,*) " Will interpolate and compute up to the"
      write (6,*) "         mth derivative; please enter m"
      write (6,*) " --------------------------------------"
      read  (5,*) m

      call fornberg (xint,n,x,m,c)

         write (6,*)
         Do k=0,m
           sum = 0.0
           Do i=1,n1
             sum = sum + c(k,n1,i)*f(i)
           End Do
           If(k.eq.0) write (6,105) sum
           If(k.eq.1) write (6,106) sum
           If(k.eq.2) write (6,107) sum
           If(k.eq.3) write (6,108) sum
           If(k.eq.4) write (6,109) sum
           If(k.eq.5) write (6,110) sum
           If(k.eq.6) write (6,111) sum
           If(k.eq.7) write (6,112) sum
           If(k.eq.8) write (6,113) sum
         End Do

c-----------
      End If
c-----------

      Go to 98     ! return to repeat

c------
c  Done
c------

 99   Continue

 100  Format (1x,i3,2(1x,f10.5))
 102  Format (20(1x,f10.5))
 105  Format ("Interpolated value =",f15.10)
 106  Format ("First   derivative =",f15.10)
 107  Format ("Second  derivative =",f15.10)
 108  Format ("Third   derivative =",f15.10)
 109  Format ("Fourth  derivative =",f15.10)
 110  Format ("Fifth   derivative =",f15.10)
 111  Format ("Sixth   derivative =",f15.10)
 112  Format ("Seventh derivative =",f15.10)
 113  Format ("Eighth  derivative =",f15.10)

      Stop
      End
