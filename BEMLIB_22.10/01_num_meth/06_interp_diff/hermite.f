      program hermite_dr

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

c------------------------------------
c Mono-variate Hermite interpolation
c for specified values of a function
c and its derivatives
c by two methods
c
c SYMBOLS:
c --------
c
c f(i,0) value of function at ith datum point
c f(i,j) value of jth derivative of the function 
c        at ith datum point
c------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(50),f(50,0:4),m(10)
      Dimension q(50,0:4)
      Dimension t(80),c(0:50)

c--------------
c read the data
c--------------

      open (8,file="hermite.dat")

        read (8,*) n1
        n = n1-1
        Do i=1,n1
          read (8,*) k,m(i),x(i),(f(i,k-1),k=1,m(i))
        End Do

      close (8)

      write (6,*)
      write (6,*) " Data points"
      write (6,*) " -----------"
      write (6,*)

      Do i=1,n1
        write (6,100) i,m(i),x(i),(f(i,k-1),k=1,m(i))
      End Do

 98   Continue

      write (6,*)
      write (6,*) " Choose the interpolation method"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " 1 for the recursive method"
      write (6,*) "   algorithms (6.7.11), (6.7.12)"
      write (6,*) " 2 for the difference-table method"
      write (6,*) " ---------------------------------"

      read  (5,*) method

      If(method.eq.0) Go to 99

      write (6,*)
      write (6,*) " Please enter the value of x for interpolation"
      write (6,*) " ---------------------------------------------"
      read  (5,*) xint

c-------------------------
      If(method.eq.1) then
c-------------------------

      call Hermite_card (n,m,x,xint,q)

      fint = 0.0D0           ! loop (6.7.9)

      Do i=1,n1
        Do k = 0,m(i)-1
          fint = fint + f(i,k)*q(i,k)
        End Do
      End Do

      write (6,*)
      write (6,105) fint
      write (6,*)

c-------------------------
      Else If(method.eq.2) then
c-------------------------

      L = m(1)

      Do i=2,n1
       L = L+m(i)
      End Do

      L = L-1
     
      call Hermite_coeff (n,m,L,x,f,xint,t,c)
 
      write (6,*)
      write (6,*) " coefficients"
      write (6,*) " ------------"
      write (6,*)

      Do i=0,L
       write (6,101) i,c(i)
      End Do

      fint = c(L)           ! Newton like summation (6.7.14)
      Do i = L,1,-1
       fint = fint*(xint-t(i))+c(i-1)
      End Do

      write (6,*)
      write (6,105) fint
      write (6,*)

c-----------
      End If
c-----------

      Go to 98

c-----
c Done
c-----

 99   Continue

 100  Format(1x,i3,1x,i3,8(1x,f10.5))
 101  Format(1x,i3,8(1x,f15.10))
 105  Format("Interpolated value =",f15.10)

      Stop
      End
