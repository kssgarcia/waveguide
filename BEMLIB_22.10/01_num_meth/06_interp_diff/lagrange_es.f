      program lagrange_es

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

c---------------------------------------------------
c Testing interpolation formulae (6.2.16) - (6.2.18)
c---------------------------------------------------

      Implicit Double Precision (a-h,o,p,r-z)

      Dimension Y(-200:200),A(-200:200)
      Integer q

 98   Continue

      write (6,*)
      write (6,*) " Please enter N"
      write (6,*) " 0 to quit"
      write (6,*) "----------" 
      read  (5,*) n

      If(n.eq.0) Go to 99

      write (6,*) " Please enter h"
      write (6,*) "---------------" 
      read  (5,*) h

      write (6,*) " Please enter x"
      write (6,*) "---------------" 
      read  (5,*) x

c---
c prepare
c---

      Na = N-1
      N1 = N+1

      index = mod(N,2)

      If(index.eq.0) then

        write (6,*)
        write (6,*) " N is even"
        write (6,*)

        Nh   = N/2
        Imin = -Nh
        Imax =  Nh

      Else

        write (6,*)
        write (6,*) " N is odd"
        write (6,*)

        Imin = -na/2
        Imax =  n1/2

      End If

c-------------------------
c Generate data points
c for a quadratic function
c-------------------------

       Do i=Imin,Imax
        xn   = i*h
        y(i) = xn**2
       End Do

c--------
c prepare
c--------

      x0 = 0.0D0
      p  = (x-x0)/h

      yint = 0.0D0

c-------------------
      If(Index.eq.0) then  ! N is even
c-------------------

      Do i=Imin,Imax

        A(i) = 1.0
        Do j=1,nh+i
         A(i) = A(i)/j
        End Do

        Do j=1,nh-i
         A(i) = A(i)/j
        End Do

        Do j=1,nh+i
         A(i) = -A(i)
        End Do

        Do q = 0,n
         A(i) = A(i)*(p-q+nh)
         if(q.gt.20) stop
        End Do

        A(i) = A(i)/(p-i)
        yint = yint + A(i)*Y(i)

      End Do

c---------
      Else    ! N is odd
c---------

       Do i=Imin,Imax

         A(i) = 1.0
         Do j=1,na/2+i
          A(i) =  A(i)/(j-1.0+1.0)
         End Do

         Do j=1,n1/2-i
          A(i) =  A(i)/(j-1.0+1.0)
         End Do

         Do j=1,n1/2+i
          A(i) = - A(i)
         End Do

         Do q = 1,n1
          A(i) = A(i)*(p-q+N1/2)
         End Do
          A(i) = A(i)/(p-i)

         yint = yint + A(i)*y(i)

       End Do

c-----------
      End If
c-----------

c---
c confirm
c---

      yexact = x**2
      write (6,100) yint,yexact

      Go to 98
c-----
c Done
c-----

  99  Continue

 100  Format (2(1x,f15.10))

      Stop
      End
